#!/usr/bin/env python3
"""Evaluate 3-D helix tracking and each roll inversion against Gazebo truth.

Ground truth is used offline only. Fix ENU yaw/translation before takeoff,
then retain that frame throughout the flight; never fit the tracked path.
"""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Register 3-D projection on Matplotlib 3.1.


def rmse(values):
    return float(np.sqrt(np.mean(np.asarray(values)**2)))


def interpolate(data, times, columns):
    return np.column_stack([np.interp(times,data[:,0],data[:,j]) for j in columns])


def quaternion_at(data, times, start_column):
    q=data[:,start_column:start_column+4].copy()
    q/=np.linalg.norm(q,axis=1)[:,None]
    for k in range(1,len(q)):
        if q[k]@q[k-1]<0:
            q[k]*=-1
    result=np.column_stack([np.interp(times,data[:,0],q[:,j]) for j in range(4)])
    return result/np.linalg.norm(result,axis=1)[:,None]


def matrices(q):
    w,x,y,z=q.T
    return np.stack([1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y),
                     2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x),
                     2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)],axis=-1).reshape(-1,3,3)


def euler_xyz(rot):
    # Display conventional roll/pitch/yaw only; SO(3) errors above do not use
    # Euler coordinates. This roll-dominant experiment stays away from pitch 90.
    return np.rad2deg(np.column_stack([np.arctan2(rot[:,2,1],rot[:,2,2]),
                                     np.arctan2(-rot[:,2,0],np.hypot(rot[:,0,0],rot[:,1,0])),
                                     np.arctan2(rot[:,1,0],rot[:,0,0])]))


def yaw(rot):
    return np.arctan2(rot[:,1,0],rot[:,0,0])


def angle_error(actual,reference):
    relative=np.einsum('nji,njk->nik',actual,reference)
    return np.rad2deg(np.arccos(np.clip((np.trace(relative,axis1=1,axis2=2)-1)/2,-1,1)))


def continuous_roll(rot):
    return np.rad2deg(np.unwrap(np.arctan2(-rot[:,1,2],rot[:,2,2])))


def analyze(output, metadata):
    s=json.loads((output/'summary.json').read_text())
    meta=json.loads(metadata.read_text())
    events=json.loads((output/'events.json').read_text())
    read=lambda name:np.loadtxt(output/name,delimiter=',',skiprows=1)
    r,g,o,c,imu=[read(name) for name in ['references.csv','ground_truth.csv','odometry.csv','commands.csv','imu.csv']]
    first_arm=next(row[0] for row in events['controller_logs'] if 'HOVER --> TAKEOFF' in row[2])
    rest=(o[:,0]>=first_arm-.5)&(o[:,0]<first_arm-.05)
    if rest.sum()<6:
        raise ValueError('Insufficient static pre-takeoff alignment data')
    rest_t=o[rest,0]
    diff=yaw(matrices(quaternion_at(g,rest_t,4)))-yaw(matrices(quaternion_at(o,rest_t,7)))
    alignment_yaw=np.arctan2(np.sin(diff).mean(),np.cos(diff).mean())
    co,si=np.cos(alignment_yaw),np.sin(alignment_yaw)
    alignment=np.array([[co,-si,0],[si,co,0],[0,0,1]])
    offset=np.median(interpolate(g,rest_t,range(1,4))-o[rest,1:4]@alignment.T,axis=0)
    active=(r[:,0]<r[0,0]+meta['duration_s'])
    r=r[active];stamps=r[:,0];t=stamps-stamps[0]
    if not (g[0,0]<=stamps[0] and g[-1,0]>=stamps[-1] and o[-1,0]>=stamps[-1]):
        raise ValueError('State logs do not cover the reference samples')
    actual_world=interpolate(g,stamps,range(1,4))
    actual_local=(actual_world-offset)@alignment
    estimated_local=interpolate(o,stamps,range(1,4))
    actual_rotation=matrices(quaternion_at(g,stamps,4))
    local_rotation=np.einsum('ij,njk->nik',alignment.T,actual_rotation)
    estimated_rotation=matrices(quaternion_at(o,stamps,7))
    reference_rotation=matrices(r[:,4:8])
    error=actual_local-r[:,1:4]
    position_error=np.linalg.norm(error,axis=1)
    estimator_position_error=np.linalg.norm(estimated_local-actual_local,axis=1)
    body_velocity=interpolate(o,stamps,range(4,7))
    estimated_world_rotation=np.einsum('ij,njk->nik',alignment,estimated_rotation)
    estimated_velocity=np.einsum('nij,nj->ni',estimated_world_rotation,body_velocity)
    velocity_error=np.linalg.norm(estimated_velocity-interpolate(g,stamps,range(8,11)),axis=1)
    feedback_angle=angle_error(local_rotation,estimated_rotation)
    reference_angle=angle_error(local_rotation,reference_rotation)
    actual_roll=continuous_roll(local_rotation)
    reference_roll=continuous_roll(reference_rotation)
    loops=[]
    for k in range(meta['loops']):
        begin=meta['entry_s']+k*meta['loop_period_s'];end=begin+meta['loop_period_s']
        if t[-1]<end-.03:
            break
        mask=(t>=begin)&(t<end)
        top=begin+meta['loop_period_s']/4
        near=abs(t-top)<.08
        turn=float(np.interp(end,t,actual_roll)-np.interp(begin,t,actual_roll))
        rzz=float(np.mean(local_rotation[near,2,2]))
        loops.append(dict(loop=k+1,roll_revolution_deg=turn,top_body_z_world_z=rzz,
                          completed=bool(abs(turn)>330 and rzz<-.7),
                          position_rmse_m=rmse(position_error[mask]),maximum_position_error_m=float(position_error[mask].max()),
                          reference_attitude_rmse_deg=rmse(reference_angle[mask])))
    s.update(trajectory_metadata=meta,ground_truth_position_rmse_m=rmse(position_error),
             ground_truth_position_max_error_m=float(position_error.max()),
             ground_truth_position_error_p95_m=float(np.percentile(position_error,95)),
             ground_truth_axis_rmse_m=np.sqrt(np.mean(error**2,axis=0)).tolist(),
             ground_truth_frame_alignment='fixed ENU yaw and translation from static pre-takeoff window',
             frame_alignment=dict(yaw_world_from_local_deg=float(np.rad2deg(alignment_yaw)),
                                  world_from_local_translation_m=offset.tolist(),samples=int(rest.sum())),
             estimator_validation=dict(position_error_rmse_m=rmse(estimator_position_error),velocity_error_rmse_m_s=rmse(velocity_error),
                                       attitude_error_rmse_deg=rmse(feedback_angle)),
             flip_validation=dict(loops=loops,all_loops_passed=len(loops)==meta['loops'] and all(x['completed'] for x in loops),
                                  attitude_error_rmse_deg=rmse(reference_angle),attitude_error_max_deg=float(reference_angle.max())))
    odom_mask=(o[:,0]>=stamps[0])&(o[:,0]<stamps[0]+meta['duration_s'])
    odom_times=o[odom_mask,0]
    s['feedback_mean_rate_hz']=float((len(odom_times)-1)/(odom_times[-1]-odom_times[0]))
    (output/'summary.json').write_text(json.dumps(s,indent=2))
    np.savetxt(output/'helix_analysis.csv',np.column_stack([t,r[:,1:4],actual_local,estimated_local,position_error,estimator_position_error,velocity_error,feedback_angle,reference_angle,reference_roll,actual_roll]),
               delimiter=',',comments='',header='ref_time,ref_x,ref_y,ref_z,true_x,true_y,true_z,est_x,est_y,est_z,position_error,estimator_position_error,estimator_velocity_error,estimator_attitude_error_deg,reference_attitude_error_deg,reference_roll_deg,actual_roll_deg')

    # Layout follows the quantities in paper Fig. 7(a). All errors are truth-based.
    fig=plt.figure(figsize=(18,9),constrained_layout=True)
    grid=fig.add_gridspec(3,5)
    path=fig.add_subplot(grid[:,0],projection='3d')
    circle=(t>=meta['entry_s'])&(t<=meta['entry_s']+meta['loops']*meta['loop_period_s'])
    path.plot(r[circle,1],r[circle,2],r[circle,3],color='tab:orange',linestyle='--',label='Reference')
    path.plot(*actual_local[circle].T,color='tab:blue',label='Gazebo truth')
    path.set(xlabel='X (m)',ylabel='Y (m)',zlabel='Z (m)',title='Helix circle segments')
    path.view_init(elev=18,azim=-55)
    if hasattr(path,'set_box_aspect'):
        path.set_box_aspect(np.ptp(r[circle,1:4],axis=0))
    path.legend(fontsize=8)
    for j,label in enumerate(['X','Y','Z']):
        a=fig.add_subplot(grid[j,1]);a.plot(t,actual_local[:,j],label='Truth');a.plot(t,r[:,j+1],'--',label='Reference')
        a.set(ylabel=label+' (m)');a.grid(True)
        if j==0:a.set_title('Position');a.legend(fontsize=8)
        if j==2:a.set_xlabel('Reference time (s)')
    pe=fig.add_subplot(grid[:,2])
    pe.plot(t,error[:,0],label='X');pe.plot(t,error[:,1],label='Y');pe.plot(t,error[:,2],label='Z')
    pe.set(xlabel='Reference time (s)',ylabel='Actual - reference (m)',title='Position errors');pe.grid(True);pe.legend()
    true_euler=euler_xyz(local_rotation)
    ref_euler=euler_xyz(reference_rotation)
    for j,label in enumerate(['Roll','Pitch','Yaw']):
        a=fig.add_subplot(grid[j,3]);a.plot(t,true_euler[:,j],label='Truth');a.plot(t,ref_euler[:,j],'--',label='Reference')
        a.set(ylabel=label+' (deg)');a.grid(True)
        if j==0:a.set_title('Attitude (wrapped)');a.legend(fontsize=8)
        if j==2:a.set_xlabel('Reference time (s)')
    inputs=fig.add_subplot(grid[:2,4]);rate=fig.add_subplot(grid[2,4])
    mapping=s.get('thrust_model')
    if mapping is None:
        manifest=output.parent/'manifest.json'
        if not manifest.exists():
            raise ValueError('Missing thrust model: use the current harness to record controller calibration')
        mapping=json.loads(manifest.read_text())['thrust_model']
    cm=(c[:,0]>=stamps[0])&(c[:,0]<stamps[0]+meta['duration_s'])
    im=(imu[:,0]>=stamps[0])&(imu[:,0]<stamps[0]+meta['duration_s'])
    thrust=mapping['gain']*(c[cm,1]+mapping['offset'])**mapping['power']/9.81
    inputs.plot(imu[im,0]-stamps[0],np.linalg.norm(imu[im,4:7],axis=1)/9.81,label='IMU specific force')
    inputs.plot(c[cm,0]-stamps[0],thrust,'--',label='Command')
    inputs.set(ylabel='Thrust acceleration (g)',title='Control inputs');inputs.legend(fontsize=8);inputs.grid(True)
    rate.plot(imu[im,0]-stamps[0],np.rad2deg(np.linalg.norm(imu[im,1:4],axis=1)),label='IMU')
    rate.plot(c[cm,0]-stamps[0],np.rad2deg(np.linalg.norm(c[cm,2:5],axis=1)),'--',label='Command')
    rate.set(xlabel='Reference time (s)',ylabel='Angular rate (deg/s)');rate.grid(True);rate.legend(fontsize=8)
    fig.suptitle('OMMPC / PX4 100 Hz synchronized feedback: five roll flips',fontsize=15)
    fig.savefig(output/'paper_style_tracking.png',dpi=170);plt.close(fig)
    fig,ax=plt.subplots(2,2,figsize=(13,8),constrained_layout=True)
    ax[0,0].plot(t,reference_roll,label='Reference');ax[0,0].plot(t,actual_roll,label='Truth')
    ax[0,0].set(xlabel='Reference time (s)',ylabel='Unwrapped roll (deg)',title='Continuous roll revolutions');ax[0,0].legend()
    ax[0,1].plot(t,position_error);ax[0,1].set(xlabel='Reference time (s)',ylabel='Position error (m)',title='Actual tracking error')
    ax[1,0].plot(t,local_rotation[:,2,2],label='Actual body Z / world Z')
    ax[1,0].plot(t,reference_rotation[:,2,2],'--',label='Reference')
    ax[1,0].set(xlabel='Reference time (s)',ylabel='Z projection',title='Downward thrust at each top');ax[1,0].legend()
    ax[1,1].plot(t,feedback_angle,label='Feedback vs truth');ax[1,1].plot(t,reference_angle,label='Truth vs reference')
    ax[1,1].set(xlabel='Reference time (s)',ylabel='Angle error (deg)',title='Feedback and control errors');ax[1,1].legend()
    for a in ax.flat:a.grid(True)
    fig.savefig(output/'flip_validation.png',dpi=160);plt.close(fig)
    print(json.dumps({k:s[k] for k in ['tracking_completed','ground_truth_position_rmse_m','ground_truth_position_max_error_m','feedback_mean_rate_hz','estimator_validation','flip_validation']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('output',type=Path)
    parser.add_argument('--metadata',type=Path,required=True)
    args=parser.parse_args()
    analyze(args.output,args.metadata)
