#!/usr/bin/env python3
"""Check each complete pitch revolution and inverted thrust at circle tops."""
import argparse, json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
parser=argparse.ArgumentParser()
parser.add_argument('output',type=Path)
parser.add_argument('--metadata',type=Path,required=True)
args=parser.parse_args(); out=args.output
meta=json.loads(args.metadata.read_text())
summary=json.loads((out/'summary.json').read_text())
r=np.loadtxt(out/'references.csv',delimiter=',',skiprows=1)
g=np.loadtxt(out/'ground_truth.csv',delimiter=',',skiprows=1)
i=np.loadtxt(out/'imu.csv',delimiter=',',skiprows=1)
c=np.loadtxt(out/'commands.csv',delimiter=',',skiprows=1)
t=r[:,0]-r[0,0]
q=np.column_stack([np.interp(r[:,0],g[:,0],g[:,j]) for j in range(4,8)])
q/=np.linalg.norm(q,axis=1)[:,None]
def axes(q):
 w,x,y,z=q.T
 return 2*(x*z+w*y),1-2*(x*x+y*y)
def pitch(q):
 x,z=axes(q)
 return np.unwrap(np.arctan2(x,z))*180/np.pi
rp=pitch(r[:,4:8]); ap=pitch(q)
att_error=np.rad2deg(2*np.arccos(np.clip(np.abs(np.sum(q*r[:,4:8],axis=1)),0,1)))
loops=[]
for k in range(meta['loops']):
 start=meta['entry_s']+k*meta['loop_period_s']; end=start+meta['loop_period_s']
 if t[-1]<end-.03: break
 rotation=float(np.interp(end,t,ap)-np.interp(start,t,ap))
 top=start+meta['loop_period_s']/4
 near=(abs(t-top)<.08)
 rzz=axes(q)[1]
 top_rzz=float(np.mean(rzz[near]))
 loops.append(dict(loop=k+1,pitch_revolution_deg=rotation,top_body_z_world_z=top_rzz,
                   completed=abs(rotation)>330 and top_rzz<-.7))
mask=t<=meta['duration_s']
active=t<=meta['duration_s']
world=np.column_stack([np.interp(r[:,0],g[:,0],g[:,j]) for j in range(1,4)])
odom=np.loadtxt(out/'odometry.csv',delimiter=',',skiprows=1)
# Estimate-based runs use an initial translation alignment; there was no full
# reference-frame orientation recording in early runs. Ideal feedback is exact.
initial_odom=np.array([np.interp(r[0,0],odom[:,0],odom[:,j]) for j in range(1,4)])
offset=world[0]-initial_odom if summary.get('feedback_topic')!='/flip/ground_truth/odom' else np.zeros(3)
world_error=np.linalg.norm(world-r[:,1:4]-offset,axis=1)
summary['ground_truth_position_rmse_m']=float(np.sqrt(np.mean(world_error[active]**2)))
summary['ground_truth_position_max_error_m']=float(world_error[active].max())
summary['ground_truth_frame_alignment']='direct' if summary.get('feedback_topic')=='/flip/ground_truth/odom' else 'initial translation, axes assumed aligned'
summary['flip_validation']={'loops':loops,'all_loops_passed':len(loops)==meta['loops'] and all(x['completed'] for x in loops),
 'attitude_error_rmse_deg':float(np.sqrt(np.mean(att_error[mask]**2))),
 'attitude_error_max_deg':float(att_error[mask].max())}
(out/'summary.json').write_text(json.dumps(summary,indent=2))
track=np.loadtxt(out/'tracking.csv',delimiter=',',skiprows=1)
fig,ax=plt.subplots(2,2,figsize=(12,8))
circle=active&(t>=meta['entry_s'])&(t<=meta['entry_s']+meta['loops']*meta['loop_period_s'])
ax[0,0].plot(r[circle,1]+offset[0],r[circle,3]+offset[2],label='Reference')
ax[0,0].plot(world[circle,0],world[circle,2],label='Actual (Gazebo)',alpha=.8)
top_indices=[np.argmin(abs(t-(meta['entry_s']+(k+.25)*meta['loop_period_s']))) for k in range(len(loops))]
if top_indices:
 bx,bz=axes(q)
 ax[0,0].quiver(world[top_indices,0],world[top_indices,2],bx[top_indices],bz[top_indices],
                angles='xy',scale_units='xy',scale=2.5,color='black',label='Thrust axis at tops')
ax[0,0].set_aspect('equal',adjustable='box')
ax[0,0].set(xlabel='X (m)',ylabel='Z (m)',title='Circle segments (XZ plane)');ax[0,0].legend(fontsize=9)
ax[0,1].plot(track[:,1],track[:,-1],label='Feedback');ax[0,1].plot(t[active],world_error[active],label='Gazebo truth',alpha=.7);ax[0,1].legend(); ax[0,1].set(xlabel='Time (s)',ylabel='Position error (m)')
ax[1,0].plot(t[active],rp[active],label='Reference');ax[1,0].plot(t[active],ap[active],label='Actual');ax[1,0].set(xlabel='Time (s)',ylabel='Continuous pitch (deg)');ax[1,0].legend()
ax[1,1].plot(i[:,0]-r[0,0],i[:,2],label='Measured pitch rate');ax[1,1].plot(c[:,0]-r[0,0],c[:,3],label='Command',alpha=.7)
ax[1,1].set(xlim=(0,meta['duration_s']),xlabel='Time (s)',ylabel='Body pitch rate (rad/s)');ax[1,1].legend()
for a in ax.flat:a.grid(True)
fig.suptitle('OMMPC inverted vertical loop SITL'); fig.tight_layout();fig.savefig(out/'flip_validation.png',dpi=160)
print(json.dumps(summary,indent=2))
