import time, json, threading, pathlib, csv, argparse
import numpy as np
import rospy
from nav_msgs.msg import Odometry
from gazebo_msgs.msg import ModelStates
from geometry_msgs.msg import PoseStamped
from mavros_msgs.msg import State, AttitudeTarget
from mavros_msgs.srv import SetMode, ParamGet, ParamSet, ParamPull, MessageInterval
from mavros_msgs.msg import ParamValue
from sensor_msgs.msg import Imu
from rosgraph_msgs.msg import Log
from dynamic_reconfigure.client import Client

parser=argparse.ArgumentParser()
parser.add_argument('--output',type=pathlib.Path,required=True)
parser.add_argument('--trajectory',type=pathlib.Path,required=True)
parser.add_argument('--agile',action='store_true')
parser.add_argument('--odometry',default='/mavros/local_position/odom')
parser.add_argument('--feedback-mode',choices=['baseline','high_rate','synchronized'],default='baseline')
parser.add_argument('--feedback-rate',type=float,default=100.)
args=parser.parse_args()
if args.feedback_rate<=0:parser.error('--feedback-rate must be positive')
OUT=args.output; OUT.mkdir(parents=True,exist_ok=True)
REF=args.trajectory
imu=[]; saved_params={}
odom_arrivals=[]
ref=np.loadtxt(REF)
references=[]; ground_truth=[]; odoms=[]; cmds=[]; logs=[]; states=[]; state=None; tracking_start=None; done=False

def reference_cb(m):
 p=m.pose.position
 q=m.pose.orientation
 references.append([m.header.stamp.to_sec(),p.x,p.y,p.z,q.w,q.x,q.y,q.z])
def odom_cb(m):
 odom_arrivals.append([rospy.get_time(),m.header.stamp.to_sec()])
 p=m.pose.pose.position; v=m.twist.twist.linear
 q=m.pose.pose.orientation
 odoms.append([m.header.stamp.to_sec(),p.x,p.y,p.z,v.x,v.y,v.z,q.w,q.x,q.y,q.z])
def ground_truth_cb(m):
 if 'iris' not in m.name: return
 index=m.name.index('iris');pose=m.pose[index];p=pose.position;q=pose.orientation
 v=m.twist[index].linear;w=m.twist[index].angular
 ground_truth.append([rospy.get_time(),p.x,p.y,p.z,q.w,q.x,q.y,q.z,v.x,v.y,v.z,w.x,w.y,w.z])
def state_cb(m):
 global state
 if state is None or (state.mode,state.armed,state.connected)!=(m.mode,m.armed,m.connected):
  states.append([rospy.get_time(),m.mode,m.armed,m.connected])
 state=m
def cmd_cb(m):
 w=m.body_rate
 cmds.append([m.header.stamp.to_sec(),m.thrust,w.x,w.y,w.z])
def log_cb(m):
 global tracking_start,done
 if 'ommpc_controller' in m.name:
  logs.append([m.header.stamp.to_sec(),m.level,m.msg])
  if 'HOVER --> POINTS' in m.msg and tracking_start is None:
   tracking_start=m.header.stamp.to_sec()
  if 'Numerical error' in m.msg and tracking_start is not None:
   done=True

def wait_until(check,seconds,label):
 deadline=time.monotonic()+seconds
 while time.monotonic()<deadline and not rospy.is_shutdown():
  if check(): return
  time.sleep(.1)
 raise RuntimeError('Timed out: '+label)

def report(msg): print(msg,flush=True)

def imu_cb(m):
 w=m.angular_velocity; a=m.linear_acceleration
 imu.append([m.header.stamp.to_sec(),w.x,w.y,w.z,a.x,a.y,a.z])
rospy.init_node('trajectory_tracking_test',anonymous=True)
rospy.Subscriber('/mavros/imu/data',Imu,imu_cb,queue_size=1000)
rospy.Subscriber(args.odometry,Odometry,odom_cb,queue_size=1000)
rospy.Subscriber('/mavros/state',State,state_cb)
rospy.Subscriber('/ommpc_controller/reference',PoseStamped,reference_cb,queue_size=1000)
rospy.Subscriber('/gazebo/model_states',ModelStates,ground_truth_cb,queue_size=1000)
rospy.Subscriber('/mavros/setpoint_raw/attitude',AttitudeTarget,cmd_cb,queue_size=1000)
rospy.Subscriber('/rosout_agg',Log,log_cb,queue_size=1000)
result={'trajectory':str(REF),'reference_rows':len(ref),'feedback_topic':args.odometry,'feedback_mode':args.feedback_mode}
client=None
try:
 wait_until(lambda: state is not None and state.connected and len(odoms)>20 and len(cmds)>20,35,'connected odometry and controller')
 client=Client('/ommpc_controller',timeout=15)
 result['thrust_model']=rospy.get_param('/ommpc_controller/thrust_model',None)
 set_mode=rospy.ServiceProxy('/mavros/set_mode',SetMode)
 report('SIM ready: '+str(states[-1]))
 if args.agile:
  rospy.wait_for_service('/mavros/param/pull',timeout=20)
  if not rospy.ServiceProxy('/mavros/param/pull',ParamPull)(True).success: raise RuntimeError('Parameter synchronization failed')
  get_param=rospy.ServiceProxy('/mavros/param/get',ParamGet)
  put_param=rospy.ServiceProxy('/mavros/param/set',ParamSet)
  for name,value in {'MC_ROLLRATE_K':0.65,'MC_PITCHRATE_K':0.65,'MC_YAWRATE_K':0.65,'MC_ROLLRATE_MAX':900.,'MC_PITCHRATE_MAX':900.}.items():
   old=get_param(name)
   if not old.success: raise RuntimeError('Cannot read '+name)
   saved_params[name]=old.value
   if not put_param(name,ParamValue(integer=0,real=value)).success: raise RuntimeError('Cannot set '+name)
  result['original_px4_parameters']={name:value.real for name,value in saved_params.items()}
  result['applied_px4_parameters']={'MC_ROLLRATE_K':.65,'MC_PITCHRATE_K':.65,'MC_YAWRATE_K':.65,'MC_ROLLRATE_MAX':900.,'MC_PITCHRATE_MAX':900.}
  delay=get_param('EKF2_GPS_DELAY')
  if not delay.success:raise RuntimeError('Cannot read GPS delay')
  result['ekf_gps_delay_ms']=delay.value.real
  report('Agile inner-loop configured; GPS delay='+str(delay.value.real)+' ms')
 if args.feedback_mode!='baseline':
  # These are MAVLink stream requests, not EKF sensor-rate changes. They persist
  # until SITL restarts, so frequency comparisons require fresh FCU instances.
  rospy.wait_for_service('/mavros/set_message_interval',timeout=10)
  set_interval=rospy.ServiceProxy('/mavros/set_message_interval',MessageInterval)
  result['message_interval_requests']={}
  for msgid in [32,31,331]:
   ok=set_interval(msgid,args.feedback_rate).success
   result['message_interval_requests'][str(msgid)]={'requested_hz':args.feedback_rate,'success':ok}
   if not ok:raise RuntimeError('Cannot change message interval '+str(msgid))
  report('Requested '+str(args.feedback_rate)+' Hz LOCAL_POSITION_NED, ATTITUDE_QUATERNION and ODOMETRY')
  time.sleep(1)
 for _ in range(10):
  r=set_mode(0,'OFFBOARD')
  time.sleep(.5)
  if state.mode=='OFFBOARD': break
 wait_until(lambda: state.mode=='OFFBOARD',10,'OFFBOARD')
 client.update_configuration({'takeoff_enabled':True})
 wait_until(lambda: any('TAKEOFF succeeded' in x[2] for x in logs),35,'automatic takeoff')
 report('TAKEOFF passed, z='+str(odoms[-1][3]))
 time.sleep(2)
 client.update_configuration({'command_or_hover':True})
 wait_until(lambda:tracking_start is not None,8,'POINTS mode')
 report('TRACKING started at simulation time '+str(tracking_start))
 deadline=time.monotonic()+180; last_print=0; bad_since=None
 while time.monotonic()<deadline:
  samples=sum(x[0]>=tracking_start for x in cmds)
  elapsed=rospy.get_time()-(references[0][0] if references else tracking_start)
  expected=np.array(references[-1][1:4]) if references else ref[0,:3]
  latest_error=float(np.linalg.norm(np.asarray(odoms[-1][1:4])-expected))
  bad_since=(bad_since if bad_since is not None else time.monotonic()) if latest_error>2.0 else None
  if bad_since is not None and time.monotonic()-bad_since>2.0:
   result['tracking_completed']=False; result['stop_reason']='position error above 2 m for 2 seconds'; result['tracking_end']=rospy.get_time(); break
  if odoms[-1][3]<0.2 and samples>100:
   result['tracking_completed']=False; result['stop_reason']='altitude below 0.2 m'; result['tracking_end']=rospy.get_time(); break
  if time.monotonic()-last_print>10:
   report('TRACKING samples='+str(samples)+'/'+str(len(ref))+' z='+str(round(odoms[-1][3],3))+' ref_time='+str(round(elapsed,3))+' error='+str(round(latest_error,3))+' mode='+state.mode)
   last_print=time.monotonic()
  if done or state.mode!='OFFBOARD' or not state.armed:
   result['tracking_completed']=False; result['stop_reason']='solver failure or mode/arming change'; result['tracking_end']=rospy.get_time(); break
  if elapsed>=(len(ref)-1)*.01:
   result['tracking_completed']=True; time.sleep(3); break
  time.sleep(.1)
 else:
  result['tracking_completed']=False; result['stop_reason']='wall-time timeout'
except Exception as e:
 result['error']=str(e);report('TEST ERROR: '+str(e))
finally:
 # Restore rate-controller parameters only after the vehicle is disarmed;
 # changing the inner loop during landing would invalidate the flight setup.
 try:
  if client and state is not None and state.armed:
   client.update_configuration({'command_or_hover':False})
   time.sleep(1)
   client.update_configuration({'land_enabled':True})
   wait_until(lambda:state is not None and not state.armed,35,'automatic landing and disarm')
   result['landing_disarmed']=True
 except Exception as e:
  result['landing_error']=str(e)
  try:
   result['fallback_mode_sent']=rospy.ServiceProxy('/mavros/set_mode',SetMode)(0,'AUTO.LAND').mode_sent
   wait_until(lambda:state is not None and not state.armed,35,'fallback AUTO.LAND disarm')
   result['fallback_disarmed']=True
  except Exception as fallback_error: result['fallback_error']=str(fallback_error)
 if state is not None and not state.armed and saved_params:
  result['px4_parameters_restored']=all(rospy.ServiceProxy('/mavros/param/set',ParamSet)(name,value).success for name,value in saved_params.items())
 np.savetxt(OUT/'odom_arrivals.csv',np.asarray(odom_arrivals),delimiter=',',header='arrival_sim_time,sample_sim_time',comments='')
 np.savetxt(OUT/'imu.csv',np.asarray(imu),delimiter=',',header='sim_time,wx,wy,wz,ax,ay,az',comments='')
 # Capture data after landing.
 np.savetxt(OUT/'ground_truth.csv',np.asarray(ground_truth),delimiter=',',header='sim_time,x_world,y_world,z_world,qw,qx,qy,qz,vx_world,vy_world,vz_world,wx_world,wy_world,wz_world',comments='')
 oa=np.asarray(odoms); ca=np.asarray(cmds)
 ra=np.asarray(references)
 np.savetxt(OUT/'references.csv',ra,delimiter=',',header='sim_time,x,y,z,qw,qx,qy,qz',comments='')
 np.savetxt(OUT/'odometry.csv',oa,delimiter=',',header='sim_time,x,y,z,vx_body,vy_body,vz_body,qw,qx,qy,qz',comments='')
 np.savetxt(OUT/'commands.csv',ca,delimiter=',',header='sim_time,thrust,wx,wy,wz',comments='')
 (OUT/'events.json').write_text(json.dumps({'controller_logs':logs,'states':states},indent=2))
 result['tracking_start']=tracking_start
 result['automatic_takeoff_passed']=any('TAKEOFF succeeded' in x[2] for x in logs)
 result['automatic_landing_passed']=any('LAND succeeded' in x[2] for x in logs)
 result['numerical_errors']=sum('Numerical error' in x[2] for x in logs)
 result['mpc_latency_messages']=[x[2] for x in logs if 'MPC Timing' in x[2]]
 if tracking_start is not None and len(oa) and len(ra) and len(ca):
  active=ca[(ca[:,0]>=ra[0,0])&(ca[:,0]<result.get('tracking_end',ra[0,0]+(len(ref)-1)*.01))]
  nearest=np.maximum(0,np.minimum(np.searchsorted(ra[:,0],active[:,0],side='right')-1,len(ra)-1))
  ref_time=ra[nearest,0]-ra[0,0]
  actual=np.column_stack([np.interp(active[:,0],oa[:,0],oa[:,j]) for j in (1,2,3)])
  desired=ra[nearest,1:4]; error=actual-desired; norm=np.linalg.norm(error,axis=1)
  result.update({'tracking_samples':len(active),'tracking_duration_sim_s':float(active[-1,0]-active[0,0]),'command_frequency_hz':float(1/np.mean(np.diff(active[:,0]))),'position_rmse_m':float(np.sqrt(np.mean(norm**2))),'axis_rmse_m':np.sqrt(np.mean(error**2,axis=0)).tolist(),'position_error_p95_m':float(np.percentile(norm,95)),'position_error_max_m':float(norm.max()),'final_position_error_m':float(norm[-1]),'thrust_saturation_fraction':float(np.mean(active[:,1]>=(.9499 if args.agile else .8999))),'bodyrate_xy_saturation_fraction':float(np.mean(np.any(np.abs(active[:,2:4])>=(13.999 if args.agile else 5.999),axis=1)))})
  segments=[]
  for start in np.arange(0,ref_time[-1],10):
   mask=(ref_time>=start)&(ref_time<start+10)
   if np.any(mask): segments.append({'reference_time_s':[float(start),float(min(start+10,ref_time[-1]))],'position_rmse_m':float(np.sqrt(np.mean(norm[mask]**2))),'max_error_m':float(norm[mask].max())})
  result['segments']=segments
  np.savetxt(OUT/'tracking.csv',np.column_stack([active[:,0],ref_time,desired,actual,error,norm]),delimiter=',',header='sim_time,reference_time,ref_x,ref_y,ref_z,x,y,z,err_x,err_y,err_z,error_norm',comments='')
  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  fig,ax=plt.subplots(2,2,figsize=(12,8))
  ax[0,0].plot(desired[:,0],desired[:,2],label='Reference');ax[0,0].plot(actual[:,0],actual[:,2],label='Actual',alpha=.8);ax[0,0].axis('equal');ax[0,0].set(xlabel='x (m)',ylabel='z (m)',title='Vertical trajectory (XZ plane)');ax[0,0].legend()
  ax[0,1].plot(ref_time,norm);ax[0,1].set(xlabel='Reference time (s)',ylabel='Position error (m)')
  ax[1,0].plot(ref_time,actual[:,2],label='Actual');ax[1,0].plot(ref_time,desired[:,2],label='Reference');ax[1,0].set(xlabel='Reference time (s)',ylabel='Altitude (m)');ax[1,0].legend()
  ax[1,1].plot(ref_time,active[:,1]);ax[1,1].set(xlabel='Reference time (s)',ylabel='Normalized thrust')
  for a in ax.flat:a.grid(True)
  fig.suptitle('OMMPC corrected vertical circle SITL test');fig.tight_layout(rect=(0,0,1,.95));fig.savefig(OUT/'tracking.png',dpi=160);plt.close(fig)
 (OUT/'summary.json').write_text(json.dumps(result,indent=2))
 report(json.dumps(result,indent=2))
