#!/usr/bin/env python3
"""Forward atomic PX4 ODOMETRY with verified frames and correct body twist.

Uses only FCU data. No Gazebo feedback, interpolation, or state extrapolation.
MAVROS supplies the FCU-synchronized timestamp, position and quaternion from
one ODOMETRY packet. Raw MAVLink supplies its frame IDs (lost by this plugin).
"""
import copy
import numpy as np
import rospy
from nav_msgs.msg import Odometry
from mavros_msgs.msg import Mavlink
from mavros.mavlink import convert_to_bytes
from pymavlink.dialects.v20 import common
from px4_odometry_math import body_velocity_from_mavros, LOCAL_NED

class SynchronizedOdometry:
    def __init__(self):
        self.frames = None
        self.decoder = common.MAVLink(None)
        self.publisher = rospy.Publisher('/ommpc_controller/px4_synced_odom',Odometry,queue_size=10)
        self.source_system = rospy.get_param('~source_system',1)
        rospy.Subscriber('/mavlink/from',Mavlink,self.raw_callback,queue_size=20,tcp_nodelay=True)
        rospy.Subscriber('/mavros/odometry/in',Odometry,self.odom_callback,queue_size=5,tcp_nodelay=True)

    def raw_callback(self,message):
        if message.msgid != 331 or message.sysid != self.source_system or message.compid != 1:
            return
        if message.framing_status != Mavlink.FRAMING_OK:
            return
        try:
            packet = self.decoder.decode(convert_to_bytes(message))
            if packet.estimator_type != common.MAV_ESTIMATOR_TYPE_AUTOPILOT:
                return
            frames = (packet.frame_id, packet.child_frame_id)
            if self.frames != frames:
                rospy.loginfo('PX4 atomic ODOMETRY frames: pose=%d velocity=%d',*frames)
            self.frames = frames
        except (ValueError, TypeError) as error:
            rospy.logerr_throttle(1,'Cannot decode PX4 odometry: '+str(error))

    def odom_callback(self,message):
        frames = self.frames
        if frames is None:
            return
        if frames[0] != LOCAL_NED:
            rospy.logerr_throttle(1,'Unsupported odometry pose frame: '+str(frames[0]))
            return
        q=message.pose.pose.orientation;v=message.twist.twist.linear
        p=message.pose.pose.position;w=message.twist.twist.angular
        if not np.isfinite([p.x,p.y,p.z,w.x,w.y,w.z]).all():
            return
        try:
            velocity,transform=body_velocity_from_mavros([q.w,q.x,q.y,q.z],[v.x,v.y,v.z],frames[1])
        except ValueError as error:
            rospy.logerr_throttle(1,str(error))
            return
        result=copy.deepcopy(message)
        result.twist.twist.linear.x,result.twist.twist.linear.y,result.twist.twist.linear.z=velocity
        norm=np.linalg.norm([q.w,q.x,q.y,q.z])
        result.pose.pose.orientation.w=q.w/norm;result.pose.pose.orientation.x=q.x/norm
        result.pose.pose.orientation.y=q.y/norm;result.pose.pose.orientation.z=q.z/norm
        # Angular rates are already body FLU. Correct covariance when known.
        covariance=np.asarray(message.twist.covariance).reshape(6,6)
        if np.isfinite(covariance).all():
            rotation=np.eye(6);rotation[:3,:3]=transform
            result.twist.covariance=(rotation@covariance@rotation.T).reshape(-1).tolist()
        self.publisher.publish(result)

if __name__=='__main__':
    rospy.init_node('px4_synchronized_odometry')
    SynchronizedOdometry()
    rospy.spin()
