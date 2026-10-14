#!/usr/bin/env python3
"""Ideal motion-capture-like feedback for the explicit agile simulation only."""
import numpy as np
import rospy
from gazebo_msgs.msg import ModelStates
from nav_msgs.msg import Odometry

rospy.init_node('flip_ground_truth_odometry')
model=rospy.get_param('~model','iris')
latest=None
publisher=rospy.Publisher('/flip/ground_truth/odom',Odometry,queue_size=10)
def callback(msg):
    global latest
    if model not in msg.name:return
    index=msg.name.index(model)
    pose=msg.pose[index]; twist=msg.twist[index];q=pose.orientation
    w,x,y,z=q.w,q.x,q.y,q.z
    rotation=np.array([[1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y)],
                       [2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x)],
                       [2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y)]])
    odom=Odometry();odom.header.stamp=rospy.Time.now();odom.header.frame_id='map';odom.child_frame_id='base_link'
    odom.pose.pose=pose
    v=rotation.T@np.array([twist.linear.x,twist.linear.y,twist.linear.z])
    rate=rotation.T@np.array([twist.angular.x,twist.angular.y,twist.angular.z])
    odom.twist.twist.linear.x,odom.twist.twist.linear.y,odom.twist.twist.linear.z=v
    odom.twist.twist.angular.x,odom.twist.twist.angular.y,odom.twist.twist.angular.z=rate
    latest=odom
rospy.Subscriber('/gazebo/model_states',ModelStates,callback,queue_size=1)
rate=rospy.Rate(100)
while not rospy.is_shutdown():
    if latest is not None:publisher.publish(latest)
    rate.sleep()
