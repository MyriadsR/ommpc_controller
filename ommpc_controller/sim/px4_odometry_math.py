"""Frame correction for MAVROS 1.20.x incoming PX4 ODOMETRY twist.

This MAVROS plugin rotates every linear twist as BODY_FRD even when PX4
labels it LOCAL_NED. Recover its source frame and produce body FLU velocity.
Pose and angular velocity already use ENU/FLU correctly; no time prediction.
"""
import numpy as np

LOCAL_NED = 1
BODY_FRD = 12

def rotation_from_wxyz(quaternion):
    q = np.asarray(quaternion, dtype=float)
    if q.shape != (4,) or not np.isfinite(q).all() or np.linalg.norm(q) < 1e-9:
        raise ValueError('Invalid quaternion')
    w, x, y, z = q / np.linalg.norm(q)
    return np.array([[1-2*(y*y+z*z), 2*(x*y-w*z), 2*(x*z+w*y)],
                     [2*(x*y+w*z), 1-2*(x*x+z*z), 2*(y*z-w*x)],
                     [2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x*x+y*y)]])

def body_velocity_from_mavros(quaternion, mapped_velocity, source_frame):
    """Return (correct body velocity, linear twist frame rotation)."""
    rotation = rotation_from_wxyz(quaternion)
    velocity = np.asarray(mapped_velocity, dtype=float)
    if velocity.shape != (3,) or not np.isfinite(velocity).all():
        raise ValueError('Invalid velocity')
    if source_frame == BODY_FRD:
        transform = np.eye(3)
    elif source_frame == LOCAL_NED:
        # Undo MAVROS's FRD->FLU mapping, apply NED->ENU, then ENU->body.
        mapped_to_enu = np.array([[0,-1,0], [1,0,0], [0,0,1]])
        transform = rotation.T @ mapped_to_enu
    else:
        raise ValueError('Unsupported PX4 velocity frame: '+str(source_frame))
    return transform @ velocity, transform
