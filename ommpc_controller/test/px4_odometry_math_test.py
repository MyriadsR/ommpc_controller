#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'sim'))
from px4_odometry_math import body_velocity_from_mavros, LOCAL_NED, BODY_FRD

class OdometryFrameTest(unittest.TestCase):
    def test_ned_world_velocity_is_correct_at_yawed_attitude(self):
        q=[2**-.5,0,0,2**-.5]
        v,_=body_velocity_from_mavros(q,[1,-2,-3],LOCAL_NED)
        np.testing.assert_allclose(v,[1,-2,-3],atol=1e-12)

    def test_large_pitch_transforms_world_velocity_to_body(self):
        v,_=body_velocity_from_mavros([2**-.5,0,2**-.5,0],[0,-1,0],LOCAL_NED)
        np.testing.assert_allclose(v,[0,0,1],atol=1e-12)

    def test_inverted_aircraft_reverses_horizontal_body_velocity(self):
        v,_=body_velocity_from_mavros([0,0,1,0],[0,-1,0],LOCAL_NED)
        np.testing.assert_allclose(v,[-1,0,0],atol=1e-12)

    def test_actual_body_frame_requires_no_extra_rotation(self):
        v,_=body_velocity_from_mavros([0,0,1,0],[1,2,3],BODY_FRD)
        np.testing.assert_array_equal(v,[1,2,3])

    def test_unknown_frames_and_invalid_states_are_rejected(self):
        for q,v,frame in [([0,0,0,0],[1,2,3],LOCAL_NED),([1,0,0,0],[np.nan,0,0],LOCAL_NED),([1,0,0,0],[1,2,3],20)]:
            with self.assertRaises(ValueError):body_velocity_from_mavros(q,v,frame)

if __name__=='__main__':unittest.main()
