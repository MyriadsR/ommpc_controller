#!/usr/bin/env python3
"""Paper Fig. 7(a)-type five-loop roll helix; parameters are not digitized data."""
import argparse
import json
from pathlib import Path
import numpy as np
from generate_vertical_flip import hermite, evaluate


def generate(output, loops=5, radius=1.2, period=1.68, axial_travel=-2.4):
    if loops < 1 or radius <= 0 or period <= 0:
        raise ValueError('Positive loops, radius and period required')
    entry_time, exit_time, height, dt = 1.5, 2.5, 4.0, .01
    angular_speed = 2*np.pi/period
    speed = radius*angular_speed
    circle_duration = loops*period
    axial_speed = axial_travel/circle_duration
    side_height = height+speed*entry_time/2
    start = np.array([axial_speed*entry_time/2, 0, side_height])
    end = start+np.array([axial_travel, 0, 0])
    finish = end+np.array([axial_speed*exit_time/2, 0, speed*exit_time/2])
    velocity = np.array([axial_speed, 0, speed])
    acceleration = np.array([0, radius*angular_speed**2, 0])
    jerk = np.array([0, 0, -radius*angular_speed**3])
    zeros = np.zeros(3)
    entry = hermite([0,0,height], zeros, zeros, zeros, start, velocity,
                    acceleration, jerk, entry_time)
    exit = hermite(end, velocity, acceleration, jerk, finish, zeros, zeros,
                   zeros, exit_time)
    total = entry_time+circle_duration+exit_time

    def sample(time):
        if time < entry_time:
            return evaluate(entry, time, entry_time)
        if time < entry_time+circle_duration:
            tau = time-entry_time
            theta = tau*angular_speed
            return [start+np.array([axial_speed*tau, radius-radius*np.cos(theta), radius*np.sin(theta)]),
                    np.array([axial_speed, speed*np.sin(theta), speed*np.cos(theta)]),
                    np.array([0, radius*angular_speed**2*np.cos(theta), -radius*angular_speed**2*np.sin(theta)]),
                    np.array([0, -radius*angular_speed**3*np.sin(theta), -radius*angular_speed**3*np.cos(theta)])]
        return evaluate(exit, time-entry_time-circle_duration, exit_time)

    # Check the transitions against the analytic circular p/v/a/j, not only
    # continuity of the sampled table. Both joins must be C3.
    for polynomial, duration, time, expected in [
        (entry, entry_time, entry_time, [start,velocity,acceleration,jerk]),
        (exit, exit_time, 0, [end,velocity,acceleration,jerk])]:
        np.testing.assert_allclose(evaluate(polynomial,time,duration),expected,atol=1e-9)
    np.testing.assert_allclose(evaluate(exit,exit_time,exit_time),[finish,zeros,zeros,zeros],atol=1e-9)

    dense = np.asarray([np.concatenate(sample(t)) for t in np.linspace(0,total,int(total/.001)+1)])
    force = dense[:,6:9]+[0,0,9.81]
    magnitude = np.linalg.norm(force,axis=1)
    rate = np.linalg.norm(np.cross(force,dense[:,9:12]),axis=1)/magnitude**2
    requirements = dict(trajectory_type='vertical_helix_roll',paper_reference='Fig. 7(a), similar geometry; not original reference data',
                        radius_m=radius,loop_period_s=period,loop_speed_m_s=speed,loops=loops,
                        axial_travel_in_circles_m=axial_travel,axial_speed_m_s=axial_speed,
                        entry_s=entry_time,exit_s=exit_time,duration_s=total,time_step_s=dt,
                        minimum_thrust_acceleration_m_s2=float(magnitude.min()),
                        maximum_thrust_acceleration_m_s2=float(magnitude.max()),
                        maximum_body_rate_rad_s=float(rate.max()),
                        minimum_height_m=float(dense[:,2].min()),maximum_height_m=float(dense[:,2].max()),
                        circle_start_position_m=start.tolist(),end_position_m=finish.tolist(),c3_joins_verified=True)
    if magnitude.min()<5 or magnitude.max()>39 or rate.max()>14:
        raise ValueError('Trajectory exceeds aggressive-example limits: '+str(requirements))
    times = np.arange(int(round(total/dt))+1)*dt
    # Ref files have fixed timestep; do not silently truncate a custom duration.
    if not np.isclose(times[-1],total,atol=1e-9):
        raise ValueError('Total duration must be a multiple of 0.01 s')
    rows = [np.r_[p,v,0.,0.,a,j] for p,v,a,j in map(sample,times)]
    output = Path(output)
    np.savetxt(output,rows,fmt='%.12f')
    output.with_suffix('.json').write_text(json.dumps(requirements,indent=2))
    print(json.dumps(requirements,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('vertical_helix.traj'))
    parser.add_argument('--loops',type=int,default=5)
    parser.add_argument('--radius',type=float,default=1.2)
    parser.add_argument('--period',type=float,default=1.68)
    parser.add_argument('--axial-travel',type=float,default=-2.4)
    args = parser.parse_args()
    generate(args.output,args.loops,args.radius,args.period,args.axial_travel)
