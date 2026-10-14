#!/usr/bin/env python3
"""C3 entry, constant-speed inverted circles, C3 exit; no near-zero-force crossing."""
from pathlib import Path
import argparse
import json
import numpy as np


def hermite(p0, v0, a0, j0, p1, v1, a1, j1, duration):
    c = np.zeros((8, 3))
    c[:4] = [p0, np.asarray(v0)*duration, np.asarray(a0)*duration**2/2,
             np.asarray(j0)*duration**3/6]
    target = np.array([p1, np.asarray(v1)*duration, np.asarray(a1)*duration**2,
                       np.asarray(j1)*duration**3])
    initial = np.array([c[:4].sum(axis=0), c[1]+2*c[2]+3*c[3],
                        2*c[2]+6*c[3], 6*c[3]])
    matrix = np.array([[1,1,1,1], [4,5,6,7], [12,20,30,42], [24,60,120,210]])
    c[4:] = np.linalg.solve(matrix, target-initial)
    return c


def evaluate(c, time, duration):
    x = time/duration
    values = []
    for derivative in range(4):
        value = np.zeros(3)
        for k in range(derivative, 8):
            factor = np.prod(np.arange(k-derivative+1, k+1)) if derivative else 1
            value += factor*c[k]*x**(k-derivative)/duration**derivative
        values.append(value)
    return values


def generate(output, loops=3):
    radius, period, entry_time, exit_time, height, dt = 1.2, 1.68, 1.5, 2.5, 4.0, .01
    speed = 2*np.pi*radius/period
    side_height = height+speed*entry_time/2
    finish_height = side_height+speed*exit_time/2
    circle_duration = loops*period
    total = entry_time+circle_duration+exit_time
    pv = [0,0,speed]; pa = [speed**2/radius,0,0]; pj = [0,0,-speed**3/radius**2]
    zeros = [0,0,0]
    entry = hermite([0,0,height],zeros,zeros,zeros,[0,0,side_height],pv,pa,pj,entry_time)
    exit = hermite([0,0,side_height],pv,pa,pj,[0,0,finish_height],zeros,zeros,zeros,exit_time)
    rows=[]
    for time in np.arange(round(total/dt)+1)*dt:
        if time<entry_time:
            p,v,a,j=evaluate(entry,time,entry_time)
        elif time<entry_time+circle_duration:
            theta=(time-entry_time)*speed/radius
            p=np.array([radius-radius*np.cos(theta),0,side_height+radius*np.sin(theta)])
            v=np.array([speed*np.sin(theta),0,speed*np.cos(theta)])
            a=np.array([speed**2/radius*np.cos(theta),0,-speed**2/radius*np.sin(theta)])
            j=np.array([-speed**3/radius**2*np.sin(theta),0,-speed**3/radius**2*np.cos(theta)])
        else:
            p,v,a,j=evaluate(exit,time-entry_time-circle_duration,exit_time)
        rows.append(np.r_[p,v,0.,0.,a,j])
    rows=np.array(rows)
    force=rows[:,8:11]+[0,0,9.81]
    magnitude=np.linalg.norm(force,axis=1)
    rate=np.linalg.norm(np.cross(force,rows[:,11:14]),axis=1)/magnitude**2
    requirements={'radius_m':radius,'loop_speed_m_s':speed,'loop_period_s':period,
                  'entry_s':entry_time,'exit_s':exit_time,'loops':loops,'duration_s':total,
                  'minimum_thrust_acceleration_m_s2':float(magnitude.min()),
                  'maximum_thrust_acceleration_m_s2':float(magnitude.max()),
                  'maximum_body_rate_rad_s':float(rate.max()),
                  'minimum_height_m':float(rows[:,2].min()),'maximum_height_m':float(rows[:,2].max())}
    if magnitude.min()<1 or magnitude.max()>39 or rate.max()>14:
        raise ValueError('Trajectory exceeds aggressive-example input limits: '+str(requirements))
    np.savetxt(output,rows,fmt='%.12f')
    Path(output).with_suffix('.json').write_text(json.dumps(requirements,indent=2))
    print(json.dumps(requirements,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=Path(__file__).with_name('vertical_flip.traj'))
    parser.add_argument('--loops',type=int,default=3)
    args=parser.parse_args()
    if args.loops<1:parser.error('--loops must be positive')
    generate(args.output,args.loops)
