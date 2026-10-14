#!/usr/bin/env python3
"""Generate an explicit test vehicle from Iris; preserve stock PX4 files."""
from pathlib import Path
import argparse
import xml.etree.ElementTree as ET
parser=argparse.ArgumentParser()
parser.add_argument('--px4',type=Path,default=Path.home()/'PX4-Autopilot')
args=parser.parse_args()
source=args.px4/'Tools/simulation/gazebo-classic/sitl_gazebo-classic/models/iris/iris.sdf'
root=ET.parse(source).getroot()
root.find('model').set('name','iris_agile')
for plugin in root.findall('.//plugin'):
    constant=plugin.find('motorConstant')
    if constant is not None:constant.text=str(float(constant.text)*2.4)
output=Path(__file__).with_name('iris_agile.sdf')
output.write_text('<!-- Experimental vehicle: stock Iris inertias/mass, motorConstant x2.4; not a hardware replica. -->\n'+ET.tostring(root,encoding='unicode'))
print(output)
