from modeller import *
import mdt
import mdt.features

env = environ()
mlib = mdt.Library(env)
mlib.angle_classes.read('${LIB}/anggrp.lib')
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
angle_type = mdt.features.AngleType(mlib)
angle = mdt.features.Angle(mlib, bins=mdt.uniform_bins(720, 0, 0.25))

m = mdt.Table(mlib, file='mdt.mdt')
m = m.reshape(features=(xray, angle_type, angle),
              offset=(0,0,0), shape=(0,-1,-1))

m = m.integrate(features=(angle_type, angle))

mdt.write_anglelib(file('angles.py', 'w'), m, density_cutoff=0.1)
