from modeller import *
from modeller import mdt

env = environ()
mlib = mdt.mdt_library(env, '${LIB}/mdt.ini', '../lib/mdt2.bin')
mlib.angle_classes.read('${LIB}/anggrp.lib')

m = mdt.mdt(mlib, file='mdt.mdt')
m = m.reshape(features=(35,111,112), offset=(0,0,0), shape=(0,-1,-1))

m = m.integrate(features=(111,112))

mdt.write_anglelib(file('angles.py', 'w'), m, density_cutoff=0.1)
