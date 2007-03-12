from modeller import *
import mdt

env = environ()
mlib = mdt.mdt_library(env, '${LIB}/mdt.ini', '../lib/mdt2.bin')
mlib.bond_classes.read('${LIB}/bndgrp.lib')

m = mdt.mdt(mlib, file='mdt.mdt')
m = m.reshape(features=(35,109,110), offset=(0,0,0), shape=(0,-1,-1))

m = m.integrate(features=(109,110))

mdt.write_bondlib(file('bonds.py', 'w'), m, density_cutoff=0.1)
