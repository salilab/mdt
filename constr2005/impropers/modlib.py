from modeller import *
from modeller import mdt

env = environ()
mlib = mdt.mdt_library(env, '${LIB}/mdt.ini', '../lib/mdt2.bin')
mlib.dihedral_classes.read('${LIB}/impgrp.lib')

m = mdt.mdt(mlib, file='mdt.mdt')
m = m.reshape(features=(35,113,114), offset=(0,0,0), shape=(0,-1,-1))

m = m.integrate(features=(113,114))

mdt.write_improperlib(file('impropers.py', 'w'), m, density_cutoff=0.1)
