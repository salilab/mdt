from modeller import *
import mdt

env = environ()
mlib = mdt.Library(env, '../lib/mdt2.bin')
mlib.dihedral_classes.read('${LIB}/impgrp.lib')

m = mdt.Table(mlib, file='mdt.mdt')
m = m.reshape(features=(35,113,114), offset=(0,0,0), shape=(1,-1,-1))

m = m.integrate(features=(113,114))

mdt.write_improperlib(file('impropers.py', 'w'), m, density_cutoff=0.1)
