from modeller import *
import mdt
import mdt.features

env = Environ()
mlib = mdt.Library(env)
mlib.bond_classes.read('${LIB}/bndgrp.lib')
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
bond_type = mdt.features.BondType(mlib)
bond_length = mdt.features.BondLength(mlib,
                                      bins=mdt.uniform_bins(400, 1.0, 0.0025))

m = mdt.Table(mlib, file='mdt.mdt')
m = m.reshape(features=(xray, bond_type, bond_length),
              offset=(0,0,0), shape=(0,-1,-1))

m = m.integrate(features=(bond_type, bond_length))

mdt.write_bondlib(file('bonds.py', 'w'), m, density_cutoff=0.1)
