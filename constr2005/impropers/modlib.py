from modeller import *
import mdt
import mdt.features

env = environ()
mlib = mdt.Library(env)
mlib.dihedral_classes.read('${LIB}/impgrp.lib')
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
impr_type = mdt.features.DihedralType(mlib)
improper = mdt.features.Dihedral(mlib, bins=mdt.uniform_bins(400, 1.0, 0.0025))

m = mdt.Table(mlib, file='mdt.mdt')
m = m.reshape(features=(xray, impr_type, improper),
              offset=(0,0,0), shape=(1,-1,-1))

m = m.integrate(features=(impr_type, improper))

mdt.write_improperlib(file('impropers.py', 'w'), m, density_cutoff=0.1)
