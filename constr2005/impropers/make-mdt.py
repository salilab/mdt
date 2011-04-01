from modeller import *
import mdt
import mdt.features

# See ../bonds/make_mdt.py for additional comments

env = environ()
log.minimal()
env.io.atom_files_directory = ['/salilab/park2/database/pdb/divided/']

mlib = mdt.Library(env)
mlib.dihedral_classes.read('${LIB}/impgrp.lib')
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
impr_type = mdt.features.DihedralType(mlib)
improper = mdt.features.Dihedral(mlib, bins=mdt.uniform_bins(400, 1.0, 0.0025))

m = mdt.Table(mlib, features=(xray, impr_type, improper))

a = alignment(env)
f = modfile.File('../cluster-PDB/pdb_60.pir', 'r')
while a.read_one(f):
    m.add_alignment(a)

m.write('mdt.mdt')
