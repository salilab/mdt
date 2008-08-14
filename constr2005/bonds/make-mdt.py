from modeller import *
import mdt
import mdt.features

env = environ()
log.minimal()
env.io.atom_files_directory = ['/park2/database/pdb/divided/']

mlib = mdt.Library(env)

# read the bond definitions in terms of the constituting atom type pairs:
mlib.bond_classes.read('${LIB}/bndgrp.lib')

# define the features:
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, '<=2.0')])
bond_type = mdt.features.BondType(mlib)
bond_length = mdt.features.BondLength(mlib,
                                      bins=mdt.uniform_bins(400, 1.0, 0.0025))

m = mdt.Table(mlib, features=(xray, bond_type, bond_length))

# make the MDT table using the pdb_60 sample chains:
a = alignment(env)
f = modfile.File('../cluster-PDB/pdb_60.pir', 'r')
while a.read_one(f):
    m.add_alignment(a)

# write out the MDT table:
m.write('mdt.mdt')
