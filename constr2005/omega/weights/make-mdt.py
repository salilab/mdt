from modeller import *
import mdt
import mdt.features

# See ../bonds/make_mdt.py for additional comments

env = environ()
log.minimal()
env.io.atom_files_directory = ['/salilab/park2/database/pdb/divided/']

mlib = mdt.Library(env)
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
restyp_1 = mdt.features.ResidueType(mlib, delta=1)
omega_class = mdt.features.OmegaClass(mlib)

# Table of the subsequent residue type relative to the omega class
m = mdt.Table(mlib, features=(xray, restyp_1, omega_class))

a = alignment(env)
f = modfile.File('../../cluster-PDB/pdb_60.pir', 'r')
while a.read_one(f):
    m.add_alignment(a)

m.write('mdt.mdt')
