from modeller import *
import mdt
import mdt.features

# See ../bonds/make_mdt.py for additional comments

env = environ()
log.minimal()
env.io.atom_files_directory = ['/park2/database/pdb/divided/']

mlib = mdt.Library(env)
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
restyp_1 = mdt.features.ResidueType(mlib, delta=1)
omega = mdt.features.OmegaDihedral(mlib, bins=mdt.uniform_bins(720, -180, 0.5))

# This table uses the subsequent residue type, relative to the omega
m = mdt.Table(mlib, features=(xray, restyp_1, omega))

a = alignment(env)
f = modfile.File('../cluster-PDB/pdb_60.pir', 'r')
while a.read_one(f):
    m.add_alignment(a)

m.write('mdt.mdt')
