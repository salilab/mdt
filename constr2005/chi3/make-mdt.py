from modeller import *
import mdt
import mdt.features

# See ../bonds/make_mdt.py for additional comments

env = environ()
log.minimal()
env.io.atom_files_directory = ['/park2/database/pdb/divided/']

mlib = mdt.Library(env)
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
restyp = mdt.features.ResidueType(mlib)
chi3 = mdt.features.Chi3Dihedral(mlib, bins=mdt.uniform_bins(144, -180, 2.5))

m = mdt.Table(mlib, features=(xray, restyp, chi3))

a = alignment(env)
f = modfile.File('../cluster-PDB/pdb_60.pir', 'r')
while a.read_one(f):
    m.add_alignment(a)

m.write('mdt.mdt')
