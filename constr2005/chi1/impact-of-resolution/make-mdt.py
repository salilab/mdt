from modeller import *
import mdt
import mdt.features

env = environ()
log.minimal()
env.io.atom_files_directory = ['/park2/database/pdb/divided/']

mlib = mdt.Library(env)
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 1.4, "under 1.4"),
                                               (1.4,  1.6, "1.4-1.6"),
                                               (1.6,  1.8, "1.6-1.8"),
                                               (1.8,  2.001, "1.8-2.0")])
restyp = mdt.features.ResidueType(mlib)
chi1 = mdt.features.Chi1Dihedral(mlib, bins=mdt.uniform_bins(144, -180, 2.5))

m = mdt.Table(mlib, features=(xray, restyp, chi1))

a = alignment(env)
f = modfile.File('../../cluster-PDB/pdb_60.pir', 'r')
while a.read_one(f):
    m.add_alignment(a)

m.write('mdt.mdt')
