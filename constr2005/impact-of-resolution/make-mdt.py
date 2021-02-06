from modeller import *
import mdt
import mdt.features

env = Environ()
log.minimal()
env.io.atom_files_directory = ['/salilab/park2/database/pdb/divided/']

mlib = mdt.Library(env)
xray = mdt.features.XRayResolution(mlib, bins=mdt.uniform_bins(20, 0, 0.2))
m = mdt.Table(mlib, features=xray)

a = Alignment(env)
f = modfile.File('../cluster-PDB/pdb_60.pir', 'r')
while a.read_one(f):
    m.add_alignment(a)

m.write('mdt2.mdt')
