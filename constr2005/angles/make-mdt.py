from modeller import *
import mdt
import mdt.features

# See ../bonds/make_mdt.py for additional comments

env = environ()
log.minimal()
env.io.atom_files_directory = ['/park2/database/pdb/divided/']

mlib = mdt.Library(env)
mlib.angle_classes.read('${LIB}/anggrp.lib')
xray = mdt.features.XRayResolution(mlib, bins=[(0.51, 2.001, 'High res(2.0A)')])
angle_type = mdt.features.AngleType(mlib)
angle = mdt.features.Angle(mlib, bins=mdt.uniform_bins(720, 0, 0.25))

m = mdt.Table(mlib, features=(xray, angle_type, angle))

a = alignment(env)
f = modfile.File('../cluster-PDB/pdb_60.pir', 'r')
while a.read_one(f):
    m.add_alignment(a)

m.write('mdt.mdt')
