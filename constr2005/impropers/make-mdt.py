from modeller import *
import mdt

# See ../bonds/make_mdt.py for additional comments

env = environ()
log.minimal()
env.io.atom_files_directory = ['/park2/database/pdb/divided/']

mlib = mdt.Library(env, '../lib/mdt2.bin')

mlib.dihedral_classes.read('${LIB}/impgrp.lib')

m = mdt.Table(mlib, features=(35,113,114))

a = alignment(env)
while (a.read_one(file='../cluster-PDB/pdb_60.pir')):
    m.add_alignment(a)

m.write('mdt.mdt')
