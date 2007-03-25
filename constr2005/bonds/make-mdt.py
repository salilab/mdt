from modeller import *
import mdt

env = environ()
log.minimal()
env.io.atom_files_directory = '/diva3/database/pdb/uncompressed_files'

mlib = mdt.mdt_library(env, '${LIB}/mdt.ini', '../lib/mdt2.bin')

# read the bond definitions in terms of the constituting atom type pairs:
mlib.bond_classes.read('${LIB}/bndgrp.lib')

# define the features: X-ray resolution, bond type, and bond length:
m = mdt.mdt(mlib, features=(35,109,110))

# make the MDT table using the pdb_60 sample chains:
a = alignment(env)
while (a.read_one(file='../lib/pdb_60.pir')):
    m.add_alignment(a)

# write out the MDT table:
m.write('mdt.mdt')