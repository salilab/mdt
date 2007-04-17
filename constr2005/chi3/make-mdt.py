from modeller import *
import mdt

# See ../bonds/make_mdt.py for additional comments

env = environ()
log.minimal()
env.io.atom_files_directory = '/diva3/database/pdb/uncompressed_files'

mlib = mdt.mdt_library(env, '${LIB}/mdt.ini', '../lib/mdt2.bin')

m = mdt.mdt(mlib, features=(35,1,53))

a = alignment(env)
while (a.read_one(file='../cluster-PDB/pdb_60.pir')):
    m.add_alignment(a)

m.write('mdt.mdt')
