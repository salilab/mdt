from modeller import *
import mdt

# See ../bonds/make_mdt.py for additional comments

env = environ()
log.minimal()
env.io.atom_files_directory = '/diva3/database/pdb/uncompressed_files'

mlib = mdt.mdt_library(env, '${LIB}/mdt.ini', '../lib/mdt2.bin', deltai=1)

# feature 66 is the subsequent residue type, relative to the omega (feature 28)
m = mdt.mdt(mlib, features=(35,66,28))

a = alignment(env)
while (a.read_one(file='../lib/pdb_60.pir')):
    m.add_alignment(a)

m.write('mdt.mdt')
