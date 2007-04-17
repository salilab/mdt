from modeller import *
import mdt

env = environ()
log.minimal()
env.io.atom_files_directory = '/diva3/database/pdb/uncompressed_files'

mlib = mdt.mdt_library(env, '${LIB}/mdt.ini', 'mdt_test_xray.bin')

m = mdt.mdt(mlib, features=35)

a = alignment(env)
while (a.read_one(file='../cluster-PDB/pdb_60.pir')):
    m.add_alignment(a)

m.write('mdt2.mdt')
