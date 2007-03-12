from modeller import *
import mdt

env = environ()
log.minimal()
env.io.atom_files_directory = '/diva3/database/pdb/uncompressed_files'

mlib = mdt.mdt_library(env, '${LIB}/mdt.ini', '../../lib/mdt1.bin')

m = mdt.mdt(mlib, features=(35,1,3))

a = alignment(env)
while (a.read_one(file='../../lib/pdb_60.pir')):
    m.add_alignment(a)

m.write('mdt.mdt')
