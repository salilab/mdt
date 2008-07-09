from modeller import *
import mdt

env = environ()
log.minimal()
env.io.atom_files_directory = ['/park2/database/pdb/divided/']

mlib = mdt.Library(env, '../../lib/mdt1.bin')

m = mdt.Table(mlib, features=(35,1,3))

a = alignment(env)
while (a.read_one(file='../../cluster-PDB/pdb_60.pir')):
    m.add_alignment(a)

m.write('mdt.mdt')
