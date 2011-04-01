import modeller
import mdt
import mdt.features

# Setup of Modeller and MDT system
env = modeller.environ()
mlib = mdt.Library(env)

# Creation of feature types
restyp = mdt.features.ResidueType(mlib)

# Create a 1D table of residue type
table = mdt.Table(mlib, features=restyp)

# Read in a PDB file and make an alignment of just this one structure
mdl = modeller.model(env, file='5fd1')
aln = modeller.alignment(env)
aln.append_model(mdl, align_codes='5fd1', atom_files='5fd1')

# Collect MDT statistics for this alignment
table.add_alignment(aln)

# Print out the MDT by treating it as a Python list
print "Distribution of residue types:"
print [bin for bin in table]
