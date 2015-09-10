[![docs](https://readthedocs.org/projects/mdt/badge/)](http://mdt.readthedocs.org/) [![nightly build](https://salilab.org/nightly/stat/?s=mdt&t=build)](https://salilab.org/nightly/mdt/) [![python coverage](https://salilab.org/nightly/stat/?s=mdt&t=python)](https://salilab.org/nightly/mdt/logs/coverage/python/) [![C coverage](https://salilab.org/nightly/stat/?s=mdt&t=c)](https://salilab.org/nightly/mdt/logs/coverage/c/)

MDT prepares a raw frequency table, given information from MODELLER alignments
and/or PDB files. More precisely, MDT uses a sample of sequences, structures,
and/or alignments to construct a table *N(a,b,c,...,d)* for features
*a, b, c, ..., d*. The sample for generating the frequencies *N* is obtained
depending on the type of features *a, b, c, ..., d*. The sample can contain
individual proteins, pairs of proteins, pairs of residues in proteins,
pairs of aligned residues, pairs of aligned pairs of residues, chemical bonds,
angles, dihedral angles, and pairs of tuples of atoms. Some features work
with triple alignments, too. All the needed features *a, b, c, ..., d*
are calculated automatically from the sequences, alignments, and/or PDB files.
The feature bins are defined by the user when each feature is created.

To build MDT, simply type `scons` in this directory. To run it, use
`scons install` to install it (you may need to set command line options to
install files in the correct locations - see `scons -h`) and then run MDT
scripts just like any other Python script. Alternatively, you can use the
`bin/mdtpy.sh` script to run MDT directly from the build directory. Simply
run it before your regular Python invocation
(e.g. "`bin/mdtpy.sh python foo.py`")

MDT is Copyright 1989-2015 Andrej Sali, and available under the terms of
version 2 of the GNU General Public License as published by the
Free Software Foundation.
