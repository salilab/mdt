.. highlight:: rest

.. currentmodule:: mdt

Introduction
============

MDT prepares a raw frequency table, given information from MODELLER
alignments and/or PDB files. It can also process the raw frequency table in
several ways (e.g., normalization with :meth:`Table.normalize`, smoothing with
:meth:`Table.smooth`, perform entropy calculations with
:meth:`Table.entropy_full`, and write out the data in various formats,
including for plotting by ASGL
(:meth:`Table.write_asgl`) and use as restraints by MODELLER.

More precisely, MDT uses a sample of sequences, structures, and/or
alignments to construct a table N(a,b,c,...,d) for features a, b, c, ..., d.
The sample for generating the frequencies N is obtained depending on the type
of features a, b, c, ..., d. The sample can contain individual proteins, pairs
of proteins, pairs of residues in proteins, pairs of aligned residues, pairs
of aligned pairs of residues, chemical bonds, angles, dihedral angles, and
pairs of tuples of atoms. Some features work with triple alignments, too.
All the needed features a, b, c, ..., d are calculated automatically from the
sequences, alignments, and/or PDB files. The feature bins are defined by the
user when each feature is created.
