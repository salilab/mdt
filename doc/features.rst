.. highlight:: rest

The :mod:`mdt.features` Python module
=====================================

.. automodule:: mdt.features

.. _protein_features:

Protein features
----------------

These features yield a single value for each protein in the alignment.
Each feature takes some common arguments:

 * *mlib*: the :class:`mdt.Library` to create the feature in.
 * *bins*: list of bins (see :ref:`binspec`).
 * *protein*: the protein index on which to evaluate the feature
   from each group of proteins (individual protein, pairs, triples)
   selected from the alignment (0 for the first, 1 for the second, 2 for
   the third). See :meth:`Table.add_alignment` for more details.

.. autoclass:: XRayResolution

.. autoclass:: RadiusOfGyration

.. autoclass:: SequenceLength

.. autoclass:: HydrogenBondSatisfaction

.. autoclass:: AlphaContent

.. _protein_pair_features:

Protein pair features
---------------------

These features yield a single value for each pair of proteins in the alignment.
Each feature takes some common arguments:

 * *mlib*: the :class:`mdt.Library` to create the feature in.
 * *bins*: list of bins (see :ref:`binspec`).
 * *protein1* and *protein2*: the indexes of proteins in each group of
   proteins selected from the alignment to evaluate the feature on; each can
   range from 0 to 2 inclusive. See :meth:`Table.add_alignment` for more
   details.

.. autoclass:: SequenceIdentity

.. _residue_features:

Residue features
----------------

These features yield a single value for each residue in each sequence
in the alignment. Each feature takes some common arguments:

 * *delta*: if non-zero, don't calculate the feature for the residue
   position returned by the residue scan - instead, offset it by
   *delta* residues in the sequence. Applied before *align_delta*.
 * *align_delta*: if non-zero, don't calculate the feature for the
   alignment position returned by the residue scan - instead, offset
   it by *align_delta* alignment positions. Applied after *delta*.
 * *pos2*: if True, force a residue pair scan, and evaluate the
   feature on the second residue in each pair.
 * *mlib*, *bins*, *protein*: see :ref:`protein_features`. Note that some
   residue features do not use the *bins* argument, because they have a
   fixed number of bins.

.. autoclass:: ResidueType

.. autoclass:: ResidueAccessibility

.. class:: Chi1Dihedral(self, mlib, protein=0, delta=0, align_delta=0, pos2=False)
           Chi2Dihedral
           Chi3Dihedral
           Chi4Dihedral
           PhiDihedral
           PsiDihedral
           OmegaDihedral
           AlphaDihedral

           Residue dihedral angle, from -180 to 180 degrees.

.. class:: Chi1Class(self, mlib, protein=0, delta=0, align_delta=0, pos2=False)
           Chi2Class
           Chi3Class
           Chi4Class
           Chi5Class
           PhiClass
           PsiClass
           OmegaClass

           Residue dihedral class. These classes are defined by MODELLER to
           group common regions of dihedral space for each residue type.

.. autoclass:: MainchainConformation

.. autoclass:: ResidueGroup

.. autoclass:: SidechainBiso

.. _residue_pair_features:

Residue pair features
---------------------

These features yield a single value for each pair of residues in each sequence
in the alignment. See :ref:`protein_features` for a description of the common
arguments.

.. autoclass:: ResidueDistance

.. autoclass:: AverageResidueAccessibility

.. autoclass:: ResidueIndexDifference

.. _aligned_residue_features:

Aligned residue features
------------------------

These features yield a single value for residues aligned between two proteins.
For each pair of proteins, every alignment position is scanned, and the
feature is evaluated for each pair of aligned residues.
See :ref:`protein_pair_features` for a description of the common arguments.

.. class:: PhiDihedralDifference(self, mlib, bins, protein1=0, protein2=1)
           PsiDihedralDifference
           OmegaDihedralDifference

           Shortest difference in dihedral angle (in degrees) between a pair of
           aligned residues.

.. autoclass:: NeighborhoodDifference

.. autoclass:: GapDistance


Aligned residue pair features
-----------------------------

These features yield a single value for each pair of residues aligned
between two proteins. For each pair of proteins, each pair of alignment
positions is scanned, and the feature is evaluated for each pair of pairs
of aligned residues.
See :ref:`protein_pair_features` for a description of the common arguments.

.. autoclass:: ResidueDistanceDifference

.. autoclass:: AverageNeighborhoodDifference

.. autoclass:: AverageGapDistance

.. _atom_features:

Atom features
-------------

These features yield a single value for each atom in the first protein
in each group of proteins selected from the alignment.
Each feature takes some common arguments:

 * *pos2*: if True, force an atom pair scan, and evaluate the
   feature on the second atom in each pair.
 * *mlib*, *bins*: see :ref:`protein_features`. Note that some
   atom features do not use the *bins* argument, because they have a
   fixed number of bins.

.. autoclass:: AtomAccessibility

.. autoclass:: FractionalAtomAccessibility

.. autoclass:: AtomType

.. autoclass:: HydrogenBondDonor

.. autoclass:: HydrogenBondAcceptor

.. autoclass:: HydrogenBondCharge

.. autoclass:: AtomTable

.. _atom_pair_features:

Atom pair features
---------------------

These features yield a single value for each pair of atoms in the first
protein in each group of proteins selected from the alignment.
See :ref:`protein_features` for a description of the common arguments.

.. autoclass:: AtomDistance
.. autoclass:: AtomBondSeparation

.. _tuple_features:

Tuple features
--------------

These features yield a single value for each tuple of atoms in the first
protein in each group of proteins selected from the alignment. (The set of
tuples must first be read into the :class:`mdt.Library`.)
Each feature takes some common arguments:

 * *mlib*: the :class:`mdt.Library` to create the feature in.
 * *pos2*: if True, force a tuple pair scan, and evaluate the
   feature on the second tuple in each pair.

.. autoclass:: TupleType

.. _tuple_pair_features:

Tuple pair features
-------------------

These features yield a single value for each pair of tuples of atoms in the
first protein in each group of proteins selected from the alignment. (The set of
tuples must first be read into the :class:`mdt.Library`.)
See :ref:`protein_features` for a description of the common arguments.

.. autoclass:: TupleType

.. autoclass:: TupleDistance

.. autoclass:: TupleAngle1

.. autoclass:: TupleAngle2

.. autoclass:: TupleDihedral1

.. autoclass:: TupleDihedral2

.. autoclass:: TupleDihedral3

.. _chemical_bond_features:

Chemical bond features
----------------------

These features yield a single value for each defined chemical bond, angle
or dihedral in the first protein in each group of proteins selected from the
alignment. (The definitions of the chemical connectivity must first be read
from a bond class file; see the *bond_clases*, *angle_classes* and
*dihedral_classes* attributes in :class:`mdt.Library`.)
See :ref:`protein_features` for a description of the common arguments.

.. autoclass:: BondType

.. autoclass:: AngleType

.. autoclass:: DihedralType

.. autoclass:: BondLength

.. autoclass:: Angle

.. autoclass:: Dihedral

Group features
--------------

These features are used to make combinations of other features.
Each feature takes some common arguments:

 * *mlib*: the :class:`mdt.Library` to create the feature in.
 * *feat1*: an existing feature object that will be included in this group.
 * *feat2*: another existing feature object to include.
 * *nbins*: the number of bins in this feature.

.. autoclass:: Cluster
