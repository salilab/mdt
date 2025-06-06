"""
   MDT features.

   Copyright 1989-2025 Andrej Sali.

   MDT is free software: you can redistribute it and/or modify
   it under the terms of version 2 of the GNU General Public License
   as published by the Free Software Foundation.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with MDT.  If not, see <http://www.gnu.org/licenses/>.
"""

__docformat__ = "restructuredtext"

import _mdt


class _Base(object):
    """Base class for all features."""

    def _setup(*args, **keys):
        raise TypeError("This feature class is abstract, and cannot "
                        "be instantiated directly.")

    def __init__(self, mlib):
        self._mlib = mlib

    def _get_ifeat(self, mlib):
        if mlib == self._mlib:
            return self._ifeat
        else:
            raise ValueError("Cannot use features from different libraries")

    def _create_bins(self, mlib, bins):
        """Set bins in the MDT library for this feature."""
        _mdt.mdt_feature_nbins_set(mlib._modpt, self._ifeat, len(bins))
        for i, (start, end, symbol) in enumerate(bins):
            _mdt.mdt_feature_bin_set(mlib._modpt, self._ifeat, i, start,
                                     end, symbol)


class Protein(_Base):
    """A feature defined on entire proteins."""
    def __init__(self, mlib, bins, protein=0):
        """
        Create the protein feature.

        :Parameters:
          - `mlib`: the `Library` to create the feature in.
          - `bins`: list of bins (list of (start, end, symbol) triples, or
            the result of `mdt.uniform_bins`).
          - `protein`: the protein index on which to evaluate the feature
            from each group of proteins selected from the alignment (0 for
            the first, 1 for the second, etc.). If 0, a scan over all
            individual proteins is requested, while 1 requests a scan over
            all pairs, 2 over all triples, etc. (However, another feature may
            request a scan over a larger group - for example, if some other
            feature requests a triples scan while protein=1, this feature will
            be evaluated on the second protein in each triple.)
        """
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, protein)
        self._create_bins(mlib, bins)


class ProteinPair(_Base):
    """A feature defined on pairs of proteins."""
    def __init__(self, mlib, bins, protein1=0, protein2=1):
        """
        Create a new protein pair feature.
        See `Protein` for an explanation of the `mlib` and `bins`
        parameters. The `protein1` and `protein2` arguments determine the
        indexes of proteins in each group of proteins selected from the
        alignment to evaluate the feature on. If set to 0 and 1 respectively,
        a scan of all protein pairs in the alignment is requested and the
        feature is evaluated on each pair in turn. If set to 0 and 2, a scan
        of all protein triples is asked for, and the feature is evaluated on
        the first and third proteins in each triple.
        """
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, protein1, protein2)
        self._create_bins(mlib, bins)


class Residue(_Base):
    """A feature defined on single residues in a protein."""
    def __init__(self, mlib, bins, protein=0, delta=0, align_delta=0,
                 pos2=False):
        """
        Create a new residue feature.
        See `Protein` for an explanation of the `mlib`, `bins` and `protein`
        parameters.

        :Parameters:
          - `delta`: if non-zero, don't calculate the feature for the residue
            position returned by the residue scan - instead, offset it by
            `delta` residues in the sequence. Applied before `align_delta`.
          - `align_delta`: if non-zero, don't calculate the feature for the
            alignment position returned by the residue scan - instead, offset
            it by `align_delta` alignment positions. Applied after `delta`.
          - `pos2`: if True, force a residue pair scan, and evaluate the
            feature on the second residue in each pair.
        """
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, protein, delta, align_delta,
                                  pos2)
        self._create_bins(mlib, bins)


class ResidueFixedBins(Residue):
    """A residue feature for which bins need not be specified."""
    def __init__(self, mlib, protein=0, delta=0, align_delta=0, pos2=False):
        """See `Residue` for a description of the arguments."""
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, protein, delta, align_delta,
                                  pos2, mlib._env.libs.modpt)


class ResiduePair(_Base):
    """A feature defined on pairs of residues in a protein."""
    def __init__(self, mlib, bins, protein=0):
        """See `Protein` for a description of the arguments."""
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, protein)
        self._create_bins(mlib, bins)


class AlignedResidue(_Base):
    """A feature defined on residues aligned between two proteins.
       For each pair of proteins, every alignment position is scanned, and the
       feature is evaluated for each pair of aligned residues."""
    def __init__(self, mlib, bins, protein1=0, protein2=1):
        """See `ResiduePair` for a description of the arguments."""
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, protein1, protein2)
        self._create_bins(mlib, bins)


class AlignedResiduePair(_Base):
    """A feature defined on pairs of residues aligned between two proteins.
       For each pair of proteins, each pair of alignment positions is scanned,
       and the feature is evaluated for each pair of pairs of aligned
       residues."""
    def __init__(self, mlib, bins, protein1=0, protein2=1):
        """See `ResiduePair` for a description of the arguments."""
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, protein1, protein2)
        self._create_bins(mlib, bins)


class Atom(_Base):
    """A feature defined on single atoms.
       This works only on the first protein in each group of proteins
       selected from the alignment."""
    def __init__(self, mlib, bins, pos2=False):
        """See `Residue` for a description of the arguments."""
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, pos2)
        self._create_bins(mlib, bins)


class AtomFixedBins(Atom):
    """An atom feature for which bins need not be specified."""
    def __init__(self, mlib, pos2=False):
        """See `Residue` for a description of the arguments."""
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, pos2)


class AtomPair(_Base):
    """A feature defined on pairs of atoms.
       This works only on the first protein in each group of proteins
       selected from the alignment."""
    def __init__(self, mlib, bins):
        """See `Protein` for a description of the arguments."""
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt)
        self._create_bins(mlib, bins)


class Tuple(_Base):
    """A feature defined on tuples of atoms.
       This works only on the first protein in each group of proteins
       selected from the alignment."""
    def __init__(self, mlib, bins, pos2=False):
        """See `Residue` for a description of the arguments."""
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, pos2)
        self._create_bins(mlib, bins)


class TupleFixedBins(Tuple):
    """A tuple feature for which bins need not be specified."""
    def __init__(self, mlib, pos2=False):
        """See `Residue` for a description of the arguments."""
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, pos2)


class TuplePair(_Base):
    """A feature defined on pairs of tuples of atoms.
       This works only on the first protein in each group of proteins
       selected from the alignment."""
    def __init__(self, mlib, bins):
        """See `Protein` for a description of the arguments."""
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt)
        self._create_bins(mlib, bins)


class ChemicalBond(_Base):
    """A feature defined on defined chemical bonds, angles, or dihedrals.
       The definitions of the chemical connectivity must first be read from
       a bond class file; see `BondClasses`. This feature works only on the
       first protein in each group of proteins selected from the alignment."""
    def __init__(self, mlib, bins):
        """See `Protein` for a description of the arguments."""
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt)
        self._create_bins(mlib, bins)


class ChemicalBondFixedBins(ChemicalBond):
    """A chemical bond feature for which bins need not be specified."""
    def __init__(self, mlib):
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt)


class XRayResolution(Protein):
    """
    Protein X-ray resolution in angstroms.
    Proteins with a resolution of -1.00 (generally NMR structures) are
    actually reported as having a resolution of *nmr*.
    This decreases the number of bins required to hold all defined
    resolutions while still separating NMR from X-ray structures.
    """
    _setup = _mdt.mdt_feature_xray_resolution

    def __init__(self, mlib, bins, protein=0, nmr=0.45):
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, protein, nmr)
        self._create_bins(mlib, bins)


class RadiusOfGyration(Protein):
    """Protein radius of gyration in angstroms. The calculation of the center
       of mass used for this feature is not mass weighted."""
    _setup = _mdt.mdt_feature_radius_of_gyration


class SequenceLength(Protein):
    """Protein sequence length (number of residues)."""
    _setup = _mdt.mdt_feature_sequence_length


class ResidueType(ResidueFixedBins):
    """Residue type (20 standard amino acids, gap, undefined)."""
    _setup = _mdt.mdt_feature_residue_type


class ResidueAccessibility(Residue):
    """Residue solvent accessibility. This is derived from the atomic solvent
       accessibility; see :class:`AtomAccessibility`."""
    _setup = _mdt.mdt_feature_residue_accessibility


class Chi1Dihedral(Residue):
    """Residue chi1 dihedral angle, from -180 to 180 degrees.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_chi1_dihedral


class Chi2Dihedral(Residue):
    """Residue chi2 dihedral angle, from -180 to 180 degrees.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_chi2_dihedral


class Chi3Dihedral(Residue):
    """Residue chi3 dihedral angle, from -180 to 180 degrees.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_chi3_dihedral


class Chi4Dihedral(Residue):
    """Residue chi4 dihedral angle, from -180 to 180 degrees.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_chi4_dihedral


class PhiDihedral(Residue):
    """Residue phi dihedral angle, from -180 to 180 degrees.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_phi_dihedral


class PsiDihedral(Residue):
    """Residue psi dihedral angle, from -180 to 180 degrees.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_psi_dihedral


class OmegaDihedral(Residue):
    """Residue omega dihedral angle, from -180 to 180 degrees.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_omega_dihedral


class AlphaDihedral(Residue):
    """Residue alpha dihedral angle, from -180 to 180 degrees.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_alpha_dihedral


class Chi1Class(ResidueFixedBins):
    """Residue chi1 dihedral class."""
    _setup = _mdt.mdt_feature_chi1_class


class Chi2Class(ResidueFixedBins):
    """Residue chi2 dihedral class."""
    _setup = _mdt.mdt_feature_chi2_class


class Chi3Class(ResidueFixedBins):
    """Residue chi3 dihedral class."""
    _setup = _mdt.mdt_feature_chi3_class


class Chi4Class(ResidueFixedBins):
    """Residue chi4 dihedral class."""
    _setup = _mdt.mdt_feature_chi4_class


class Chi5Class(ResidueFixedBins):
    """Residue chi5 dihedral class."""
    _setup = _mdt.mdt_feature_chi5_class


class PhiClass(ResidueFixedBins):
    """Residue phi dihedral class."""
    _setup = _mdt.mdt_feature_phi_class


class PsiClass(ResidueFixedBins):
    """Residue psi dihedral class."""
    _setup = _mdt.mdt_feature_psi_class


class OmegaClass(ResidueFixedBins):
    """Residue omega dihedral class."""
    _setup = _mdt.mdt_feature_omega_class


class MainchainConformation(ResidueFixedBins):
    """Residue mainchain conformation (Ramachandran) class.
       This is a classification of the residue's phi/psi angles into classes
       as defined in Modeller's modlib/af_mnchdef.lib file and described in
       Sali and Blundell, JMB (1993) 234, p785. The default classes are
       A (right-handed alpha-helix), P (poly-proline conformation),
       B (idealized beta-strand), L (left-handed alpha-helix), and
       E (extended conformation)."""
    _setup = _mdt.mdt_feature_mainchain_conformation


class ResidueGroup(ResidueFixedBins):
    """Residue group."""
    _setup = _mdt.mdt_feature_residue_group
    def __init__(self, mlib, protein=0, delta=0, align_delta=0, pos2=False,
                 residue_grouping=0):
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, protein, delta, align_delta,
                                  pos2, residue_grouping, mlib._env.libs.modpt)


class SidechainBiso(Residue):
    r"""
    Residue average sidechain |Biso|. A zero average |Biso| is treated as
    undefined. If the average of these values over the whole protein is less
    than 2, each residue's value is multiplied by 4 |pi| :sup:`2`.

    .. |pi| unicode:: U+03C0
    .. |Biso| replace:: B\ :sub:`iso`
    """
    _setup = _mdt.mdt_feature_sidechain_biso


class ResidueDistance(ResiduePair):
    """Distance between a pair of residues. This is defined as the distance
       between the 'special' atoms in each residue. The type of this special
       atom can be specified by the distance_atoms argument when creating a
       :class:`mdt.Library` object.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_residue_distance


class AverageResidueAccessibility(ResiduePair):
    """Average solvent accessibility of a pair of residues.
       See :class:`ResidueAccessibility`."""
    _setup = _mdt.mdt_feature_average_residue_accessibility


class ResidueIndexDifference(ResiduePair):
    """Difference in sequence index between a pair of residues. This can
       either be the simple difference (if *absolute* is False) in which case
       the feature is asymmetric, or the absolute value (if *absolute* is True)
       which gives a symmetric feature."""
    def __init__(self, mlib, bins, protein=0, absolute=False):
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, protein, absolute)
        self._create_bins(mlib, bins)
    _setup = _mdt.mdt_feature_residue_index_difference


class PhiDihedralDifference(AlignedResidue):
    """Difference in phi dihedral between a pair of aligned residues.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_phi_dihedral_difference


class PsiDihedralDifference(AlignedResidue):
    """Difference in psi dihedral between a pair of aligned residues.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_psi_dihedral_difference


class OmegaDihedralDifference(AlignedResidue):
    """Difference in omega dihedral between a pair of aligned residues.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_omega_dihedral_difference


class NeighborhoodDifference(AlignedResidue):
    """Residue neighborhood difference. This is the average of the distance
       scores (from a residue-residue scoring matrix) of all aligned residues
       where the residue in the first sequence is within a cutoff distance
       of the scanned residue. (This cutoff is set by the distngh argument to
       :meth:`mdt.Table.add_alignment`.)"""
    _setup = _mdt.mdt_feature_neighborhood_difference


class GapDistance(AlignedResidue):
    """Distance, in alignment positions, to the nearest gap. Note that
       positions which are gapped in both sequences are ignored for the
       purposes of this calculation (a 'gap' is defined as a gap in one
       sequence aligned with a residue in the other)."""
    _setup = _mdt.mdt_feature_gap_distance


class ResidueDistanceDifference(AlignedResiduePair):
    """Distance between two residues in the second protein, minus the distance
       between the equivalent residues in the first protein.
       See :class:`ResidueDistance`.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_residue_distance_difference


class AverageNeighborhoodDifference(AlignedResiduePair):
    """Average residue neighborhood difference for a pair of alignment
       positions. See :class:`NeighborhoodDifference`."""
    _setup = _mdt.mdt_feature_average_neighborhood_difference


class AverageGapDistance(AlignedResiduePair):
    """Average distance to a gap from a pair of alignment positions.
       See :class:`GapDistance`."""
    _setup = _mdt.mdt_feature_average_gap_distance


class AtomAccessibility(Atom):
    """Atom solvent accessibility. This is calculated by the PSA algorithm,
       and controlled by the surftyp and accessibility_type arguments to
       :meth:`mdt.Table.add_alignment`.
       The feature is considered undefined if the atom's Cartesian
       coordinates are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_atom_accessibility


class AtomZCoordinate(Atom):
    """Atom Z-coordinate. No orientation of the structure is performed.
       The feature is considered undefined if the coordinate is equal
       to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_z_coordinate


class FractionalAtomAccessibility(Atom):
    """Fractional atom solvent accessibility, from 0 to 1. This is the atom
       solvent accessibility (see :class:`AtomAccessibility`) divided by
       the volume of the atom, derived from its van der Waals radius.
       The feature is considered undefined if the atom's Cartesian
       coordinates are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_fractional_atom_accessibility


class AtomType(AtomFixedBins):
    """Type of an atom, as classified by the atom class file.
       See :attr:`mdt.Library.atom_classes`."""
    _setup = _mdt.mdt_feature_atom_type


class HydrogenBondDonor(Atom):
    """Number of hydrogen bond donors. It is defined as the sum, over all atoms
       within hbond_cutoff (see :class:`mdt.Library`) of the atom, of their
       donor valencies as defined in the hydrogen bond file
       (see :attr:`mdt.Library.hbond_classes`).
       The feature is considered undefined if the atom's Cartesian
       coordinates are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_hydrogen_bond_donor


class HydrogenBondAcceptor(Atom):
    """Number of hydrogen bond acceptors. It is defined as the sum, over all
       atoms within hbond_cutoff (see :class:`mdt.Library`) of the atom,
       of their acceptor valencies as defined in the hydrogen bond file
       (see :attr:`mdt.Library.hbond_classes`).
       The feature is considered undefined if the atom's Cartesian
       coordinates are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_hydrogen_bond_acceptor


class HydrogenBondCharge(Atom):
    """Hydrogen bond charge. It is defined as the sum, over all
       atoms within hbond_cutoff (see :class:`mdt.Library`) of the atom,
       of their charges as defined in the hydrogen bond file (see
       :attr:`mdt.Library.hbond_classes`)."""
    _setup = _mdt.mdt_feature_hydrogen_bond_charge


class AtomTable(Atom):
    """A tabulated atom feature. The feature is simply a table of N
       floating-point numbers, where N is the number of atoms in the system.
       This table is provided by a Python function, so can be used to implement
       user-defined features or to pass in features from other software.
       A simple example to use the x coordinate as a feature::

          def func(aln, struc, mlib, libs):
              return [a.x for a in struc.atoms]
          f = mdt.features.AtomTable(mlib, bins, "x coordinate", func)
    """
    _setup = _mdt.mdt_feature_atom_table

    class _SizeCheck(object):
        """Make sure the function returns a sequence of the right length"""
        def __init__(self, func):
            self.func = func

        def __call__(self, aln, iseq, mlib, libs):
            prop = self.func(aln, aln[iseq], mlib, libs)
            if len(prop) != len(aln[iseq].atoms):
                raise ValueError("Should return a sequence of length %d "
                                 "(number of atoms in the structure)"
                                 % len(aln[iseq].atoms))
            return prop

    def __init__(self, mlib, bins, table_name, func, pos2=False):
        """
        :Parameters:
          - `table_name`: the name of the feature that is tabulated
          - `func`: A Python function or other callable, which is expected
                    to return a sequence of floats, one per atom

        See :class:`Atom` for a description of the other arguments.
        """
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, pos2, table_name,
                                  self._SizeCheck(func))
        self._create_bins(mlib, bins)


class AtomDistance(AtomPair):
    """Distance in angstroms between a pair of atoms.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_atom_distance


class AtomBondSeparation(AtomPair):
    """Number of bonds between a pair of atoms.
       For example, two atoms that are directly bonded return '1', while two
       at opposite ends of an angle return '2'. The bonds between atoms in
       each standard amino acid are derived from the bond class file, so this
       must be read in first (see :attr:`mdt.Library.bond_classes`). For
       atoms in different residues, the residues are assumed to be linked by
       a peptide backbone, and the number of bonds is calculated accordingly.
       Atoms in different chains, or atoms of types not referenced in the bond
       class file, are not connected. If disulfide is set to True, disulfide
       bridges are also considered (if two residues have SG atoms within 2.5
       angstroms, they are counted as bonded). If disulfide is set to False
       (the default) any disulfide bridges are ignored. Either way, no account
       is taken of patches and other modifications such as terminal oxygens
       (unless bonds to OXT are explicitly listed in the bond class file).
       If a pair of atoms is not connected it is placed in the 'undefined'
       bin."""
    _setup = _mdt.mdt_feature_atom_bond_separation

    def __init__(self, mlib, bins, disulfide=False):
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, disulfide)
        self._create_bins(mlib, bins)


class HydrogenBondSatisfaction(Protein):
    """Hydrogen bond satisfaction index for a protein. This is the average
       difference, over all atoms in the protein, between the HydrogenBondDonor
       value and the atom's donor valency plus the same for the acceptor,
       as defined in the hydrogen bond file (see
       :attr:`mdt.Library.hbond_classes`)."""
    _setup = _mdt.mdt_feature_hydrogen_bond_satisfaction


class AlphaContent(Protein):
    """Alpha content of the protein. This is simply the fraction, between 0
       and 1, of residues in the first mainchain conformation class
       (see :class:`MainchainConformation`)."""
    _setup = _mdt.mdt_feature_alpha_content


class SequenceIdentity(ProteinPair):
    """Fractional sequence identity, between 0 and 1, between two sequences.
       This is the number of identical aligned residues divided by the length
       of the shorter sequence."""
    _setup = _mdt.mdt_feature_sequence_identity


class TupleType(TupleFixedBins):
    """Type of an atom tuple, as classified by the tuple class file.
       See :attr:`mdt.Library.tuple_classes`."""
    _setup = _mdt.mdt_feature_tuple_type


class TupleDistance(TuplePair):
    """Distance in angstroms between the first atom in each of two tuples
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_tuple_distance


class TupleAngle1(TuplePair):
    """Angle (0-180) between the first atom in the first tuple, the first atom
       in the second tuple, and the second atom in the second tuple.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_tuple_angle1


class TupleAngle2(TuplePair):
    """Angle (0-180) between the second atom in the first tuple, the first atom
       in the first tuple, and the first atom in the second tuple.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_tuple_angle2


class TupleDihedral1(TuplePair):
    """Dihedral (-180-180) between the second atom in the first tuple, the
       first atom in the first tuple, the first atom in the second tuple, and
       the second atom in the second tuple.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_tuple_dihedral1


class TupleDihedral2(TuplePair):
    """Dihedral (-180-180) between the third atom in the first tuple, the
       second atom in the first tuple, the first atom in the first tuple, and
       the first atom in the second tuple. Only works for atom triplets.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_tuple_dihedral2


class TupleDihedral3(TuplePair):
    """Dihedral (-180-180) between the first atom in the first tuple, the
       first atom in the second tuple, the second atom in the second tuple, and
       the third atom in the second tuple. Only works for atom triplets.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_tuple_dihedral3


class BondType(ChemicalBondFixedBins):
    """Type of a bond, as classified by the bond class file.
       See :attr:`mdt.Library.bond_classes`."""
    _setup = _mdt.mdt_feature_bond_type


class AngleType(ChemicalBondFixedBins):
    """Type of an angle, as classified by the angle class file.
       See :attr:`mdt.Library.angle_classes`."""
    _setup = _mdt.mdt_feature_angle_type


class DihedralType(ChemicalBondFixedBins):
    """Type of a dihedral, as classified by the dihedral class file.
       See :attr:`mdt.Library.dihedral_classes`."""
    _setup = _mdt.mdt_feature_dihedral_type


class BondLength(ChemicalBond):
    """Length of a bond in angstroms. See :attr:`mdt.Library.bond_classes`.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_bond_length


class Angle(ChemicalBond):
    """Angle (0-180). See :attr:`mdt.Library.angle_classes`.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_angle


class Dihedral(ChemicalBond):
    """Dihedral angle (-180-180). See :attr:`mdt.Library.dihedral_classes`.
       The feature is considered undefined if any of the atom coordinates
       are equal to the Modeller 'undefined' value (-999.0)."""
    _setup = _mdt.mdt_feature_dihedral


class Group(_Base):
    """A feature that groups other features."""
    def __init__(self, mlib, feat1, feat2, nbins):
        """
        Create the Group feature.

        :Parameters:
          - `mlib`: the `Library` to create the feature in.
          - `feat1`: an existing feature object that will be included
                     in this group.
          - `feat2`: another existing feature object to include.
          - `nbins`: the number of bins in this feature.
        """
        _Base.__init__(self, mlib)
        self._ifeat = self._setup(mlib._modpt, feat1._get_ifeat(mlib),
                                  feat2._get_ifeat(mlib), nbins)


class Cluster(Group):
    """Cluster feature. When evaluated, it evaluates the two other features
       grouped in this feature, and converts the pair of bin indices for
       those features into a single bin index, which is returned. Use the
       :meth:`add` method to control this conversion."""
    _setup = _mdt.mdt_feature_cluster

    def add(self, child_bins, bin_index):
        """Add a single mapping from a pair of child feature bin indices into
           this feature's bin index (all indexes start at 0). For example,
           calling `add((1,2), 3)` would cause this Cluster feature to return
           bin index 3 if the child features were in bins 1 and 2 respectively.
           This method can be called multiple times (even for the same
           `bin_index`) to add additional mappings from child bin indices
           to bin index. If no mapping from a given pair of child indices is
           present, the undefined bin index is returned."""
        _mdt.mdt_cluster_add(self._mlib._modpt, self._ifeat, child_bins[0],
                             child_bins[1], bin_index)
