"""MDT features."""

__docformat__ = "restructuredtext"

import _mdt

class _Base(object):
    """Base class for all features."""
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
        self._ifeat = self._setup(mlib._modpt, protein, delta, align_delta,
                                  pos2)
        self._create_bins(mlib, bins)


class ResidueFixedBins(Residue):
    """A residue feature for which bins need not be specified."""
    def __init__(self, mlib, protein=0, delta=0, align_delta=0, pos2=False):
        """See `Residue` for a description of the arguments."""
        self._ifeat = self._setup(mlib._modpt, protein, delta, align_delta,
                                  pos2, mlib._env.libs.modpt)


class ResiduePair(_Base):
    """A feature defined on pairs of residues in a protein."""
    def __init__(self, mlib, bins, protein=0):
        """See `Protein` for a description of the arguments."""
        self._ifeat = self._setup(mlib._modpt, protein)
        self._create_bins(mlib, bins)


class AlignedResidue(_Base):
    """A feature defined on residues aligned between two proteins.
       For each pair of proteins, every alignment position is scanned, and the
       feature is evaluated for each pair of aligned residues."""
    def __init__(self, mlib, bins, protein1=0, protein2=1):
        """See `ResiduePair` for a description of the arguments."""
        self._ifeat = self._setup(mlib._modpt, protein1, protein2)
        self._create_bins(mlib, bins)


class AlignedResiduePair(_Base):
    """A feature defined on pairs of residues aligned between two proteins.
       For each pair of proteins, each pair of alignment positions is scanned,
       and the feature is evaluated for each pair of pairs of aligned
       residues."""
    def __init__(self, mlib, bins, protein1=0, protein2=1):
        """See `ResiduePair` for a description of the arguments."""
        self._ifeat = self._setup(mlib._modpt, protein1, protein2)
        self._create_bins(mlib, bins)


class Atom(_Base):
    """A feature defined on single atoms.
       This works only on the first protein in each group of proteins
       selected from the alignment."""
    def __init__(self, mlib, bins, pos2=False):
        """See `Residue` for a description of the arguments."""
        self._ifeat = self._setup(mlib._modpt, pos2)
        self._create_bins(mlib, bins)


class AtomFixedBins(Atom):
    """An atom feature for which bins need not be specified."""
    def __init__(self, mlib, pos2=False):
        """See `Residue` for a description of the arguments."""
        self._ifeat = self._setup(mlib._modpt, pos2)


class AtomPair(_Base):
    """A feature defined on pairs of atoms.
       This works only on the first protein in each group of proteins
       selected from the alignment."""
    def __init__(self, mlib, bins):
        """See `Protein` for a description of the arguments."""
        self._ifeat = self._setup(mlib._modpt)
        self._create_bins(mlib, bins)


class Tuple(_Base):
    """A feature defined on tuples of atoms.
       This works only on the first protein in each group of proteins
       selected from the alignment."""
    def __init__(self, mlib, bins, pos2=False):
        """See `Residue` for a description of the arguments."""
        self._ifeat = self._setup(mlib._modpt, pos2)
        self._create_bins(mlib, bins)


class TupleFixedBins(Tuple):
    """A tuple feature for which bins need not be specified."""
    def __init__(self, mlib, pos2=False):
        """See `Residue` for a description of the arguments."""
        self._ifeat = self._setup(mlib._modpt, pos2)


class TuplePair(_Base):
    """A feature defined on pairs of tuples of atoms.
       This works only on the first protein in each group of proteins
       selected from the alignment."""
    def __init__(self, mlib, bins):
        """See `Protein` for a description of the arguments."""
        self._ifeat = self._setup(mlib._modpt)
        self._create_bins(mlib, bins)


class ChemicalBond(_Base):
    """A feature defined on defined chemical bonds, angles, or dihedrals.
       The definitions of the chemical connectivity must first be read from
       a bond class file; see `BondClasses`. This feature works only on the
       first protein in each group of proteins selected from the alignment."""
    def __init__(self, mlib, bins):
        """See `Protein` for a description of the arguments."""
        self._ifeat = self._setup(mlib._modpt)
        self._create_bins(mlib, bins)


class ChemicalBondFixedBins(ChemicalBond):
    """A chemical bond feature for which bins need not be specified."""
    def __init__(self, mlib):
        self._ifeat = self._setup(mlib._modpt)


class XRayResolution(Protein):
    """
    Protein X-ray resolution in angstroms.
    Note that proteins with a resolution of -1.00 (generally NMR structures) are
    actually treated as 0.45 by default. This decreases the number of bins
    required to hold all defined resolutions while still separating NMR from
    X-ray structures.
    """
    _setup = _mdt.mdt_feature_xray_resolution
    def __init__(self, mlib, bins, protein=0, nmr=0.45):
        """Create the feature. `nmr` is the resolution to artificially map
           proteins with resolution -1.00 to. See `Protein` for the other
           arguments."""
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
       accessibility; see `AtomAccessibility`."""
    _setup = _mdt.mdt_feature_residue_accessibility

class Chi1Dihedral(Residue):
    """Residue chi1 dihedral angle, from -180 to 180 degrees."""
    _setup = _mdt.mdt_feature_chi1_dihedral

class Chi2Dihedral(Residue):
    """Residue chi2 dihedral angle, from -180 to 180 degrees."""
    _setup = _mdt.mdt_feature_chi2_dihedral

class Chi3Dihedral(Residue):
    """Residue chi3 dihedral angle, from -180 to 180 degrees."""
    _setup = _mdt.mdt_feature_chi3_dihedral

class Chi4Dihedral(Residue):
    """Residue chi4 dihedral angle, from -180 to 180 degrees."""
    _setup = _mdt.mdt_feature_chi4_dihedral

class PhiDihedral(Residue):
    """Residue phi dihedral angle, from -180 to 180 degrees."""
    _setup = _mdt.mdt_feature_phi_dihedral

class PsiDihedral(Residue):
    """Residue psi dihedral angle, from -180 to 180 degrees."""
    _setup = _mdt.mdt_feature_psi_dihedral

class OmegaDihedral(Residue):
    """Residue omega dihedral angle, from -180 to 180 degrees."""
    _setup = _mdt.mdt_feature_omega_dihedral

class AlphaDihedral(Residue):
    """Residue alpha dihedral angle, from -180 to 180 degrees."""
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
    """Residue mainchain conformation (Ramachandran) class."""
    _setup = _mdt.mdt_feature_mainchain_conformation

class ResidueGroup(ResidueFixedBins):
    """Residue group."""
    _setup = _mdt.mdt_feature_residue_group
    def __init__(self, mlib, protein=0, delta=0, align_delta=0, pos2=False,
                 residue_grouping=0):
        self._ifeat = self._setup(mlib._modpt, protein, delta, align_delta,
                                  pos2, residue_grouping, mlib._env.libs.modpt)

class SidechainBiso(Residue):
    """
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
       `Library` object."""
    _setup = _mdt.mdt_feature_residue_distance

class AverageResidueAccessibility(ResiduePair):
    """Average solvent accessibility of a pair of residues.
       See `ResidueAccessibility`."""
    _setup = _mdt.mdt_feature_average_residue_accessibility

class ResidueIndexDifference(ResiduePair):
    """Difference in sequence index between a pair of residues. This can
       either be the simple difference, or the absolute value."""
    def __init__(self, mlib, bins, protein=0, absolute=False):
        """Create the feature. If `absolute` is True, the absolute value of
           the index difference is used (and thus the feature is symmetric).
           See `ResiduePair` for a description of the other arguments."""
        self._ifeat = self._setup(mlib._modpt, protein, absolute)
        self._create_bins(mlib, bins)
    _setup = _mdt.mdt_feature_residue_index_difference

class PhiDihedralDifference(AlignedResidue):
    """Difference in phi dihedral between a pair of aligned residues."""
    _setup = _mdt.mdt_feature_phi_dihedral_difference

class PsiDihedralDifference(AlignedResidue):
    """Difference in psi dihedral between a pair of aligned residues."""
    _setup = _mdt.mdt_feature_psi_dihedral_difference

class OmegaDihedralDifference(AlignedResidue):
    """Difference in omega dihedral between a pair of aligned residues."""
    _setup = _mdt.mdt_feature_omega_dihedral_difference

class NeighborhoodDifference(AlignedResidue):
    """Residue neighborhood difference. This is the average of the distance
       scores (from a residue-residue scoring matrix) of all aligned residues
       where the residue in the first sequence is within a cutoff distance
       of the scanned residue. (This cutoff is set by the distngh argument to
       `Table.add_alignment`.)"""
    _setup = _mdt.mdt_feature_neighborhood_difference

class GapDistance(AlignedResidue):
    """Distance, in alignment positions, to the nearest gap. Note that positions       which are gapped in both sequences are ignored for the purposes of this
       calculation (a 'gap' is defined as a gap in one sequence aligned with a
       residue in the other)."""
    _setup = _mdt.mdt_feature_gap_distance

class ResidueDistanceDifference(AlignedResiduePair):
    """Distance between two residues in the second protein, minus the distance
       between the equivalent residues in the first protein.
       See `ResidueDistance`."""
    _setup = _mdt.mdt_feature_residue_distance_difference

class AverageNeighborhoodDifference(AlignedResiduePair):
    """Average residue neighborhood difference for a pair of alignment
       positions. See `NeighborhoodDifference`."""
    _setup = _mdt.mdt_feature_average_neighborhood_difference

class AverageGapDistance(AlignedResiduePair):
    """Average distance to a gap from a pair of alignment positions.
       See `GapDistance`."""
    _setup = _mdt.mdt_feature_average_gap_distance

class AtomAccessibility(Atom):
    """Atom solvent accessibility. This is calculated by the PSA algorithm,
       and controlled by the surftyp and accessibility_type arguments to
       `Table.add_alignment`."""
    _setup = _mdt.mdt_feature_atom_accessibility

class FractionalAtomAccessibility(Atom):
    """Fractional atom solvent accessibility, from 0 to 1. This is the atom
       solvent accessibility (see `AtomAccessibility`) divided by the volume
       of the atom, derived from its van der Waals radius."""
    _setup = _mdt.mdt_feature_fractional_atom_accessibility

class AtomType(AtomFixedBins):
    """Type of an atom, as classified by the atom class file.
       See `Library.atom_classes`."""
    _setup = _mdt.mdt_feature_atom_type

class HydrogenBondDonor(Atom):
    """Number of hydrogen bond donors. It is defined as the sum, over all atoms
       within hbond_cutoff (see `Library`) of the atom, of their donor
       valencies as defined in the hydrogen bond file
       (see `Library.hbond_classes`)."""
    _setup = _mdt.mdt_feature_hydrogen_bond_donor

class HydrogenBondAcceptor(Atom):
    """Number of hydrogen bond acceptors. It is defined as the sum, over all
       atoms within hbond_cutoff (see `Library`) of the atom, of their acceptor
       valencies as defined in the hydrogen bond file
       (see `Library.hbond_classes`)."""
    _setup = _mdt.mdt_feature_hydrogen_bond_acceptor

class HydrogenBondCharge(Atom):
    """Hydrogen bond charge. It is defined as the sum, over all
       atoms within hbond_cutoff (see `Library`) of the atom, of their charges
       as defined in the hydrogen bond file (see `Library.hbond_classes`)."""
    _setup = _mdt.mdt_feature_hydrogen_bond_charge

class AtomDistance(AtomPair):
    """Distance in angstroms between a pair of atoms."""
    _setup = _mdt.mdt_feature_atom_distance

class HydrogenBondSatisfaction(Protein):
    """Hydrogen bond satisfaction index for a protein. This is the average
       difference, over all atoms in the protein, between the HydrogenBondDonor
       value and the atom's donor valency plus the same for the acceptor,
       as defined in the hydrogen bond file (see `Library.hbond_classes`)."""
    _setup = _mdt.mdt_feature_hydrogen_bond_satisfaction

class AlphaContent(Protein):
    """Alpha content of the protein. This is simply the fraction, between 0
       and 1, of residues in the first mainchain conformation class
       (see `MainchainConformation`)."""
    _setup = _mdt.mdt_feature_alpha_content

class SequenceIdentity(ProteinPair):
    """Fractional sequence identity, between 0 and 1, between two sequences.
       This is the number of identical aligned residues divided by the length
       of the shorter sequence."""
    _setup = _mdt.mdt_feature_sequence_identity

class TupleType(TupleFixedBins):
    """Type of an atom tuple, as classified by the tuple class file.
       See `Library.tuple_classes`."""
    _setup = _mdt.mdt_feature_tuple_type

class TupleDistance(TuplePair):
    """Distance in angstroms between the first atom in each of two tuples."""
    _setup = _mdt.mdt_feature_tuple_distance

class TupleAngle1(TuplePair):
    """Angle (0-180) between the first atom in the first tuple, the first atom
       in the second tuple, and the second atom in the second tuple."""
    _setup = _mdt.mdt_feature_tuple_angle1

class TupleAngle2(TuplePair):
    """Angle (0-180) between the second atom in the first tuple, the first atom
       in the first tuple, and the first atom in the second tuple."""
    _setup = _mdt.mdt_feature_tuple_angle2

class TupleDihedral1(TuplePair):
    """Dihedral (-180-180) between the second atom in the first tuple, the
       first atom in the first tuple, the first atom in the second tuple, and
       the second atom in the second tuple."""
    _setup = _mdt.mdt_feature_tuple_dihedral1

class TupleDihedral2(TuplePair):
    """Dihedral (-180-180) between the third atom in the first tuple, the
       second atom in the second tuple, the first atom in the first tuple, and
       the first atom in the second tuple. Only works for atom triplets."""
    _setup = _mdt.mdt_feature_tuple_dihedral2

class TupleDihedral3(TuplePair):
    """Dihedral (-180-180) between the first atom in the first tuple, the
       second atom in the first tuple, the second atom in the first tuple, and
       the third atom in the second tuple. Only works for atom triplets."""
    _setup = _mdt.mdt_feature_tuple_dihedral3

class BondType(ChemicalBondFixedBins):
    """Type of a bond, as classified by the bond class file.
       See `Library.bond_classes`."""
    _setup = _mdt.mdt_feature_bond_type

class AngleType(ChemicalBondFixedBins):
    """Type of an angle, as classified by the angle class file.
       See `Library.angle_classes`."""
    _setup = _mdt.mdt_feature_angle_type

class DihedralType(ChemicalBondFixedBins):
    """Type of a dihedral, as classified by the dihedral class file.
       See `Library.dihedral_classes`."""
    _setup = _mdt.mdt_feature_dihedral_type

class BondLength(ChemicalBond):
    """Length of a bond in angstroms. See `Library.bond_classes`."""
    _setup = _mdt.mdt_feature_bond_length

class Angle(ChemicalBond):
    """Angle (0-180). See `Library.angle_classes`."""
    _setup = _mdt.mdt_feature_angle

class Dihedral(ChemicalBond):
    """Dihedral angle (-180-180). See `Library.dihedral_classes`."""
    _setup = _mdt.mdt_feature_dihedral
