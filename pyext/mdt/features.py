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
        self._ifeat = self._setup(mlib._modpt, protein)
        self._create_bins(mlib, bins)


class ProteinPair(_Base):
    """A feature defined on pairs of proteins."""
    def __init__(self, mlib, bins, protein1=0, protein2=1):
        self._ifeat = self._setup(mlib._modpt, protein1, protein2)
        self._create_bins(mlib, bins)


class Residue(_Base):
    """A feature defined on single residues in a protein."""
    def __init__(self, mlib, bins, protein=0, delta=0, align_delta=0,
                 pos2=False):
        self._ifeat = self._setup(mlib._modpt, protein, delta, align_delta,
                                  pos2)
        self._create_bins(mlib, bins)


class ResidueFixedBins(_Base):
    def __init__(self, mlib, protein=0, delta=0, align_delta=0, pos2=False):
        self._ifeat = self._setup(mlib._modpt, protein, delta, align_delta,
                                  pos2, mlib._env.libs.modpt)


class ResiduePair(_Base):
    """A feature defined on pairs of residues in a protein."""
    def __init__(self, mlib, bins, protein=0):
        self._ifeat = self._setup(mlib._modpt, protein)
        self._create_bins(mlib, bins)


class AlignedResidue(_Base):
    """A feature defined on residues aligned between two proteins."""
    def __init__(self, mlib, bins, protein1=0, protein2=1):
        self._ifeat = self._setup(mlib._modpt, protein1, protein2)
        self._create_bins(mlib, bins)


class AlignedResiduePair(_Base):
    """A feature defined on pairs of residues aligned between two proteins."""
    def __init__(self, mlib, bins, protein1=0, protein2=1):
        self._ifeat = self._setup(mlib._modpt, protein1, protein2)
        self._create_bins(mlib, bins)


class Atom(_Base):
    """A feature defined on single atoms.
       This works only on the first protein in the alignment."""
    def __init__(self, mlib, bins, pos2=False):
        self._ifeat = self._setup(mlib._modpt, pos2)
        self._create_bins(mlib, bins)


class AtomFixedBins(Atom):
    def __init__(self, mlib, pos2=False):
        self._ifeat = self._setup(mlib._modpt, pos2)


class AtomPair(_Base):
    """A feature defined on pairs of atoms.
       This works only on the first protein in the alignment."""
    def __init__(self, mlib, bins):
        self._ifeat = self._setup(mlib._modpt)
        self._create_bins(mlib, bins)


class Tuple(_Base):
    """A feature defined on tuples of atoms.
       This works only on the first protein in the alignment."""
    def __init__(self, mlib, bins, pos2=False):
        self._ifeat = self._setup(mlib._modpt, pos2)
        self._create_bins(mlib, bins)


class TupleFixedBins(Tuple):
    def __init__(self, mlib, pos2=False):
        self._ifeat = self._setup(mlib._modpt, pos2)


class TuplePair(_Base):
    """A feature defined on pairs of tuples of atoms.
       This works only on the first protein in the alignment."""
    def __init__(self, mlib, bins):
        self._ifeat = self._setup(mlib._modpt)
        self._create_bins(mlib, bins)


class ChemicalBond(_Base):
    """A feature defined on chemical bonds.
       This works only on the first protein in the alignment."""
    def __init__(self, mlib, bins):
        self._ifeat = self._setup(mlib._modpt)
        self._create_bins(mlib, bins)


class ChemicalBondFixedBins(_Base):
    def __init__(self, mlib):
        self._ifeat = self._setup(mlib._modpt)


class XRayResolution(Protein):
    """Protein X-ray resolution."""
    _setup = _mdt.mdt_feature_xray_resolution

class RadiusOfGyration(Protein):
    """Protein radius of gyration."""
    _setup = _mdt.mdt_feature_radius_of_gyration

class SequenceLength(Protein):
    """Protein sequence length (number of residues)."""
    _setup = _mdt.mdt_feature_sequence_length

class ResidueType(ResidueFixedBins):
    """Residue type (20 standard amino acids, gap, undefined)."""
    _setup = _mdt.mdt_feature_residue_type

class ResidueAccessibility(Residue):
    """Residue solvent accessibility."""
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
    """Residue average sidechain Biso."""
    _setup = _mdt.mdt_feature_sidechain_biso

class ResidueDistance(ResiduePair):
    """Distance between a pair of residues."""
    _setup = _mdt.mdt_feature_residue_distance

class AverageResidueAccessibility(ResiduePair):
    """Average solvent accessibility of a pair of residues."""
    _setup = _mdt.mdt_feature_average_residue_accessibility

class ResidueIndexDifference(ResiduePair):
    """Difference in sequence index between a pair of residues.
       Note that this can be positive or negative."""
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
    _setup = _mdt.mdt_feature_neighborhood_difference

class GapDistance(AlignedResidue):
    _setup = _mdt.mdt_feature_gap_distance

class ResidueDistanceDifference(AlignedResiduePair):
    _setup = _mdt.mdt_feature_residue_distance_difference

class AverageNeighborhoodDifference(AlignedResiduePair):
    _setup = _mdt.mdt_feature_average_neighborhood_difference

class AverageGapDistance(AlignedResiduePair):
    _setup = _mdt.mdt_feature_average_gap_distance

class AtomAccessibility(Atom):
    _setup = _mdt.mdt_feature_atom_accessibility

class FractionalAtomAccessibility(Atom):
    _setup = _mdt.mdt_feature_fractional_atom_accessibility

class AtomType(AtomFixedBins):
    _setup = _mdt.mdt_feature_atom_type

class HydrogenBondDonor(Atom):
    _setup = _mdt.mdt_feature_hydrogen_bond_donor

class HydrogenBondAcceptor(Atom):
    _setup = _mdt.mdt_feature_hydrogen_bond_acceptor

class HydrogenBondCharge(Atom):
    _setup = _mdt.mdt_feature_hydrogen_bond_charge

class AtomDistance(AtomPair):
    _setup = _mdt.mdt_feature_atom_distance

class HydrogenBondSatisfaction(Protein):
    _setup = _mdt.mdt_feature_hydrogen_bond_satisfaction

class AlphaContent(Protein):
    _setup = _mdt.mdt_feature_alpha_content

class SequenceIdentity(ProteinPair):
    _setup = _mdt.mdt_feature_sequence_identity

class TupleType(TupleFixedBins):
    _setup = _mdt.mdt_feature_tuple_type

class TupleDistance(TuplePair):
    _setup = _mdt.mdt_feature_tuple_distance

class TupleAngle1(TuplePair):
    _setup = _mdt.mdt_feature_tuple_angle1

class TupleAngle2(TuplePair):
    _setup = _mdt.mdt_feature_tuple_angle2

class TupleDihedral1(TuplePair):
    _setup = _mdt.mdt_feature_tuple_dihedral1

class TupleDihedral2(TuplePair):
    _setup = _mdt.mdt_feature_tuple_dihedral2

class TupleDihedral3(TuplePair):
    _setup = _mdt.mdt_feature_tuple_dihedral3

class BondType(ChemicalBondFixedBins):
    _setup = _mdt.mdt_feature_bond_type

class AngleType(ChemicalBondFixedBins):
    _setup = _mdt.mdt_feature_angle_type

class DihedralType(ChemicalBondFixedBins):
    _setup = _mdt.mdt_feature_dihedral_type

class BondLength(ChemicalBond):
    _setup = _mdt.mdt_feature_bond_length

class Angle(ChemicalBond):
    _setup = _mdt.mdt_feature_angle

class Dihedral(ChemicalBond):
    _setup = _mdt.mdt_feature_dihedral
