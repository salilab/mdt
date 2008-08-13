"""MDT features."""

__docformat__ = "restructuredtext"

import _mdt

class _Base(object):
    def _create_bins(self, mlib, bins):
        _mdt.mdt_feature_nbins_set(mlib._modpt, self._ifeat, len(bins))
        for i, (start, end, symbol) in enumerate(bins):
            _mdt.mdt_feature_bin_set(mlib._modpt, self._ifeat, i, start,
                                     end, symbol)


class Protein(_Base):
    def __init__(self, mlib, bins, protein=0):
        self._ifeat = self._setup(mlib._modpt, protein)
        self._create_bins(mlib, bins)


class Residue(_Base):
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
    def __init__(self, mlib, bins, protein=0):
        self._ifeat = self._setup(mlib._modpt, protein)
        self._create_bins(mlib, bins)


class AlignedResidue(_Base):
    def __init__(self, mlib, bins, protein1=0, protein2=1):
        self._ifeat = self._setup(mlib._modpt, protein1, protein2)
        self._create_bins(mlib, bins)


class AlignedResiduePair(_Base):
    def __init__(self, mlib, bins, protein1=0, protein2=1):
        self._ifeat = self._setup(mlib._modpt, protein1, protein2)
        self._create_bins(mlib, bins)


class Atom(_Base):
    def __init__(self, mlib, bins, pos2=False):
        self._ifeat = self._setup(mlib._modpt, pos2)
        self._create_bins(mlib, bins)


class AtomFixedBins(Atom):
    def __init__(self, mlib, pos2=False):
        self._ifeat = self._setup(mlib._modpt, pos2)


class AtomPair(_Base):
    def __init__(self, mlib, bins):
        self._ifeat = self._setup(mlib._modpt)
        self._create_bins(mlib, bins)


class Tuple(_Base):
    def __init__(self, mlib, bins, pos2=False):
        self._ifeat = self._setup(mlib._modpt, pos2)
        self._create_bins(mlib, bins)


class TupleFixedBins(Tuple):
    def __init__(self, mlib, pos2=False):
        self._ifeat = self._setup(mlib._modpt, pos2)


class TuplePair(_Base):
    def __init__(self, mlib, bins):
        self._ifeat = self._setup(mlib._modpt)
        self._create_bins(mlib, bins)


class ChemicalBond(_Base):
    def __init__(self, mlib, bins):
        self._ifeat = self._setup(mlib._modpt)
        self._create_bins(mlib, bins)


class ChemicalBondFixedBins(_Base):
    def __init__(self, mlib):
        self._ifeat = self._setup(mlib._modpt)


class XRayResolution(Protein):
    _setup = _mdt.mdt_feature_xray_resolution

class RadiusOfGyration(Protein):
    _setup = _mdt.mdt_feature_radius_of_gyration

class SequenceLength(Protein):
    _setup = _mdt.mdt_feature_sequence_length

class ResidueType(ResidueFixedBins):
    _setup = _mdt.mdt_feature_residue_type

class ResidueAccessibility(Residue):
    _setup = _mdt.mdt_feature_residue_accessibility

class Chi1Dihedral(Residue):
    _setup = _mdt.mdt_feature_chi1_dihedral

class Chi2Dihedral(Residue):
    _setup = _mdt.mdt_feature_chi2_dihedral

class Chi3Dihedral(Residue):
    _setup = _mdt.mdt_feature_chi3_dihedral

class Chi4Dihedral(Residue):
    _setup = _mdt.mdt_feature_chi4_dihedral

class PhiDihedral(Residue):
    _setup = _mdt.mdt_feature_phi_dihedral

class PsiDihedral(Residue):
    _setup = _mdt.mdt_feature_psi_dihedral

class OmegaDihedral(Residue):
    _setup = _mdt.mdt_feature_omega_dihedral

class AlphaDihedral(Residue):
    _setup = _mdt.mdt_feature_alpha_dihedral

class Chi1Class(ResidueFixedBins):
    _setup = _mdt.mdt_feature_chi1_class

class Chi2Class(ResidueFixedBins):
    _setup = _mdt.mdt_feature_chi2_class

class Chi3Class(ResidueFixedBins):
    _setup = _mdt.mdt_feature_chi3_class

class Chi4Class(ResidueFixedBins):
    _setup = _mdt.mdt_feature_chi4_class

class Chi5Class(ResidueFixedBins):
    _setup = _mdt.mdt_feature_chi5_class

class PhiClass(ResidueFixedBins):
    _setup = _mdt.mdt_feature_phi_class

class PsiClass(ResidueFixedBins):
    _setup = _mdt.mdt_feature_psi_class

class OmegaClass(ResidueFixedBins):
    _setup = _mdt.mdt_feature_omega_class

class MainchainConformation(ResidueFixedBins):
    _setup = _mdt.mdt_feature_mainchain_conformation

class ResidueDistance(ResiduePair):
    _setup = _mdt.mdt_feature_residue_distance

class AverageResidueAccessibility(ResiduePair):
    _setup = _mdt.mdt_feature_average_residue_accessibility

class ResidueIndexDifference(ResiduePair):
    _setup = _mdt.mdt_feature_residue_index_difference

class PhiDihedralDifference(AlignedResidue):
    _setup = _mdt.mdt_feature_phi_dihedral_difference

class PsiDihedralDifference(AlignedResidue):
    _setup = _mdt.mdt_feature_psi_dihedral_difference

class OmegaDihedralDifference(AlignedResidue):
    _setup = _mdt.mdt_feature_omega_dihedral_difference

class ResidueDistanceDifference(AlignedResiduePair):
    _setup = _mdt.mdt_feature_residue_distance_difference

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
