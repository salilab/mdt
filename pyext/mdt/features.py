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
    def __init__(self, mlib, bins, protein=0, delta=0, pos2=False):
        self._ifeat = self._setup(mlib._modpt, protein, delta, pos2)
        self._create_bins(mlib, bins)


class Atom(_Base):
    def __init__(self, mlib, bins, pos2=False):
        self._ifeat = self._setup(mlib._modpt, pos2)
        self._create_bins(mlib, bins)


class AtomFixedBins(Atom):
    def __init__(self, mlib, pos2=False):
        self._ifeat = self._setup(mlib._modpt, pos2)


class XRayResolution(Protein):
    _setup = _mdt.mdt_feature_xray_resolution

class RadiusOfGyration(Protein):
    _setup = _mdt.mdt_feature_radius_of_gyration

class ResidueAccessibility(Residue):
    _setup = _mdt.mdt_feature_residue_accessibility

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

class HydrogenBondSatisfaction(Protein):
    _setup = _mdt.mdt_feature_hydrogen_bond_satisfaction
