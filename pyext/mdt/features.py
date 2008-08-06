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


class XRayResolution(Protein):
    _setup = _mdt.mdt_feature_xray_resolution

class RadiusOfGyration(Protein):
    _setup = _mdt.mdt_feature_radius_of_gyration
