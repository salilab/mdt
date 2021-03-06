import unittest
from mdt_test import MDTTest
import mdt
import mdt.features


class SuperSmoothTests(MDTTest):

    def _get_weights(self, weight, nbins, norm):
        """Calculate weights for combining apriori PDF with data"""
        w1 = 1. / (1. + norm / (weight * float(nbins)))
        w2 = 1. - w1
        return w1, w2

    def test_empty_2d(self):
        """Super-smoothed empty 2D MDT should be the even distribution"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = mdt.Table(mlib, features=(restyp0, restyp1))
        # A super-smoothed empty MDT should be the even distribution,
        # i.e. weight/nbin, for any value of the prior weight:
        for weight in (1.0, 2.0, 0.0):
            m2 = m.super_smooth(1, weight, True)
            inds = []
            while self.roll_inds(inds, m.shape, m.offset):
                self.assertAlmostEqual(1./22., m2[inds], places=3)

    def test_dimensions(self):
        """Check acceptable values for dimensions for super_smooth"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        m = mdt.Table(mlib, features=(restyp0, restyp1, chi1))
        for dimensions in (1, 2):
            _ = m.super_smooth(dimensions, 1.0, True)
        for dimensions in (-1, 0, 3, 4):
            self.assertRaises(ValueError, m.super_smooth, dimensions,
                              1.0, True)

    def test_3d(self):
        """Super-smoothed 3D MDT should not crash"""
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        chi1class = mdt.features.Chi1Class(mlib)
        m = self.get_test_mdt(mlib, features=(restyp, chi1, chi1class))
        m1 = m.super_smooth(1, 0.5, False)
        # Every 1D section should be normalized:
        for sec in m1:
            for sec2 in sec:
                self.assertSectionNormalized(sec2)
        m1 = m.super_smooth(2, 0.5, False)
        # Every 2D section should be normalized:
        for sec in m1:
            self.assertSectionNormalized(sec)

    def test_only_data(self):
        """Only data should be used if prior_weight=0"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = mdt.Table(mlib, features=(restyp0, restyp1))
        m[0, 0] = m[1, 1] = 1.0
        m2 = m.super_smooth(dimensions=1, prior_weight=0.0,
                            entropy_weighing=False)
        # Every section should be normalized:
        for sec in m2:
            self.assertSectionNormalized(sec)
        # Where we have data, m2==m:
        for i in (0, 1):
            for (a, b) in zip(m[i], m2[i]):
                self.assertEqual(a, b)
        # Where no data, m2==even distribution:
        for i in range(2, 22):
            for a in m2[i]:
                self.assertAlmostEqual(a, 1.0 / 22.0, delta=1e-5)

    def test_simple_smooth(self):
        """Test smoothing of a simple input distribution"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = mdt.Table(mlib, features=(restyp0, restyp1))
        m[0, 0] = m[1, 1] = m[1, 2] = 1.0
        m2 = m.super_smooth(dimensions=1, prior_weight=1.0,
                            entropy_weighing=False)
        # Get first level distribution from all data (1,1,1,0,0,0...):
        w1, w2 = self._get_weights(1.0, 22, 3.0)
        lev1 = [w1 / 22. + w2 / 3.] * 3 + [w1 / 22.] * 19
        # First row should be a combination of the apriori (level1)
        # distribution and the data (1,0,0,0...):
        w1, w2 = self._get_weights(1.0, 22, 1.0)
        self.assertAlmostEqual(m2[0, 0], w1 * lev1[0] + w2 * 1.0, delta=1e-5)
        for i in range(1, 22):
            self.assertAlmostEqual(m2[0, i], w1 * lev1[i], delta=1e-5)
        # Same deal for second row, using data (0,1,1,0,0,0...):
        w1, w2 = self._get_weights(1.0, 22, 2.0)
        for i in (1, 2):
            self.assertAlmostEqual(m2[1, i], w1 * lev1[i] + w2 * 0.5,
                                   delta=1e-5)
        for i in list(range(3, 22)) + [0]:
            self.assertAlmostEqual(m2[1, i], w1 * lev1[i], delta=1e-5)
        # Every other row is just the level1 data:
        for i in range(2, 22):
            for j in range(22):
                self.assertAlmostEqual(m2[i, j], lev1[j], delta=1e-5)


if __name__ == '__main__':
    unittest.main()
