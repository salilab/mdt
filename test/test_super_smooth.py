import unittest
from mdt_test import MDTTest
import mdt

class SuperSmoothTests(MDTTest):

    def test_empty_2d(self):
        """Super-smoothed empty 2D MDT should be the even distribution"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m = mdt.Table(mlib, features=(1,2))
        # A super-smoothed empty MDT should be the even distribution,
        # i.e. weight/nbin, for any value of the prior weight:
        for weight in (1.0, 2.0, 0.0):
            m2 = m.super_smooth(1, weight, True)
            inds = []
            while self.roll_inds(inds, m.shape, m.offset):
                self.assertAlmostEqual(1./22., m2[inds], places=3)

    def test_dimensions(self):
        """Check acceptable values for dimensions for super_smooth"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m = mdt.Table(mlib, features=(1,2,3))
        for dimensions in (1, 2):
            m2 = m.super_smooth(dimensions, 1.0, True)
        for dimensions in (-1, 0, 3, 4):
            self.assertRaises(ValueError, m.super_smooth, dimensions, 1.0, True)

    def test_3d(self):
        """Super-smoothed 3D MDT should not crash"""
        m = self.get_test_mdt(features=(1,3,18))
        m1 = m.super_smooth(1, 0.5, False)
        m1 = m.super_smooth(2, 0.5, False)

if __name__ == '__main__':
    unittest.main()
