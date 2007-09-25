import unittest
from mdt_test import MDTTest
import mdt

class SuperSmoothTests(MDTTest):

    def test_empty_2d(self):
        """Super-smoothed empty 2D MDT should be the prior"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m = mdt.Table(mlib, features=(1,2))
        # A super-smoothed empty MDT should be the prior, i.e. 1/nbin:
        m2 = m.super_smooth(1.0, True)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            self.assertAlmostEqual(1./22., m2[inds], places=3)

    def test_3d(self):
        """Super-smoothed 3D MDT should not crash"""
        m = self.get_test_mdt(features=(1,3,18))
        m1 = m.super_smooth(0.5, False)

if __name__ == '__main__':
    unittest.main()
