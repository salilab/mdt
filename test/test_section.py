import unittest
from mdt_test import MDTTest
import mdt

class SectionTests(MDTTest):

    def test_section(self):
        """Test access to sections of MDTs"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        m = mdt.Table(mlib, features=(1,chi1,18))
        m[0,1,2] = 1.0
        self.assertEqual(m[0][1][2], 1.0)
        m[1][-2][3] = 4.0
        self.assertEqual(m[1,-2,3], 4.0)
        sec = m[1]
        self.assertEqual([f.ifeat for f in sec.features],
                         [f.ifeat for f in m.features[1:]])
        self.assertEqual(sec.shape, m.shape[1:])
        self.assertEqual(sec.offset, m.offset[1:])
        sec = m[2][4]
        self.assertEqual([f.ifeat for f in sec.features],
                         [f.ifeat for f in m.features[2:]])
        self.assertEqual(sec.shape, m.shape[2:])
        self.assertEqual(sec.offset, m.offset[2:])

if __name__ == '__main__':
    unittest.main()
