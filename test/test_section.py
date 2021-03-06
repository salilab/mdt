import unittest
from mdt_test import MDTTest
import mdt
import mdt.features


class SectionTests(MDTTest):

    def test_section(self):
        """Test access to sections of MDTs"""
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        chi1class = mdt.features.Chi1Class(mlib)
        m = mdt.Table(mlib, features=(restyp, chi1, chi1class))
        m[0, 1, 2] = 1.0
        self.assertEqual(m[0][1][2], 1.0)
        m[1][-2][3] = 4.0
        self.assertEqual(m[1, -2, 3], 4.0)
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

    def test_section_bins(self):
        """Test checks of section bin indices"""
        import _mdt
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        m = mdt.Table(mlib, features=(restyp, chi1))
        self.assertRaises(ValueError, _mdt.mdt_section_sum, m._basept,
                          [0, 0, 0])
        self.assertRaises(IndexError, _mdt.mdt_section_sum, m._basept, [22])
        self.assertRaises(IndexError, _mdt.mdt_section_entropy, m._basept,
                          [22])
        self.assertRaises(IndexError, _mdt.mdt_section_meanstdev,
                          m._basept, mlib._modpt, [22])


if __name__ == '__main__':
    unittest.main()
