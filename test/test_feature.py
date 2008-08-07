import unittest
from mdt_test import MDTTest
import mdt
import mdt.features
import modeller

class FeatureTests(MDTTest):

    def test_deltaij(self):
        """Test residue type at deltai,j features"""
        env = self.get_environ()
        aln = modeller.alignment(env)
        aln.append_sequence("AFVVTDNCIKXCKYTDCVEVCPVDCFYEG")
        aln.append_sequence("DNCIKXCCYCDCVEPCPVDCFGEGAFVVT")
        # When deltai=j=0, features 66,67 (residue type in A,B at deltai)
        # and 77,78 (residue type in A,B at deltaj) should match 1,2:
        mlib = self.get_mdt_library(deltai=0, deltaj=0)
        m1 = mdt.Table(mlib, features=(1,2))
        m1.add_alignment(aln)
        m2 = mdt.Table(mlib, features=(66,67))
        m2.add_alignment(aln)
        m3 = mdt.Table(mlib, features=(77,78))
        m3.add_alignment(aln)
        self.assertMDTDataEqual(m1, m2)
        self.assertMDTDataEqual(m1, m3)
        # When deltai=j != 0, 66,67 should still match 77,78, but not the
        # original MDT (m3)
        mlib = self.get_mdt_library(deltai=3, deltaj=3)
        m1 = mdt.Table(mlib, features=(66,67))
        m1.add_alignment(aln)
        m2 = mdt.Table(mlib, features=(77,78))
        m2.add_alignment(aln)
        self.assertMDTDataEqual(m1, m2)
        self.assertInTolerance(m2[0,2], 0.0, 0.0005)
        self.assertInTolerance(m3[0,2], 1.0, 0.0005)

    def test_feature_iatta(self):
        """Check for atom type features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.atom_classes.read('${LIB}/atmcls-melo.lib')
        attyp = mdt.features.AtomType(mlib)
        attyp2 = mdt.features.AtomType(mlib, pos2=True)
        aln = modeller.alignment(env, file='test/data/tiny.ali')
        m = mdt.Table(mlib, features=attyp)
        m.add_alignment(aln)
        self.assertInTolerance(m[0], 6.0, 0.0005)
        self.assertInTolerance(m[1], 0.0, 0.0005)
        self.assertInTolerance(m[2], 6.0, 0.0005)
        self.assertEqual(m.shape, (41,))
        # Using the pos2 feature should force a scan of all atom pairs:
        m = mdt.Table(mlib, features=attyp2)
        m.add_alignment(aln)
        self.assertEqual(m.shape, (41,))
        self.assertInTolerance(m[0], 74.0, 0.0005)
        self.assertInTolerance(m[1], 0.0, 0.0005)

    def test_feature_hbond(self):
        """Check hydrogen bond features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.hbond_classes.read('data/atmcls-hbda.lib')
        m = mdt.Table(mlib, features=84)
        m2 = mdt.Table(mlib, features=85)
        m3 = mdt.Table(mlib, features=86)
        aln = modeller.alignment(env, file='test/data/alignment.ali')
        m.add_alignment(aln)
        m2.add_alignment(aln)
        m3.add_alignment(aln)
        self.assertInTolerance(m[0], 295.0, 0.0005)
        self.assertInTolerance(m[1], 139.0, 0.0005)
        self.assertEqual(m[-1], 349.0)
        self.assertInTolerance(m2[0], 236.0, 0.0005)
        self.assertInTolerance(m2[1], 223.0, 0.0005)
        self.assertEqual(m2[-1], 168.0)
        self.assertInTolerance(m3[0], 1.0, 0.0005)
        self.assertInTolerance(m3[1], 0.0, 0.0005)
        self.assertEqual(m3[-1], 0.0)

    def test_feature_radius_gyration(self):
        """Check radius of gyration feature"""
        mlib = self.get_mdt_library()
        rgyr = mdt.features.RadiusOfGyration(mlib,
                                             mdt.uniform_bins(45, 5.0, 1.0))
        m = self.get_test_mdt(mlib, features=rgyr)
        self.assertEqual(m[7], 1.0)
        for (n, bin) in enumerate(m):
            if n != 7:
                self.assertEqual(bin, 0.0)

    def test_feature_iresol(self):
        """Check resolution features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(3, -1.0, 1.5)
        xray0 = mdt.features.XRayResolution(mlib, bins, protein=0)
        xray1 = mdt.features.XRayResolution(mlib, bins, protein=1)
        m = self.get_test_mdt(mlib, features=xray0)
        m2 = self.get_test_mdt(mlib, features=xray1)
        self.assertEqual(m.shape, (4,))
        self.assertEqual([b for b in m], [0., 1., 1., 0.])
        self.assertMDTDataEqual(m, m2)

        for (code, bin) in (('bin0', 0), ('bin1', 1), ('bin2', 2),
                            ('undef1', 3), ('undef2', 3)):
            m = mdt.Table(mlib, features=xray0)
            aln = modeller.alignment(env, file='test/data/resol.ali',
                                     align_codes=code)
            m.add_alignment(aln)
            self.assertEqual(m[bin], 1.0)

    def test_feature_atmacc(self):
        """Check atom accessibility features"""
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(120, 0.0, 0.20)
        atmacc = mdt.features.AtomAccessibility(mlib, bins)
        bins = mdt.uniform_bins(60, 0.0, 0.01)
        fatmacc = mdt.features.FractionalAtomAccessibility(mlib, bins)
        m = self.get_test_mdt(mlib, features=atmacc)
        self.assertEqual(m.shape, (121,))
        self.assertInTolerance(m[0], 425.0, 1.0005)
        self.assertInTolerance(m[1], 35.0, 2.0005)
        self.assertInTolerance(m[2], 17.0, 0.0005)
        self.assertEqual(m[-1], 0.0)
        m = self.get_test_mdt(mlib, features=fatmacc)
        self.assertEqual(m.shape, (61,))
        self.assertInTolerance(m[0], 457.0, 1.0005)
        self.assertInTolerance(m[1], 39.0, 1.0005)
        self.assertInTolerance(m[2], 35.0, 2.0005)
        self.assertEqual(m[-1], 0.0)

    def test_feature_resacc(self):
        """Check residue accessibility features"""
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(30, 0.0, 5.0)
        resacc = mdt.features.ResidueAccessibility(mlib, bins)
        m = self.get_test_mdt(mlib, features=resacc)
        self.assertEqual(m.shape, (31,))
        self.assertInTolerance(m[0], 24.0, 1.0005)
        self.assertInTolerance(m[1], 10.0, 2.0005)
        self.assertInTolerance(m[2], 4.0, 1.0005)
        self.assertEqual(m[-1], 0.0)

    def test_feature_bond_type(self):
        """Check bond type features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.bond_classes.read('data/bndgrp.lib')
        m = mdt.Table(mlib, features=109)
        m2 = mdt.Table(mlib, features=110)
        aln = modeller.alignment(env, file='test/data/alignment.ali')
        m.add_alignment(aln)
        m2.add_alignment(aln)
        self.assertInTolerance(m[0], 7.0, 0.0005)
        self.assertInTolerance(m[5], 9.0, 0.0005)
        self.assertInTolerance(m[19], 14.0, 0.0005)
        self.assertEqual(m[-1], 0.0)
        self.assertEqual(m.shape, (174,))
        self.assertInTolerance(m2[0], 0.0, 0.0005)
        self.assertInTolerance(m2[43], 3.0, 0.0005)
        self.assertInTolerance(m2[44], 10.0, 0.0005)
        self.assertEqual(m2[-1], 0.0)
        self.assertEqual(m2.shape, (201,))

    def test_feature_angle_type(self):
        """Check angle type features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.angle_classes.read('data/anggrp.lib')
        m = mdt.Table(mlib, features=111)
        m2 = mdt.Table(mlib, features=112)
        aln = modeller.alignment(env, file='test/data/alignment.ali')
        m.add_alignment(aln)
        m2.add_alignment(aln)
        self.assertInTolerance(m[0], 7.0, 0.0005)
        self.assertInTolerance(m[7], 9.0, 0.0005)
        self.assertInTolerance(m[15], 11.0, 0.0005)
        self.assertEqual(m.shape, (236,))
        self.assertEqual(m[-1], 0.0)
        self.assertInTolerance(m2[176], 48.0, 1.0005)
        self.assertInTolerance(m2[177], 42.0, 0.0005)
        self.assertInTolerance(m2[178], 38.0, 0.0005)
        self.assertEqual(m2.shape, (289,))
        self.assertEqual(m2[-1], 0.0)

    def test_feature_dihedral_type(self):
        """Check dihedral type features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.dihedral_classes.read('data/impgrp.lib')
        m = mdt.Table(mlib, features=113)
        m2 = mdt.Table(mlib, features=114)
        aln = modeller.alignment(env, file='test/data/alignment.ali')
        m.add_alignment(aln)
        m2.add_alignment(aln)
        self.assertInTolerance(m[0], 7.0, 0.0005)
        self.assertInTolerance(m[2], 9.0, 0.0005)
        self.assertInTolerance(m[4], 11.0, 0.0005)
        self.assertEqual(m.shape, (79,))
        self.assertEqual(m[-1], 0.0)
        self.assertInTolerance(m2[143], 60.0, 1.0005)
        self.assertInTolerance(m2[144], 53.0, 1.0005)
        self.assertInTolerance(m2[145], 24.0, 0.0005)
        self.assertEqual(m2.shape, (289,))
        self.assertEqual(m2[-1], 0.0)

    def test_feature_doublet_type(self):
        """Check doublet type features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.tuple_classes.read('data/dblcls.lib')
        m1 = mdt.Table(mlib, features=105)
        m2 = mdt.Table(mlib, features=106)
        aln = modeller.alignment(env, file='test/data/tiny.ali')
        for m in (m1, m2):
            m.add_alignment(aln, residue_span_range=(-9999, 0, 0, 9999))
        for f in (107, 108):
            m = mdt.Table(mlib, features=f)
            self.assertRaises(mdt.MDTError, m.add_alignment, aln)
        self.assertEqual(m1.shape, (7,))
        self.assertEqual(m2.shape, (7,))
        self.assertInTolerance(m1[0], 311.0, 0.0005)
        self.assertInTolerance(m2[0], 302.0, 0.0005)

    def test_feature_triplet_type(self):
        """Check triplet type features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.tuple_classes.read('data/trpcls.lib')
        m1 = mdt.Table(mlib, features=101)
        m2 = mdt.Table(mlib, features=102)
        m3 = mdt.Table(mlib, features=103)
        m4 = mdt.Table(mlib, features=104)
        m5 = mdt.Table(mlib, features=106)
        aln = modeller.alignment(env, file='test/data/tiny.ali')
        for m in (m1, m2, m3, m4, m5):
            m.add_alignment(aln, residue_span_range=(-9999, 0, 0, 9999))
        self.assertInTolerance(m1[0], 1.0, 0.0005)
        self.assertInTolerance(m1[1], 0.0, 0.0005)
        self.assertInTolerance(m1[2], 1.0, 0.0005)
        self.assertEqual(m1.shape, (236,))
        self.assertEqual(m1[-1], 0.0)
        self.assertInTolerance(m2[0], 60.0, 0.0005)
        self.assertInTolerance(m2[1], 0.0, 0.0005)
        self.assertInTolerance(m2[2], 60.0, 0.0005)
        self.assertEqual(m2.shape, (236,))
        self.assertEqual(m2[-1], 0.0)
        self.assertInTolerance(m3[0], 0.0, 0.0005)
        self.assertInTolerance(m3[1], 82.0, 0.0005)
        self.assertInTolerance(m3[2], 226.0, 0.0005)
        self.assertEqual(m3.shape, (10,))
        self.assertInTolerance(m3[-1], 3018.0, 0.0005)
        self.assertInTolerance(m4[0], 479.0, 0.0005)
        self.assertInTolerance(m4[1], 806.0, 0.0005)
        self.assertInTolerance(m4[2], 471.0, 0.0005)
        self.assertEqual(m4.shape, (7,))
        self.assertEqual(m4[-1], 0.0)
        self.assertInTolerance(m5[0], 556.0, 0.0005)
        self.assertInTolerance(m5[1], 642.0, 0.0005)
        self.assertInTolerance(m5[2], 470.0, 6.0005)
        self.assertEqual(m5.shape, (7,))
        self.assertInTolerance(m5[-1], 180.0, 0.0005)

    def test_feature_alpha_dihedral(self):
        """Check alpha dihedral feature"""
        mlib = self.get_mdt_library()
        m = self.get_test_mdt(mlib, features=88)
        self.assertEqual(m.shape, (37,))
        self.assertInTolerance(m[30], 3.0, 0.0005)
        self.assertInTolerance(m[31], 4.0, 0.0005)
        self.assertInTolerance(m[32], 2.0, 0.0005)

    def test_feature_distance(self):
        """Check atom-atom distance feature"""
        mlib = self.get_mdt_library()
        m = self.get_test_mdt(mlib, features=82)
        self.assertEqual(m.shape, (61,))
        self.assertInTolerance(m[30], 10057.0, 0.0005)
        self.assertInTolerance(m[31], 10214.0, 0.0005)
        self.assertInTolerance(m[32], 10095.0, 0.0005)
        self.assertInTolerance(m[-1], 4892.0, 0.0005)

    def test_symmetric(self):
        """Test symmetric/asymmetric residue pair features"""
        sym_features = (24, 25, 39, 48)
        asym_features = (16, 17, 40, 49, 51, 52)
        env = self.get_environ()
        mlib = self.get_mdt_library()
        for a in sym_features:
            m = mdt.Table(mlib, features=a)
            self.assertEqual(m.symmetric, True)
        for a in asym_features:
            m = mdt.Table(mlib, features=a)
            self.assertEqual(m.symmetric, False)
        # Combinations are only symmetric if all features are symmetric:
        for (a, b, symm) in ((sym_features[0], sym_features[1], True),
                             (asym_features[0], asym_features[1], False),
                             (sym_features[0], asym_features[0], False)):
            m = mdt.Table(mlib, features=(a,b))
            self.assertEqual(m.symmetric, symm)

if __name__ == '__main__':
    unittest.main()
