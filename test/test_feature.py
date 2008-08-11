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
        donor = mdt.features.HydrogenBondDonor(mlib,
                                               mdt.uniform_bins(7, 1., 1.))
        accep = mdt.features.HydrogenBondAcceptor(mlib,
                                                  mdt.uniform_bins(7, 1., 1.))
        totchg = mdt.features.HydrogenBondCharge(mlib,
                                                 mdt.uniform_bins(9, 1., 1.))
        satisf = mdt.features.HydrogenBondSatisfaction(mlib,
                                                 mdt.uniform_bins(100, 0., 10.))
        m = mdt.Table(mlib, features=donor)
        m2 = mdt.Table(mlib, features=accep)
        m3 = mdt.Table(mlib, features=satisf)
        m4 = mdt.Table(mlib, features=totchg)
        aln = modeller.alignment(env, file='test/data/alignment.ali')
        m.add_alignment(aln)
        m2.add_alignment(aln)
        m3.add_alignment(aln)
        m4.add_alignment(aln)
        self.assertInTolerance(m[0], 295.0, 0.0005)
        self.assertInTolerance(m[1], 139.0, 0.0005)
        self.assertEqual(m[-1], 349.0)
        self.assertInTolerance(m2[0], 236.0, 0.0005)
        self.assertInTolerance(m2[1], 223.0, 0.0005)
        self.assertEqual(m2[-1], 168.0)
        self.assertInTolerance(m3[0], 1.0, 0.0005)
        self.assertInTolerance(m3[1], 0.0, 0.0005)
        self.assertEqual(m3[-1], 0.0)
        self.assertInTolerance(m4[0], 78.0, 0.0005)
        self.assertInTolerance(m4[1], 24.0, 0.0005)
        self.assertEqual(m4[-1], 739.0)

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

    def test_feature_seqlen(self):
        """Check sequence length feature"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        seqlen = mdt.features.SequenceLength(mlib,
                                             bins=mdt.uniform_bins(49, 0, 5))
        m = self.get_test_mdt(mlib, features=seqlen)
        self.assertEqual(m.shape, (50,))
        self.assertEqual(m[10], 1.0)
        self.assertEqual(m[21], 1.0)
        self.assertEqual(sum([b for b in m]), 2.0)

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
        bondtype = mdt.features.BondType(mlib)
        bondlen = mdt.features.BondLength(mlib,
                                      bins=mdt.uniform_bins(200, 1.0, 0.005))
        m = mdt.Table(mlib, features=bondtype)
        m2 = mdt.Table(mlib, features=bondlen)
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
        angletype = mdt.features.AngleType(mlib)
        angle = mdt.features.Angle(mlib,
                                   bins=mdt.uniform_bins(288, 0.0, 0.625))
        m = mdt.Table(mlib, features=angletype)
        m2 = mdt.Table(mlib, features=angle)
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
        dihedtype = mdt.features.DihedralType(mlib)
        dihedral = mdt.features.Dihedral(mlib,
                                         bins=mdt.uniform_bins(288, -180, 1.25))
        m = mdt.Table(mlib, features=dihedtype)
        m2 = mdt.Table(mlib, features=dihedral)
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
        tuple_angle2 = mdt.features.TupleAngle2(mlib,
                                            bins=mdt.uniform_bins(6, 0, 30.0))
        tuple_dihed1 = mdt.features.TupleDihedral1(mlib,
                                          bins=mdt.uniform_bins(6, -180, 60.0))
        # These features only work on atom triplets:
        for f in (mdt.features.TupleDihedral2, mdt.features.TupleDihedral3):
            self.assertRaises(mdt.MDTError, f, mlib,
                              bins=mdt.uniform_bins(6, -180, 60.0))
        m1 = mdt.Table(mlib, features=tuple_angle2)
        m2 = mdt.Table(mlib, features=tuple_dihed1)
        aln = modeller.alignment(env, file='test/data/tiny.ali')
        for m in (m1, m2):
            m.add_alignment(aln, residue_span_range=(-9999, 0, 0, 9999))
        self.assertEqual(m1.shape, (7,))
        self.assertEqual(m2.shape, (7,))
        self.assertInTolerance(m1[0], 311.0, 0.0005)
        self.assertInTolerance(m2[0], 302.0, 0.0005)

    def test_feature_triplet_type(self):
        """Check triplet type features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.tuple_classes.read('data/trpcls.lib')
        tuple_type = mdt.features.TupleType(mlib)
        tuple_type2 = mdt.features.TupleType(mlib, pos2=True)
        tuple_dist = mdt.features.TupleDistance(mlib,
                                            bins=mdt.uniform_bins(9, 2.0, 0.2))
        tuple_angle1 = mdt.features.TupleAngle1(mlib,
                                            bins=mdt.uniform_bins(6, 0, 30.0))
        tuple_dihed1 = mdt.features.TupleDihedral1(mlib,
                                          bins=mdt.uniform_bins(6, -180, 60.0))
        m1 = mdt.Table(mlib, features=tuple_type)
        m2 = mdt.Table(mlib, features=tuple_type2)
        m3 = mdt.Table(mlib, features=tuple_dist)
        m4 = mdt.Table(mlib, features=tuple_angle1)
        m5 = mdt.Table(mlib, features=tuple_dihed1)
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

    def test_feature_chi1_dihedral(self):
        """Check chi1 dihedral and dihedral class features"""
        mlib = self.get_mdt_library()
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        m = self.get_test_mdt(mlib, features=chi1)
        self.assertEqual(m.shape, (37,))
        self.assertInTolerance(m[0], 8.0, 0.0005)
        self.assertInTolerance(m[1], 7.0, 0.0005)
        self.assertInTolerance(m[2], 2.0, 0.0005)
        m = self.get_test_mdt(mlib, features=18)
        self.assertEqual([b for b in m], [50, 31, 15, 10])

    def test_feature_chi2_dihedral(self):
        """Check chi2 dihedral and dihedral class features"""
        mlib = self.get_mdt_library()
        chi2 = mdt.features.Chi2Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        m = self.get_test_mdt(mlib, features=chi2)
        self.assertEqual(m.shape, (37,))
        self.assertInTolerance(m[0], 6.0, 0.0005)
        self.assertInTolerance(m[1], 4.0, 0.0005)
        self.assertInTolerance(m[2], 3.0, 0.0005)
        m = self.get_test_mdt(mlib, features=20)
        self.assertEqual([b for b in m], [47, 23, 4, 32])

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
        dist = mdt.features.AtomDistance(mlib,
                                         bins=mdt.uniform_bins(60, 0, 0.5))
        m = self.get_test_mdt(mlib, features=dist)
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
