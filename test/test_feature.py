import unittest
from mdt_test import MDTTest
import mdt
import mdt.features
import modeller
import os
import math


class FeatureTests(MDTTest):

    def build_mdt_from_model(self, mlib, features, mdl, **keys):
        """Build a simple test MDT for a given model"""
        env = self.get_environ()
        m = mdt.Table(mlib, features=features)
        a = modeller.Alignment(env)
        a.append_model(mdl, atom_files='test', align_codes='test')
        m.add_alignment(a, **keys)
        return m

    def build_test_model(self):
        """Build a simple test model"""
        env = self.get_environ()
        m = modeller.Model(env)
        m.build_sequence('C')
        return m

    def test_feature_atom_table(self):
        """Test AtomTable feature"""
        class UniqueError(Exception):
            pass

        def build_mdt(mlib, func, aln):
            bins = mdt.uniform_bins(5, 0, 1.0)
            f = mdt.features.AtomTable(mlib, bins, "test data", func)
            m = mdt.Table(mlib, features=f)
            m.add_alignment(aln)
            return m

        def ex_func(*args):
            raise UniqueError()

        def bad_type_func(*args):
            return 42

        def bad_val_func(*args):
            return ['foo']*7

        def bad_len_func(*args):
            return [42]*3

        def func(aln, struc, mlib, libs):
            return [a.x for a in struc.atoms]
        mlib = self.get_mdt_library()
        env = self.get_environ()
        mdl = self.build_test_model()
        aln = modeller.Alignment(env)
        aln.append_model(mdl, atom_files='test', align_codes='test')

        # Python exceptions raised in the function should be propagated
        self.assertRaises(UniqueError, build_mdt, mlib, ex_func, aln)
        # Catch return values that aren't sequences
        self.assertRaises(TypeError, build_mdt, mlib, bad_type_func, aln)
        # Catch return values that aren't sequences of floats
        self.assertRaises(ValueError, build_mdt, mlib, bad_val_func, aln)
        # Catch return values that are the wrong length
        self.assertRaises(ValueError, build_mdt, mlib, bad_len_func, aln)
        t = build_mdt(mlib, func, aln)
        self.assertEqual([b for b in t], [2.0, 2.0, 1.0, 2.0, 0.0, 0.0])

    def test_feature_gap_distance(self):
        """Check distance from a gap features"""
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(9, 0, 1)
        gapdist = mdt.features.GapDistance(mlib, bins)
        avgapdist = mdt.features.AverageGapDistance(mlib, bins)
        m = self.get_test_mdt(mlib, features=gapdist)
        self.assertEqual(m.shape, (10,))
        self.assertEqual(m.sum(), 212)
        self.assertEqual(m[0], 104)
        self.assertEqual(m[1], 12)
        self.assertEqual(m[-1], 28)
        m = self.get_test_mdt(mlib, features=avgapdist)
        self.assertEqual(m.shape, (10,))
        self.assertEqual(m.sum(), 10920)
        self.assertEqual(m[0], 3168)
        self.assertEqual(m[1], 1338)
        self.assertEqual(m[-1], 406)

    def test_feature_alpha_content(self):
        """Check alpha content feature"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        self.assertRaises(ValueError, mdt.features.AlphaContent, mlib,
                          bins=mdt.uniform_bins(10, 0, 0.1), protein=3)
        alpha = mdt.features.AlphaContent(mlib,
                                          bins=mdt.uniform_bins(10, 0, 0.1))
        for (alnfile, bin) in (('tiny.ali', 0), ('alignment.ali', 5)):
            m = mdt.Table(mlib, features=alpha)
            a = modeller.Alignment(env,
                                   file=os.path.join('test', 'data', alnfile))
            m.add_alignment(a)
            self.assertEqual(m.shape, (11,))
            self.assertEqual(m.sum(), 1)
            self.assertEqual(m[bin], 1)

    def test_feature_sidechain_biso(self):
        """Check average sidechain Biso feature"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        self.assertRaises(ValueError, mdt.features.SidechainBiso, mlib,
                          bins=mdt.uniform_bins(5, 0, 10), protein=3)
        sidechain_biso = mdt.features.SidechainBiso(
            mlib, bins=mdt.uniform_bins(5, 0, 10))
        mdl = modeller.Model(env)
        mdl.build_sequence('A')
        aln = modeller.Alignment(env)
        aln.append_model(mdl, align_codes='test')
        chain = aln[0].chains[0]
        # Mainchain atom Biso should be ignored:
        for mainchain in ('N:1', 'C:1', 'O:1', 'OXT:1', 'CA:1'):
            chain.atoms[mainchain].biso = 1000
        for (biso, bin) in ((22, 2), (32, 3),  # Map regular values to bins
                            (0, -1),  # Zero Biso should be "undefined"
                            (1, 3)):  # Biso < 2 is multiplied by 4pi^2
            chain.atoms['CB:1'].biso = biso
            m = mdt.Table(mlib, features=sidechain_biso)
            m.add_alignment(aln)
            self.assertEqual(m.shape, (6,))
            self.assertEqual(m.sum(), 1)
            self.assertEqual(m[bin], 1)

    def test_feature_sequence_identity(self):
        """Check sequence identity feature"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        # Put into 25% bins
        sid = mdt.features.SequenceIdentity(mlib,
                                            bins=mdt.uniform_bins(5, 0, 0.250))
        for (seq, id) in (('GGG', 0), ('AFV', 100), ('A--', 100),
                          ('AV-', 50)):
            aln = modeller.Alignment(env)
            aln.append_sequence('AFV')
            aln.append_sequence(seq)
            m = mdt.Table(mlib, features=sid)
            m.add_alignment(aln)
            self.assertEqual(m.shape, (6,))
            self.assertEqual(m.sum(), 2.0)
            self.assertEqual(m[int(id / 25)], 2.0)

    def test_delta(self):
        """Test residue type at delta features"""
        env = self.get_environ()
        aln = modeller.Alignment(env)
        aln.append_sequence("AFVVTDNCIKXCKYTDCVEVCPVDCFYEG")
        aln.append_sequence("DNCIKXCCYCDCVEPCPVDCFGEGAFVVT")
        mlib = self.get_mdt_library()
        self.assertRaises(ValueError, mdt.features.ResidueType, mlib,
                          protein=3)
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        restyp0_del3 = mdt.features.ResidueType(mlib, protein=0, delta=3)
        restyp1_del3 = mdt.features.ResidueType(mlib, protein=1, delta=3)

        m1 = mdt.Table(mlib, features=(restyp0, restyp1))
        m1.add_alignment(aln)

        # When deltai=j != 0, offset MDTs should not match the original:
        m2 = mdt.Table(mlib, features=(restyp0_del3, restyp1_del3))
        m2.add_alignment(aln)

        self.assertAlmostEqual(m2[0, 2], 0.0, delta=0.0005)

    def test_feature_residue_distance(self):
        """Check residue-residue distance feature"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        dist = mdt.features.ResidueDistance(mlib,
                                            bins=mdt.uniform_bins(7, 0, 2.0))
        aln = modeller.Alignment(env, file='test/data/tiny.ali')
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln)
        self.assertEqual([b for b in m], [0, 0, 0, 8, 2, 4, 4, 2])

    def test_feature_residue_distance_difference(self):
        """Check residue-residue distance difference feature"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        ddist = mdt.features.ResidueDistanceDifference(
            mlib, bins=mdt.uniform_bins(20, -10, 1))
        aln = modeller.Alignment(env, file='test/data/struc-struc.ali')
        m = mdt.Table(mlib, features=ddist)
        m.add_alignment(aln)
        self.assertEqual(m[9], 20)
        self.assertEqual(m[10], 20)
        self.assertEqual(sum([b for b in m]), 40)
        self.assertEqual(m[-1], 0)

        # Undefined (-999) coordinates in either structure should put
        # features in the undefined bin
        oldx = aln[0].residues[0].atoms['CA'].x
        aln[0].residues[0].atoms['CA'].x = -999
        m = mdt.Table(mlib, features=ddist)
        m.add_alignment(aln)
        self.assertEqual(m[-1], 16)

        aln[0].residues[0].atoms['CA'].x = oldx
        aln[1].residues[0].atoms['CA'].x = -999
        m = mdt.Table(mlib, features=ddist)
        m.add_alignment(aln)
        self.assertEqual(m[-1], 16)

    def test_feature_neighborhood_difference(self):
        """Check residue neighborhood difference features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(9, 0, 0.2)
        ndif = mdt.features.NeighborhoodDifference(mlib, bins)
        avndif = mdt.features.AverageNeighborhoodDifference(mlib, bins)
        aln = modeller.Alignment(env, file='test/data/struc-struc.ali')
        m = mdt.Table(mlib, features=ndif)
        m.add_alignment(aln)
        self.assertEqual([b for b in m], [4, 6, 2] + [0]*7)
        m = mdt.Table(mlib, features=avndif)
        m.add_alignment(aln)
        self.assertEqual([b for b in m], [6, 12, 2] + [0]*7)

    def test_feature_rama(self):
        """Check Ramachandran mainchain conformation class feature"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        self.assertRaises(ValueError, mdt.features.MainchainConformation,
                          mlib, protein=3)
        conf = mdt.features.MainchainConformation(mlib)
        aln = modeller.Alignment(env, file='test/data/tiny.ali')
        m = mdt.Table(mlib, features=conf)
        m.add_alignment(aln)
        self.assertEqual([b.symbol for b in m.features[0].bins],
                         ['A', 'P', 'B', 'L', 'E', 'U'])
        self.assertEqual([b for b in m], [0, 2, 2, 0, 0, 2])

    def test_feature_resgrp(self):
        """Check residue group feature"""
        mlib = self.get_mdt_library()
        mnch = mdt.features.ResidueGroup(mlib, residue_grouping=0)
        hydro = mdt.features.ResidueGroup(mlib, residue_grouping=1)
        for out_of_range in (-1, 2):
            self.assertRaises(ValueError, mdt.features.ResidueGroup,
                              mlib, residue_grouping=out_of_range)
        self.assertRaises(ValueError, mdt.features.ResidueGroup,
                          mlib, protein=3)
        m = self.get_test_mdt(mlib, features=mnch)
        self.assertEqual([b for b in m], [139, 7, 14, 0])
        m = self.get_test_mdt(mlib, features=hydro)
        self.assertEqual([b for b in m], [97, 47, 16])

    def test_feature_iatta_special(self):
        """Check atom type feature with disulfide/termini special handling"""
        env = self.get_environ()
        mlib = mdt.Library(env, special_atoms=True)
        mlib.atom_classes.read('${LIB}/atmcls-melo.lib')
        attyp = mdt.features.AtomType(mlib)
        aln = modeller.Alignment(env, file='test/data/tiny.ali')
        m = mdt.Table(mlib, features=attyp)
        m.add_alignment(aln)
        self.assertAlmostEqual(m[0], 6.0, delta=0.0005)
        self.assertAlmostEqual(m[1], 0.0, delta=0.0005)
        self.assertAlmostEqual(m[2], 5.0, delta=0.0005)

    def test_feature_iatta(self):
        """Check for atom type features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.atom_classes.read('${LIB}/atmcls-melo.lib')
        attyp = mdt.features.AtomType(mlib)
        attyp2 = mdt.features.AtomType(mlib, pos2=True)
        self.assertRaises(mdt.MDTError, mlib.atom_classes.read,
                          '${LIB}/atmcls-melo.lib')
        aln = modeller.Alignment(env, file='test/data/tiny.ali')
        m = mdt.Table(mlib, features=attyp)
        m.add_alignment(aln)
        self.assertAlmostEqual(m[0], 6.0, delta=0.0005)
        self.assertAlmostEqual(m[1], 0.0, delta=0.0005)
        self.assertAlmostEqual(m[2], 6.0, delta=0.0005)
        self.assertEqual(m.shape, (41,))
        # Using the pos2 feature should force a scan of all atom pairs:
        m = mdt.Table(mlib, features=attyp2)
        m.add_alignment(aln)
        self.assertEqual(m.shape, (41,))
        self.assertAlmostEqual(m[0], 74.0, delta=0.0005)
        self.assertAlmostEqual(m[1], 0.0, delta=0.0005)

    def test_bad_hbond_class_file(self):
        """Check reading of bad hbond class file"""
        mlib = self.get_mdt_library()
        with open('badhbond', 'w') as fh:
            fh.write("ATMGRP 'N'    not-a-float")
        self.assertRaises(mdt.MDTError, mlib.hbond_classes.read, 'badhbond')
        os.unlink('badhbond')

    def test_bad_atom_class_file(self):
        """Check reading of bad atom class file"""
        mlib = self.get_mdt_library()

        self.assertRaises(IOError, mlib.atom_classes.read,
                          '/does/not/exist/test.lib')

        # Bad fields after ATMGRP
        with open('badatmcls', 'w') as fh:
            fh.write("ATMGRP not-a-quoted-string")
        self.assertRaises(mdt.MDTError, mlib.atom_classes.read, 'badatmcls')

        # ATOM outside of ATMGRP
        with open('badatmcls', 'w') as fh:
            fh.write("ATOM 'foo' 'bar'")
        self.assertRaises(mdt.MDTError, mlib.atom_classes.read, 'badatmcls')

        # Bad fields after ATOM
        with open('badatmcls', 'w') as fh:
            fh.write("ATMGRP 'foo'\nATOM 'ok' not-a-string")
        self.assertRaises(mdt.MDTError, mlib.atom_classes.read, 'badatmcls')

        # Line starting with something not ATOM or ATMGRP
        with open('badatmcls', 'w') as fh:
            fh.write("DBLGRP 'foo'")
        self.assertRaises(mdt.MDTError, mlib.atom_classes.read, 'badatmcls')
        os.unlink('badatmcls')

    def test_feature_hbond_undef(self):
        """Check hydrogen bond features undefined bin"""
        mdl = self.build_test_model()
        modeller.Selection(mdl).unbuild()
        mlib = self.get_mdt_library()
        mlib.hbond_classes.read('data/atmcls-hbda.lib')
        bins = mdt.uniform_bins(1, -20000, 40000)

        # All features should go in the undefined bin, even though
        # the raw value of the feature should fit in the first bin
        for f in [mdt.features.HydrogenBondDonor(mlib, bins),
                  mdt.features.HydrogenBondAcceptor(mlib, bins),
                  mdt.features.HydrogenBondCharge(mlib, bins)]:
            m = self.build_mdt_from_model(mlib, f, mdl)
            self.assertEqual(m[0], 0)
            self.assertEqual(m[-1], 7)

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
        satisf = mdt.features.HydrogenBondSatisfaction(
            mlib, mdt.uniform_bins(100, 0., 10.))
        self.assertRaises(mdt.MDTError, mlib.hbond_classes.read,
                          'data/atmcls-hbda.lib')
        m = mdt.Table(mlib, features=donor)
        m2 = mdt.Table(mlib, features=accep)
        m3 = mdt.Table(mlib, features=satisf)
        m4 = mdt.Table(mlib, features=totchg)
        aln = modeller.Alignment(env, file='test/data/alignment.ali')
        m.add_alignment(aln)
        m2.add_alignment(aln)
        m3.add_alignment(aln)
        m4.add_alignment(aln)
        self.assertAlmostEqual(m[0], 295.0, delta=0.0005)
        self.assertAlmostEqual(m[1], 139.0, delta=0.0005)
        self.assertEqual(m[-1], 349.0)
        self.assertAlmostEqual(m2[0], 236.0, delta=0.0005)
        self.assertAlmostEqual(m2[1], 223.0, delta=0.0005)
        self.assertEqual(m2[-1], 168.0)
        self.assertAlmostEqual(m3[0], 1.0, delta=0.0005)
        self.assertAlmostEqual(m3[1], 0.0, delta=0.0005)
        self.assertEqual(m3[-1], 0.0)
        self.assertAlmostEqual(m4[0], 78.0, delta=0.0005)
        self.assertAlmostEqual(m4[1], 24.0, delta=0.0005)
        self.assertEqual(m4[-1], 739.0)
        # Exercise writing of hbond information to HDF5 files:
        for t in (m, m2, m3, m4):
            t.write_hdf5('test.hdf5')
            os.unlink('test.hdf5')

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
        xray0_nmr = mdt.features.XRayResolution(mlib, bins, protein=0, nmr=1.0)
        xray1 = mdt.features.XRayResolution(mlib, bins, protein=1)
        _ = mdt.features.XRayResolution(mlib, bins, protein=2)
        # Check valid range for protein argument
        for p in (-1, 3):
            self.assertRaises(ValueError, mdt.features.XRayResolution,
                              mlib, bins, protein=p)
        m = self.get_test_mdt(mlib, features=xray0)
        m2 = self.get_test_mdt(mlib, features=xray1)
        self.assertEqual(m.shape, (4,))
        self.assertEqual([b for b in m], [0., 1., 1., 0.])
        self.assertMDTDataEqual(m, m2)

        for (code, feat, bin) in (('bin0', xray0, 0), ('bin0', xray0_nmr, 1),
                                  ('bin1', xray0, 1), ('bin2', xray0, 2),
                                  ('undef1', xray0, 3), ('undef2', xray0, 3)):
            m = mdt.Table(mlib, features=feat)
            aln = modeller.Alignment(env, file='test/data/resol.ali',
                                     align_codes=code)
            m.add_alignment(aln)
            self.assertEqual(m[bin], 1.0)

    def test_feature_atmacc_undef(self):
        """Check atom accessibility features undefined bin"""
        mdl = self.build_test_model()
        modeller.Selection(mdl).unbuild()
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(1, -20000, 40000)

        # All features should go in the undefined bin, even though
        # the raw value of the feature should fit in the first bin
        for f in [mdt.features.AtomAccessibility(mlib, bins),
                  mdt.features.FractionalAtomAccessibility(mlib, bins)]:
            m = self.build_mdt_from_model(mlib, f, mdl)
            self.assertEqual(m[0], 0)
            self.assertEqual(m[-1], 7)

    def test_feature_atmacc(self):
        """Check atom accessibility features"""
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(120, 0.0, 0.20)
        atmacc = mdt.features.AtomAccessibility(mlib, bins)
        bins = mdt.uniform_bins(60, 0.0, 0.01)
        fatmacc = mdt.features.FractionalAtomAccessibility(mlib, bins)
        m = self.get_test_mdt(mlib, features=atmacc)
        self.assertEqual(m.shape, (121,))
        self.assertAlmostEqual(m[0], 425.0, delta=1.0005)
        self.assertAlmostEqual(m[1], 35.0, delta=2.0005)
        self.assertAlmostEqual(m[2], 17.0, delta=0.0005)
        self.assertEqual(m[-1], 0.0)
        m = self.get_test_mdt(mlib, features=fatmacc)
        self.assertEqual(m.shape, (61,))
        self.assertAlmostEqual(m[0], 457.0, delta=1.0005)
        self.assertAlmostEqual(m[1], 39.0, delta=1.0005)
        self.assertAlmostEqual(m[2], 35.0, delta=2.0005)
        self.assertEqual(m[-1], 0.0)

    def test_feature_z_coordinate(self):
        """Check atom Z-coordinate feature"""
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(40, 0.0, 1.0)
        z = mdt.features.AtomZCoordinate(mlib, bins)
        m = self.get_test_mdt(mlib, features=z)
        self.assertEqual(m.shape, (41,))
        self.assertEqual(m[0], 0.0)
        self.assertEqual(m[1], 1.0)
        self.assertEqual(m[2], 3.0)
        self.assertEqual(m[3], 6.0)
        self.assertEqual(m[4], 10.0)
        self.assertEqual(m[5], 31.0)
        self.assertEqual(m[40], 0.0)

    def test_feature_z_coordinate_undefined(self):
        """Check atom Z-coordinate feature undefined bin"""
        mdl = self.build_test_model()
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(1, -20000, 40000)
        z = mdt.features.AtomZCoordinate(mlib, bins)
        m = self.build_mdt_from_model(mlib, z, mdl)
        self.assertEqual(m.shape, (2,))
        self.assertEqual(m[0], 7.0)
        self.assertEqual(m[1], 0.0)
        # Make one z-coordinate undefined
        modeller.Selection(mdl.atoms[0]).unbuild()
        m = self.build_mdt_from_model(mlib, z, mdl)
        self.assertEqual(m.shape, (2,))
        self.assertEqual(m[0], 6.0)
        # undefined bin should contain one count even though actual value
        # should fall within first bin
        self.assertEqual(m[1], 1.0)

    def test_feature_resacc(self):
        """Check residue accessibility features"""
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(30, 0.0, 5.0)
        resacc = mdt.features.ResidueAccessibility(mlib, bins)
        avresacc = mdt.features.AverageResidueAccessibility(
            mlib, bins=mdt.uniform_bins(10, 0, 15))
        m = self.get_test_mdt(mlib, features=resacc)
        self.assertEqual(m.shape, (31,))
        self.assertAlmostEqual(m[0], 24.0, delta=1.0005)
        self.assertAlmostEqual(m[1], 10.0, delta=2.0005)
        self.assertAlmostEqual(m[2], 4.0, delta=1.0005)
        self.assertEqual(m[-1], 0.0)
        m = self.get_test_mdt(mlib, features=avresacc)
        self.assertEqual(m.shape, (11,))
        self.assertAlmostEqual(m[0], 988.0, delta=10.0005)
        self.assertAlmostEqual(m[1], 1379.0, delta=10.0005)
        self.assertAlmostEqual(m[2], 1317.0, delta=10.0005)

    def test_feature_resind_diff(self):
        """Test the residue index difference feature"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        diff = mdt.features.ResidueIndexDifference(
            mlib, bins=mdt.uniform_bins(21, -10, 1))
        absdiff = mdt.features.ResidueIndexDifference(
            mlib, absolute=True, bins=mdt.uniform_bins(21, -10, 1))
        aln = modeller.Alignment(env, file='test/data/alignment.ali',
                                 align_codes='5fd1')
        m1 = mdt.Table(mlib, features=diff)
        m2 = mdt.Table(mlib, features=absdiff)
        self.assertEqual(m1.symmetric, False)
        self.assertEqual(m2.symmetric, True)
        m1.add_alignment(aln, residue_span_range=(-999, -2, 2, 999))
        m2.add_alignment(aln, residue_span_range=(-999, -2, 2, 999))
        self.assertEqual(m1.sum(), 10920)
        self.assertEqual(m2.sum(), 5460)
        # span range should result in 0, +/- 1 bins being zero:
        for m in (m1, m2):
            self.assertEqual(m[9], 0.)
            self.assertEqual(m[10], 0.)
            self.assertEqual(m[11], 0.)
        # Non-absolute feature should have other bins
        # symmetrically distributed:
        for i in range(9):
            self.assertEqual(m1[i], m[-2 - i])
        # Absolute feature should have no negative values:
        for i in range(9):
            self.assertEqual(m2[i], 0.)

    def test_feature_bond_type(self):
        """Check bond type features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.bond_classes.read('data/bndgrp.lib')
        bondtype = mdt.features.BondType(mlib)
        bondlen = mdt.features.BondLength(
            mlib, bins=mdt.uniform_bins(200, 1.0, 0.005))
        self.assertRaises(mdt.MDTError, mlib.bond_classes.read,
                          'data/bndgrp.lib')
        m = mdt.Table(mlib, features=bondtype)
        m2 = mdt.Table(mlib, features=bondlen)
        aln = modeller.Alignment(env, file='test/data/alignment.ali')
        m.add_alignment(aln)
        m2.add_alignment(aln)
        self.assertAlmostEqual(m[0], 7.0, delta=0.0005)
        self.assertAlmostEqual(m[5], 9.0, delta=0.0005)
        self.assertAlmostEqual(m[19], 14.0, delta=0.0005)
        self.assertEqual(m[-1], 0.0)
        self.assertEqual(m.shape, (174,))
        self.assertAlmostEqual(m2[0], 0.0, delta=0.0005)
        self.assertAlmostEqual(m2[43], 3.0, delta=0.0005)
        self.assertAlmostEqual(m2[44], 10.0, delta=0.0005)
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
        self.assertRaises(mdt.MDTError, mlib.angle_classes.read,
                          'data/anggrp.lib')
        m = mdt.Table(mlib, features=angletype)
        m2 = mdt.Table(mlib, features=angle)
        aln = modeller.Alignment(env, file='test/data/alignment.ali')
        m.add_alignment(aln)
        m2.add_alignment(aln)
        self.assertAlmostEqual(m[0], 7.0, delta=0.0005)
        self.assertAlmostEqual(m[7], 9.0, delta=0.0005)
        self.assertAlmostEqual(m[15], 11.0, delta=0.0005)
        self.assertEqual(m.shape, (236,))
        self.assertEqual(m[-1], 0.0)
        self.assertAlmostEqual(m2[176], 48.0, delta=1.0005)
        self.assertAlmostEqual(m2[177], 42.0, delta=0.0005)
        self.assertAlmostEqual(m2[178], 38.0, delta=0.0005)
        self.assertEqual(m2.shape, (289,))
        self.assertEqual(m2[-1], 0.0)
        # Exercise writing of angle class information to HDF5 files:
        m.write_hdf5('test.hdf5')
        os.unlink('test.hdf5')

    def test_feature_distance_undefined(self):
        """Check atom-atom distance feature undefined bin"""
        mlib = self.get_mdt_library()
        dist = mdt.features.AtomDistance(mlib,
                                         bins=mdt.uniform_bins(1, -1000, 2000))
        mdl = self.build_test_model()
        m = self.build_mdt_from_model(mlib, dist, mdl,
                                      residue_span_range=(-99999, 0, 0, 99999))
        self.assertEqual(m.shape, (2,))
        self.assertEqual(m[0], 21.0)
        self.assertEqual(m[1], 0.0)
        # If any coordinate is undefined, the distance is
        modeller.Selection(mdl.atoms[0]).unbuild()
        m = self.build_mdt_from_model(mlib, dist, mdl,
                                      residue_span_range=(-99999, 0, 0, 99999))
        self.assertEqual(m.shape, (2,))
        self.assertEqual(m[0], 15.0)
        self.assertEqual(m[1], 6.0)

    def test_feature_angle_undefined(self):
        """Check angle feature undefined bin"""
        mlib = self.get_mdt_library()
        mlib.angle_classes.read('data/anggrp.lib')
        angle = mdt.features.Angle(mlib,
                                   bins=mdt.uniform_bins(1, -500, 1000))
        mdl = self.build_test_model()
        m = self.build_mdt_from_model(mlib, angle, mdl,
                                      residue_span_range=(-99999, 0, 0, 99999))
        self.assertEqual(m.shape, (2,))
        self.assertEqual(m[0], 5.0)
        self.assertEqual(m[1], 0.0)
        # If any coordinate is undefined, the angle is
        modeller.Selection(mdl.atoms[0]).unbuild()
        m = self.build_mdt_from_model(mlib, angle, mdl,
                                      residue_span_range=(-99999, 0, 0, 99999))
        self.assertEqual(m.shape, (2,))
        self.assertEqual(m[0], 3.0)
        self.assertEqual(m[1], 2.0)

    def test_feature_dihedral_undefined(self):
        """Check dihedral feature undefined bin"""
        mlib = self.get_mdt_library()
        mlib.dihedral_classes.read('data/impgrp.lib')
        dih = mdt.features.Dihedral(mlib,
                                    bins=mdt.uniform_bins(1, -500, 1000))
        mdl = self.build_test_model()
        m = self.build_mdt_from_model(mlib, dih, mdl,
                                      residue_span_range=(-99999, 0, 0, 99999))
        self.assertEqual(m.shape, (2,))
        self.assertEqual(m[0], 1.0)
        self.assertEqual(m[1], 0.0)
        # If any coordinate is undefined, the dihedral is
        modeller.Selection(mdl.atoms[0]).unbuild()
        m = self.build_mdt_from_model(mlib, dih, mdl,
                                      residue_span_range=(-99999, 0, 0, 99999))
        self.assertEqual(m.shape, (2,))
        self.assertEqual(m[0], 0.0)
        self.assertEqual(m[1], 1.0)

    def test_feature_dihedral_type(self):
        """Check dihedral type features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.dihedral_classes.read('data/impgrp.lib')
        dihedtype = mdt.features.DihedralType(mlib)
        dihedral = mdt.features.Dihedral(
            mlib, bins=mdt.uniform_bins(288, -180, 1.25))
        self.assertRaises(mdt.MDTError, mlib.dihedral_classes.read,
                          'data/impgrp.lib')
        m = mdt.Table(mlib, features=dihedtype)
        m2 = mdt.Table(mlib, features=dihedral)
        aln = modeller.Alignment(env, file='test/data/alignment.ali')
        m.add_alignment(aln)
        m2.add_alignment(aln)
        self.assertAlmostEqual(m[0], 7.0, delta=0.0005)
        self.assertAlmostEqual(m[2], 9.0, delta=0.0005)
        self.assertAlmostEqual(m[4], 11.0, delta=0.0005)
        self.assertEqual(m.shape, (79,))
        self.assertEqual(m[-1], 0.0)
        self.assertAlmostEqual(m2[143], 60.0, delta=1.0005)
        self.assertAlmostEqual(m2[144], 53.0, delta=1.0005)
        self.assertAlmostEqual(m2[145], 24.0, delta=0.0005)
        self.assertEqual(m2.shape, (289,))
        self.assertEqual(m2[-1], 0.0)
        # Exercise writing of dihedral class information to HDF5 files:
        m.write_hdf5('test.hdf5')
        os.unlink('test.hdf5')

    def test_feature_doublet_type(self):
        """Check doublet type features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.tuple_classes.read('data/dblcls.lib')
        tuple_angle2 = mdt.features.TupleAngle2(
            mlib, bins=mdt.uniform_bins(6, 0, 30.0))
        tuple_dihed1 = mdt.features.TupleDihedral1(
            mlib, bins=mdt.uniform_bins(6, -180, 60.0))
        self.assertRaises(mdt.MDTError, mlib.tuple_classes.read,
                          'data/dblcls.lib')
        # These features only work on atom triplets:
        for f in (mdt.features.TupleDihedral2, mdt.features.TupleDihedral3):
            self.assertRaises(mdt.MDTError, f, mlib,
                              bins=mdt.uniform_bins(6, -180, 60.0))
        m1 = mdt.Table(mlib, features=tuple_angle2)
        m2 = mdt.Table(mlib, features=tuple_dihed1)
        aln = modeller.Alignment(env, file='test/data/tiny.ali')
        for m in (m1, m2):
            m.add_alignment(aln, residue_span_range=(-9999, 0, 0, 9999))
        self.assertEqual(m1.shape, (7,))
        self.assertEqual(m2.shape, (7,))
        self.assertAlmostEqual(m1[0], 311.0, delta=0.0005)
        self.assertAlmostEqual(m2[0], 302.0, delta=0.0005)

    def test_feature_triplet_residue(self):
        """Check triplet features with residue qualifier"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.tuple_classes.read('test/data/trpcls-residue.lib')
        feat = mdt.features.TupleType(mlib)
        m = mdt.Table(mlib, features=feat)

        mdl = modeller.Model(env)
        mdl.build_sequence('AAACAAACSAA')
        a = modeller.Alignment(env)
        a.append_model(mdl, align_codes='test')

        m.add_alignment(a)
        self.assertEqual([x for x in m], [6.0, 2.0, 1.0, 1.0, 0.0, 0.0])

    def test_feature_cluster(self):
        """Check Cluster feature"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.tuple_classes.read('test/data/trpcls-residue.lib')
        t1 = mdt.features.TupleType(mlib)
        t2 = mdt.features.TupleType(mlib, pos2=True)
        c = mdt.features.Cluster(mlib, t1, t2, nbins=5)
        c.add((0, 0), 0)
        c.add((1, 1), 0)
        c.add((1, 2), 1)
        self.assertRaises(ValueError, c.add, (19, 2), 1)
        self.assertRaises(ValueError, c.add, (1, 29), 1)
        self.assertRaises(ValueError, c.add, (1, 2), 6)

        mdl = modeller.Model(env)
        mdl.build_sequence('AAACAAACSAA')
        a = modeller.Alignment(env)
        a.append_model(mdl, align_codes='test')

        m1 = mdt.Table(mlib, features=(t1, t2))
        m1.add_alignment(a)
        self.assertAlmostEqual(m1[0, 0], 24.0, delta=1e-5)
        self.assertAlmostEqual(m1[1, 1], 2.0, delta=1e-5)
        self.assertAlmostEqual(m1[1, 2], 2.0, delta=1e-5)
        self.assertAlmostEqual(m1.sample_size, 70.0, delta=1e-5)

        m2 = mdt.Table(mlib, features=c)
        m2.add_alignment(a)
        m2.write_hdf5('testcluster.hdf5')
        os.unlink('testcluster.hdf5')
        # m2[0] should be the sum of m1[0,0] and m1[1,1]
        self.assertAlmostEqual(m2[0], 26.0, delta=1e-5)
        # m2[1] should be equal to m1[1,2]
        self.assertAlmostEqual(m2[1], 2.0, delta=1e-5)
        # Everything else should be undefined
        self.assertAlmostEqual(m2[5], 42.0, delta=1e-5)
        # Sample size should be the same
        self.assertAlmostEqual(m2.sample_size, 70.0, delta=1e-5)
        self.assertEqual(m1.symmetric, m2.symmetric)

    def test_feature_triplet_type(self):
        """Check triplet type features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.tuple_classes.read('data/trpcls.lib')
        tuple_type = mdt.features.TupleType(mlib)
        tuple_type2 = mdt.features.TupleType(mlib, pos2=True)
        tuple_dist = mdt.features.TupleDistance(
            mlib, bins=mdt.uniform_bins(9, 2.0, 0.2))
        tuple_angle1 = mdt.features.TupleAngle1(
            mlib, bins=mdt.uniform_bins(6, 0, 30.0))
        tuple_dihed1 = mdt.features.TupleDihedral1(
            mlib, bins=mdt.uniform_bins(6, -180, 60.0))
        tuple_dihed2 = mdt.features.TupleDihedral2(
            mlib, bins=mdt.uniform_bins(6, -180, 60.0))
        tuple_dihed3 = mdt.features.TupleDihedral3(
            mlib, bins=mdt.uniform_bins(6, -180, 60.0))
        self.assertRaises(mdt.MDTError, mlib.tuple_classes.read,
                          'data/trpcls.lib')
        m1 = mdt.Table(mlib, features=tuple_type)
        m2 = mdt.Table(mlib, features=tuple_type2)
        m3 = mdt.Table(mlib, features=tuple_dist)
        m4 = mdt.Table(mlib, features=tuple_angle1)
        m5 = mdt.Table(mlib, features=tuple_dihed1)
        m6 = mdt.Table(mlib, features=tuple_dihed2)
        m7 = mdt.Table(mlib, features=tuple_dihed3)
        aln = modeller.Alignment(env, file='test/data/tiny.ali')
        for m in (m1, m2, m3, m4, m5, m6, m7):
            m.add_alignment(aln, residue_span_range=(-9999, 0, 0, 9999))
        self.assertAlmostEqual(m1[0], 1.0, delta=0.0005)
        self.assertAlmostEqual(m1[1], 0.0, delta=0.0005)
        self.assertAlmostEqual(m1[2], 1.0, delta=0.0005)
        self.assertEqual(m1.shape, (236,))
        self.assertEqual(m1[-1], 0.0)
        self.assertAlmostEqual(m2[0], 60.0, delta=0.0005)
        self.assertAlmostEqual(m2[1], 0.0, delta=0.0005)
        self.assertAlmostEqual(m2[2], 60.0, delta=0.0005)
        self.assertEqual(m2.shape, (236,))
        self.assertEqual(m2[-1], 0.0)
        self.assertAlmostEqual(m3[0], 0.0, delta=0.0005)
        self.assertAlmostEqual(m3[1], 82.0, delta=0.0005)
        self.assertAlmostEqual(m3[2], 226.0, delta=0.0005)
        self.assertEqual(m3.shape, (10,))
        self.assertAlmostEqual(m3[-1], 3018.0, delta=0.0005)
        self.assertAlmostEqual(m4[0], 479.0, delta=0.0005)
        self.assertAlmostEqual(m4[1], 806.0, delta=0.0005)
        self.assertAlmostEqual(m4[2], 471.0, delta=0.0005)
        self.assertEqual(m4.shape, (7,))
        self.assertEqual(m4[-1], 0.0)
        self.assertAlmostEqual(m5[0], 556.0, delta=0.0005)
        self.assertAlmostEqual(m5[1], 642.0, delta=0.0005)
        self.assertAlmostEqual(m5[2], 470.0, delta=6.0005)
        self.assertEqual(m5.shape, (7,))
        self.assertAlmostEqual(m5[-1], 180.0, delta=0.0005)
        self.assertAlmostEqual(m6[0], 661.0, delta=0.0005)
        self.assertAlmostEqual(m6[1], 520.0, delta=0.0005)
        self.assertAlmostEqual(m6[2], 545.0, delta=6.0005)
        self.assertEqual(m6.shape, (7,))
        self.assertAlmostEqual(m6[-1], 112.0, delta=0.0005)
        self.assertAlmostEqual(m7[0], 661.0, delta=0.0005)
        self.assertAlmostEqual(m7[1], 520.0, delta=0.0005)
        self.assertAlmostEqual(m7[2], 545.0, delta=6.0005)
        self.assertEqual(m7.shape, (7,))
        self.assertAlmostEqual(m7[-1], 112.0, delta=0.0005)

    def test_feature_chi1_dihedral(self):
        """Check chi1 dihedral and dihedral class features"""
        mlib = self.get_mdt_library()
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        self.assertRaises(ValueError, mdt.features.Chi1Class, mlib, protein=3)
        chi1class = mdt.features.Chi1Class(mlib)
        m = self.get_test_mdt(mlib, features=chi1)
        self.assertEqual(m.features[0].periodic, True)
        self.assertEqual(m.shape, (37,))
        self.assertAlmostEqual(m[0], 8.0, delta=0.0005)
        self.assertAlmostEqual(m[1], 7.0, delta=0.0005)
        self.assertAlmostEqual(m[2], 2.0, delta=0.0005)
        m = self.get_test_mdt(mlib, features=chi1class)
        self.assertEqual(m.features[0].periodic, False)
        self.assertEqual([b for b in m], [50, 31, 15, 10])

    def test_feature_chi2_dihedral(self):
        """Check chi2 dihedral and dihedral class features"""
        mlib = self.get_mdt_library()
        chi2 = mdt.features.Chi2Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        chi2class = mdt.features.Chi2Class(mlib)
        m = self.get_test_mdt(mlib, features=chi2)
        self.assertEqual(m.shape, (37,))
        self.assertAlmostEqual(m[0], 6.0, delta=0.0005)
        self.assertAlmostEqual(m[1], 4.0, delta=0.0005)
        self.assertAlmostEqual(m[2], 3.0, delta=0.0005)
        m = self.get_test_mdt(mlib, features=chi2class)
        self.assertEqual([b for b in m], [47, 23, 4, 32])

    def test_feature_chi3_dihedral(self):
        """Check chi3 dihedral and dihedral class features"""
        mlib = self.get_mdt_library()
        chi3 = mdt.features.Chi3Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        chi3class = mdt.features.Chi3Class(mlib)
        m = self.get_test_mdt(mlib, features=chi3)
        self.assertEqual(m.shape, (37,))
        self.assertAlmostEqual(m[11], 6.0, delta=0.0005)
        self.assertAlmostEqual(m[12], 4.0, delta=0.0005)
        self.assertAlmostEqual(m[13], 4.0, delta=0.0005)
        m = self.get_test_mdt(mlib, features=chi3class)
        self.assertEqual([b for b in m], [30, 5, 0, 71])

    def test_feature_chi4_dihedral(self):
        """Check chi4 dihedral and dihedral class features"""
        mlib = self.get_mdt_library()
        chi4 = mdt.features.Chi4Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        chi4class = mdt.features.Chi4Class(mlib)
        m = self.get_test_mdt(mlib, features=chi4)
        self.assertEqual(m.shape, (37,))
        self.assertAlmostEqual(m[0], 1.0, delta=0.0005)
        self.assertAlmostEqual(m[1], 1.0, delta=0.0005)
        self.assertAlmostEqual(m[2], 0.0, delta=0.0005)
        m = self.get_test_mdt(mlib, features=chi4class)
        self.assertEqual([b for b in m], [5, 1, 1, 99])

    def test_feature_chi5_dihedral(self):
        """Check chi5 dihedral class feature"""
        mlib = self.get_mdt_library()
        chi5class = mdt.features.Chi5Class(mlib)
        m = self.get_test_mdt(mlib, features=chi5class)
        self.assertEqual([b for b in m], [0, 0, 0, 106])

    def test_feature_phi_dihedral(self):
        """Check phi dihedral and dihedral class features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        phi = mdt.features.PhiDihedral(mlib,
                                       mdt.uniform_bins(36, -180, 10))
        phidiff = mdt.features.PhiDihedralDifference(
            mlib, mdt.uniform_bins(36, -180, 10))
        phiclass = mdt.features.PhiClass(mlib)
        m = self.get_test_mdt(mlib, features=phi)
        self.assertEqual(m.shape, (37,))
        self.assertAlmostEqual(m[5], 4.0, delta=0.0005)
        self.assertAlmostEqual(m[6], 10.0, delta=0.0005)
        self.assertAlmostEqual(m[7], 6.0, delta=0.0005)
        m = self.get_test_mdt(mlib, features=phiclass)
        self.assertEqual([b for b in m], [62, 9, 34, 1])
        m = mdt.Table(mlib, features=phidiff)
        a = modeller.Alignment(env, file='test/data/struc-struc.ali')
        m.add_alignment(a)
        self.assertEqual(m[18], 2)
        self.assertEqual(m[19], 3)
        self.assertEqual(m[20], 0)

    def test_feature_psi_dihedral(self):
        """Check psi dihedral and dihedral class features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        psi = mdt.features.PsiDihedral(mlib,
                                       mdt.uniform_bins(36, -180, 10))
        psidiff = mdt.features.PsiDihedralDifference(
            mlib, mdt.uniform_bins(36, -180, 10))
        psiclass = mdt.features.PsiClass(mlib)
        m = self.get_test_mdt(mlib, features=psi)
        self.assertEqual(m.shape, (37,))
        self.assertAlmostEqual(m[14], 15.0, delta=0.0005)
        self.assertAlmostEqual(m[15], 13.0, delta=0.0005)
        self.assertAlmostEqual(m[16], 6.0, delta=0.0005)
        m = self.get_test_mdt(mlib, features=psiclass)
        self.assertEqual([b for b in m], [47, 23, 35, 1])
        m = mdt.Table(mlib, features=psidiff)
        a = modeller.Alignment(env, file='test/data/struc-struc.ali')
        m.add_alignment(a)
        self.assertEqual(m[19], 2)
        self.assertEqual(m[20], 1)
        self.assertEqual(m[21], 0)

    def test_dihedral_diff_periodic(self):
        """Make sure that dihedral difference features are periodic"""
        def set_omega(mdl, angle):
            chain = mdl.chains[0]
            ca = chain.atoms['CA:1']
            c = chain.atoms['C:1']
            n2 = chain.atoms['N:2']
            ca2 = chain.atoms['CA:2']
            n2.x = n2.y = n2.z = 0.0
            c.x = -2.0
            c.y = c.z = 0.0
            ca.x = -2.0
            ca.y = 2.0
            ca.z = 0.0
            ca2.x = 0.0
            ca2.y = 2.0 * math.cos(math.pi * angle / 180.0)
            ca2.z = 2.0 * math.sin(math.pi * angle / 180.0)
        env = self.get_environ()
        mlib = self.get_mdt_library()
        # Make bins start at slightly less than -180, to allow for floating
        # point rounding
        omegadiff = mdt.features.OmegaDihedralDifference(
            mlib, mdt.uniform_bins(36, -180.01, 10))
        # Note that difference must be shortest around the circle, so
        # 100.0 - (-100.0) is not 200 degrees but -160 degrees
        for dih1, dih2, expected in ((80.0, 80.0, 0.0),
                                     (80.0, -80.0, -160.0),
                                     (-80.0, 80.0, 160.0),
                                     (-100.0, 100.0, -160.0),
                                     (100.0, -100.0, 160.0)):
            m = mdt.Table(mlib, features=omegadiff)
            a = modeller.Alignment(env)
            for d in dih1, dih2:
                mdl = modeller.Model(env)
                mdl.build_sequence('CC')
                set_omega(mdl, d)
                a.append_model(mdl, atom_files='test', align_codes='test')
            m.add_alignment(a, sympairs=True)
            # 2 data points, 1 for each residue
            self.assertAlmostEqual(m.sample_size, 2.0, delta=1e-5)
            # Last residue has no omega, so is always undefined
            self.assertAlmostEqual(m[-1], 1.0, delta=1e-5)
            expected_bin = int((expected + 180.0) / 10.0)
            self.assertAlmostEqual(m[expected_bin], 1.0, delta=1e-5)

    def test_feature_omega_dihedral(self):
        """Check omega dihedral and dihedral class features"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        omega = mdt.features.OmegaDihedral(mlib,
                                           mdt.uniform_bins(36, -180, 10))
        omegadiff = mdt.features.OmegaDihedralDifference(
            mlib, mdt.uniform_bins(36, -180, 10))
        omegaclass = mdt.features.OmegaClass(mlib)
        m = self.get_test_mdt(mlib, features=omega)
        self.assertEqual(m.shape, (37,))
        self.assertAlmostEqual(m[0], 44.0, delta=0.0005)
        self.assertAlmostEqual(m[1], 5.0, delta=0.0005)
        self.assertAlmostEqual(m[2], 0.0, delta=0.0005)
        m = self.get_test_mdt(mlib, features=omegaclass)
        self.assertEqual([b for b in m], [105, 0, 0, 1])
        m = mdt.Table(mlib, features=omegadiff)
        a = modeller.Alignment(env, file='test/data/struc-struc.ali')
        m.add_alignment(a)
        self.assertEqual(m[17], 5)
        self.assertEqual(m[18], 5)
        self.assertEqual(m[36], 2)

    def test_feature_alpha_dihedral(self):
        """Check alpha dihedral feature"""
        mlib = self.get_mdt_library()
        alpha = mdt.features.AlphaDihedral(mlib,
                                           mdt.uniform_bins(36, -180, 10))
        m = self.get_test_mdt(mlib, features=alpha)
        self.assertEqual(m.shape, (37,))
        self.assertAlmostEqual(m[30], 3.0, delta=0.0005)
        self.assertAlmostEqual(m[31], 4.0, delta=0.0005)
        self.assertAlmostEqual(m[32], 2.0, delta=0.0005)

    def test_feature_distance(self):
        """Check atom-atom distance feature"""
        mlib = self.get_mdt_library()
        dist = mdt.features.AtomDistance(mlib,
                                         bins=mdt.uniform_bins(60, 0, 0.5))
        m = self.get_test_mdt(mlib, features=dist)
        self.assertEqual(m.shape, (61,))
        self.assertAlmostEqual(m[30], 10057.0, delta=0.0005)
        self.assertAlmostEqual(m[31], 10214.0, delta=0.0005)
        self.assertAlmostEqual(m[32], 10095.0, delta=0.0005)
        self.assertAlmostEqual(m[-1], 4892.0, delta=0.0005)

    def test_symmetric(self):
        """Test symmetric/asymmetric residue pair features"""
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(10, 0, 1.0)
        dist = mdt.features.ResidueDistance(mlib, bins)
        avresacc = mdt.features.AverageResidueAccessibility(mlib, bins)
        avndif = mdt.features.AverageNeighborhoodDifference(mlib, bins)
        diff = mdt.features.ResidueIndexDifference(mlib, bins)
        ddist = mdt.features.ResidueDistanceDifference(mlib, bins)
        avgapdist = mdt.features.AverageGapDistance(mlib, bins)
        sym_features = (avndif, avresacc, avgapdist)
        asym_features = (dist, ddist, diff)
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
            m = mdt.Table(mlib, features=(a, b))
            self.assertEqual(m.symmetric, symm)

    def test_abstract(self):
        """Should not be able to instantiate abstract features"""
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(10, 0, 1.0)
        for feat in [mdt.features.Protein, mdt.features.ProteinPair,
                     mdt.features.Residue, mdt.features.ResiduePair,
                     mdt.features.AlignedResidue,
                     mdt.features.AlignedResiduePair, mdt.features.Atom,
                     mdt.features.AtomPair, mdt.features.Tuple,
                     mdt.features.TuplePair, mdt.features.ChemicalBond]:
            self.assertRaises(TypeError, feat, mlib, bins)
        for feat in [mdt.features.ResidueFixedBins, mdt.features.AtomFixedBins,
                     mdt.features.TupleFixedBins,
                     mdt.features.ChemicalBondFixedBins]:
            self.assertRaises(TypeError, feat, mlib)

    def test_tuple_base(self):
        """Test otherwise unexercised Tuple constructor"""
        class DummyTuple(mdt.features.Tuple):
            def _setup(self, mlib, pos2):
                return "dummy ifeat"

            def _create_bins(self, mlib, bins):
                self.bins_created = bins
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(10, 0, 1.0)
        f = DummyTuple(mlib, bins)
        self.assertEqual(f._ifeat, 'dummy ifeat')
        self.assertEqual(f.bins_created, bins)


if __name__ == '__main__':
    unittest.main()
