import unittest
from modeller import *
from mdt_test import MDTTest
import mdt
import mdt.features
import math
import os
import sys

class TableTests(MDTTest):

    restyps = "ACDEFGHIKLMNPQRSTVWY"

    def __mdt_restyp(self, a):
        """Convert a one-letter residue code into the MDT bin number"""
        ind = self.restyps.find(a)
        # Map to 'undefined' bin:
        if ind == -1:
            ind = 21
        return ind

    def test_version(self):
        """Check version number"""
        self.assertGreater(len(mdt.version), 0)
        self.assertGreater(len(mdt.__version__), 0)
        self.assertGreater(len(mdt.version_info), 0)

    def test_add(self):
        """Check adding MDTs"""
        mlib = self.get_mdt_library()
        bins = mdt.uniform_bins(3, -1.0, 1.5)
        xray0 = mdt.features.XRayResolution(mlib, bins, protein=0)
        xray1 = mdt.features.XRayResolution(mlib, bins, protein=1)
        m1 = mdt.Table(mlib, features=xray0)
        for (n, val) in enumerate((1,2,3,4)):
            m1[n] = val
        m2 = mdt.Table(mlib, features=xray0)
        for (n, val) in enumerate((10,20,30,40)):
            m2[n] = val
        m3 = m1 + m2
        m1 += m2
        self.assertMDTsEqual(m1, m3)
        for (n, val) in enumerate((11,22,33,44)):
            self.assertEqual(m3[n], val)
        # Cannot add if numbers of features are different
        badmdt = mdt.Table(mlib, features=(xray0, xray1))
        self.assertRaises(ValueError, m1.__add__, badmdt)

        # Cannot add if feature types are different
        badmdt = mdt.Table(mlib, features=xray1)
        self.assertRaises(ValueError, m1.__add__, badmdt)

        # Cannot add if starts are different
        badmdt = m2.reshape(features=xray0, offset=1, shape=0)
        self.assertRaises(ValueError, m1.__add__, badmdt)

        # Cannot add if nbins are different
        badmdt = m2.reshape(features=xray0, offset=0, shape=-1)
        self.assertRaises(ValueError, m1.__add__, badmdt)

    def test_features_same_library(self):
        """Features must all come from the same Library"""
        mlib1 = self.get_mdt_library()
        mlib2 = self.get_mdt_library()
        restyp1 = mdt.features.ResidueType(mlib1, protein=0)
        restyp2 = mdt.features.ResidueType(mlib2, protein=1)
        m = mdt.Table(mlib1, features=restyp1)
        m = mdt.Table(mlib2, features=restyp2)
        for mlib in (mlib1, mlib2):
            self.assertRaises(ValueError, mdt.Table, mlib,
                              features=(restyp1,restyp2))

    def test_make_outrange(self):
        """Check handling of out-of-range features in mdt.make()"""
        class OutOfRangeFeature(object):
            def _get_ifeat(self, mlib):
                return 0
        mlib = self.get_mdt_library()
        f1 = OutOfRangeFeature()
        self.assertRaises(IndexError, mdt.Table, mlib, features=f1)

    def test_make_shape(self):
        """Test Table.make() with shape"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        for shape in [[0,0], [-22], [23]]:
            self.assertRaises(ValueError, mdt.Table, mlib, features=restyp,
                              shape=shape)
        for shape in [[], [0], [22]]:
            t = mdt.Table(mlib, features=restyp, shape=shape)
            self.assertEqual(t.shape, (22,))
        t = mdt.Table(mlib, features=restyp, shape=[10])
        self.assertEqual(t.shape, (10,))
        t = mdt.Table(mlib, features=restyp, shape=[-4])
        self.assertEqual(t.shape, (18,))

    def test_clear(self):
        """Check Table.clear()"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        t = mdt.Table(mlib, features=restyp)
        t[0] = 5
        t.clear()
        self.assertEqual(t[0], 0)

    def test_bad_read(self):
        """Check error handling in Table.read()"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        feat = """
  #  FEATURE NBINS NAME
  1        1    22 Residue type of protein 0

"""

        with open('bad.mdt', 'w') as fh:
            fh.write("%-37s: garbage\n" % "Number of alignments")
        self.assertRaises(mdt.FileFormatError, mdt.Table, mlib, file='bad.mdt')

        with open('bad.mdt', 'w') as fh:
            fh.write("%-37s: garbage\n" % "Sample size")
        self.assertRaises(mdt.FileFormatError, mdt.Table, mlib, file='bad.mdt')

        with open('bad.mdt', 'w') as fh:
            fh.write("  #  FEATURE NBINS NAME\ngarbage\n")
        self.assertRaises(mdt.FileFormatError, mdt.Table, mlib, file='bad.mdt')

        with open('bad.mdt', 'w') as fh:
            fh.write(feat + "  # ISTART   IEND\n")
        self.assertRaises(mdt.FileFormatError, mdt.Table, mlib, file='bad.mdt')

        with open('bad.mdt', 'w') as fh:
            fh.write(feat + "  # ISTART   IEND\ngarbage\n")
        self.assertRaises(mdt.FileFormatError, mdt.Table, mlib, file='bad.mdt')

        with open('bad.mdt', 'w') as fh:
            fh.write(feat + "MDT TABLE START: garbage\n")
        self.assertRaises(mdt.FileFormatError, mdt.Table, mlib, file='bad.mdt')

        with open('bad.mdt', 'w') as fh:
            fh.write(feat + "MDT TABLE START:       1\n")
        self.assertRaises(mdt.FileFormatError, mdt.Table, mlib, file='bad.mdt')

        with open('bad.mdt', 'w') as fh:
            fh.write(feat + "MDT TABLE START:       1\ngarbage\n")
        self.assertRaises(mdt.FileFormatError, mdt.Table, mlib, file='bad.mdt')

        with open('bad.mdt', 'w') as fh:
            fh.write(feat + "MDT TABLE START:       1\n1.0\n")
        self.assertRaises(mdt.FileFormatError, mdt.Table, mlib, file='bad.mdt')

        with open('bad.mdt', 'w') as fh:
            fh.write(feat + "MDT TABLE START:       1\n1.0\nf\n")
        self.assertRaises(mdt.FileFormatError, mdt.Table, mlib, file='bad.mdt')

        os.unlink('bad.mdt')

    def test_mdt_formats(self):
        """Make sure we can read and write MDT files"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = mdt.Table(mlib, features=(restyp0,restyp1))
        aln = alignment(env)
        aln.append_sequence('AFVVTDNCIKCKYTDCVEVXPVDCFYEG')
        aln.append_sequence('CVEVCPVDCFYEGAFVVTDNCIKCKYTX')
        m.add_alignment(aln)
        m.write('test.mdt')
        self.assertRaises(IOError, m.write, '/does/not/exist/foo.mdt')
        m2 = m.copy()
        self.assertMDTsEqual(m, m2)
        m2 = mdt.Table(mlib, file='test.mdt')
        self.assertMDTsEqual(m, m2)
        m.write_hdf5('test.hdf5')
        m.write_hdf5('test.hdf5', gzip=True)
        m.write_hdf5('test.hdf5', gzip=4)
        self.assertRaises(ValueError, m.write_hdf5, 'test.hdf5',
                          gzip=True, chunk_size=(2,3,4))
        m.write_hdf5('test.hdf5', gzip=9, chunk_size=(4,4))
        self.assertRaises(mdt.MDTError, m.write_hdf5,
                          '/does/not/exist/foo.hdf5')
        m2 = mdt.Table(mlib, file='test.hdf5')
        self.assertMDTsEqual(m, m2)
        m2 = mdt.Table(mlib)
        # Trying to read an HDF5 file in text format should fail gracefully:
        self.assertRaises(mdt.MDTError, m2.read, 'test.hdf5')
        # Same problem should occur starting with a non-empty MDT:
        m2 = mdt.Table(mlib, features=(restyp0,restyp1))
        self.assertRaises(mdt.MDTError, m2.read, 'test.hdf5')
        os.unlink('test.hdf5')

    def test_feature_check(self):
        """When rereading MDTs, features should be the same"""
        mlib = self.get_mdt_library()
        feat = mdt.features.AtomDistance(mlib, mdt.uniform_bins(10, 0, 1))
        m = mdt.Table(mlib, features=feat)
        m.write(file='test.mdt')
        m.write_hdf5(file='test.hdf5')
        # Different feature type, same bins should fail:
        mlib = self.get_mdt_library()
        feat = mdt.features.ResidueDistance(mlib, mdt.uniform_bins(10, 0, 1))
        for f in ('test.mdt', 'test.hdf5'):
            self.assertRaises(mdt.MDTError, mdt.Table, mlib, file=f)
        # Same feature type, different number of bins should fail:
        mlib = self.get_mdt_library()
        feat = mdt.features.AtomDistance(mlib, mdt.uniform_bins(20, 0, 1))
        for f in ('test.mdt', 'test.hdf5'):
            self.assertRaises(mdt.MDTError, mdt.Table, mlib, file=f)
        os.unlink('test.mdt')
        os.unlink('test.hdf5')

    def test_guess_chunk_size(self):
        """Test guess_chunk_size method"""
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        m = self.get_test_mdt(mlib, features=restyp)
        self.assertEqual(m._guess_chunk_size([12,6,6,12,306,306], 1024*1024*10),
                         [6,3,3,6,161,161])

    def test_bin_info(self):
        """Test query of bin symbol and range"""
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        m = self.get_test_mdt(mlib, features=restyp)
        for (n, bin) in enumerate(m.features[0].bins):
            self.assertEqual(bin.range, (n, n+1))
        symbols = [bin.symbol for bin in m.features[0].bins]
        self.assertEqual(symbols, ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY',
                                   'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN',
                                   'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL',
                                   'TRP', 'TYR', 'GAP', 'u'])
        m2 = m.reshape(features=restyp, offset=1, shape=20)
        for (n, bin) in enumerate(m2.features[0].bins):
            self.assertEqual(bin.symbol, m.features[0].bins[n+1].symbol)

    def test_set(self):
        """Test set of MDT data"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = mdt.Table(mlib, features=(restyp0,restyp1))
        # Make sure that index checks work:
        self.assertRaises(ValueError, m.__setitem__, [0,22,10], 0.0)
        self.assertRaises(IndexError, m.__setitem__, [0,22], 0.0)
        self.assertRaises(IndexError, m.__setitem__, [22,0], 0.0)
        # Cannot set a section
        self.assertRaises(ValueError, m.__setitem__, [0], [0.0])
        for i in range(0, 22):
            for j in range(0, 22):
                val = 1.0 + i * 40.0 + j   # Some value different for each i,j
                self.assertEqual(m[i,j], 0.0)
                m[i,j] = val
                self.assertEqual(m[i,j], val)

    def test_empty_aln(self):
        """Make sure that adding an empty alignment works OK"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        bondtype = mdt.features.BondType(mlib)
        m = mdt.Table(mlib, features=(restyp,bondtype))
        aln = alignment(env)
        self.assertRaises(mdt.MDTError, m.add_alignment, aln)

    def test_1d_mdt(self):
        """Make sure a 1D MDT matches known residue data"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        m = mdt.Table(mlib, features=restyp)
        aln = alignment(env)
        seq = "AFVVTDNCIKXCKYTDCVEVCPVDCFYEG"
        aln.append_sequence(seq)
        # Get known counts of all residue types (including undefined)
        seq_types = [self.__mdt_restyp(a) for a in seq]
        known_dist = [seq_types.count(n) for n in range(22)]

        # Now compare with that obtained by mdt.add_alignment()
        m.add_alignment(aln)
        for i in range(22):
            self.assertEqual(known_dist[i], int(m[i]))

        # Now do the same test on individually-queried indices:
        source = m.open_alignment(aln)
        bins = [ source.index(restyp, 0, n, 0, n, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0) - 1 for n in range(len(seq)) ]
        bin_dist = [bins.count(n) for n in range(22)]
        self.assertEqual(bin_dist, known_dist)

        # source.sum() should return each sample nsample times, so it should be
        # equal to the sum of the squares:
        self.assertEqual(source.sum(), sum([i * i for i in known_dist]))

    def test_residue_span_range(self):
        """Test residue_span_range argument"""
        env = self.get_environ()
        aln = alignment(env, file='test/data/tiny.ali')
        mlib = self.get_mdt_library()
        mlib.tuple_classes.read('data/dblcls.lib')
        tuple_dist = mdt.features.TupleDistance(mlib,
                                          bins=mdt.uniform_bins(49, 2.0, 0.2))
        # All residue-residue pairs should be out of range, so this MDT should
        # end up empty:
        m = mdt.Table(mlib, features=tuple_dist)
        m.add_alignment(aln, residue_span_range=(-999, -999, 999, 999))
        self.assertEqual(m.sum(), 0.0)
        m1 = mdt.Table(mlib, features=tuple_dist)
        m2 = mdt.Table(mlib, features=tuple_dist)
        # When excluding only intra-residue interactions, the short-range
        # bins should differ but the long-range behavior should be the same:
        m1.add_alignment(aln, residue_span_range=(-999, 0, 0, 999))
        m2.add_alignment(aln, residue_span_range=(-999, -1, 1, 999))
        self.assertEqual(m1.shape, m2.shape)
        for i in range(1, 5):
            self.assertNotEqual(m1[i], m2[i])
        for i in range(38, 49):
            self.assertEqual(m1[i], m2[i])

    def test_atom_pair_exclusions(self):
        """Test exclusion of atom pairs from atom pair features"""
        env = self.get_environ()
        mdl = model(env)
        mdl.build_sequence('G')
        aln = alignment(env)
        aln.append_model(mdl, align_codes='test')
        mlib = self.get_mdt_library()
        mlib.bond_classes.read('data/bndgrp.lib')
        mlib.angle_classes.read('data/anggrp.lib')
        mlib.dihedral_classes.read('data/impgrp.lib')
        dist = mdt.features.AtomDistance(mlib,
                                         bins=mdt.uniform_bins(49, 0.0, 0.2))
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, residue_span_range=(-9999,0,0,9999),
                        exclude_bonds=False)
        self.assertEqual(m.sample_size, 10)
        # 3 bonds (N-CA, O-C, C-CA) should be excluded in Gly
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, residue_span_range=(-9999,0,0,9999),
                        exclude_bonds=True)
        self.assertEqual(m.sample_size, 7)
        # A further 2 angles (CA-C-O, N-CA-C) should be excluded in Gly
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, residue_span_range=(-9999,0,0,9999),
                        exclude_bonds=True, exclude_angles=True)
        self.assertEqual(m.sample_size, 5)
        # No improper dihedrals
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, residue_span_range=(-9999,0,0,9999),
                        exclude_bonds=True, exclude_angles=True,
                        exclude_dihedrals=True)
        self.assertEqual(m.sample_size, 5)

        mdl = model(env)
        mdl.build_sequence('GG')
        aln = alignment(env)
        aln.append_model(mdl, align_codes='test')
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, residue_span_range=(-9999,0,0,9999))
        self.assertEqual(m.sample_size, 36)
        # One dihedral (C:CA:+N:O) in Gly-Gly should be excluded
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, residue_span_range=(-9999,0,0,9999),
                        exclude_dihedrals=True)
        self.assertEqual(m.sample_size, 35)

    def test_atom_tuple_pair_exclusions(self):
        """Test exclusion of atom pairs from atom tuple pair features"""
        env = self.get_environ()
        mdl = model(env)
        mdl.build_sequence('GG')
        aln = alignment(env)
        aln.append_model(mdl, align_codes='test')
        mlib = self.get_mdt_library()
        mlib.tuple_classes.read('data/trpcls.lib')
        mlib.bond_classes.read('data/bndgrp.lib')
        mlib.angle_classes.read('data/anggrp.lib')
        mlib.dihedral_classes.read('data/impgrp.lib')
        dist = mdt.features.TupleDistance(mlib,
                                         bins=mdt.uniform_bins(49, 0.0, 0.2))
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, residue_span_range=(-9999,0,0,9999),
                        exclude_bonds=False)
        self.assertEqual(m.sample_size, 38)
        # Exclusions should cut number of sample points
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, residue_span_range=(-9999,0,0,9999),
                        exclude_bonds=True, exclude_angles=True)
        self.assertEqual(m.sample_size, 20)

    def test_chain_span_range(self):
        """Test chain_span_range argument"""
        env = self.get_environ()
        mdl = model(env)
        mdl.build_sequence('G/G')
        aln = alignment(env)
        aln.append_model(mdl, align_codes='test')
        mlib = self.get_mdt_library()
        mlib.tuple_classes.read('data/dblcls.lib')
        tuple_dist = mdt.features.TupleDistance(mlib,
                                          bins=mdt.uniform_bins(49, 2.0, 0.2))
        # All chain differences should be out of range, so this MDT should
        # end up empty:
        m = mdt.Table(mlib, features=tuple_dist)
        m.add_alignment(aln, chain_span_range=(-999, -999, 999, 999),
                        residue_span_range=(-999, 0, 0, 999))
        self.assertEqual(m.sum(), 0.0)
        # Default chain separation should allow intra-chain interactions, so
        # should yield more (56) than only allowing inter-chain
        # interactions (32)
        m = mdt.Table(mlib, features=tuple_dist)
        m.add_alignment(aln, residue_span_range=(-999, 0, 0, 999))
        self.assertEqual(m.sum(), 56.0)
        m = mdt.Table(mlib, features=tuple_dist)
        m.add_alignment(aln, chain_span_range=(-999, -1, 1, 999),
                        residue_span_range=(-999, 0, 0, 999))
        self.assertEqual(m.sum(), 32.0)

    def test_tuple_pair_bond_span_range(self):
        """Test bond_span_range with tuple pair scan"""
        env = self.get_environ()
        mdl = model(env)
        mdl.build_sequence('A')
        aln = alignment(env)
        aln.append_model(mdl, align_codes='test')
        mlib = self.get_mdt_library()
        mlib.bond_classes.read('data/bndgrp.lib')
        mlib.tuple_classes.read('data/trpcls.lib')
        typ = mdt.features.TupleType(mlib)
        typ2 = mdt.features.TupleType(mlib, pos2=True)
        dist = mdt.features.TupleDistance(mlib,
                                          bins=mdt.uniform_bins(9, 2.0, 0.2))

        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, residue_span_range=(0,0,0,0))
        self.assertEqual(m.sample_size, 10.0)

        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, bond_span_range=(1,1),
                        residue_span_range=(0,0,0,0))
        # Bond span should restrict interactions to 6
        # (C:CA:CB-CA:C:O, CA:C:O-N:CA:C, CA:C:O-N:CA:CB, and the reverse)
        self.assertEqual(m.sample_size, 6.0)

    def test_bond_span_range(self):
        """Test bond_span_range argument"""
        env = self.get_environ()
        mdl = model(env)
        mdl.build_sequence('A')
        aln = alignment(env)
        aln.append_model(mdl, align_codes='test')
        mlib = self.get_mdt_library()
        mlib.bond_classes.read('data/bndgrp.lib')
        dist = mdt.features.AtomDistance(mlib,
                                         bins=mdt.uniform_bins(60, 0, 0.5))

        # Only 4 direct chemical bonds (N-CA, CA-CB, CA-C, C-O) in ALA; note
        # that bond library does not include OXT so C-OXT interaction is
        # excluded
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, bond_span_range=(1,1),
                        residue_span_range=(0,0,0,0))
        self.assertEqual(m.sample_size, 4.0)

        # Only 2 dihedrals (N-CA-C-O, O-C-CA-CB)
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, bond_span_range=(3,3),
                        residue_span_range=(0,0,0,0))
        self.assertEqual(m.sample_size, 2.0)

        # 4 bonds, 4 angles and 2 dihedrals: 10 in total
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, bond_span_range=(1,3),
                        residue_span_range=(0,0,0,0))
        self.assertEqual(m.sample_size, 10.0)

        # Check for bonds between residues (just the N-C bond here)
        mdl = model(env)
        mdl.build_sequence('AA')
        aln = alignment(env)
        aln.append_model(mdl, align_codes='test')

        # Force a non-symmetric scan (to check handling of bond separation
        # regardless of which order atom indices are in)
        diff = mdt.features.ResidueIndexDifference(mlib,
                                      bins=mdt.uniform_bins(21, -10, 1))
        m = mdt.Table(mlib, features=(dist,diff))
        m.add_alignment(aln, bond_span_range=(0,1),
                        residue_span_range=(-10,-1,1,10))
        self.assertEqual(m.sample_size, 2.0)

        # Bonds never span chains
        mdl = model(env)
        mdl.build_sequence('A/A')
        aln = alignment(env)
        aln.append_model(mdl, align_codes='test')
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, bond_span_range=(0,99999),
                        residue_span_range=(-10,-1,1,10))
        self.assertEqual(m.sample_size, 0.0)

    def test_sources(self):
        """Make sure that alignments and models both work as sources"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        dist = mdt.features.AtomDistance(mlib,
                                         bins=mdt.uniform_bins(60, 0, 0.5))
        m1 = mdt.Table(mlib, features=dist)
        m2 = mdt.Table(mlib, features=dist)
        a1 = alignment(env, file='test/data/tiny.ali', align_codes='5fd1')
        m1.add_alignment(a1)
        mdl = model(env, file='test/data/5fd1.atm', model_segment=('1:', '6:'))
        a2 = alignment(env)
        # Atom file 'foo' does not exist; all data should be taken from mdl
        a2.append_model(mdl, align_codes='foo', atom_files='foo')
        m2.add_alignment(a2)
        self.assertMDTsEqual(m1, m2)

    def test_feature_combination(self):
        """Check that invalid feature combinations are rejected"""
        mlib = self.get_mdt_library()
        atmdist = mdt.features.AtomDistance(mlib,
                                            bins=mdt.uniform_bins(60, 0, 0.5))
        resdist = mdt.features.ResidueDistance(mlib, protein=1,
                                               bins=mdt.uniform_bins(7, 0, 2.0))
        self.assertRaises(ValueError, self.get_test_mdt, mlib,
                          features=(resdist,atmdist))

    def test_integrate(self):
        """Make sure MDT integration works"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        m1 = mdt.Table(mlib, features=restyp0)
        m2 = mdt.Table(mlib, features=restyp1)
        mboth = mdt.Table(mlib, features=(restyp0,restyp1))
        seq1 = "AFVVTDNCIK"
        seq2 = "DCVEVCPVDC"
        aln = alignment(env)
        aln.append_sequence(seq1)
        aln.append_sequence(seq2)
        for m in (m1, m2, mboth):
            m.add_alignment(aln)

        # Number of features must be correct:
        for features in ((), (restyp0,restyp1), (restyp0,restyp1,chi1)):
            self.assertRaises(ValueError, mboth.integrate, features=features)
        # Features must exist in input MDT:
        self.assertRaises(ValueError, mboth.integrate, features=chi1)
        m1int = mboth.integrate(features=restyp0)
        self.assertMDTsEqual(m1, m1int)
        m2int = mboth.integrate(features=restyp1)
        self.assertMDTsEqual(m2, m2int)

    def test_entropy_hx(self):
        """Check for expected dependent entropy value"""
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        m = self.get_test_mdt(mlib, features=(restyp,chi1))
        # Check against ENTROPY_MDT_HX value for this system in MDT:
        self.assertAlmostEqual(m.entropy_hx(), 2.7048, places=3)

        # Empty table should give an error
        empty = mdt.Table(mlib, features=(restyp,chi1))
        self.assertRaises(mdt.MDTError, empty.entropy_hx)

    def test_entropy_full(self):
        """Test entropy_full()"""
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        m = self.get_test_mdt(mlib, features=(restyp,chi1))
        m.entropy_full()

        # Empty table should give an error
        empty = mdt.Table(mlib, features=(restyp,chi1))
        self.assertRaises(mdt.MDTError, empty.entropy_full)

    def test_exp_transform(self):
        """Check for correctness of exp transform"""
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        phi = mdt.features.PhiDihedral(mlib,
                                       mdt.uniform_bins(36, -180, 10))
        m = self.get_test_mdt(mlib, features=(restyp,phi))
        offset = 1.
        expoffset = 0.1
        multiplier = 0.8
        power = 0.7
        m2 = m.exp_transform(offset, expoffset, multiplier, power)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            y = offset + math.exp(expoffset + multiplier * m[inds] ** power)
            self.assertAlmostEqual(y, m2[inds], places=3)

    def test_log_transform(self):
        """Check for correctness of log transform"""
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        phi = mdt.features.PhiDihedral(mlib,
                                       mdt.uniform_bins(36, -180, 10))
        m = self.get_test_mdt(mlib, features=(restyp,phi))
        offset = 0.0
        multiplier = 0.7
        undefined = 500.0
        m2 = m.log_transform(offset, multiplier, undefined)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            y = offset + multiplier * m[inds]
            if y < 1e-10:
                y = undefined
            else:
                y = math.log(y)
            self.assertAlmostEqual(y, m2[inds], places=3)

    def test_linear_transform(self):
        """Check for correctness of linear transform"""
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        m = self.get_test_mdt(mlib, features=(restyp,chi1))
        offset = -1.3
        multiplier = 0.8
        m2 = m.linear_transform(offset, multiplier)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            y = offset + multiplier * m[inds]
            self.assertAlmostEqual(y, m2[inds], places=3)

    def test_inverse_transform(self):
        """Check for correctness of inverse transform"""
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        phi = mdt.features.PhiDihedral(mlib,
                                       mdt.uniform_bins(36, -180, 10))
        m = self.get_test_mdt(mlib, features=(restyp,phi))
        offset = 0.0
        multiplier = 0.8
        undefined = 300.0
        m2 = m.inverse_transform(offset, multiplier, undefined)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            if abs(m[inds]) < 1e-15:
                y = undefined
            else:
                y = offset + multiplier / m[inds]
            self.assertAlmostEqual(y, m2[inds], places=3)

    def test_offset(self):
        """Check for correctness of offset transforms"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = self.get_test_mdt(mlib, features=(restyp0,restyp1))
        # Get a better range of values in our starting MDT
        m = m.log_transform(0.0, 0.7, 10.0)
        # Only 1D or 2D offsets should be allowed
        for dim in (0, 3):
            self.assertRaises(ValueError, m.offset_min, dimensions=dim)
        m2 = m.offset_min(dimensions=1)
        nx, ny = m.shape
        minminval = 1e9
        for x in range(nx):
            minval = 1e9
            for y in range(ny):
                minval = min(minval, m[x,y])
            minminval = min(minminval, minval)
            for y in range(ny):
                self.assertAlmostEqual(m[x,y] - minval, m2[x,y], places=3)
        m2 = m.offset_min(dimensions=2)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            self.assertAlmostEqual(m[inds] - minminval, m2[inds], places=3)

    def test_close(self):
        """Check for correctness of spline closing"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = self.get_test_mdt(mlib, features=(restyp0,restyp1))
        # Get a better range of values in our starting MDT
        m = m.log_transform(0.0, 0.7, 10.0)
        # Only 1D or 2D should be allowed
        for dim in (0, 3):
            self.assertRaises(ValueError, m.close, dimensions=dim)
        nx, ny = m.shape
        m2 = m.close(dimensions=1)
        for x in range(nx):
            self.assertEqual(m2[x,0], m2[x,ny-1])
        m2 = m.close(dimensions=2)
        self.assertEqual(m2[0,0], m2[0,ny-1])
        self.assertEqual(m2[0,0], m2[nx-1,0])
        self.assertEqual(m2[0,0], m2[nx-1,ny-1])
        for x in range(nx):
            self.assertEqual(m2[x,0], m2[x,ny-1])
        for y in range(ny):
            self.assertEqual(m2[0,y], m2[nx-1,y])

    def test_smooth(self):
        """Check for correctness of smoothing"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = mdt.Table(mlib, features=(restyp0,restyp1))
        # Dimensions must be 1 or 2:
        for dim in (0, 3):
            self.assertRaises(ValueError, m.smooth, dim, 1.0)
        # A smoothed empty MDT should be the prior, i.e. 1/nbin:
        m2 = m.smooth(2, 1.0)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            self.assertAlmostEqual(1./484., m2[inds], places=3)
        m2 = m.smooth(1, 1.0)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            self.assertAlmostEqual(1./22., m2[inds], places=3)

    def test_save_reshape(self):
        """Check that we can correctly load and save reshaped MDTs"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = self.get_test_mdt(mlib, features=(restyp0,restyp1))
        m = m.reshape(features=(restyp0,restyp1), offset=(4,2), shape=(11,10))
        m.write('test.mdt')
        m.write_hdf5('test.hdf5')
        m2 = mdt.Table(mlib, 'test.mdt')
        self.assertMDTsEqual(m, m2)
        m2 = mdt.Table(mlib, 'test.hdf5')
        self.assertMDTsEqual(m, m2)
        os.unlink('test.mdt')
        os.unlink('test.hdf5')

    def test_reshape(self):
        """Check that reshaping works correctly"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        m = self.get_test_mdt(mlib, features=(restyp0,restyp1))
        m2 = self.get_test_mdt(mlib, features=(restyp1,restyp0))
        # New features must be a subset of the old
        for features in ((chi1,restyp0), (restyp1,restyp0,restyp0)):
            self.assertRaises(ValueError, m.reshape, features=features,
                              offset=m.offset, shape=(22,22))
        m3 = m.reshape(features=(restyp1,restyp0), offset=m.offset,
                       shape=(22,22))
        self.assertMDTsEqual(m2, m3)
        m3 = m.reshape(features=(restyp0,restyp1), offset=(0,0), shape=m.shape)
        self.assertMDTsEqual(m, m3)
        m3 = m.reshape(features=(restyp0,restyp1), offset=(4,2), shape=(11,10))
        # Reshaping to same offset and shape should be a no-op
        m4 = m3.reshape(features=(restyp0,restyp1), offset=(4,2), shape=(11,10))
        self.assertEqual(m3.shape, (11,10))
        self.assertEqual(m3.offset, (4,2))
        inds = []
        while self.roll_inds(inds, m3.shape, m3.offset):
            self.assertAlmostEqual(m[inds], m3[inds], places=3)
            self.assertAlmostEqual(m[inds], m4[inds], places=3)
        # Offset cannot be less than old
        for offset in [(3,2), (4,1)]:
            self.assertRaises(IndexError, m3.reshape,
                              features=(restyp0,restyp1), offset=offset,
                              shape=(1,1))
        # Negative shape (old end+shape) must result in at least size 1
        for shape in [(-11,10), (11,-10)]:
            self.assertRaises(IndexError, m3.reshape,
                              features=(restyp0,restyp1), offset=(4,2),
                              shape=shape)
        # Shape cannot be larger than old
        for shape in [(12,10), (11,11)]:
            self.assertRaises(IndexError, m3.reshape,
                              features=(restyp0,restyp1), offset=(4,2),
                              shape=shape)

    def test_reshape_distancemdt(self):
        """Check that reshaping distance mdt produces the table correctly"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        dist = mdt.features.AtomDistance(mlib,
                                         bins=mdt.uniform_bins(30, 0, 0.5))
        aln = alignment(env, file='test/data/alignment.ali', align_codes='5fd1')
        m1 = mdt.Table(mlib, features=dist)
        m1.reshape(dist,[0],[-1])
        m1.add_alignment(aln, residue_span_range=(-999, -1, 1, 999))
        self.assertEqual(m1[29],10014.0)
        self.assertAlmostEqual(m1[28], 9758.5, delta=0.6)

    def test_sum(self):
        """Check that sum of each row sums to that of the whole table"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = self.get_test_mdt(mlib, features=(restyp0,restyp1))
        sum = m.sum()
        rowsum = 0.
        for i in range(len(m.features[0].bins)):
            rowsum += m[i].sum()
        self.assertAlmostEqual(sum, rowsum, places=3)

    def test_entropy(self):
        """Check the entropy of each row"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = self.get_test_mdt(mlib, features=(restyp0,restyp1))
        nbins = len(m.features[0].bins)
        rowentr = [m[i].entropy() for i in range(nbins)]
        known_rowentr = [1.908535, 0.578325, 1.229918, 1.559581,
                         1.332179, 1.747868, 0.693147, 1.671595,
                         0.796311, 0.796311, 0.0,      1.475076,
                         1.569152, 0.950270, 0.0,      1.560710,
                         1.098612, 1.671595, 0.0,      1.386294,
                         2.642088, 3.091042]
        for (row, known) in zip(rowentr, known_rowentr):
            self.assertAlmostEqual(row, known, places=3)

    def test_mean_stdev(self):
        """Check the mean and entropy of each row"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = self.get_test_mdt(mlib, features=(restyp0,restyp1))
        mean_stdev = [m[i].mean_stdev() for i in range(3)]
        known = [[9.357142, 8.131571],
                 [3.264705, 4.492207],
                 [8.8125,   8.228902]]
        for (row, known) in zip(mean_stdev, known):
            for (a, b) in zip(row, known):
                self.assertAlmostEqual(a, b, places=3)

    def test_normalize(self):
        """Check that normalize works"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = mdt.Table(mlib, features=(restyp0,restyp1))
        # Dimensions must be 1 or 2:
        for dim in (0, 3):
            self.assertRaises(ValueError, m.normalize, dimensions=dim,
                              dx_dy=(1,)*dim, to_zero=False, to_pdf=True)

        # Size of dx_xy must match dimensions
        self.assertRaises(ValueError, m.normalize, dimensions=2, dx_dy=1,
                          to_zero=False, to_pdf=True)

        # A normalized empty MDT should be the prior, i.e. 1/nbin:
        m2 = m.normalize(dimensions=2, dx_dy=(1,1), to_zero=False, to_pdf=True)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            self.assertAlmostEqual(1./484., m2[inds], places=3)
        m2 = m.normalize(dimensions=1, dx_dy=1, to_zero=False, to_pdf=True)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            self.assertAlmostEqual(1./22., m2[inds], places=3)

        # A normalized empty MDT should be empty if to_zero=True
        m2 = m.normalize(dimensions=2, dx_dy=(1,1), to_zero=True, to_pdf=False)
        self.assertEqual(m.pdf, 0)
        self.assertEqual(m2.pdf, 1)
        self.assertMDTsEqual(m, m2, check_pdf=False)

    def test_write_asgl(self):
        """Check that ASGL output works"""
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        m = self.get_test_mdt(mlib, features=(restyp0,restyp1))
        root = 'asgl1-a'
        for dim in (0, 3):
            self.assertRaises(ValueError, m.write_asgl, dimensions=dim,
                              asglroot=root, plots_per_page=8, plot_position=1,
                              every_x_numbered=999, text="test text",
                              x_decimal=0)
        m.write_asgl(asglroot=root, plots_per_page=8, dimensions=1,
                     plot_position=1, every_x_numbered=999, text="test text",
                     x_decimal=0)
        with open(root+'.top', 'r') as f:
            self.assertEqual(len(f.readlines()), 486)
        os.unlink(root+'.top')
        for n in range(1,23):
            fname = '%s.%d' % (root, n)
            with open(fname, 'r') as f:
                self.assertEqual(len(f.readlines()), 22)
            os.unlink(fname)

    def test_bin_type(self):
        """Check building tables with different bin storage types"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        restyp0 = mdt.features.ResidueType(mlib, protein=0)
        restyp1 = mdt.features.ResidueType(mlib, protein=1)
        features = (restyp0, restyp1)
        aln = alignment(env, file='test/data/alignment.ali')
        # bin_type must be a valid MDT BinType
        self.assertRaises(TypeError, mdt.Table, mlib, features=features,
                          bin_type='garbage')
        m1 = mdt.Table(mlib, features=features, bin_type=mdt.Double)
        m1.add_alignment(aln)
        for bin_type in (mdt.Int8, mdt.Int16, mdt.Int32, mdt.UnsignedInt8,
                         mdt.UnsignedInt16, mdt.UnsignedInt32, mdt.Float,
                         mdt.Double):
            m2 = mdt.Table(mlib, features=features, bin_type=bin_type)
            m2.add_alignment(aln)
            m2.write_hdf5('test.hdf5')
            self.assertMDTsEqual(m1, m2)
            m3 = m1.copy(bin_type=bin_type)
            self.assertMDTsEqual(m1, m3)
        os.unlink('test.hdf5')

    def test_copy_bad_bin_type(self):
        """Table.copy() should complain if given a bad bin type"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m = mdt.Table(mlib, features=mdt.features.ResidueType(mlib))
        self.assertRaises(TypeError, m.copy, bin_type='garbage')

    def test_features_to_ifeat(self):
        """Test Table._features_to_ifeat method"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        f = mdt.features.ResidueType(mlib)
        m = mdt.Table(mlib, features=f)
        for arg in [f, [f], (f,)]:
            ifeat = m._features_to_ifeat(arg)
            self.assertEqual(ifeat, [1])
        self.assertRaises(TypeError, m._features_to_ifeat, 'garbage')

    def test_get_splinerange(self):
        """Test _get_splinerange utility function"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        angle = mdt.features.Angle(mlib,
                                   bins=mdt.uniform_bins(288, 0.0, 0.625))
        dih = mdt.features.Dihedral(mlib,
                                    bins=mdt.uniform_bins(1, -500, 1000))
        t = mdt.Table(mlib, features=[angle,dih])
        t = t.reshape(features=[angle,dih], offset=[0,0], shape=[-1,-1])

        periodic, dx, x1, x2 = mdt._get_splinerange(t.features[0])
        self.assertEqual(periodic, 0)
        self.assertAlmostEqual(dx, 0.011, delta=1e-3)
        self.assertAlmostEqual(x1, 0.005, delta=1e-3)
        self.assertAlmostEqual(x2, 3.136, delta=1e-3)

        periodic, dx, x1, x2 = mdt._get_splinerange(t.features[1])
        self.assertEqual(periodic, 1)
        self.assertAlmostEqual(dx, 17.453, delta=1e-3)
        self.assertAlmostEqual(x1, 0.0, delta=1e-3)
        self.assertAlmostEqual(x2, 17.453, delta=1e-3)

    def test_pass_cutoffs(self):
        """Test _pass_cutoffs utility function"""
        class DummySection(object):
            def __init__(self, sum, entropy):
                self._sum = sum
                self._entropy = entropy
            def sum(self):
                return self._sum
            def entropy(self):
                return self._entropy
        class DummyBin(object):
            symbol = 'symbol'
        dummy_table = [DummySection(10.0, 20.0)]
        dummy_bin = DummyBin()
        self.assertEqual(mdt._pass_cutoffs(dummy_table, 0, dummy_bin,
                              density_cutoff=None, entropy_cutoff=None), True)
        self.assertEqual(mdt._pass_cutoffs(dummy_table, 0, dummy_bin,
                              density_cutoff=10.1, entropy_cutoff=None), False)
        self.assertEqual(mdt._pass_cutoffs(dummy_table, 0, dummy_bin,
                              density_cutoff=None, entropy_cutoff=19.9), False)
        self.assertEqual(mdt._pass_cutoffs(dummy_table, 0, dummy_bin,
                              density_cutoff=9.9, entropy_cutoff=20.1), True)

    def test_mdt_witherr(self):
        """Test the calculation of mdt with error"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        mlib.tuple_classes.read('data/dblcls.lib')
        dbldist = mdt.features.TupleDistance(mlib,
                                             bins=mdt.uniform_bins(30, 0, 0.5))
        aln = alignment(env, file='test/data/tiny.ali',
                        align_codes='5fd1')
        m1 = mdt.Table(mlib, features=dbldist)
        m1.add_alignment_witherr(aln, residue_span_range=(-999, -1, 1, 999),
                                 errorscale=0.1)
        self.assertAlmostEqual(m1.sum(), 1475.91, delta=0.01)
        self.assertAlmostEqual(m1[29], 22.00, delta=0.5)
        self.assertAlmostEqual(m1[25], 38.00, delta=0.5)
        self.assertAlmostEqual(m1[0], 0.00, delta=1.0)

    def test_triple_protein_scan(self):
        """Test scan of triples of proteins"""
        env = self.get_environ()
        mlib = self.get_mdt_library()

        a = alignment(env, file='test/data/resol.ali',
                      align_codes=('bin1', 'bin2', 'undef2'))

        f1 = mdt.features.XRayResolution(mlib, mdt.uniform_bins(5, 0, 1.0),
                                         protein=0)
        f2 = mdt.features.XRayResolution(mlib, mdt.uniform_bins(5, 0, 1.0),
                                         protein=1)
        f3 = mdt.features.XRayResolution(mlib, mdt.uniform_bins(5, 0, 1.0),
                                         protein=2)

        # When we have features that cover proteins 0,1,2, triple scan should
        # be forced; different number of samples for symmetric/non-symmetric
        for (sym,size) in ((False, 6.0), (True, 3.0)):
            t = mdt.Table(mlib, features=(f1,f2,f3))
            t.add_alignment(a, symtriples=sym)
            self.assertAlmostEqual(t.sample_size, size, delta=1e-6)

        # When feature only covers some of the proteins, triplet scan does
        # not occur.
        for sym in (True, False):
            t = mdt.Table(mlib, features=f3)
            t.add_alignment(a, symtriples=sym)
            self.assertAlmostEqual(t.sample_size, 3.0, delta=1e-6)

        # Exercise residue pair features in combination with triplet scans
        ri = mdt.features.ResidueIndexDifference(mlib,
                                                 mdt.uniform_bins(5, 0, 1.0),
                                                 protein=2)
        t = mdt.Table(mlib, features=(f1,f2,ri))
        t.add_alignment(a, symtriples=sym)
        self.assertAlmostEqual(t.sample_size, 12, delta=1e-6)

    def test_bond_span_range_disulfide(self):
        """Test bond_span_range argument with disulfides"""
        env = self.get_environ()
        mdl = model(env)
        mdl.read('1HEL.pdb')
        aln = alignment(env)
        aln.append_model(mdl, align_codes='test')
        mlib = self.get_mdt_library()
        mlib.bond_classes.read('data/bndgrp.lib')
        dist = mdt.features.AtomDistance(mlib,
                                         bins=mdt.uniform_bins(60, 0, 0.5))
        # Four disulfide bond in this structure
        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, bond_span_range=(1,1),
                        residue_span_range=(-9999,0,0,9999))

        m2 = mdt.Table(mlib, features=dist)
        m2.add_alignment(aln, bond_span_range=(1,1),
                        residue_span_range=(-9999,0,0,9999), disulfide=True)
        self.assertEqual(m2.sample_size-m.sample_size, 4.0)

        m = mdt.Table(mlib, features=dist)
        m.add_alignment(aln, bond_span_range=(3,3),
                        residue_span_range=(-9999,0,0,9999))

        m2 = mdt.Table(mlib, features=dist)
        m2.add_alignment(aln, bond_span_range=(3,3),
                        residue_span_range=(-9999,0,0,9999), disulfide=True)
        self.assertEqual(m2.sample_size-m.sample_size, 12.0)

if __name__ == '__main__':
    unittest.main()
