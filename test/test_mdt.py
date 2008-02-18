import unittest
from modeller import *
from mdt_test import MDTTest
import mdt
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

    def test_bad_bins(self):
        """Check that bad bin files raise an error"""
        env = self.get_environ()
        mlib = mdt.Library(env, 'test/data/bad.bin')
        self.assertRaises(mdt.MDTError, mdt.Table, mlib, features=3)
        self.assertRaises(mdt.MDTError, mdt.Table, mlib, features=30)

    def test_mdt_formats(self):
        """Make sure we can read and write MDT files"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m = mdt.Table(mlib, features=(1,2))
        aln = alignment(env)
        aln.append_sequence('AFVVTDNCIKCKYTDCVEVXPVDCFYEG')
        aln.append_sequence('CVEVCPVDCFYEGAFVVTDNCIKCKYTX')
        m.add_alignment(aln)
        m.write('test.mdt')
        m2 = m.copy()
        self.assertMDTsEqual(m, m2)
        m2 = mdt.Table(mlib, file='test.mdt')
        self.assertMDTsEqual(m, m2)
        m.write_hdf5('test.hdf5')
        self.assertRaises(mdt.MDTError, m.write_hdf5,
                          '/does/not/exist/foo.hdf5')
        m2 = mdt.Table(mlib, file='test.hdf5')
        self.assertMDTsEqual(m, m2)
        m2 = mdt.Table(mlib)
        # Trying to read an HDF5 file in text format should fail gracefully:
        self.assertRaises(mdt.MDTError, m2.read, 'test.hdf5')
        # Same problem should occur starting with a non-empty MDT:
        m2 = mdt.Table(mlib, features=(1,2))
        self.assertRaises(mdt.MDTError, m2.read, 'test.hdf5')
        os.unlink('test.hdf5')

    def test_bin_info(self):
        """Test query of bin symbol and range"""
        m = self.get_test_mdt(features=66)
        for (n, bin) in enumerate(m.features[0].bins):
            self.assertEqual(bin.range, (n, n+1))
        symbols = [bin.symbol for bin in m.features[0].bins]
        self.assertEqual(symbols, ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY',
                                   'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN',
                                   'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL',
                                   'TRP', 'TYR', 'GAP', 'u'])
        m2 = m.reshape(features=66, offset=1, shape=20)
        for (n, bin) in enumerate(m2.features[0].bins):
            self.assertEqual(bin.symbol, m.features[0].bins[n+1].symbol)

    def test_set(self):
        """Test set of MDT data"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m = mdt.Table(mlib, features=(1,2))
        # Make sure that index checks work:
        self.assertRaises(ValueError, m.__setitem__, [0,22,10], 0.0)
        self.assertRaises(IndexError, m.__setitem__, [0,22], 0.0)
        self.assertRaises(IndexError, m.__setitem__, [22,0], 0.0)
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
        m = mdt.Table(mlib, features=(1,109))
        aln = alignment(env)
        self.assertRaises(mdt.MDTError, m.add_alignment, aln)

    def test_feature_resind_diff(self):
        """Test the residue index difference feature"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        aln = alignment(env, file='test/data/alignment.ali', align_codes='5fd1')
        m = mdt.Table(mlib, features=51)
        m.add_alignment(aln, residue_span_range=(-999, -2, 2, 999))
        # span range should result in 0, +/- 1 bins being zero:
        self.assertEqual(m[9], 0.)
        self.assertEqual(m[10], 0.)
        self.assertEqual(m[11], 0.)
        # other bins should be symmetrically distributed:
        for i in range(9):
            self.assertEqual(m[i], m[-2 - i])

    def test_1d_mdt(self):
        """Make sure a 1D MDT matches known residue data"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m = mdt.Table(mlib, features=1)
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
        bins = [ source.index(1, 0, n, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
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
        mlib = mdt.Library(env, 'test/data/dist.bin')
        mlib.tuple_classes.read('data/dblcls.lib')
        # All residue-residue pairs should be out of range, so this MDT should
        # end up empty:
        m = mdt.Table(mlib, features=103)
        m.add_alignment(aln, residue_span_range=(-999, -999, 999, 999))
        self.assertEqual(m.sum(), 0.0)
        m1 = mdt.Table(mlib, features=103)
        m2 = mdt.Table(mlib, features=103)
        # When excluding only intra-residue interactions, the short-range
        # bins should differ but the long-range behavior should be the same:
        m1.add_alignment(aln, residue_span_range=(-999, 0, 0, 999))
        m2.add_alignment(aln, residue_span_range=(-999, -1, 1, 999))
        self.assertEqual(m1.shape, m2.shape)
        for i in range(1, 5):
            self.assertNotEqual(m1[i], m2[i])
        for i in range(38, 49):
            self.assertEqual(m1[i], m2[i])

    def test_sources(self):
        """Make sure that alignments and models both work as sources"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m1 = mdt.Table(mlib, features=82)
        m2 = mdt.Table(mlib, features=82)
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
        self.assertRaises(ValueError, self.get_test_mdt, features=(17,82))

    def test_integrate(self):
        """Make sure MDT integration works"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m1 = mdt.Table(mlib, features=1)
        m2 = mdt.Table(mlib, features=2)
        mboth = mdt.Table(mlib, features=(1,2))
        seq1 = "AFVVTDNCIK"
        seq2 = "DCVEVCPVDC"
        aln = alignment(env)
        aln.append_sequence(seq1)
        aln.append_sequence(seq2)
        for m in (m1, m2, mboth):
            m.add_alignment(aln)

        # Number of features must be correct:
        for features in ((), (1,2), (1,2,3)):
            self.assertRaises(ValueError, mboth.integrate, features=features)
        # Features must exist in input MDT:
        self.assertRaises(ValueError, mboth.integrate, features=3)
        m1int = mboth.integrate(features=1)
        self.assertMDTsEqual(m1, m1int)
        m2int = mboth.integrate(features=2)
        self.assertMDTsEqual(m2, m2int)

    def test_entropy_hx(self):
        """Check for expected dependent entropy value"""
        m = self.get_test_mdt(features=(1,3))
        # Check against ENTROPY_MDT_HX value for this system in MDT:
        self.assertAlmostEqual(m.entropy_hx(), 2.7048, places=3)

    def test_exp_transform(self):
        """Check for correctness of exp transform"""
        m = self.get_test_mdt(features=(1,7))
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
        m = self.get_test_mdt(features=(1,7))
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
        m = self.get_test_mdt(features=(1,3))
        offset = -1.3
        multiplier = 0.8
        m2 = m.linear_transform(offset, multiplier)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            y = offset + multiplier * m[inds]
            self.assertAlmostEqual(y, m2[inds], places=3)

    def test_inverse_transform(self):
        """Check for correctness of inverse transform"""
        m = self.get_test_mdt(features=(1,7))
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
        m = self.get_test_mdt(features=(1,2))
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
        m = self.get_test_mdt(features=(1,2))
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
        m = mdt.Table(mlib, features=(1,2))
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
        m = self.get_test_mdt(features=(1,2))
        m = m.reshape(features=(1,2), offset=(4,2), shape=(11,10))
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
        m = self.get_test_mdt(features=(1,2))
        m2 = self.get_test_mdt(features=(2,1))
        for features in ((3,1), (2,1,1)):
            self.assertRaises(ValueError, m.reshape, features=features,
                              offset=m.offset, shape=(22,22))
        m3 = m.reshape(features=(2,1), offset=m.offset, shape=(22,22))
        self.assertMDTsEqual(m2, m3)
        m3 = m.reshape(features=(1,2), offset=(0,0), shape=m.shape)
        self.assertMDTsEqual(m, m3)
        m3 = m.reshape(features=(1,2), offset=(4,2), shape=(11,10))
        self.assertEqual(m3.shape, (11,10))
        self.assertEqual(m3.offset, (4,2))
        inds = []
        while self.roll_inds(inds, m3.shape, m3.offset):
            self.assertAlmostEqual(m[inds], m3[inds], places=3)

    def test_sum(self):
        """Check that sum of each row sums to that of the whole table"""
        m = self.get_test_mdt(features=(1,2))
        sum = m.sum()
        rowsum = 0.
        for i in range(len(m.features[0].bins)):
            rowsum += m[i].sum()
        self.assertAlmostEqual(sum, rowsum, places=3)

    def test_entropy(self):
        """Check the entropy of each row"""
        m = self.get_test_mdt(features=(1,2))
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
        m = self.get_test_mdt(features=(1,2))
        mean_stdev = [m[i].mean_stdev() for i in range(3)]
        known = [[9.357142, 8.131571],
                 [3.264705, 4.492207],
                 [8.8125,   8.228902]]
        for (row, known) in zip(mean_stdev, known):
            for (a, b) in zip(row, known):
                self.assertAlmostEqual(a, b, places=3)

    def test_normalize(self):
        """Check that normalize works"""
        m = self.get_test_mdt(features=(1,2))
        mlib = self.get_mdt_library()
        m = mdt.Table(mlib, features=(1,2))
        # Dimensions must be 1 or 2:
        for dim in (0, 3):
            self.assertRaises(ValueError, m.normalize, dimensions=dim,
                              dx_dy=(1,)*dim, to_zero=False, to_pdf=True)
                               
        # A normalized empty MDT should be the prior, i.e. 1/nbin:
        m2 = m.normalize(dimensions=2, dx_dy=(1,1), to_zero=False, to_pdf=True)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            self.assertAlmostEqual(1./484., m2[inds], places=3)
        m2 = m.normalize(dimensions=1, dx_dy=1, to_zero=False, to_pdf=True)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            self.assertAlmostEqual(1./22., m2[inds], places=3)

    def test_write_asgl(self):
        """Check that ASGL output works"""
        m = self.get_test_mdt(features=(1,2))
        root = 'asgl1-a'
        for dim in (0, 3):
            self.assertRaises(ValueError, m.write_asgl, dimensions=dim,
                              asglroot=root, plots_per_page=8, plot_position=1,
                              every_x_numbered=999, text="test text",
                              x_decimal=0)
        m.write_asgl(asglroot=root, plots_per_page=8, dimensions=1,
                     plot_position=1, every_x_numbered=999, text="test text",
                     x_decimal=0)
        f = file(root+'.top', 'r')
        self.assertEqual(len(f.readlines()), 486)
        del f
        os.unlink(root+'.top')
        for n in range(1,23):
            fname = '%s.%d' % (root, n)
            f = file(fname, 'r')
            self.assertEqual(len(f.readlines()), 22)
            del f
            os.unlink(fname)

    def test_bin_type(self):
        """Check building tables with different bin storage types"""
        env = self.get_environ()
        features = (1, 2)
        mlib = self.get_mdt_library()
        aln = alignment(env, file='test/data/alignment.ali')
        m1 = mdt.Table(mlib, features=features, bin_type=mdt.Double)
        m1.add_alignment(aln)
        for bin_type in (mdt.Int8, mdt.Int16, mdt.Int32, mdt.UnsignedInt8,
                         mdt.UnsignedInt16, mdt.UnsignedInt32, mdt.Float,
                         mdt.Double):
            m2 = mdt.Table(mlib, features=features, bin_type=bin_type)
            m2.add_alignment(aln)
            self.assertMDTsEqual(m1, m2)
            m3 = m1.copy(bin_type=bin_type)
            self.assertMDTsEqual(m1, m3)

if __name__ == '__main__':
    unittest.main()
