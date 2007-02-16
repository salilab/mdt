import unittest
from modeller import *
from modeller.test import ModellerTest
import mdt
import math

class MDTTests(ModellerTest):

    restyps = "ACDEFGHIKLMNPQRSTVWY"

    def __mdt_restyp(self, a):
        """Convert a one-letter residue code into the MDT bin number"""
        ind = self.restyps.find(a)
        # Map to 'undefined' bin:
        if ind == -1:
            ind = 21
        return ind

    def get_mdt_library(self):
        """Read in MDT library and bin definitions"""
        env = self.get_environ()
        return mdt.mdt_library(env, '${LIB}/mdt.ini', '${LIB}/mdt.bin')

    def get_test_mdt(self, features):
        """Build a simple test MDT"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m = mdt.mdt(mlib, features=features)
        aln = alignment(env, file='test/data/alignment.ali')
        m.add_alignment(aln)
        return m

    def roll_inds(self, inds, shape, offset):
        """Return the next set of indices within an array of given shape"""
        i = len(shape) - 1
        if len(inds) == 0:
            inds.extend(offset)
            return sum(shape) != 0
        while (i >= 0):
            if inds[i] < shape[i] - 1:
                inds[i] += 1
                return True
            elif i == 0:
                return False
            else:
                inds[i] = offset[i]
                i -= 1
        return False

    def assertMDTDataEqual(self, mdt1, mdt2):
        """Make sure that the actual data points in two MDTs are equal"""
        self.assertEqual(mdt1.shape, mdt2.shape)
        shape = mdt1.shape
        inds = []
        npoints = 0
        while self.roll_inds(inds, mdt1.shape, mdt1.offset):
            npoints += 1
            self.assertAlmostEqual(mdt1[inds], mdt2[inds], places=3)
        self.assertEqual(npoints, reduce(lambda x,y: x*y, shape))

    def assertMDTsEqual(self, mdt1, mdt2):
        """Make sure that two MDTs are equal"""
        self.assertEqual(len(mdt1.features), len(mdt2.features))
        self.assertEqual(mdt1.n_proteins, mdt2.n_proteins)
        self.assertEqual(mdt1.n_protein_pairs, mdt2.n_protein_pairs)
        self.assertEqual(mdt1.sample_size, mdt2.sample_size)
        self.assertEqual(mdt1.pdf, mdt2.pdf)
        for (f1, f2) in zip(mdt1.features, mdt2.features):
            self.assertEqual(len(f1.bins), len(f2.bins))
            self.assertEqual(f1.ifeat, f2.ifeat)
        self.assertMDTDataEqual(mdt1, mdt2)

    def test_mdt_formats(self):
        """Make sure we can read and write MDT files"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m = mdt.mdt(mlib, features=(1,2))
        aln = alignment(env)
        aln.append_sequence('AFVVTDNCIKCKYTDCVEVXPVDCFYEG')
        aln.append_sequence('CVEVCPVDCFYEGAFVVTDNCIKCKYTX')
        m.add_alignment(aln)
        m.write('test.mdt')
        m2 = m.copy()
        self.assertMDTsEqual(m, m2)

        m2 = mdt.mdt(mlib, file='test.mdt')
        self.assertMDTsEqual(m, m2)

    def test_1d_mdt(self):
        """Make sure a 1D MDT matches known residue data"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m = mdt.mdt(mlib, features=1)
        aln = alignment(env)
        seq = "AFVVTDNCIKXCKYTDCVEVCPVDCFYEG"
        seq_types = [self.__mdt_restyp(a) for a in seq]
        known_dist = [seq_types.count(n) for n in range(22)]
        aln.append_sequence(seq)
        m.add_alignment(aln)
        for i in range(22):
            self.assertEqual(known_dist[i], int(m[i]))

    def test_integrate(self):
        """Make sure MDT integration works"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m1 = mdt.mdt(mlib, features=1)
        m2 = mdt.mdt(mlib, features=2)
        mboth = mdt.mdt(mlib, features=(1,2))
        seq1 = "AFVVTDNCIK"
        seq2 = "DCVEVCPVDC"
        aln = alignment(env)
        aln.append_sequence(seq1)
        aln.append_sequence(seq2)
        for m in (m1, m2, mboth):
            m.add_alignment(aln)

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

    def test_super_smooth(self):
        """Check for correctness of super smoothing"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m = mdt.mdt(mlib, features=(1,2))
        # A super-smoothed empty MDT should be the prior, i.e. 1/nbin:
        m2 = m.super_smooth(1.0, True)
        inds = []
        while self.roll_inds(inds, m.shape, m.offset):
            self.assertAlmostEqual(1./22., m2[inds], places=3)

    def test_smooth(self):
        """Check for correctness of smoothing"""
        env = self.get_environ()
        mlib = self.get_mdt_library()
        m = mdt.mdt(mlib, features=(1,2))
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

    def test_reshape(self):
        """Check that reshaping works correctly"""
        m = self.get_test_mdt(features=(1,2))
        m2 = self.get_test_mdt(features=(2,1))
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

    def test_normalize(self):
        """Check that normalize works"""
        m = self.get_test_mdt(features=(1,2))
        mlib = self.get_mdt_library()
        m = mdt.mdt(mlib, features=(1,2))
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

if __name__ == '__main__':
    unittest.main()
