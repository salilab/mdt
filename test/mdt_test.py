from modeller.test import ModellerTest
import modeller
import mdt
import sys

if sys.version_info[0] >= 3:
    from functools import reduce


class MDTTest(ModellerTest):

    def get_mdt_library(self, **vars):
        """Read in MDT library and bin definitions"""
        env = self.get_environ()
        return mdt.Library(env, **vars)

    def get_test_mdt(self, mlib, features):
        """Build a simple test MDT"""
        env = self.get_environ()
        m = mdt.Table(mlib, features=features)
        aln = modeller.alignment(env, file='test/data/alignment.ali')
        m.add_alignment(aln)
        return m

    def roll_inds(self, inds, shape, offset):
        """Return the next set of indices within an array of given shape"""
        i = len(shape) - 1
        if len(inds) == 0:
            inds.extend(offset)
            return sum(shape) != 0
        while (i >= 0):
            if inds[i] < offset[i] + shape[i] - 1:
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
        self.assertEqual(npoints, reduce(lambda x, y: x*y, shape))

    def assertMDTsEqual(self, mdt1, mdt2, check_pdf=True):
        """Make sure that two MDTs are equal"""
        self.assertEqual(len(mdt1.features), len(mdt2.features))
        self.assertEqual(mdt1.n_proteins, mdt2.n_proteins)
        self.assertEqual(mdt1.n_protein_pairs, mdt2.n_protein_pairs)
        self.assertEqual(mdt1.sample_size, mdt2.sample_size)
        if check_pdf:
            self.assertEqual(mdt1.pdf, mdt2.pdf)
        for (f1, f2) in zip(mdt1.features, mdt2.features):
            self.assertEqual(len(f1.bins), len(f2.bins))
            self.assertEqual(f1.ifeat, f2.ifeat)
        self.assertMDTDataEqual(mdt1, mdt2)

    def assertSectionNormalized(self, section):
        """Make sure that a table section is normalized"""
        sum = section.sum()
        self.assertAlmostEqual(sum, 1.0, delta=1e-5)
