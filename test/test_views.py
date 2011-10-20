import unittest
from mdt_test import MDTTest
import mdt
import mdt.features
import weakref
import os
import sys

class ViewTests(MDTTest):

    def get_mdt_and_view(self):
        mlib = self.get_mdt_library()
        restyp = mdt.features.ResidueType(mlib)
        chi1 = mdt.features.Chi1Dihedral(mlib,
                                         mdt.uniform_bins(36, -180, 10))
        m = self.get_test_mdt(mlib, features=(restyp,chi1))
        try:
            a = m.get_array_view()
        except NotImplementedError:
            sys.stderr.write("No NumPy support; skipping test; ")
            a = None
        return m, a, restyp

    def test_ref_count(self):
        """Test refcounting of views"""
        m, a, restyp = self.get_mdt_and_view()
        if a is None:
            return
        b = a.view()
        mr = weakref.ref(m)
        ar = weakref.ref(a)
        br = weakref.ref(b)
        # m should stay around until all views are destroyed
        del m
        self.assertNotEqual(mr(), None)
        del a
        self.assertNotEqual(mr(), None)
        del b
        self.assertEqual(mr(), None)

    def test_view_shared_data(self):
        """Check that views share data with tables"""
        m, a, restyp = self.get_mdt_and_view()
        if a is None:
            return
        self.assertEqual(m.shape, a.shape)
        m[3][2] = 1000.
        self.assertInTolerance(a[3][2], 1000., 1e-2)
        # numpy array view should be writeable
        self.assertEqual(a.flags.writeable, True)
        a[4][1] = 2000.
        self.assertInTolerance(m[4][1], 2000., 1e-2)

    def test_view_no_own(self):
        """Check that views do not own data"""
        m, a, restyp = self.get_mdt_and_view()
        if a is None:
            return
        self.assertEqual(a.flags.owndata, False)
        self.assertRaises(ValueError, a.resize, 2, 2)

    def test_view_locks(self):
        """Check that an active view locks the table"""
        m, a, restyp = self.get_mdt_and_view()
        if a is None:
            return
        m.write('test.mdt')
        m.write_hdf5('test.hdf5')
        # Cannot destroy data while a view exists
        self.assertRaises(ValueError, m.read, 'test.mdt')
        self.assertRaises(ValueError, m.read_hdf5, 'test.hdf5')
        self.assertRaises(ValueError, m.make, features=restyp)
        # OK once the view is destroyed
        del a
        m.read('test.mdt')
        m.read_hdf5('test.hdf5')
        m.make(features=restyp)
        os.unlink('test.mdt')
        os.unlink('test.hdf5')

if __name__ == '__main__':
    unittest.main()
