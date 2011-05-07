import unittest
import modeller
from mdt_test import MDTTest
import mdt
import mdt.features
import os

class BondLibTests(MDTTest):

    def get_test_mdt(self, mlib, features):
        env = self.get_environ()
        mdl = modeller.model(env)
        mdl.build_sequence('C')

        m = mdt.Table(mlib, features=features)
        a = modeller.alignment(env)
        a.append_model(mdl, atom_files='test', align_codes='test')
        m.add_alignment(a)
        m = m.reshape(features, [0] * len(features), [-1] * len(features))
        return m

    def test_write_bondlib(self):
        """Test the write_bondlib function"""
        mlib = self.get_mdt_library()
        mlib.bond_classes.read('data/bndgrp.lib')

        bondtype = mdt.features.BondType(mlib)
        bondlen = mdt.features.BondLength(mlib,
                                      bins=mdt.uniform_bins(200, 1.0, 0.005))
        m = self.get_test_mdt(mlib, [bondtype, bondlen])
        mdt.write_bondlib(open('test.out', 'w'), m)

        # Make sure that valid Python code was produced
        code = compile(open('test.out').read(), 'test.out', 'exec')
        os.unlink('test.out')

    def test_write_anglelib(self):
        """Test the write_anglelib function"""
        mlib = self.get_mdt_library()
        mlib.angle_classes.read('data/anggrp.lib')

        angtype = mdt.features.AngleType(mlib)
        ang = mdt.features.Angle(mlib, bins=mdt.uniform_bins(180, 0.0, 2.0))
        m = self.get_test_mdt(mlib, [angtype, ang])
        mdt.write_anglelib(open('test.out', 'w'), m)

        # Make sure that valid Python code was produced
        code = compile(open('test.out').read(), 'test.out', 'exec')
        os.unlink('test.out')

    def test_write_improperlib(self):
        """Test the write_improperlib function"""
        mlib = self.get_mdt_library()
        mlib.dihedral_classes.read('data/impgrp.lib')

        angtype = mdt.features.DihedralType(mlib)
        ang = mdt.features.Dihedral(mlib, bins=mdt.uniform_bins(180, 0.0, 2.0))
        m = self.get_test_mdt(mlib, [angtype, ang])
        mdt.write_improperlib(open('test.out', 'w'), m)

        # Make sure that valid Python code was produced
        code = compile(open('test.out').read(), 'test.out', 'exec')
        os.unlink('test.out')

    def test_write_splinelib(self):
        """Test the write_splinelib function"""
        mlib = self.get_mdt_library()

        restype = mdt.features.ResidueType(mlib)
        ang = mdt.features.Chi1Dihedral(mlib,
                                        bins=mdt.uniform_bins(180, -180.0, 4.0))
        m = self.get_test_mdt(mlib, [restype, ang])
        mdt.write_splinelib(open('test.out', 'w'), m, 'chi1')

        # Make sure that valid Python code was produced
        code = compile(open('test.out').read(), 'test.out', 'exec')
        os.unlink('test.out')

    def test_write_2dsplinelib(self):
        """Test the write_2dsplinelib function"""
        mlib = self.get_mdt_library()

        restype = mdt.features.ResidueType(mlib)
        phi = mdt.features.PhiDihedral(mlib,
                                       bins=mdt.uniform_bins(18, -180.0, 40.0))
        psi = mdt.features.PsiDihedral(mlib,
                                       bins=mdt.uniform_bins(18, -180.0, 40.0))
        m = self.get_test_mdt(mlib, [restype, phi, psi])
        mdt.write_2dsplinelib(open('test.out', 'w'), m)

        # Make sure that valid Python code was produced
        code = compile(open('test.out').read(), 'test.out', 'exec')
        os.unlink('test.out')

if __name__ == '__main__':
    unittest.main()
