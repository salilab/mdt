import unittest
from mdt_test import MDTTest
import mdt
import mdt.features
import modeller
import os

class BondSeparationFeatureTests(MDTTest):

    def build_mdt_from_sequence(self, mlib, features, seq, **keys):
        """Build a simple test MDT for a given sequence"""
        env = self.get_environ()
        mdl = modeller.model(env)
        mdl.build_sequence(seq)

        m = mdt.Table(mlib, features=features)
        a = modeller.alignment(env)
        a.append_model(mdl, atom_files='test', align_codes='test')
        m.add_alignment(a, **keys)
        return m

    def get_all_libraries(self):
        """Create an MDT library and read in bond and atom classes"""
        mlib = self.get_mdt_library()
        mlib.atom_classes.read('${LIB}/atmcls-mf.lib')
        mlib.bond_classes.read('data/bndgrp.lib')
        return mlib

    def test_chain_span(self):
        """Atom pairs spanning chains should not be connected"""
        mlib = self.get_all_libraries()
        bsep = mdt.features.AtomBondSeparation(mlib,
                                        bins=mdt.uniform_bins(20, 0, 1.0))
        m = self.build_mdt_from_sequence(mlib, bsep, 'C/C',
                                         residue_span_range=(-999,-1,1,999))
        # All atom pairs should be undefined, since chains are not connected
        self.assertEqual(m.sample_size, 49.0)
        self.assertEqual(m[-1], 49.0)

    def test_unknown_atom_type(self):
        """Atoms of unknown type should count as undefined"""
        mlib = self.get_all_libraries()
        bsep = mdt.features.AtomBondSeparation(mlib,
                                 bins=mdt.uniform_bins(10, -1000, 1000.0))
        m = self.build_mdt_from_sequence(mlib, bsep, 'C',
                                         residue_span_range=(0,0,0,0))
        # The OXT atom is not referenced in the bond library, so bond separation
        # to any of the other 6 atoms should count as undefined (even though
        # the raw "distance" of -1 would otherwise fall in the first bin)
        self.assertEqual(m[-1], 6.0)

    def test_unconnected_internal(self):
        """Separation between unconnected internal atoms is undefined"""
        mlib = self.get_mdt_library()
        mlib.atom_classes.read('${LIB}/atmcls-mf.lib')
        # Make dummy bond class library that does not connect C and CA atoms
        open('dummy.lib', 'w').write("""
BNDGRP 'ALA:CB:CA'
  BOND 'ALA' 'CB' 'CA'
BNDGRP 'ALA:N:CA'
  BOND 'ALA' 'N' 'CA'
BNDGRP 'ALA:O:C'
  BOND 'ALA' 'O' 'C'
""")
        mlib.bond_classes.read('dummy.lib')

        attyp1 = mdt.features.AtomType(mlib)
        attyp2 = mdt.features.AtomType(mlib, pos2=True)
        bsep = mdt.features.AtomBondSeparation(mlib,
                                        bins=mdt.uniform_bins(20, 0, 1.0))
        m = self.build_mdt_from_sequence(mlib, [attyp1, attyp2, bsep],
                        'A', residue_span_range=(0,0,0,0))
        atom_types = {}
        for n, b in enumerate(m.features[0].bins):
            atom_types[b.symbol] = n

        def assertBondSep(at1, at2, numbond, sep):
            bins = [b for b in m[atom_types[at1]][atom_types[at2]]]
            self.assertEqual(sum(bins), numbond)
            self.assertEqual(bins[sep], numbond)

        # All bonds that pass through C-CA should be undefined, others are OK
        assertBondSep('AN', 'ACA', numbond=1, sep=1)
        assertBondSep('AN', 'AC', numbond=1, sep=-1)
        assertBondSep('AN', 'AO', numbond=1, sep=-1)
        assertBondSep('AN', 'ACB', numbond=1, sep=2)
        os.unlink('dummy.lib')

    def test_unconnected_external(self):
        """Separation between unconnected external atoms is undefined"""
        mlib = self.get_mdt_library()
        mlib.atom_classes.read('${LIB}/atmcls-mf.lib')
        # Make dummy bond class library that does not connect C and CA atoms
        open('dummy.lib', 'w').write("""
BNDGRP 'ALA:CB:CA'
  BOND 'ALA' 'CB' 'CA'
BNDGRP 'ALA:N:CA'
  BOND 'ALA' 'N' 'CA'
BNDGRP 'ALA:O:C'
  BOND 'ALA' 'O' 'C'
""")
        mlib.bond_classes.read('dummy.lib')

        attyp1 = mdt.features.AtomType(mlib)
        attyp2 = mdt.features.AtomType(mlib, pos2=True)
        bsep = mdt.features.AtomBondSeparation(mlib,
                                        bins=mdt.uniform_bins(20, 0, 1.0))
        m = self.build_mdt_from_sequence(mlib, [attyp1, attyp2, bsep],
                        'AA', residue_span_range=(-1,-1,1,1))
        atom_types = {}
        for n, b in enumerate(m.features[0].bins):
            atom_types[b.symbol] = n

        def assertBondSep(at1, at2, numbond, sep):
            bins = [b for b in m[atom_types[at1]][atom_types[at2]]]
            self.assertEqual(sum(bins), numbond)
            self.assertEqual(bins[sep], numbond)

        # All bonds that pass through C-CA should be undefined, others are OK
        assertBondSep('AN', 'ACA', numbond=1, sep=-1)
        assertBondSep('AN', 'AC', numbond=1, sep=-1)
        assertBondSep('AN', 'AO', numbond=1, sep=-1)
        assertBondSep('AN', 'ACB', numbond=1, sep=-1)

        assertBondSep('AC', 'AN', numbond=1, sep=1)
        assertBondSep('AO', 'ACB', numbond=1, sep=4)
        os.unlink('dummy.lib')

    def test_internal(self):
        """Check bond separation feature within residues"""
        mlib = self.get_all_libraries()
        attyp1 = mdt.features.AtomType(mlib)
        attyp2 = mdt.features.AtomType(mlib, pos2=True)
        bsep = mdt.features.AtomBondSeparation(mlib,
                                        bins=mdt.uniform_bins(20, 0, 1.0))
        m = self.build_mdt_from_sequence(mlib, [attyp1, attyp2, bsep],
                        'ARNDCQEHFWY', residue_span_range=(0,0,0,0))
        atom_types = {}
        for n, b in enumerate(m.features[0].bins):
            atom_types[b.symbol] = n

        def assertBondSep(at1, at2, numbond, sep):
            bins = [b for b in m[atom_types[at1]][atom_types[at2]]]
            self.assertEqual(sum(bins), numbond)
            self.assertEqual(bins[sep], numbond)

        # Check known ALA bond separations
        assertBondSep('AN', 'ACA', numbond=1, sep=1)
        assertBondSep('AN', 'AC', numbond=1, sep=2)
        assertBondSep('AN', 'AO', numbond=1, sep=3)
        assertBondSep('AN', 'ACB', numbond=1, sep=2)

        # Check known ARG bond separations
        assertBondSep('RCA', 'RCZ', numbond=1, sep=5)
        assertBondSep('RCA', 'RNH', numbond=2, sep=6)

        # Check known HIS bond separations
        assertBondSep('HCD2', 'HC', numbond=1, sep=4)
        assertBondSep('HNE2', 'HC', numbond=1, sep=5)
        assertBondSep('HCE1', 'HC', numbond=1, sep=5)
        assertBondSep('HND1', 'HC', numbond=1, sep=4)

        # Check known PHE bond separations
        assertBondSep('FN', 'FCZ', numbond=1, sep=6)

        # Check known TRP bond separations
        assertBondSep('WN', 'WCE3', numbond=1, sep=5)
        assertBondSep('WN', 'WNE1', numbond=1, sep=5)
        assertBondSep('WN', 'WCZ2', numbond=1, sep=6)
        assertBondSep('WN', 'WCZ3', numbond=1, sep=6)

        # Check known TYR bond separations
        assertBondSep('YN', 'YOH', numbond=1, sep=7)

    def test_external(self):
        """Check bond separation feature between residues"""
        mlib = self.get_all_libraries()
        attyp1 = mdt.features.AtomType(mlib)
        attyp2 = mdt.features.AtomType(mlib, pos2=True)
        bsep = mdt.features.AtomBondSeparation(mlib,
                                        bins=mdt.uniform_bins(20, 0, 1.0))
        m = self.build_mdt_from_sequence(mlib, [attyp1, attyp2, bsep],
                        'ARN', residue_span_range=(-999,-1,1,999))
        atom_types = {}
        for n, b in enumerate(m.features[0].bins):
            atom_types[b.symbol] = n

        def assertBondSep(at1, at2, numbond, sep):
            bins = [b for b in m[atom_types[at1]][atom_types[at2]]]
            self.assertEqual(sum(bins), numbond)
            self.assertEqual(bins[sep], numbond)

        # Check known bond separations between adjacent residues
        assertBondSep('AN', 'RN', numbond=1, sep=3)
        assertBondSep('AC', 'RN', numbond=1, sep=1)

        # Check known bond separations between non-adjacent residues
        assertBondSep('AN', 'NN', numbond=1, sep=6)
        assertBondSep('AN', 'NC', numbond=1, sep=8)
        assertBondSep('ACB', 'NCG', numbond=1, sep=9)

    def test_disulfide(self):
        """Test handling of disulfide bonds"""
        mlib = self.get_all_libraries()
        bsep = mdt.features.AtomBondSeparation(mlib,
                                        bins=mdt.uniform_bins(20, 0, 1.0))
        bsep_ss = mdt.features.AtomBondSeparation(mlib,
                                        bins=mdt.uniform_bins(20, 0, 1.0),
                                        disulfide=True)
        env = self.get_environ()
        mdl = modeller.model(env)
        mdl.build_sequence('CC')
        # When SG-SG distance is small enough, an extra bond
        # (separation feature = 1) should be detected, but only with
        # disulfide=True
        for (dist, num) in [(2.6, 11.0), (2.4, 12.0)]:
            sg1 = mdl.residues[0].atoms['SG']
            sg2 = mdl.residues[1].atoms['SG']
            sg1.x = sg1.y = sg1.z = 0.
            sg2.x = sg2.y = 0.
            sg2.z = dist
            a = modeller.alignment(env)
            a.append_model(mdl, atom_files='test', align_codes='test')
            m = mdt.Table(mlib, features=bsep)
            m.add_alignment(a, residue_span_range=(-999,0,0,999))
            self.assertEqual(m[1], 11.0)
            m2 = mdt.Table(mlib, features=bsep_ss)
            m2.add_alignment(a, residue_span_range=(-999,0,0,999))
            self.assertEqual(m2[1], num)

if __name__ == '__main__':
    unittest.main()
