import unittest
import modeller

class ModellerTest(unittest.TestCase):
    _env = None

    def get_environ(self):
        """Set up the Modeller environment, and keep a reference to it so that
           we only need to do it once for all tests"""
        if not self._env:
            env = modeller.environ()
            env.io.atom_files_directory = 'test/data'
            env.edat.dynamic_sphere = False
            env.libs.topology.read(file='${LIB}/top_heav.lib')
            env.libs.parameters.read(file='${LIB}/par.lib')
            ModellerTest._env = env
        return self._env

    def assertModelsEqual(self, mdl1, mdl2):
        """Check to make sure that models |mdl1| and |mdl2| are the same"""
        self.assertEqual(len(mdl1.atoms), len(mdl2.atoms))
        for (at1, at2) in zip(mdl1.atoms, mdl2.atoms):
            self.assertAtomsEqual(at1, at2)

    def assertAtomsEqual(self, at1, at2):
        """Check to make sure atoms |at1| and |at2| are the same"""
        self.assertEqual(at1.name, at2.name)
        self.assertAlmostEqual(at1.x, at2.x, places=3)

    def assertDistance(self, o1, o2, mean, tolerance):
        """Check to make sure the distance between o1 and o2 is as expected"""
        distsq = (o1.x-o2.x)**2 + (o1.y-o2.y)**2 + (o1.z-o2.z)**2
        dist = distsq ** 0.5
        msg = "Distance between %s and %s is %f - expected %f" % (o1, o2,
                                                                  dist, mean)
        self.assertInTolerance(dist, mean, tolerance, msg)

    def assertAlignmentsEqual(self, refaln, aln):
        """Check to make sure that alignments |refaln| and |aln| are the same"""
        self.assertEqual(len(refaln), len(aln),
                         "Inconsistent number of sequences (%d vs. %d)" \
                         % (len(refaln), len(aln)))
        for (seq1, seq2) in zip(refaln, aln):
            self.assertAlnSequencesEqual(seq1, seq2)
         
    def assertAlnSequencesEqual(self, seq1, seq2):
        """Check to make sure that sequences |seq1| and |seq2| are the same"""
        self.assertEqual(len(seq1), len(seq2),
                         "Inconsistent number of residues (%d vs. %d)" \
                         % (len(seq1), len(seq2)))
        for (res1, res2) in zip(seq1.residues, seq2.residues):
            self.assertResiduesEqual(res1, res2)

    def assertResiduesEqual(self, res1, res2):
        """Check to make sure that residues |res1| and |res2| are the same"""
        self.assertEqual(res1.name, res2.name,
                         "Inconsistent residue names (%s vs. %s)" \
                         % (res1.name, res2.name))

    def assertInTolerance(self, num1, num2, tolerance, msg=None):
        """Assert that the difference between num1 and num2 is less than
           tolerance"""
        diff = abs(num1 - num2)
        if msg is None:
            msg = "%f != %f within %g" % (num1, num2, tolerance)
        self.assert_(diff < tolerance, msg)
