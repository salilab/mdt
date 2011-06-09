import unittest, sys, os, re
import glob

# Our bundled version of coverage only works with Python 2
if sys.version_info[0] < 3:
    import coverage
else:
    coverage = None

try:
    # Custom TestResult class that does not duplicate test name and docstring
    # when using unittest2
    from unittest import TextTestResult

    class _TestRunner(unittest.TextTestRunner):
        class _TestResult(TextTestResult):
            def getDescription(self, test):
                doc_first_line = test.shortDescription()
                if self.descriptions and doc_first_line:
                    return doc_first_line
                else:
                    return str(test)
        def _makeResult(self):
            return self._TestResult(self.stream, self.descriptions,
                                    self.verbosity)
except ImportError:
    # Older Pythons with unittest1 don't export TextTestResult and also
    # don't duplicate test names
    _TestRunner = unittest.TextTestRunner


class RunAllTests(unittest.TestProgram):
    """Custom main program that also displays a final coverage report"""
    def __init__(self, *args, **keys):
        if coverage:
            # Start coverage testing now before we import any modules
            coverage.start()

        # Run the tests
        unittest.TestProgram.__init__(self, *args, **keys)

    def runTests(self):
        self.testRunner = _TestRunner(verbosity=self.verbosity)
        result = self.testRunner.run(self.test)

        if coverage:
            coverage.stop()
            coverage.the_coverage.collect()
            coverage.use_cache(False)
            print >> sys.stderr, "\nPython coverage report\n"

            # Don't show full paths to modules in coverage output
            cwd = os.path.dirname(sys.argv[0])
            topdir = os.path.abspath(os.path.join(cwd, '..', 'pyext')) + '/'
            coverage.the_coverage.relative_dir = topdir

            mods = [topdir + 'mdt/*.py']
            coverage.report(mods, file=sys.stderr)
            for cov in glob.glob('.coverage.*'):
                os.unlink(cov)
        sys.exit(not result.wasSuccessful())


def regressionTest():
    path = os.path.abspath(os.path.dirname(sys.argv[0]))
    files = os.listdir(path)
    test = re.compile("^test_.*\.py$", re.IGNORECASE)
    files = filter(test.search, files)
    modnames = [os.path.splitext(f)[0] for f in files]

    modobjs = [__import__(m) for m in modnames]
    tests = [unittest.defaultTestLoader.loadTestsFromModule(o) for o in modobjs]
    return unittest.TestSuite(tests)

if __name__ == "__main__":
    RunAllTests(defaultTest="regressionTest")
