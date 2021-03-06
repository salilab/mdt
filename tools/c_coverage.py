import sys
import os
import glob
import tempfile
import shutil

class _TempDir(object):
    """Simple RAII-style class to make a temporary directory"""
    def __init__(self):
        self._tmpdir = tempfile.mkdtemp()
    def tempfile(self, fname):
        return os.path.join(self._tmpdir, fname)
    def __del__(self):
        shutil.rmtree(self._tmpdir, ignore_errors=True)


class _CCoverageTester(object):
    def __init__(self, env):
        self._env = env
        self._files = []
        self._coverage = self._env.get('coverage', False)
        self._html_coverage = self._env.get('html_coverage', None)

    def add(self, directory, pattern, report):
        self._files.append([directory, pattern, report])

    def Execute(self, *args, **keys):
        self._cleanup_coverage_files()
        ret = self._env.Execute(*args, **keys)
        if self._coverage:
            self._report()
            if self._html_coverage:
                self._report_html(self._html_coverage)
        self._cleanup_coverage_files()
        return ret

    def _cleanup_coverage_files(self):
        for dir, pattern, report in self._files:
            for f in glob.glob(os.path.join(dir, '*.gcda')):
                os.unlink(f)

    def _report_html(self, html_coverage):
        import subprocess
        def call(args):
            r = subprocess.call(args)
            if r != 0:
                raise OSError("%s failed with exit code %d" % (args[0], r))
        topdir = os.getcwd()
        d = _TempDir()
        call(['lcov', '-c', '-b', '.', '-d', 'src',
              '-o', d.tempfile('all.info')])
        call(['lcov', '-e', d.tempfile('all.info'),
              os.path.join(topdir, 'src', '*'), '-o', d.tempfile('mdt.info')])
        call(['genhtml', '--legend', '--title', 'MDT C Coverage ',
              '-o', os.path.join(html_coverage, 'c'), d.tempfile('mdt.info')])

    def _report(self):
        print >> sys.stderr, "\n\nC coverage report\n"
        print >> sys.stderr, "%-41s Stmts   Exec  Cover   Missing" % "Name"
        divider = "-" * 71
        print >> sys.stderr, divider
        total_statements = total_executed = 0
        for dir, pattern, report in self._files:
            if report:
                self._run_gcov(dir, pattern)
                covs = glob.glob("*.gcov")
                covs.sort()
                for cov in covs:
                    source, statements, executed, missing \
                                 = self._parse_gcov_file(cov)
                    # Skip system and non-MDT headers
                    if not source.startswith('/'):
                        self._report_gcov_file(source, statements, executed,
                                               missing)
                        total_statements += statements
                        total_executed += executed
                    os.unlink(cov)
        print >> sys.stderr, divider
        self._report_gcov_file('TOTAL', total_statements, total_executed, [])

    def _get_missing_ranges(self, missing):
        ranges = []
        start_range = None
        end_range = None
        def add_range():
            if start_range is not None:
                if end_range == start_range:
                    ranges.append('%d' % start_range)
                else:
                    ranges.append('%d-%d' % (start_range, end_range))
        for line in missing:
            if end_range is not None and end_range == line - 1:
                end_range = line
            else:
                add_range()
                start_range = line
                end_range = line
        add_range()
        return ", ".join(ranges)

    def _parse_gcov_file(self, cov):
        executable_statements = 0
        missing = []
        for line in open(cov):
            spl = line.split(':', 2)
            calls = spl[0].strip()
            line_number = int(spl[1])
            if calls == '#####':
                missing.append(line_number)
            if line_number == 0:
                if spl[2].startswith('Source:'):
                    source = spl[2][7:].strip()
            elif calls != '-':
                executable_statements += 1
        return (source, executable_statements,
                executable_statements - len(missing), missing)

    def _report_gcov_file(self, source, statements, executed, missing):
        cover = float(executed) * 100. / float(statements)

        if len(source) > 40:
            source = "[..]" + source[-36:]

        print >> sys.stderr, "%-40s %6d %6d %5d%%   %s" \
              % (source, statements, executed, cover,
                 self._get_missing_ranges(missing))

    def _run_gcov(self, dir, pattern):
        import subprocess
        cmd = 'gcov -o %s %s' % (dir, os.path.join(dir, pattern))
        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        out = p.stdout.read()
        err = p.stderr.read()
        ret = p.wait()
        if ret != 0:
            raise OSError("'%s' failed with code %d, error %s" \
                          % (cmd, ret, err))

def CCoverageTester(env):
    return _CCoverageTester(env)
