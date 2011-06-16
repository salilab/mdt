import sys
import os
import glob

class _CCoverageTester(object):
    def __init__(self, env):
        self._env = env
        self._files = []
        self._coverage = self._env.get('coverage', False)

    def add(self, directory, pattern, report):
        self._files.append([directory, pattern, report])

    def Execute(self, *args, **keys):
        self._cleanup_coverage_files()
        ret = self._env.Execute(*args, **keys)
        if self._coverage:
            self._report()
        self._cleanup_coverage_files()
        return ret

    def _cleanup_coverage_files(self):
        for dir, pattern, report in self._files:
            for f in glob.glob(os.path.join(dir, '*.gcda')):
                os.unlink(f)

    def _report(self):
        for dir, pattern, report in self._files:
            if report:
                self._run_gcov(dir, pattern)
        print >> sys.stderr, "\n\nC coverage report\n"
        print >> sys.stderr, "%-41s Stmts   Exec  Cover   Missing" % "Name"
        divider = "-" * 71
        print >> sys.stderr, divider
        total_statements = total_executed = 0
        for cov in glob.glob("*.gcov"):
            source, statements, executed, missing = self._parse_gcov_file(cov)
            self._report_gcov_file(source, statements, executed, missing)
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
        return source, line_number, line_number - len(missing), missing

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
