"""SCons.Tool.swig

Tool-specific initialization for swig.

There normally shouldn't be any need to import this module directly.
It will usually be imported through the generic SCons.Tool.Tool()
selection method.

This is a modified version of the SWIG tool, to work with older versions of
scons and still populate SWIGVERSION.

"""

from __future__ import print_function

#
# Copyright (c) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 The SCons Foundation
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
# KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

__revision__ = "src/engine/SCons/Tool/swig.py 3266 2008/08/12 07:31:01 knight"

import os.path
import re
import sys

import SCons.Action
import SCons.Defaults
import SCons.Scanner
import SCons.Tool
import SCons.Util
try:
    import subprocess
except ImportError:
    subprocess = None


def _get_swig_version(env):
    """Run the SWIG command line tool to get and return the version number"""
    if not env['SWIG']:
        return ""
    if subprocess is not None and hasattr(SCons.Action, '_subproc'):
        pipe = SCons.Action._subproc(env, [env['SWIG'], '-version'],
                                     stdin='devnull',
                                     stderr='devnull',
                                     stdout=subprocess.PIPE)
        if pipe.wait() != 0:
            return ""
        out = pipe.stdout.read()
    else:
        out = os.popen(env['SWIG'] + ' -version').read()
    if sys.version_info[0] >= 3 and isinstance(out, bytes):
        out = out.decode()
    match = re.search(r'SWIG Version\s+(\S+)\s*$', out, re.MULTILINE)
    if match:
        return match.group(1)
    else:
        return ""


SwigAction = SCons.Action.Action('$SWIGCOM', '$SWIGCOMSTR')


def swigSuffixEmitter(env, source):
    if '-c++' in SCons.Util.CLVar(env.subst("$SWIGFLAGS", source=source)):
        return '$SWIGCXXFILESUFFIX'
    else:
        return '$SWIGCFILESUFFIX'


# Match '%module test', as well as '%module(directors="1") test'
_reModule = re.compile(r'%module(?:\s*\(.*\))?\s+(.+)')


def _swigEmitter(target, source, env):
    swigflags = env.subst("$SWIGFLAGS", target=target, source=source)
    flags = SCons.Util.CLVar(swigflags)
    for src in source:
        src = str(src.rfile())
        mnames = None
        if "-python" in flags and "-noproxy" not in flags:
            if mnames is None:
                mnames = _reModule.findall(open(src).read())
            target.extend(map(lambda m, d=target[0].dir:
                              d.File(m + ".py"), mnames))
        if "-java" in flags:
            if mnames is None:
                mnames = _reModule.findall(open(src).read())
            java_files = map(lambda m: [m + ".java", m + "JNI.java"], mnames)
            java_files = SCons.Util.flatten(java_files)
            outdir = env.subst('$SWIGOUTDIR', target=target, source=source)
            if outdir:
                java_files = map(lambda j, o=outdir: os.path.join(o, j), java_files)
            java_files = map(env.fs.File, java_files)
            for jf in java_files:
                def t_from_s(t, p, s, x):
                    return t.dir
                SCons.Util.AddMethod(jf, t_from_s, 'target_from_source')
            target.extend(java_files)
    return (target, source)


def generate(env):
    """Add Builders and construction variables for swig to an Environment."""
    c_file, cxx_file = SCons.Tool.createCFileBuilders(env)

    c_file.suffix['.i'] = swigSuffixEmitter
    cxx_file.suffix['.i'] = swigSuffixEmitter

    c_file.add_action('.i', SwigAction)
    c_file.add_emitter('.i', _swigEmitter)
    cxx_file.add_action('.i', SwigAction)
    cxx_file.add_emitter('.i', _swigEmitter)

    java_file = SCons.Tool.CreateJavaFileBuilder(env)

    java_file.suffix['.i'] = swigSuffixEmitter

    java_file.add_action('.i', SwigAction)
    java_file.add_emitter('.i', _swigEmitter)

    env['SWIG'] = env.WhereIs('swig')
    env['SWIGVERSION'] = _get_swig_version(env)
    env['SWIGFLAGS'] = SCons.Util.CLVar('')
    env['SWIGCFILESUFFIX'] = '_wrap$CFILESUFFIX'
    env['SWIGCXXFILESUFFIX'] = '_wrap$CXXFILESUFFIX'
    env['_SWIGOUTDIR'] = '${"-outdir " + str(SWIGOUTDIR)}'
    env['SWIGPATH'] = []
    env['SWIGINCPREFIX'] = '-I'
    env['SWIGINCSUFFIX'] = ''
    env['_SWIGINCFLAGS'] = '$( ${_concat(SWIGINCPREFIX, SWIGPATH, SWIGINCSUFFIX, __env__, RDirs, TARGET, SOURCE)} $)'
    env['SWIGCOM'] = '$SWIG -o $TARGET ${_SWIGOUTDIR} ${_SWIGINCFLAGS} $SWIGFLAGS $SOURCES'

    expr = '^[ \t]*%[ \t]*(?:include|import|extern)[ \t]*(<|"?)([^>\\s"]+)(?:>|"?)'
    scanner = SCons.Scanner.ClassicCPP("SWIGScan", ".i", "SWIGPATH", expr)

    env.Append(SCANNERS=scanner)


def exists(env):
    return env.Detect(['swig'])
