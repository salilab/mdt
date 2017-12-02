"""Utility functions used by all IMP modules"""

from __future__ import print_function
import os.path
import re
import sys
import subst
import c_coverage
import sizeof_check
from SCons.Script import *

__all__ = ["add_common_variables", "MyEnvironment", "get_pyext_environment",
           "get_sharedlib_environment", "invalidate_environment", "embed"]

import SCons
_SWIGScanner = SCons.Scanner.ClassicCPP(
    "SWIGScan",
    ".i",
    "CPPPATH",
    '^[ \t]*[%,#][ \t]*(?:include|import)[ \t]*(<|")([^>"]+)(>|")'
)

# On older Pythons that don't have subprocess, fall back to os.popen3
try:
    import subprocess
    def MyPopen(cmd):
        return subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
except ImportError:
    class MyPopen(object):
        def __init__(self, cmd):
            (self.stdin, self.stdout, self.stderr) \
                = os.popen3(cmd, 't', -1)
        def wait(self):
            return 0

class WineEnvironment(Environment):
    """Environment to build Windows binaries under Linux, by running the
       MSVC compiler (cl) and linker (link) through wine, using the w32cc
       and w32link shell scripts"""
    def __init__(self, x64=False, platform='win32', CC=None, LINK=None,
                 path=None, **kw):
        self.x64 = x64
        if x64:
            bitprefix = 'w64'
        else:
            bitprefix = 'w32'
        if CC is None:
            CC = '%scc' % bitprefix
        if LINK is None:
            LINK = '%slink' % bitprefix
        if sys.platform != 'linux2':
            print("ERROR: Wine is supported only on Linux systems")
            Exit(1)
        self._fix_scons_msvc_detect()

        posix_env = Environment(platform='posix')
        if path is not None:
            posix_env['ENV']['PATH'] = path
        Environment.__init__(self, platform=platform, CC=CC, LINK=LINK,
                             ENV=posix_env['ENV'], **kw)
        self['PSPAWN'] = posix_env['PSPAWN']
        self['SPAWN'] = posix_env['SPAWN']
        self['SHELL'] = posix_env['SHELL']
        if x64:
            self['PYTHON'] = self.get('python', 'w64python')
        else:
            self['PYTHON'] = self.get('python', 'w32python')
        self['PATHSEP'] = ';'
        # Use / rather than \ path separator:
        try:
            self['LINKCOM'] = self['LINKCOM'].replace('.windows', '')
        except AttributeError:
            self['LINKCOM'] = '$LINK $LINKFLAGS /OUT:$TARGET $_LIBDIRFLAGS $_LIBFLAGS $_PDB $SOURCES'

    def _fix_scons_msvc_detect(self):
        """Ensure that MSVC auto-detection finds tools on Wine builds"""
        def _wine_read_reg32(value, hkroot=None):
            return '/usr/lib/w32comp/Program Files/' + \
                   'Microsoft Visual Studio 10.0'
        def _wine_read_reg64(value, hkroot=None):
            return '/usr/lib/w64comp/Program Files/' + \
                   'Microsoft Visual Studio 11.0'
        try:
            import SCons.Tool.MSCommon.common
        except ImportError:
            return # Older versions of scons don't have this module
        if self.x64:
            SCons.Tool.MSCommon.common.read_reg = _wine_read_reg64
        else:
            SCons.Tool.MSCommon.common.read_reg = _wine_read_reg32

def _get_python_include(env):
    """Get the directory containing Python.h"""
    if env['pythoninclude']:
        return env['pythoninclude']
    elif env['wine64']:
        return '/usr/lib/w64comp/w64python/2.7/include/'
    elif env['wine']:
        return '/usr/lib/w32comp/w32python/2.7/include/'
    else:
        import distutils.sysconfig
        return distutils.sysconfig.get_python_inc()

def _add_release_flags(env):
    """Add compiler flags for release builds, if requested"""
    if env.get('release', False):
        env.Append(CPPDEFINES='NDEBUG')

def CheckGCCVisibility(context):
    """Check if the compiler supports setting visibility of symbols"""
    context.Message('Checking whether compiler supports -fvisibility...')
    lastCCFLAGS = context.env['CCFLAGS']
    context.env.Append(CCFLAGS="-fvisibility=hidden")
    text = """
__attribute__ ((visibility("default")))
int main(void)
{ return 0; }
"""
    res = context.TryLink(text, '.c')
    context.env.Replace(CCFLAGS=lastCCFLAGS)
    if not res:
        context.Result("no")
    else:
        context.env.Append(VIS_CPPDEFINES=['GCC_VISIBILITY'],
                           VIS_CCFLAGS="-fvisibility=hidden")
        context.Result("yes")
    return res

def check_pkgconfig(context, pkgconfig_name, human_name, env_key):
    context.Message('Checking for %s using pkg-config...' % human_name)
    try:
        flags = context.env.ParseFlags('!pkg-config --cflags --libs ' \
                                       + pkgconfig_name)
    except OSError, detail:
        context.Result("failed: %s" % str(detail))
        return False
    context.env[env_key] = flags
    context.Result("found, headers in %s" % flags['CPPPATH'][0])
    return True

def CheckGlib2(context):
    # On Windows, assume glib-2.0 library is in the path, since pkg-config
    # is not usually present
    if context.env['PLATFORM'] == 'win32':
        context.Result("Assuming glib2 can be linked by adding glib-2.0.lib")
        context.env['GLIB'] = {'SHLINKFLAGS':['glib-2.0.lib']}
        return True
    elif check_pkgconfig(context, pkgconfig_name='glib-2.0',
                         human_name='GLib2', env_key='GLIB'):
        return True
    else:
        print("GLib2 and pkg-config are required to install this software.")
        print("Install the glib2-devel and pkg-config (or similar) packages,")
        print("or download the GLib2 sourcecode from www.gtk.org")
        Exit(1)

def check_modeller_python(context):
    context.Message("Checking for MODELLER in Python path...")
    try:
        import modeller
    except ImportError, e:
        context.Result("not found (specify installation path with 'modeller' "
                       "scons option): %s" % str(e))
        return False
    try:
        exetype = modeller.info.exe_type
        bindir = modeller.info.bindir
    except AttributeError, e:
        context.Result("'import modeller' succeeded, but the package does "
                       "not appear to be MODELLER; perhaps a 'modeller.py' "
                       "file in the current directory?")
        return False

    modeller = os.path.normpath(os.path.join(bindir, '..'))
    if not parse_modeller_dirs(context, modeller, exetype):
        return False
    else:
        context.Result("found in %s" % modeller)
        return True

def parse_modeller_dirs(context, modeller, exetype):
    include = ['%s/src/include' % modeller,
               '%s/src/include/%s' % (modeller, exetype)]
    platform = context.env['PLATFORM']
    if exetype in ('i386-w32', 'x86_64-w64'):
        libpath = ['%s/src/main' % modeller]
        if not os.path.exists(libpath[0]):
            libpath = ['%s/bin' % modeller]
        if platform != 'win32':
            context.Result("MODELLER is built for Windows, but this is not " + \
                           "a Windows scons run (tip: can run on Linux " + \
                           "using Wine with 'scons wine=true'")
            return False
        # Since we don't install a symlink for
        # libmodeller.lib -> libmodeller${SOVER}.lib, glob to find a suitable
        # Modeller import library
        libs = []
        for d in os.listdir(libpath[0]):
            if d.startswith('libmodeller') and d.endswith('.lib'):
                libs.append(d[3:-4])
                break
        libs.append('saxs')
    else:
        libpath = ['%s/lib/%s' % (modeller, exetype)]
        if platform == 'win32':
            context.Result("this is a Windows scons run, but this is not a " + \
                           "Windows MODELLER binary")
            return False
        libs = ["modeller", "saxs"]
    if exetype in ('mac10v4-xlf', 'mac10v4-gnu'):
        libs += ["hdf5", "hdf5_hl"]
    elif exetype in ('mac10v4-intel', 'mac10v4-intel64'):
        libs += ["hdf5", "hdf5_hl", "imf", "svml", "ifcore", "irc"]
    modpy = "%s/bin/modpy.sh" % modeller
    # If the modpy.sh script doesn't exist, assume that Modeller will work
    # without it (e.g. on Macs, using the binary .dmg install):
    if not os.path.exists(modpy):
        modpy = ''
    context.env['MODELLER_MODPY'] = modpy
    context.env['MODELLER'] = {'CPPPATH':include, 'LIBPATH':libpath,
                               'LIBS':libs}
    return True

def CheckModeller(context):
    """Find Modeller include and library directories"""
    modeller = context.env['modeller']
    if (modeller is False or modeller is 0) \
       and (check_pkgconfig(context, pkgconfig_name='modeller',
                            human_name='MODELLER', env_key='MODELLER') \
            or check_modeller_python(context)):
        # Assume that modpy.sh script is not necessary
        context.env['MODELLER_MODPY'] = ''
        return True
    context.Message("Checking for MODELLER using 'modeller' scons option...")
    if modeller is False or modeller is 0:
        context.Result("not found")
        return False
    # Find MODELLER script
    moddir = "%s/bin" % modeller
    try:
        files = os.listdir(moddir)
    except OSError, e:
        context.Result("could not find MODELLER directory %s: %s" % (moddir, e))
        return False
    files.sort()
    r = re.compile('mod(SVN|\d+[v.]\d+)(\.exe)?$')
    files = [f for f in files if r.match(f)]
    if len(files) == 0:
        context.Result("could not find MODELLER script in %s" % moddir)
        return False
    # Last matching entry is probably the latest version:
    modbin = os.path.join(moddir, files[-1])
    try:
        p = MyPopen(modbin + " -")
        print("print 'EXE type: ', info.exe_type", file=p.stdin)
        p.stdin.close()
    except IOError, e:
        context.Result("could not run MODELLER script %s: %s" % (modbin, e))
        return False
    err = p.stderr.read()
    exetype = None
    for line in p.stdout:
        if line.startswith("EXE type"):
            exetype=line[11:].rstrip('\r\n')
    ret = p.wait()
    if exetype is None:
        if err or ret != 0:
            context.Result("could not run MODELLER script %s; %d, %s" \
                           % (modbin, ret, err))
        else:
            context.Result("unknown error running MODELLER script %s" % modbin)
        return False
    if not parse_modeller_dirs(context, modeller, exetype):
        return False
    context.Result(modeller)
    return True

def _modeller_check_failed(require_modeller):
    """Print an informative message if the Modeller check failed"""
    msg = "  Use the modeller command line option (or options file) to\n" + \
          "  set the directory where Modeller is installed\n" + \
          "  (run 'scons -h' for help.)"

    print()
    if require_modeller:
        print("ERROR: MODELLER is required to build this package\n\n" + msg)
        Exit(1)
    else:
        print("  MODELLER was not found: build will continue but some")
        print("  functionality will be missing.\n\n" + msg)


def MyEnvironment(variables=None, require_modeller=True, *args, **kw):
    """Create an environment suitable for building IMP modules"""
    import platform
    # First make a dummy environment in order to evaluate all variables, since
    # env['wine'] will tell us which 'real' environment to create:
    env = Environment(tools=[], variables=variables)
    path = env['ENV']['PATH']
    if env.get('path') is not None:
        path = env['path'] + os.path.pathsep + path

    if env['wine64']:
        env = WineEnvironment(x64=True,
                              variables=variables, path=path, *args, **kw)
    elif env['wine']:
        env = WineEnvironment(x64=False,
                              variables=variables, path=path, *args, **kw)
    else:
        env = Environment(variables=variables,
                          ENV = {'PATH':path}, *args, **kw)
        env['PYTHON'] = env.get('python', 'python')
        env['PATHSEP'] = os.path.pathsep
    try:
        env['SHLINKFLAGS'] = [x for x in env['SHLINKFLAGS'] \
                              if not '-no_archive' in x]
    except ValueError:
        pass
    env.Prepend(SCANNERS = _SWIGScanner)
    subst.TOOL_SUBST(env)
    env.AddMethod(c_coverage.CCoverageTester)

    if env['PLATFORM'] == 'win32' or env['wine'] or env['wine64']:
        env['SHLIBPREFIX'] = env['LIBLINKPREFIX'] = env['LIBPREFIX'] = 'lib'
        env['WINDOWSEXPPREFIX'] = 'lib'
        env['LIBSUFFIX'] = '.lib'
        # Make sure we get the same Windows C/C++ library as Modeller, and
        # enable C++ exception handling
        env.Append(CFLAGS="/MD")
        env.Append(CXXFLAGS="/MD /GR /GX")

    if env['PLATFORM'] == 'win32' and not env['wine'] and not env['wine64']:
        env['DEVNULL'] = 'nul'
        # 'python' typically isn't in the path on Windows systems
        env['PYTHON'] = sys.executable
    else:
        env['DEVNULL'] = '/dev/null'

    # Make sure destdir is relative to the toplevel directory,
    # if a relative path
    if env['destdir'] != '':
        env['destdir'] = Dir(env['destdir']).abspath

    if env['CC'] == 'gcc':
        if env.get('coverage', False):
            env.Append(CCFLAGS="-Wall -Werror -g -O0 -fprofile-arcs "
                               "-ftest-coverage")
            env.Append(LINKFLAGS="-fprofile-arcs -ftest-coverage")
        else:
            env.Append(CCFLAGS="-Wall -Werror -g -O3")
    elif env.get('coverage', False):
        print("ERROR: C coverage testing currently only works with gcc")
        Exit(1)
    _add_release_flags(env)

    if env.get('includepath', None):
        env['includepath'] = [os.path.abspath(x) for x in \
                              env['includepath'].split(os.path.pathsep)]
        env.Append(CPPPATH=env['includepath'])
    if env.get('libpath', None):
        env['libpath'] = [os.path.abspath(x) for x in \
                          env['libpath'].split(os.path.pathsep)]
        env.Prepend(LIBPATH=env['libpath'])

    system = platform.system()
    if system == 'SunOS':
        # Find locally-installed libraries in /usr/local (e.g. for SWIG)
        env['ENV']['LD_LIBRARY_PATH'] = '/usr/local/lib'
    # Make Modeller exetype variable available:
    if os.environ.has_key('EXECUTABLE_TYPESVN'):
        env['ENV']['EXECUTABLE_TYPESVN'] = os.environ['EXECUTABLE_TYPESVN']
    # Set empty variables in case checks fail or are not run (e.g. clean)
    env['MODELLER_MODPY'] = ''
    env['MODELLER'] = {}
    env['GLIB'] = {}
    if not env.GetOption('clean') and not env.GetOption('help'):
        custom_tests = {'CheckGCCVisibility': CheckGCCVisibility,
                        'CheckModeller': CheckModeller,
                        'CheckGlib2': CheckGlib2,
                        'CheckSizeof': sizeof_check._check}
        conf = env.Configure(custom_tests = custom_tests)
        env['SIZEOF_POINTER'] = conf.CheckSizeof()
        conf.CheckGCCVisibility()
        conf.CheckGlib2()
        # Look for the hdf5.h header file. todo: make sure its version
        # matches that of the HDF5 libraries included with Modeller.
        if not conf.CheckCHeader('hdf5.h'):
            print("HDF5 is required to build this package. Make sure you have")
            print("the same version that MODELLER is built with (check the")
            print("Modeller changelog). If it is installed already, set the")
            print("scons 'includepath' option to the directory containing")
            print("the hdf5.h header file.")
            Exit(1)

        # Check explicitly for False, since all checks will return Null if
        # configure has been disabled
        if conf.CheckModeller() is False:
            _modeller_check_failed(require_modeller)
        conf.Finish()
    else:
        env['SIZEOF_POINTER'] = ''
    return env

def _fix_aix_cpp_link(env, cplusplus, linkflags):
    """On AIX things get confused if AIX C but not AIX C++ is installed - AIX C
       options get passed to g++ - so hard code GNU link flags"""
    if cplusplus and 'aixcc' in env['TOOLS'] and 'aixc++' not in env['TOOLS'] \
       and 'g++' in env['TOOLS']:
        slflags = str(env[linkflags])
        env[linkflags] = slflags.replace('-qmkshrobj -qsuppress=1501-218',
                                         '-shared')

def get_sharedlib_environment(env, cppdefine, cplusplus=False):
    """Get a modified environment suitable for building shared libraries
       (i.e. using gcc ELF visibility macros or MSVC dllexport/dllimport macros
       to mark dynamic symbols as exported or private). `cppdefine` should be
       the name of a cpp symbol to define to tell MSVC that we are building the
       library (by convention something of the form FOO_EXPORTS).
       If `cplusplus` is True, additional configuration suitable for a C++
       shared library is done."""
    e = env.Clone()
    e.Append(CPPDEFINES=[cppdefine, '${VIS_CPPDEFINES}'],
             CCFLAGS='${VIS_CCFLAGS}')

    _fix_aix_cpp_link(e, cplusplus, 'SHLINKFLAGS')
    return e

# Workaround for SWIG bug #1863647: Ensure that the PySwigIterator class is
# renamed with a module-specific prefix, to avoid collisions when using
# multiple modules
class _swig_postprocess(object):
    def __init__(self, modprefix):
        self.modprefix = modprefix
    def builder(self, source, target, env):
        wrap_c = target[0].path
        lines = file(wrap_c, 'r').readlines()
        repl = '"swig::%s_PySwigIterator *"' % self.modprefix
        fh = file(wrap_c, 'w')
        for line in lines:
            fh.write(line.replace('"swig::PySwigIterator *"', repl))
        fh.close()
        return 0

def _check_swig_version(env, needversion):
    """Make sure that a new enough version of SWIG was found."""
    if env.GetOption('clean') or env.GetOption('help'):
        return

    needversion_str = ".".join([str(x) for x in needversion])
    failmsg = """
SWIG version %s or later must be installed to build MDT, but
%s.
Please make sure 'swig' is found in the path passed to scons.

In particular, if you have SWIG installed in a non-standard location,
please use the 'path' option to add this location to the search path.
For example, if you have SWIG installed in /opt/local/bin/, edit (or create)
config.py and add the line

path='/opt/local/bin'
"""
    version = env.get('SWIGVERSION', None)
    try:
        v = [int(x) for x in version.split(".")]
    except ValueError:
        print(failmsg % (needversion_str,
                         "it could not be found on your system"))
        Exit(1)
    if v >= needversion:
        print("Checking for SWIG >= %s: %s found" \
              % (needversion_str, ".".join([str(x) for x in v])))
        return
    else:
        print(failmsg % (needversion_str,
                         "only an older version (%s) " % version + \
                         "was found on your system"))
        Exit(1)

def get_pyext_environment(env, mod_prefix, cplusplus=False):
    """Get a modified environment for building a Python extension.
       `mod_prefix` should be a unique prefix for this module.
       If `cplusplus` is True, additional configuration suitable for a C++
       extension is done."""
    from platform import system
    e = env.Clone()
    _check_swig_version(e, [1,3,39])

    if cplusplus and isinstance(e['SWIGCOM'], str):
        # See _swig_postprocess class comments:
        repl = '$SWIG -DPySwigIterator=%s_PySwigIterator ' % mod_prefix
        e['SWIGCOM'] = e['SWIGCOM'].replace('$SWIG ', repl)
        e['SWIGCOM'] = [e['SWIGCOM'], _swig_postprocess(mod_prefix).builder]

    e['LDMODULEPREFIX'] = ''
    # We're not going to link against the extension, so don't need a Windows
    # import library (.lib file):
    e['no_import_lib'] = 1
    platform = e['PLATFORM']
    if e['wine'] or e['wine64']:
        # Have to set SHLIBSUFFIX and PREFIX on Windows otherwise the
        # mslink tool complains
        e['SHLIBPREFIX'] = ''
        e['LDMODULESUFFIX'] = e['SHLIBSUFFIX'] = '.pyd'
        # Directory containing python27.lib:
        if e['wine64']:
            e.Append(LIBPATH=['/usr/lib/w64comp/w64python/2.7/lib/'])
        else:
            e.Append(LIBPATH=['/usr/lib/w32comp/w32python/2.7/lib/'])
    else:
        if platform == 'aix':
            # Make sure compilers are in the PATH, so that Python's script for
            # building AIX extension modules can find them:
            e['ENV']['PATH'] += ':/usr/vac/bin'
        from distutils.sysconfig import get_config_vars
        vars = get_config_vars('CC', 'CXX', 'OPT', 'BASECFLAGS', 'LDSHARED',
                               'SO')
        (cc, cxx, opt, basecflags, ldshared, so) = vars
        # distutils on AIX can get confused if AIX C but GNU C++ is installed:
        if platform == 'aix' and cxx == '':
            cxx = 'g++'
        # Don't require stack protector stuff on Linux, as this adds a
        # requirement for glibc-2.4:
        if opt:
            opt = opt.replace("-fstack-protector-strong", "")
            opt = opt.replace("-fstack-protector", "")
        # Remove options that don't work with C++ code:
        if cplusplus and opt:
            opt = opt.replace("-Wstrict-prototypes", "")
        if cc:
            e.Replace(CC=cc)
        if cxx:
            e.Replace(CXX=cxx)
        if so:
            e.Replace(LDMODULESUFFIX=so)
        if basecflags or opt:
            e.Replace(CPPFLAGS=basecflags.split() + opt.split())

        # Remove NDEBUG preprocessor stuff if defined (we do it ourselves for
        # release builds)
        try:
            e['CPPFLAGS'].remove('-DNDEBUG')
        except (ValueError, KeyError):
            pass

        # Some gcc versions don't like the code that SWIG generates - but let
        # that go, because we have no control over it:
        try:
            e['CCFLAGS'].remove('-Werror')
        except ValueError:
            pass
        # AIX tries to use the C compiler rather than g++, so hardcode it here:
        if platform == 'aix' and cplusplus:
            ldshared = ldshared.replace(' cc_r', ' g++')
        # Default link flags on OS X don't work for us:
        if platform == 'darwin':
            e.Replace(LDMODULEFLAGS= \
                      '$LINKFLAGS -bundle -flat_namespace -undefined suppress')
        # Don't set link flags on Linux, as that would negate our GNU_HASH check
        elif system() != "Linux" and ldshared:
            e['LDMODULEFLAGS'] = []
            e['SHLINK'] = e['LDMODULE'] = ldshared
    e.Append(CPPPATH=[_get_python_include(e)])
    _fix_aix_cpp_link(e, cplusplus, 'SHLINKFLAGS')
    return e

def invalidate_environment(env, fail_builder):
    """'Break' an environment, so that any builds with it use the fail_builder
       function (which should be an Action which terminates the build)"""
    for var in ('SHLINKCOM', 'CCCOM', 'CXXCOM', 'SHCCCOM', 'SHCXXCOM',
                'SWIGCOM'):
        env[var] = fail_builder

def add_common_variables(vars, package):
    """Add common variables to an SCons Variables object."""
    libdir = '${prefix}/lib'
    if hasattr(os, 'uname') and sys.platform == 'linux2' \
       and os.uname()[-1] == 'x86_64':
        # Install in /usr/lib64 rather than /usr/lib on x86_64 Linux boxes
        libdir += '64'

    vars.Add(PathVariable('prefix', 'Top-level installation directory', '/usr',
                          PathVariable.PathAccept))
    # Note that destdir should not affect any compiled-in paths; see
    # http://www.gnu.org/prep/standards/html_node/DESTDIR.html
    vars.Add(PathVariable('destdir',
                          'String to prepend to every installed filename',
                          '', PathVariable.PathAccept))
    vars.Add(PathVariable('datadir', 'Data file installation directory',
                          '${prefix}/share/%s' % package,
                          PathVariable.PathAccept))
    vars.Add(PathVariable('libdir', 'Shared library installation directory',
                          libdir, PathVariable.PathAccept))
    vars.Add(PathVariable('includedir', 'Include file installation directory',
                          '${prefix}/include', PathVariable.PathAccept))
    vars.Add(PathVariable('pythondir', 'Python module installation directory',
                          '${libdir}/python%d.%d/site-packages' \
                          % sys.version_info[0:2], PathVariable.PathAccept))
    vars.Add(PathVariable('pyextdir',
                          'Python extension module installation directory',
                          '${pythondir}', PathVariable.PathAccept))
    vars.Add(PathVariable('docdir', 'Documentation installation directory',
                          '${prefix}/share/doc/%s' % package,
                          PathVariable.PathAccept))
    vars.Add(PackageVariable('pythoninclude',
                             'Directory holding Python include files ' + \
                             '(if unspecified, distutils location is used)',
                             'no'))
    vars.Add('python', "Python interpreter to use to run unit tests "
             "(if unspecified, 'python' is used on regular builds, "
             "'w64python' if wine64=True, 'w32python' if wine=True)", None)
    vars.Add(PackageVariable('modeller', 'Location of the MODELLER package',
                             'no'))
    vars.Add(BoolVariable('wine',
                          'Build using 32bit Windows tools via Wine emulation',
                          False))
    vars.Add(BoolVariable('wine64',
                          'Build using 64bit Windows tools via Wine emulation',
                          False))
    vars.Add(BoolVariable('release',
                          'Disable most runtime checks (e.g. for releases)',
                          False))
    vars.Add(BoolVariable('coverage',
                          'Enable coverage testing of C code (Python code '
                          'is always tested if at least version 3 of the '
                          'coverage module is available) in unit tests '
                          '(gcc only)', False))
    vars.Add(PathVariable('html_coverage',
                          'Directory to output HTML coverage reports into '
                          '(requires coverage=true and the lcov package for '
                          'C output, and the Python coverage module for '
                          'Python output).',
                          None, PathVariable.PathIsDirCreate))
    vars.Add(PathVariable('includepath', 'Include search path ' + \
                          '(e.g. "/usr/local/include:/opt/local/include")',
                          None, PathVariable.PathAccept))
    vars.Add(PathVariable('libpath', 'Library search path ' + \
                          '(e.g. "/usr/local/lib:/opt/local/lib")',
                          None, PathVariable.PathAccept))
    vars.Add(PathVariable('path', 'Path to search for build tools ' + \
                          '(e.g. "/usr/local/bin:/opt/local/bin")',
                          None, PathVariable.PathAccept))
