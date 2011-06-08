# Include build utility functions:
from tools import *

# We need scons 0.98 or later
EnsureSConsVersion(0, 98)

# Set up build environment:
vars = Variables('config.py', ARGUMENTS)
add_common_variables(vars, 'mdt')
env = MyEnvironment(variables=vars, require_modeller=True,
                    tools=["default", "swig", "sphinx"],
                    toolpath=["tools"])

# Version number
env['MDT_VERSION'] = 'SVN'

Help(vars.GenerateHelpText(env))

# Make these objects available to SConscript files:
Export('env', 'get_pyext_environment', 'get_sharedlib_environment')

# Subdirectories to build:
bin = SConscript('bin/SConscript')
Export('bin')
test = SConscript('test/SConscript')
pyso, pyext = SConscript('pyext/SConscript')
src = SConscript('src/SConscript')
data = SConscript('data/SConscript')
Export('pyext', 'pyso')
doc = SConscript('doc/SConscript')
rpm = SConscript('tools/rpm/SConscript')

# bin script first require Python extensions to be built:
env.Depends(bin, [pyext, pyso])

# Build the binaries by default:
env.Default(bin)
