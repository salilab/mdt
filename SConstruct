# Include IMP build utility functions:
import sys
sys.path.append('../')
from tools import *
from tools import sizeof_check

# We need scons 0.98 or later
EnsureSConsVersion(0, 98)

# Set up build environment:
vars = Variables('config.py', ARGUMENTS)
add_common_variables(vars, 'mdt')
env = MyEnvironment(variables=vars, require_modeller=True,
                    tools=["default", "sphinx"],
                    toolpath=["../tools"])
Help(vars.GenerateHelpText(env))

sizeof_check.configure_check(env)

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

# bin script first require Python extensions to be built:
env.Depends(bin, [pyext, pyso])

# Build the binaries by default:
env.Default(bin)
