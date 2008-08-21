# Include IMP build utility functions:
import sys
sys.path.append('../')
from tools import *

# Set up build environment:
opts = Options('config.py', ARGUMENTS)
add_common_options(opts, 'mdt')
env = MyEnvironment(options=opts, require_modeller=True,
                    tools=["default", "doxygen", "docbook", "epydoc"],
                    toolpath=["../tools"])
Help(opts.GenerateHelpText(env))

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
