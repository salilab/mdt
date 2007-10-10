# Include IMP build utility functions:
execfile('../SCons.include')

# Set up build environment:
opts = Options()
add_common_options(opts, 'mdt')
env = MyEnvironment(options=opts, require_modeller=True,
                    tools=["default", "doxygen"], toolpath=["../tools"])
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
