# Include IMP build utility functions:
execfile('../SCons.include')

# Get Modeller location, and set up a build environment:
modconfig = get_modeller_config()
env = MyEnvironment(modconfig, tools=["default", "doxygen"],
                    toolpath=["../tools"])

# Make these objects available to SConscript files:
Export('env', 'modconfig', 'get_pyext_environment', 'get_sharedlib_environment',
       'is_wine_platform')

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
