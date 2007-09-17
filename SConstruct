# Include NMP build utility functions:
execfile('../SCons.include')

# Get Modeller and Python locations, and set up a build environment:
modconfig = get_modeller_config()
pythoninc = get_python_include(modconfig)
env = MyEnvironment(modconfig, tools=["default", "doxygen"],
                    toolpath=["../tools"])

# Make these objects available to SConscript files:
Export('env', 'modconfig', 'pythoninc', 'configure_for_pyext',
       'is_wine_platform')

# Subdirectories to build:
test = SConscript('test/SConscript')
pyso, pyext = SConscript('pyext/SConscript')
src = SConscript('src/SConscript')
data = SConscript('data/SConscript')
Export('pyext', 'pyso')
doc = SConscript('doc/SConscript')

# testcases first require Python extensions to be built:
env.Depends(test, pyext)

# Build the C library (src directory) and Python extension by default:
env.Default(src, pyext)
