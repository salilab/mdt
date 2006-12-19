# Include TNT build utility functions:
execfile('../SCons.include')

modconfig = get_modeller_config()
pythoninc = get_python_include(modconfig)
env = MyEnvironment(modconfig)

Export('env', 'modconfig', 'pythoninc', 'configure_for_pyext')

SConscript(['src/SConscript',
            'pyext/SConscript',
            'test/SConscript'])
