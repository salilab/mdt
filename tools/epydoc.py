# Simple Epydoc tool and builder.

import os
from SCons.Script import *

# Build epydoc documentation:
def _action_epydoc(target, source, env):
    outdir = os.path.dirname(target[0].path)
    bin = source[0].path
    app = "%s %s %s -o %s %s" % (bin, env['EPYDOC'], env['EPYDOC_OPTS'], outdir,
                                 ' '.join([x.path for x in source[1:]]))
    return env.Execute(app)

def generate(env):
    """Add builders and construction variables for the epydoc tool."""
    import SCons.Builder
    builder = SCons.Builder.Builder(action=_action_epydoc)
    # Use Unix 'install' rather than env.InstallAs(), due to scons bug #1751
    install = SCons.Builder.Builder(action="install -d ${TARGET.dir} && " + \
                                    "install ${SOURCE.dir}/* ${TARGET.dir}")
    env.Append(BUILDERS = {'Epydoc': builder, 'EpydocInstall':install})

    env.AppendUnique(EPYDOC='/usr/bin/epydoc')
    env.AppendUnique(EPYDOC_OPTS="--inheritance grouped --no-private " + \
                                 "--no-frames --html")

def exists(env):
    """Make sure epydoc tools exist."""
    return env.Detect("epydoc")
