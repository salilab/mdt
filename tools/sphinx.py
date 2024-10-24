# Simple Sphinx tool and builder.

from __future__ import print_function
import os


# Build sphinx documentation:
def _action_sphinx(target, source, env):
    sourcedir = os.path.dirname(source[1].path)
    outdir = os.path.dirname(target[0].path)
    bin = source[0].path
    app = "%s %s %s %s %s" % (bin, env['SPHINX_BUILD'], env['SPHINX_OPTS'],
                              sourcedir, outdir)
    ret = env.Execute(app)
    if not ret:
        print("Build finished. The HTML pages are in " + outdir)
    return ret


def generate(env):
    """Add builders and construction variables for the sphinx tool."""
    import SCons.Builder
    builder = SCons.Builder.Builder(action=_action_sphinx)
    # Use Unix 'install' rather than env.InstallAs(), due to scons bug #1751
    install = SCons.Builder.Builder(
        action="install -d ${TARGET.dir} && " +
               "install -d ${TARGET.dir}/_static && " +
               "install -d ${TARGET.dir}/_sources && " +
               "install ${SOURCE.dir}/*.html ${TARGET.dir} && " +
               "install ${SOURCE.dir}/*.js ${TARGET.dir} && " +
               "install ${SOURCE.dir}/_sources/* ${TARGET.dir}/_sources && " +
               "install ${SOURCE.dir}/_static/* ${TARGET.dir}/_static")
    env.Append(BUILDERS={'Sphinx': builder, 'SphinxInstall': install})

    env.AppendUnique(SPHINX_BUILD='/usr/bin/sphinx-build')
    env.AppendUnique(SPHINX_OPTS='-a -E -b html')


def exists(env):
    """Make sure sphinx tools exist."""
    return env.Detect("sphinx")
