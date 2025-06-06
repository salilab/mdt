# vim: set et sw=4 tw=0 fo=awqorc ft=python:
#
# Astxx, the Asterisk C++ API and Utility Library.
# Copyright (C) 2005, 2006  Matthew A. Nicholson
# Copyright (C) 2006  Tim Blechmann
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License version 2.1 as published by the Free Software Foundation.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import os
import os.path
import glob
from fnmatch import fnmatch

def DoxyfileParse(file_contents):
    """
    Parse a Doxygen source file and return a dictionary of all the values.
    Values will be strings and lists of strings.
    """
    data = {}

    import shlex
    lex = shlex.shlex(instream = file_contents, posix = True)
    lex.wordchars += "*+./-:"
    lex.whitespace = lex.whitespace.replace("\n", "")
    lex.escape = ""

    lineno = lex.lineno
    token = lex.get_token()
    key = token   # the first token should be a key
    last_token = ""
    key_token = False
    next_key = False
    new_data = True

    def append_data(data, key, new_data, token):
        if new_data or len(data[key]) == 0:
            data[key].append(token)
        else:
            data[key][-1] += token

    while token:
        if token in ['\n']:
            if last_token not in ['\\']:
                key_token = True
        elif token in ['\\']:
            pass
        elif key_token:
            key = token
            key_token = False
        else:
            if token == "+=":
                if not data.has_key(key):
                    data[key] = list()
            elif token == "=":
                if key == "TAGFILES" and data.has_key(key):
                    append_data( data, key, False, "=" )
                    new_data=False
                else:
                    data[key] = list()
            else:
                append_data( data, key, new_data, token )
                new_data = True

        last_token = token
        token = lex.get_token()

        if last_token == '\\' and token != '\n':
            new_data = False
            append_data( data, key, new_data, '\\' )

    # compress lists of len 1 into single strings
    for (k, v) in data.items():
        if len(v) == 0:
            data.pop(k)

        # items in the following list will be kept as lists and not converted to strings
        if k in ["INPUT", "FILE_PATTERNS", "EXCLUDE_PATTERNS", "TAGFILES"]:
            continue

        if len(v) == 1:
            data[k] = v[0]

    return data

def DoxySourceScan(node, env, path):
    """
    Doxygen Doxyfile source scanner.  This should scan the Doxygen file and add
    any files used to generate docs to the list of source files.
    """
    default_file_patterns = [
       '*.c', '*.cc', '*.cxx', '*.cpp', '*.c++', '*.java', '*.ii', '*.ixx',
       '*.ipp', '*.i++', '*.inl', '*.h', '*.hh ', '*.hxx', '*.hpp', '*.h++',
       '*.idl', '*.odl', '*.cs', '*.php', '*.php3', '*.inc', '*.m', '*.mm',
       '*.py',
    ]

    default_exclude_patterns = [
       '*~',
    ]

    sources = []

    data = DoxyfileParse(node.get_contents())

    if data.get("RECURSIVE", "NO") == "YES":
        recursive = True
    else:
        recursive = False

    file_patterns = data.get("FILE_PATTERNS", default_file_patterns)
    exclude_patterns = data.get("EXCLUDE_PATTERNS", default_exclude_patterns)

    # We're running in the top-level directory, but the doxygen
    # configuration file is in the same directory as node; this means
    # that relative pathnames in node must be adjusted before they can
    # go onto the sources list
    conf_dir = os.path.dirname(str(node))

    for node in data.get("INPUT", []):
        if not os.path.isabs(node):
            node = os.path.join(conf_dir, node)
        if os.path.isfile(node):
            sources.append(node)
        elif os.path.isdir(node):
            if recursive:
                for root, dirs, files in os.walk(node):
                    for f in files:
                        filename = os.path.join(root, f)

                        pattern_check = reduce(lambda x, y: x or bool(fnmatch(filename, y)), file_patterns, False)
                        exclude_check = reduce(lambda x, y: x and fnmatch(filename, y), exclude_patterns, True)

                        if pattern_check and not exclude_check:
                            sources.append(filename)
            else:
                for pattern in file_patterns:
                    sources.extend(glob.glob("/".join([node, pattern])))

    # Add tagfiles to the list of source files:
    for node in data.get("TAGFILES", []):
        file = node.split("=")[0]
        if not os.path.isabs(file):
            file = os.path.join(conf_dir, file)
        sources.append(file)

    # Add additional files to the list of source files:
    def append_additional_source(option):
        file = data.get(option, "")
        if file != "":
            if not os.path.isabs(file):
                file = os.path.join(conf_dir, file)
            if os.path.isfile(file):
                sources.append(file)

    append_additional_source("HTML_STYLESHEET")
    append_additional_source("HTML_HEADER")
    append_additional_source("HTML_FOOTER")

    sources = map( lambda path: env.File(path), sources )
    return sources


def DoxySourceScanCheck(node, env):
    """Check if we should scan this file"""
    return os.path.isfile(node.path)

def DoxyEmitter(source, target, env):
    """Doxygen Doxyfile emitter"""
    # possible output formats and their default values and output locations
    output_formats = {
       "HTML": ("YES", "html", "index.html"),
       "LATEX": ("YES", "latex", None),
       "RTF": ("NO", "rtf", None),
       "MAN": ("YES", "man", None),
       "XML": ("NO", "xml", None),
    }

    data = DoxyfileParse(source[0].get_contents())

    targets = []
    out_dir = data.get("OUTPUT_DIRECTORY", ".")
    if not os.path.isabs(out_dir):
        conf_dir = os.path.dirname(str(source[0]))
        out_dir = os.path.join(conf_dir, out_dir)

    # add our output locations
    for (k, v) in output_formats.items():
        if data.get("GENERATE_" + k, v[0]) == "YES":
            dirname = os.path.join(out_dir, data.get(k + "_OUTPUT", v[1]))
            dir = env.Dir(dirname)
            # Since scons currently doesn't handle directory targets well,
            # work around this by defining a single known "always installed"
            # file (e.g. index.html) as the target:
            if v[2]:
                f = env.File(os.path.join(dirname, v[2]))
                targets.append(f)
                env.Clean(f, dir)
            else:
                targets.append(dir)

    # add the tag file if necessary:
    tagfile = data.get("GENERATE_TAGFILE", "")
    if tagfile != "":
        if not os.path.isabs(tagfile):
            conf_dir = os.path.dirname(str(source[0]))
            tagfile = os.path.join(conf_dir, tagfile)
        targets.append(env.File(tagfile))

    # don't clobber targets
    for node in targets:
        env.Precious(node)

    return (targets, source)

def generate(env):
    """
    Add builders and construction variables for the
    Doxygen tool.  This is currently for Doxygen 1.4.6.
    """
    doxyfile_scanner = env.Scanner(
       DoxySourceScan,
       "DoxySourceScan",
       scan_check = DoxySourceScanCheck,
    )

    import SCons.Builder
    doxyfile_builder = SCons.Builder.Builder(
       action="cd ${SOURCE.dir}  &&  ${DOXYGEN} ${SOURCE.file}",
       emitter=DoxyEmitter,
       target_factory=env.fs.Entry,
       single_source=True,
       source_scanner=doxyfile_scanner,
    )

    env.Append(BUILDERS={'Doxygen': doxyfile_builder,})

    env.AppendUnique(DOXYGEN='doxygen',)


def exists(env):
    """
    Make sure doxygen exists.
    """
    return env.Detect("doxygen")
