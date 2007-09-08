#!/bin/sh
prefix=$1
shift
LD_LIBRARY_PATH=../src/ ${prefix}modpy.sh ~/bin/bin/epydoc -o html/epydoc --inheritance grouped --no-private --no-frames --html $* ../pyext/mdt.py
