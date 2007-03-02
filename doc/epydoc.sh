#!/bin/sh
LD_LIBRARY_PATH=../src/ modpy.sh ~/bin/bin/epydoc --inheritance grouped --no-private --no-frames --html $* ../pyext/mdt.py
