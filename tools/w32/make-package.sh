#!/bin/sh

if [ ! -d tools/w32 ]; then
  echo "Run this script in the top-level mdt directory"
  exit 1
fi

if [ $# -ne 1 ]; then
  echo "Usage: tools/w32/make-package.sh modeller_dir"
  echo
  echo "modeller_dir should be the full path to a Windows Modeller install"
  exit 1
fi
modeller=$1
ROOT=w32-inst

rm -rf ${ROOT}

# Install MDT and support for all Python versions
for PYVER in 2.3 2.4 2.5 2.6 2.7 3.0 3.1 3.2; do
  cmd="scons destdir=${ROOT} modeller=$modeller wine=true libdir=/bin \
             datadir=/data pythondir=/python pyextdir=/python/python${PYVER} \
             pythoninclude=/usr/lib/w32comp/w32python/${PYVER}/include \
             libpath=/usr/lib/w32comp/w32python/${PYVER}/lib install tools/w32"
  echo $cmd
  $cmd || exit 1
  # Needed to build Python 3.2 extensions successfully with older SWIG
  if ! grep -q USE_CAPSULE pyext/mdt_wrap.c; then
    patch -f -d pyext -p1 < tools/w32/swig-pycapsule.patch || exit 1
  fi
done

# Note: no need to add MSVC runtime since Modeller already has it

# Modify Python search path to find version-specific modules
patch -f -d ${ROOT}/python/mdt -p1 < tools/w32/python-search-path.patch || exit 1

# Add README
cp tools/w32/README.txt ${ROOT}

tools/w32/gen-w32instlist ${ROOT} > w32files.tmp || exit 1
sed -e '/\.pyc"$/d' < w32files.tmp > w32files.install || exit 1
tac w32files.tmp | sed -e 's/File "w32-inst\\/Delete "$INSTDIR\\/' -e 's/^SetOutPath/RMDir/' > w32files.uninstall || exit 1
makensis -NOCD tools/w32/w32-install.nsi || exit 1

rm -rf w32files.tmp w32files.install w32files.uninstall ${ROOT}
