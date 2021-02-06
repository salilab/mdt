#!/bin/bash -e

# Install dependencies for Travis or GitHub Actions.

if [ $# -ne 1 ]; then
  echo "Usage: $0 python_version"
  exit 1
fi

python_version=$1
MODELLER_VERSION=10.0
HDF5_VERSION=1.10.6
HDF5_SOVER=103
HDF5_HL_SOVER=100

sudo apt-get -qq update
sudo apt-get install -y swig bc
pip install coverage pytest-cov pytest-flake8

# Install Modeller
modeller_url=https://salilab.org/modeller/${MODELLER_VERSION}
wget ${modeller_url}/modeller_${MODELLER_VERSION}-1_amd64.deb
sudo --preserve-env=KEY_MODELLER dpkg -i modeller_${MODELLER_VERSION}-1_amd64.deb

# Modeller installs for system Python, so link its Python packages into the
# virtualenv Python path
export PYTHON=`pip show coverage |grep Location|cut -b11-`
ln -sf /usr/lib/python${python_version}/dist-packages/_modeller.so ${PYTHON}
ln -sf /usr/lib/python${python_version}/dist-packages/modeller ${PYTHON}

# Need to build MDT against the same version of HDF5 as Modeller, but Modeller
# only includes the HDF5 libraries, so we need to get the corresponding
# headers. These were built on an Ubuntu machine with
# ./configure --enable-shared --enable-hl --prefix=/tmp/hdf5 \
#        && make -j8 && make install
# and then the /tmp/hdf5/include directory archived.
wget https://salilab.org/modeller/downloads/hdf5-${HDF5_VERSION}-headers-64bit-ubuntu.tar.gz
tar -xzf hdf5-${HDF5_VERSION}-headers-64bit-ubuntu.tar.gz && mv hdf5-${HDF5_VERSION}-headers-64bit-ubuntu hdf5-inc

# Make .so symlinks for each .so.{ver} HDF5 library
mkdir hdf5-libs
ln -sf /usr/lib/modeller${MODELLER_VERSION}/lib/x86_64-intel8/libhdf5_hl.so.${HDF5_HL_SOVER} hdf5-libs/libhdf5_hl.so
ln -sf /usr/lib/modeller${MODELLER_VERSION}/lib/x86_64-intel8/libhdf5.so.${HDF5_SOVER} hdf5-libs/libhdf5.so
