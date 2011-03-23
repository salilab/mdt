#!/bin/sh

# Script to get SVN HEAD and make a .tar.gz on a shared disk, for use by
# autobuild scripts on our build machines.
#
# Should be run on the SVN server (or some other system that can access SVN
# readonly without a password) from a crontab, e.g.
#
# 10 1 * * * /cowbell1/home/ben/impmod/tools/auto-build.sh
#
# In this case, the build user is in the apache group, and so has readonly
# access to the repositories.

VER=SVN
IMPSVNDIR=file:///cowbell1/svn/impmod/trunk/

TMPDIR=/var/tmp/modeller-build-$$
MODINSTALL=/salilab/diva1/home/modeller/.${VER}-new
IMPSRCTGZ=${MODINSTALL}/build/sources/private/impmod.tar.gz

rm -rf ${TMPDIR}
mkdir ${TMPDIR}
cd ${TMPDIR}

# Get top-most revision number (must be a nicer way of doing this?)
rev=$(svn log -q --limit 1 ${IMPSVNDIR} |grep '^r' | cut -f 1 -d' ')

# Get IMP module code from SVN
svn export -q -${rev} ${IMPSVNDIR} impmod

# Write out a version file
verfile="${MODINSTALL}/build/impmod-version"
echo "${rev}" > $verfile

# Write out a tarball:
tar -czf ${IMPSRCTGZ} impmod

# Cleanup
cd /
rm -rf ${TMPDIR}
