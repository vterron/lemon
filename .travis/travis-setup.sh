#! /bin/bash

# Install LEMON dependencies in the Travis CI environment OS
# Author: Victor Terron (c) 2013
# License: GNU GPLv3

set -e # exit if any statement returns a non-true return value

REGEXP=".*64$"
HARDWARE_NAME=`uname -m`

echo -n "Machine architecture: "
if [[ $HARDWARE_NAME =~ $REGEXP ]]; then
    ARCH_64_BITS=1;
    echo "64-bit"
else
    ARCH_64_BITS=0;
    echo "32-bit"
fi

# Use apt-get, instead of pip, in order to install as many Python dependencies
# as possible: some packages, like NumPY and SciPy, contain a lot of C code,
# which would otherwise have to be compiled every time the Travis tests are run.
# The binary packages allow us to avoid this, resulting in an enormous speed-up.

apt-get install python-dev python-pip python-numpy python-scipy \
                python-matplotlib python-lxml python-mock libfreetype6-dev \
                libpng12-dev csh libx11-dev libplplot11 alien realpath

pip install astropy>=0.2.4
pip install d2to1>=0.2.10
pip install APLpy>=0.9.9
pip install pyfits>=3.1.2
pip install pyraf>=2.1.1
pip install uncertainties>=2.4.1

CWD=$(pwd)

cd ~

########### Install IRAF ################

DOWNLOADS_SERVER="http://iaa.es/lemon/travis/"
if [[ $ARCH_64_BITS == 1 ]]; then
    IRAF_TAR="iraf.lnux.x86_64.tar.gz"
else
    IRAF_TAR="iraf.lnux.x86.tar.gz"
fi
IRAF_URL=$DOWNLOADS_SERVER$IRAF_TAR

IRAF_DIR="/iraf/iraf"
mkdir -p $IRAF_DIR
IRAF_DIR=`realpath $IRAF_DIR`
cd $IRAF_DIR

wget $IRAF_URL
tar xfz $IRAF_TAR

iraf=$IRAF_DIR
export $iraf
yes "" | $iraf/unix/hlib/install
rm $IRAF_TAR

########### Install CDSClient ###########

CDSCLIENT_TAR="cdsclient.tar.gz"

cd ~
wget $DOWNLOADS_SERVER$CDSCLIENT_TAR
tar xfz $CDSCLIENT_TAR

cd cdsclient-* # e.g., cdsclient-3.72/
./configure
make
make install

########### Install SExtractor, SCAMP and SWarp ###########

cd ~

if [[ $ARCH_64_BITS == 1 ]]; then
    SEXTRACTOR_RPM="sextractor-2.8.6-1.x86_64.rpm"
    SCAMP_RPM="scamp-1.7.0-1.x86_64.rpm"
    SWARP_RPM="swarp-2.19.1-1.x86_64.rpm"
else
    SEXTRACTOR_RPM="sextractor-2.8.6-1.i386.rpm"
    SCAMP_RPM="scamp-1.7.0-1.i386.rpm"
    SWARP_RPM="swarp-2.19.1-1.i386.rpm"
fi

SEXTRACTOR_URL=$DOWNLOADS_SERVER$SEXTRACTOR_RPM
SCAMP_URL=$DOWNLOADS_SERVER$SCAMP_RPM
SWARP_URL=$DOWNLOADS_SERVER$SWARP_RPM

wget $SEXTRACTOR_URL
alien -i $SEXTRACTOR_RPM
rm $SEXTRACTOR_RPM

wget $SCAMP_URL
alien -i $SCAMP_RPM
rm $SCAMP_RPM

wget $SWARP_URL
alien -i $SWARP_RPM
rm $SWARP_RPM

cd $PWD

exit 0
