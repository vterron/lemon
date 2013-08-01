#! /bin/bash

# Install LEMON dependencies in the Travis CI environment OS
# Author: Victor Terron (c) 2013
# License: GNU GPLv3

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

apt-get install python-dev python-pip libfreetype6-dev libpng12-dev \
                csh libx11-dev libplplot11 alien realpath

pip install -I numpy>=1.7.1
pip install -Ir pre-requirements.txt
pip install -Ir requirements.txt

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
