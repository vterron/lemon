#! /bin/bash

# Install LEMON dependencies in the Travis CI environment OS
# Author: Victor Terron (c) 2013
# License: GNU GPLv3

CWD=$(pwd)

cd ~

apt-get install git python-dev python-pip python-scipy python-pyfits \
        python-lxml python-uncertainties python-mock libfreetype6-dev \
        libpng12-dev csh libx11-dev libplplot11 alien realpath -qq

pip install numpy
pip install matplotlib
pip install astropy
pip install aplpy
pip install d2to1
pip install pyraf

########### Install IRAF ################

DOWNLOADS_SERVER="http://iaa.es/lemon/travis/"
IRAF_TAR="iraf.lnux.x86_64.tar.gz"
IRAF_URL=$DOWNLOADS_SERVER$IRAF_TAR

IRAF_DIR="/iraf/iraf"
mkdir -p $IRAF_DIR
IRAF_DIR=`realpath $IRAF_DIR`
cd $IRAF_DIR

wget $IRAF_URL
tar xvfz $IRAF_TAR

iraf=$IRAF_DIR
export $iraf
yes "" | $iraf/unix/hlib/install
rm $IRAF_TAR

########### Install CDSClient ###########

CDSCLIENT_TAR="cdsclient.tar.gz"

cd ~
wget $DOWNLOADS_SERVER$CDSCLIENT_TAR
tar xvfz $CDSCLIENT_TAR

cd cdsclient-* # e.g., cdsclient-3.72/
./configure
make
make install

########### Install SExtractor, SCAMP and SWarp ###########

cd ~

SEXTRACTOR_RPM="sextractor-2.8.6-1.x86_64.rpm"
SCAMP_RPM="scamp-1.7.0-1.x86_64.rpm"
SWARP_RPM="swarp-2.19.1-1.x86_64.rpm"

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
