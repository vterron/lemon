#! /bin/bash

# Install LEMON dependencies in the Travis CI environment OS
# Author: Victor Terron (c) 2013
# License: GNU GPLv3

set -e # exit if any statement returns a non-true return value
set -u # any attempt to use an undefined variable is an error
set -x # print commands and their arguments as they are executed

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

apt-get install git python-pip csh realpath
apt-get build-dep python-matplotlib python-scipy
easy_install -U distribute

pip install "numpy>=1.7.1"
pip install -r pre-requirements.txt
pip install -r requirements.txt

CWD=$(pwd)

cd ~

########### Install IRAF ################

if [[ $ARCH_64_BITS == 1 ]]; then
    IRAF_TAR="iraf.lnux.x86_64.tar.gz"
else
    IRAF_TAR="iraf.lnux.x86.tar.gz"
fi
IRAF_SERVER="ftp://iraf.noao.edu/iraf/v216/PCIX/"
IRAF_URL=$IRAF_SERVER$IRAF_TAR

IRAF_DIR="/iraf/iraf/"
mkdir -p $IRAF_DIR
IRAF_DIR=`realpath $IRAF_DIR`
cd $IRAF_DIR

wget $IRAF_URL
tar xfz $IRAF_TAR

iraf=$IRAF_DIR
export iraf
yes "" | $iraf/unix/hlib/install
rm $IRAF_TAR

########### Install SExtractor ###########

cd ~

if [[ $ARCH_64_BITS == 1 ]]; then
    SEXTRACTOR_RPM="sextractor-2.19.5-1.x86_64.rpm"
else
    SEXTRACTOR_RMP="sextractor-2.19.5-1.i386.rpm"
fi

SEXTRACTOR_SERVER="http://www.astromatic.net/download/sextractor/"
SEXTRACTOR_URL=$SEXTRACTOR_SERVER$SEXTRACTOR_RPM
wget $SEXTRACTOR_URL
alien -i $SEXTRACTOR_RPM
rm $SEXTRACTOR_RPM

cd $CWD # back to the LEMON directory

# The unit tests use several FITS images that are downloaded from the
# STScI Digitized Sky Survey to test/test_data/fits/. Be considerate
# and, instead of downloading them every time the tests are run, keep
# a copy on our server.

TEST_FITS_DIR="test/test_data/fits/"
DSS_IMAGES_URL="https://github.com/lemon/test-data/raw/master/DSS/"

DSS_FILENAMES=(
 "Barnard's_Star.fits"
 "IC_5070.fits"
 "IC_5146.fits"
 "Messier_92.fits"
 "NGC_2264.fits"
 "Orion.fits"
 "RMC_136.fits"
 "Serpens.fits"
 "Trapezium.fits"
 "Trumpler_37.fits"
)

mkdir -p $TEST_FITS_DIR
cd $TEST_FITS_DIR

echo "Downloading test FITS images to $(pwd)"
for filename in "${DSS_FILENAMES[@]}"; do
    wget $DSS_IMAGES_URL$filename -O $filename;
done;

cd $CWD

exit 0
