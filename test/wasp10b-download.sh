#! /bin/bash

set -Eeuo pipefail

# Takes the path to a directory as argument and downloads the WASP10b test data
# if any of the expected files is missing. In this manner, the .xz file is only
# downloaded from the remote server if necessary.

# Name of the file stored in $WASP10_URL.
XZ_FILE="WASP10b-2011-08-23.tar.xz"

pushd $1

if sha1sum -c SHA1SUMS; then
    # Cache archives expire after 45 days for repositories on https://travis-ci.com.
    # [https://docs.travis-ci.com/user/caching/#caches-expiration] Write the current
    # date and time to a file so that (because of this modification inside the cache
    # directory) the expiration delay is reset.
    date > last-access.txt

    popd
    exit 0
fi

# Don't use rm *; explicitly list the files to remove.
rm -vrf *.fits SHA1SUMS $XZ_FILE
curl --remote-header-name --remote-name $WASP10_URL
tar -xvf $XZ_FILE
rm -fv $XZ_FILE
chmod a-w -v *.fits

popd
