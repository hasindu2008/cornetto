#!/bin/bash

die() {
    echo "$1" >&2
    echo
    exit 1
}

tar zcvf cornetto-realdata.tar.gz test/real/ || die "Creating tarball failed"
mv cornetto-realdata.tar.gz "/mnt/c/Users/hasin/OneDrive - UNSW/filehost/cornetto_data/" || die "Moving tarball failed"