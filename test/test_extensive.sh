#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

test/test.sh || die "test/test.sh failed"

test/test.sh mem || die "test/test.sh mem failed"

test/realtest.sh || die "test/realtest.sh failed"

# TODO not all are tested here. Need to add some small representative test cases to test.sh
test/realtest.sh  mem || die "test/realtest.sh mem failed"

