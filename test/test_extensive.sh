#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

test/test.sh

test/test.sh mem

test/realtest.sh

# TODO not all are tested here. Need to add some small representative test cases to test.sh
test/realtest.sh  mem


echo "Tests passed"