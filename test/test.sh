#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

if [ "$1" = 'mem' ]; then
    mem=1
else
    mem=0
fi

ex() {
    if [ $mem -eq 1 ]; then
        valgrind --leak-check=full --error-exitcode=1 "$@"
    else
        "$@"
    fi
}

echo "Test 1"
ex  ./cornetto boringbits test/cov-total.bg -q test/cov-mq20.bg -m 10000 -e 1000 > test/tmp.txt  || die "Running the tool failed"
diff -q test/example_boring_t1.exp test/tmp.txt || die "diff failed"

echo "Test 2"
ex  ./cornetto funbits -H 2.5 -L 0.5 -Q 0.5 test/cov-total.bg -q test/cov-mq20.bg -m 10000 -e 1000 > test/tmp.txt  || die "Running the tool failed"
diff -q test/example_fun_t2.exp test/tmp.txt || die "diff failed"

echo "Tests passed"