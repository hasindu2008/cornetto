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
ex  ./cornetto boringbits test/cov-total.bg -q test/cov-mq20.bg -m 10000 -e 1000 -L 0.6 -Q 0.6 -H 1.6 > test/tmp.txt  || die "Running the tool failed"
diff -q test/example_boring_t1.exp test/tmp.txt || die "diff failed"

echo "Test 2"
ex  ./cornetto noboringbits -H 2.5 -L 0.5 -Q 0.5 test/cov-total.bg -q test/cov-mq20.bg -m 10000 -e 1000 > test/tmp.txt  || die "Running the tool failed"
diff -q test/example_fun_t2.exp test/tmp.txt || die "diff failed"

echo "bigenough test"
./cornetto bigenough test/bigenough/hg002-cornetto-E_3/chroms.bed test/bigenough/hg002-cornetto-E_3/in.boringbits.bed -r a.txt > a.bed || die "Running the tool failed"
diff -q test/bigenough/hg002-cornetto-E_3/out.boringbits.bed a.bed || die "diff failed"
diff -q test/bigenough/hg002-cornetto-E_3/out.boringbits.csv a.txt || die "diff failed"

./cornetto bigenough test/bigenough/hg002-cornetto-E_3/chroms.bed test/bigenough/hg002-cornetto-E_3/in_dip.boringbits.bed -r a.txt > a.bed || die "Running the tool failed"
diff -q test/bigenough/hg002-cornetto-E_3/out_dip.boringbits.bed a.bed || die "diff failed"
diff -q test/bigenough/hg002-cornetto-E_3/out_dip.boringbits.csv a.txt || die "diff failed"

echo "Tests passed"