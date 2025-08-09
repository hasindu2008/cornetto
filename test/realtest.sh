#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

echo "fixasm test"

./cornetto fixasm test/real/E_2/in/hg002-cornetto-E_2.fasta test/real/E_2/in/hg002-cornetto-E_2.fasta.tmp.paf -r b.tsv -m b.txt -w b.paf  > b.fasta || die "fixasm failed running"

diff -q b.tsv test/real/E_2/exp/fixasm/report.tsv || die "fixasm output mismatch in TSV report"
diff -q b.txt test/real/E_2/exp/fixasm/missing.txt || die "fixasm output mismatch in missing contigs"
diff -q b.paf test/real/E_2/exp/fixasm/fixed.paf || die "fixasm output mismatch in paf file"
diff -q b.fasta test/real/E_2/exp/fixasm/fixed.fasta || die "fixasm output mismatch in fasta file"

echo "Tests passed"