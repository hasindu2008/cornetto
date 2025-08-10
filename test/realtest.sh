#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

PREFIX=hg002-cornetto-E_2
export CORNETTO=./cornetto

echo "fixasm test"

${CORNETTO} fixasm test/real/E_2/${PREFIX}.fasta test/real/E_2/${PREFIX}.fasta.tmp.paf -r b.tsv -m b.txt -w b.paf  > b.fasta || die "fixasm failed running"
diff -q b.tsv test/real/E_2/fixasm/report.tsv || die "fixasm output mismatch in TSV report"
diff -q b.txt test/real/E_2/fixasm/missing.txt || die "fixasm output mismatch in missing contigs"
diff -q b.paf test/real/E_2/fixasm/fixed.paf || die "fixasm output mismatch in paf file"
diff -q b.fasta test/real/E_2/fixasm/fixed.fasta || die "fixasm output mismatch in fasta file"

echo "telostats test"
test -d ${PREFIX} && rm -r ${PREFIX}
scripts/telostats.sh test/real/E_2/${PREFIX}.fasta > telostats.txt || die "telostats failed"
diff -q ${PREFIX}/${PREFIX}.telomere test/real/E_2/telostats/${PREFIX}.telomere || die "telostats output mismatch in telomere file"
diff -q ${PREFIX}/${PREFIX}.lens test/real/E_2/telostats/${PREFIX}.lens || die "telostats output mismatch in lens file"
diff -q ${PREFIX}/${PREFIX}.windows.0.4 test/real/E_2/telostats/${PREFIX}.windows.0.4 || die "telostats output mismatch in windows.0.4 file"
diff -q ${PREFIX}/${PREFIX}.windows.0.4.50kb.ends.bed test/real/E_2/telostats/${PREFIX}.windows.0.4.50kb.ends.bed || die "telostats output mismatch in windows.0.4.50kb.ends.bed file"
diff -q telostats.txt test/real/E_2/telostats/telostats.txt || die "telostats output mismatch in telostat.txt file"

echo "extra telo tests"
${CORNETTO} sdust test/real/E_2/${PREFIX}.fasta > $PREFIX/$PREFIX.sdust || die "cornetto sdust failed"
${CORNETTO} telowin $PREFIX/$PREFIX.telomere 99.9 0.1 > $PREFIX/$PREFIX.windows || die "cornetto telowin failed"
${CORNETTO} telobreaks $PREFIX/$PREFIX.lens $PREFIX/$PREFIX.sdust $PREFIX/$PREFIX.telomere > $PREFIX/$PREFIX.breaks || die "cornetto telobreak failed"
diff -q ${PREFIX}/${PREFIX}.windows test/real/E_2/telostats/${PREFIX}.windows || die "telostats output mismatch in windows file"
diff -q ${PREFIX}/${PREFIX}.breaks test/real/E_2/telostats/${PREFIX}.breaks || die "telostats output mismatch in breaks file"

echo "asmstats test"
./cornetto asmstats test/real/E_2/${PREFIX}.fasta.tmp.paf  test/real/E_2/telostats/${PREFIX}.windows.0.4.50kb.ends.bed -r test/real/E_2/fixasm/report.tsv

echo "Tests passed"