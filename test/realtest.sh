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


PREFIX=hg002-cornetto-E_2
export CORNETTO=./cornetto

echo "fixasm test"

ex ${CORNETTO} fixasm test/real/E_2/${PREFIX}.fasta test/real/E_2/${PREFIX}.fasta.tmp.paf -r ${PREFIX}.report.tsv -m b.txt -w b.paf --trim-pat-mat > b.fasta || die "fixasm failed running"
diff -q ${PREFIX}.report.tsv test/real/E_2/fixasm/report.tsv || die "fixasm output mismatch in TSV report"
diff -q b.txt test/real/E_2/fixasm/missing.txt || die "fixasm output mismatch in missing contigs"
diff -q b.paf test/real/E_2/fixasm/fixed.paf || die "fixasm output mismatch in paf file"
diff -q b.fasta test/real/E_2/fixasm/fixed.fasta || die "fixasm output mismatch in fasta file"
rm ${PREFIX}.report.tsv b.txt b.paf b.fasta

echo "telostats test"
test -d ${PREFIX} && rm -r ${PREFIX}
scripts/telostats.sh test/real/E_2/${PREFIX}.fasta > telostats.txt || die "telostats failed"
TEMPDIR=tmp_${PREFIX}_telostats
diff -q ${TEMPDIR}/${PREFIX}.telomere test/real/E_2/telostats/${PREFIX}.telomere || die "telostats output mismatch in telomere file"
diff -q ${TEMPDIR}/${PREFIX}.lens test/real/E_2/telostats/${PREFIX}.lens || die "telostats output mismatch in lens file"
diff -q ${TEMPDIR}/${PREFIX}.windows.0.4 test/real/E_2/telostats/${PREFIX}.windows.0.4 || die "telostats output mismatch in windows.0.4 file"
diff -q ${PREFIX}.windows.0.4.50kb.ends.bed test/real/E_2/telostats/${PREFIX}.windows.0.4.50kb.ends.bed || die "telostats output mismatch in windows.0.4.50kb.ends.bed file"
diff -q telostats.txt test/real/E_2/telostats/telostats.txt || die "telostats output mismatch in telostat.txt file"
rm telostats.txt

echo "extra telo tests"
${CORNETTO} sdust test/real/E_2/${PREFIX}.fasta > $TEMPDIR/$PREFIX.sdust || die "cornetto sdust failed"
${CORNETTO} telowin $TEMPDIR/$PREFIX.telomere 99.9 0.1 > $TEMPDIR/$PREFIX.windows || die "cornetto telowin failed"
${CORNETTO} telobreaks $TEMPDIR/$PREFIX.lens $TEMPDIR/$PREFIX.sdust $TEMPDIR/$PREFIX.telomere > $TEMPDIR/$PREFIX.breaks || die "cornetto telobreak failed"
diff -q ${TEMPDIR}/${PREFIX}.windows test/real/E_2/telostats/${PREFIX}.windows || die "telostats output mismatch in windows file"
diff -q ${TEMPDIR}/${PREFIX}.breaks test/real/E_2/telostats/${PREFIX}.breaks || die "telostats output mismatch in breaks file"
rm -r ${TEMPDIR}

echo "asmstats test"
ex ./cornetto asmstats test/real/E_2/${PREFIX}.fasta.tmp.paf  test/real/E_2/telostats/${PREFIX}.windows.0.4.50kb.ends.bed -r test/real/E_2/fixasm/report.tsv -s human1 --trim-pat-mat > asmstats.txt
diff -q asmstats.txt test/real/E_2/asmstats/asmstats.txt || die "asmstats output mismatch"
rm asmstats.txt

echo "Tests passed"