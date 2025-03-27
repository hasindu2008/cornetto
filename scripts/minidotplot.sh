#!/bin/bash

set -o pipefail

die() {
    echo "$1"
    exit 1
}

[ $# -ne 2 ] && die "Usage: $0 <reference> <myassembly>"

MINIMAP2=minimap2
SAMTOOLS=samtools
CORNETTO=cornetto
test -z ${MINIDOT} && MINIDOT=/install/miniasm/minidot

REF=$1
ASM=$2

PREFIX=$(basename ${ASM})

[ ! -f $REF ] && die "Reference $REF not found"
[ ! -f $ASM ] && die "Assembly $ASM not found"

$MINIMAP2 --version > /dev/null 2>&1 || die "minimap2 not found"
$SAMTOOLS --version > /dev/null 2>&1 || die "samtools not found"
test -e $MINIDOT  > /dev/null 2>&1 || die "minidot not found"

${MINIMAP2} -t16 --eqx -cx asm5 $REF $ASM > ${PREFIX}.tmp.paf || die "minimap2 failed"
cut -f 1 ${PREFIX}.tmp.paf | sort -u > ${PREFIX}.tmp.ctg.list

${CORNETTO} fixdir ${ASM} ${PREFIX}.tmp.fix.fasta || die "cornetto failed" # TODO:direct fixdir output to stdout instead of corrected_contigs.fasta

grep '^>' corrected_contigs.fasta | sed 's/^>//' | while read p;
do
echo -n -e "${p}\t"
grep -w $p ${PREFIX}.tmp.paf | cut -f 6 | sort | uniq -c | sort -k1,1 -rn | head -1 | awk '{print $2}' | awk -F'_' '{print $1}'
done  > ${PREFIX}.chr.annotate.txt || die "grep failed"

awk '{count[$2]++; if (count[$2] > 1) {print $1"\t"$2 "_" count[$2]-1} else {print $1"\t"$2 "_0"}}' ${PREFIX}.chr.annotate.txt > ${PREFIX}.chr.rename.txt || die "awk failed"

awk '{print "s/"$1"/"$2"/g"}' ${PREFIX}.chr.rename.txt | sed -f - corrected_contigs.fasta > ${PREFIX}.tmp.renamed.fasta || die "sed failed"

$MINIMAP2 -t16 --eqx -cx asm5 $REF ${PREFIX}.tmp.renamed.fasta > ${PREFIX}.fix.tmp.paf || die "minimap2 failed"

${MINIDOT} ${PREFIX}.fix.tmp.paf -f 2  > ${PREFIX}.eps || die "minidot failed"

echo "yey, all done"
