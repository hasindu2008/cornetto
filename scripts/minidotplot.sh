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

REF=$1
ASM=$2

PREFIX=$(basename ${ASM})

[ ! -f $REF ] && die "Reference $REF not found"
[ ! -f $ASM ] && die "Assembly $ASM not found"

test -z ${CORNETTO} && CORNETTO=cornetto
${CORNETTO} --version > /dev/null 2>&1 || die "cornetto executable not found! Either put cornetto under path or set CORNETTO variable, e.g.,export CORNETTO=/path/to/cornetto"

test -z ${MINIMAP2} && MINIMAP2=minimap2
$MINIMAP2 --version > /dev/null 2>&1 || die "minimap2 not found!. Either put minimap2 under path or set MINIMAP2 variable, e.g.,export MINIMAP2=/path/to/minimap2"
test -z ${SAMTOOLS} && SAMTOOLS=samtools
$SAMTOOLS --version > /dev/null 2>&1 || die "samtools not found!. Either put samtools under path or set SAMTOOLS variable, e.g.,export SAMTOOLS=/path/to/samtools"

${MINIMAP2} -t16 --eqx -cx asm5 $REF $ASM > ${PREFIX}.tmp.paf || die "minimap2 failed"
cut -f 1 ${PREFIX}.tmp.paf | sort -u > ${PREFIX}.tmp.ctg.list

${CORNETTO} fixdir ${ASM} ${PREFIX}.tmp.paf > ${PREFIX}.tmp.fix.fasta 2> ${PREFIX}.missing_sequences.log || die "cornetto failed"

grep '^>' ${PREFIX}.tmp.fix.fasta | sed 's/^>//' | while read p;
do
echo -n -e "${p}\t"
grep -w $p ${PREFIX}.tmp.paf | cut -f 6 | sort | uniq -c | sort -k1,1 -rn | head -1 | awk '{print $2}' | awk -F'_' '{print $1}'
done  > ${PREFIX}.chr.annotate.txt || die "grep failed"

awk '{count[$2]++; if (count[$2] > 1) {print $1"\t"$2 "_" count[$2]-1} else {print $1"\t"$2 "_0"}}' ${PREFIX}.chr.annotate.txt > ${PREFIX}.chr.rename.txt || die "awk failed"

awk '{print "s/"$1"/"$2"/g"}' ${PREFIX}.chr.rename.txt | sed -f - ${PREFIX}.tmp.fix.fasta > ${PREFIX}.tmp.renamed.fasta || die "sed failed"

$MINIMAP2 -t16 --eqx -cx asm5 $REF ${PREFIX}.tmp.renamed.fasta > ${PREFIX}.fix.tmp.paf || die "minimap2 failed"

${CORNETTO} minidot ${PREFIX}.fix.tmp.paf -f 2  > ${PREFIX}.eps || die "minidot failed"

echo "yey, all done"
