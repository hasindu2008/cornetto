#!/bin/bash

die() {
    echo "$1"
    exit 1
}

[ $# -ne 2 ] && die "Usage: $0 <reference> <myassembly>"

MINIMAP2=minimap2
SAMTOOLS=samtools
MINIDOT=/install/miniasm/minidot

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

while read p;
do
    grep $p ${PREFIX}.tmp.paf | awk 'BEGIN{sump=0;sumn=0} {if($5=="-"){sumn+=($9-$8)}else{sump+=($9-$8)}} END{if(sump>sumn){print $1"\t+"}else{print $1"\t-"}}'
done < ${PREFIX}.tmp.ctg.list > ${PREFIX}.tmp.dir.txt

grep "+" ${PREFIX}.tmp.dir.txt | cut -f 1 > ${PREFIX}.tmp.ctg_plus.txt || die "grep failed"
grep "-" ${PREFIX}.tmp.dir.txt | cut -f 1 > ${PREFIX}.tmp.ctg_mins.txt || die "grep failed"

$SAMTOOLS faidx  ${PREFIX} -r ${PREFIX}.tmp.ctg_plus.txt > ${PREFIX}.tmp.fix.fasta || die "samtools failed"
$SAMTOOLS faidx  ${PREFIX} -r ${PREFIX}.tmp.ctg_mins.txt  -i >> ${PREFIX}.tmp.fix.fasta || die "samtools failed"

$MINIMAP2 -t16 --eqx -cx asm5 $REF ${PREFIX}.tmp.fix.fasta > ${PREFIX}.fix.tmp.paf || die "minimap2 failed"

${MINIDOT} ${PREFIX}.fix.tmp.paf -f 4  > ${PREFIX}.eps || die "minidot failed"

echo "yey, all done"
