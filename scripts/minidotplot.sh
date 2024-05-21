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

cat ${PREFIX}.tmp.ctg_plus.txt ${PREFIX}.tmp.ctg_mins.txt  | while read p;
do
echo -n -e "${p}\t"
grep $p ${PREFIX}.tmp.paf | cut -f 6 | sort | uniq -c | sort -k1,1 -rn | head -1 | awk '{print $2}' | awk -F'_' '{print $1}'
done  > ${PREFIX}.chr.annotate.txt || die "grep failed"

awk '{count[$2]++; if (count[$2] > 1) {print $1, $2 "_" count[$2]-1} else {print $1, $2}}' ${PREFIX}.chr.annotate.txt > ${PREFIX}.chr.rename.txt || die "awk failed"

awk '{print "s/"$1"/"$2"/g"}' ${PREFIX}.chr.rename.txt | sed -f - ${PREFIX}.tmp.fix.fasta > ${PREFIX}.tmp.renamed.fasta || die "sed failed"

$MINIMAP2 -t16 --eqx -cx asm5 $REF ${PREFIX}.tmp.renamed.fasta > ${PREFIX}.fix.tmp.paf || die "minimap2 failed"

${MINIDOT} ${PREFIX}.fix.tmp.paf -f 4  > ${PREFIX}.eps || die "minidot failed"

echo "yey, all done"
