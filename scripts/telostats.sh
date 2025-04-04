#!/bin/bash

set -o pipefail

die(){
    echo $1
    exit 1
}

[ $# -ne 1 ] && die "Usage: $0 <file>"

test -z ${CORNETTO} && CORNETTO=cornetto
${CORNETTO} --version || die "cornetto executable not found! Either put cornetto under path or set CORNETTO variable, e.g.,export CORNETTO=/path/to/cornetto"
test -z ${BEDTOOLS} && BEDTOOLS=bedtools
$BEDTOOLS --version > /dev/null 2>&1 || die "bedtools not found!. Either put bedtools under path or set BEDTOOLS variable, e.g.,export BEDTOOLS=/path/to/bedtools"

FILE=$1

test -f $FILE || die "File $FILE not found"
PREFIX=$(basename $FILE .fa)
PREFIX=$(basename $PREFIX .fasta)
BED=${PREFIX}.windows.0.4.50kb.ends.bed

THRESHOLD=0.4
ENDS=50000

mkdir -p $PREFIX
cd $PREFIX

echo "genome: $PREFIX"
echo "THRESHOLD: $THRESHOLD"
echo "ends: $ENDS"
echo "asm: $FILE"

ln -s $FILE 2> /dev/null

${CORNETTO} teltelofind $FILE | awk '{print $1"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' - > $PREFIX.telomere || die "cornetto teltelofind failed"
${CORNETTO} sdust $FILE > $PREFIX.sdust || die "cornetto sdust failed"

test -f ${FILE}.fai || samtools faidx ${FILE} || die "samtools faidx on ${FILE} failed"
awk '{print $1"\t"$2}' ${FILE}.fai > ${PREFIX}.lens  || die "awk failed"

${CORNETTO} telowin --windows $PREFIX.telomere 99.9 0.1 > $PREFIX.windows || die "cornetto telowin failed"
${CORNETTO} telobreak --breaks $PREFIX.lens $PREFIX.sdust $PREFIX.telomere > $PREFIX.breaks || die "cornetto telobreak failed"
${CORNETTO} telowin --windows $PREFIX.telomere 99.9 $THRESHOLD > $PREFIX.windows.$THRESHOLD || die "cornetto telowin failed"

echo "Merge telomere motifs in 100bp"
cat $PREFIX.windows.$THRESHOLD | awk '{print $2"\t"$(NF-2)"\t"$(NF-1)}' | sed 's/>//g' | ${BEDTOOLS} merge -d 100  > $PREFIX.windows.$THRESHOLD.bed || die "bedtools merge failed"
echo

echo "Find those at end of scaffolds, within < $ENDS"
cat $PREFIX.lens | awk -v ends=$ENDS '{if ($2>(ends*2)) {print $1"\t0\t"ends"\n"$1"\t"($NF-ends)"\t"$NF} else {print $1"\t0\t"$NF}}' > asm.ends.bed || die "awk failed"

ENDS=`echo $ENDS | awk '{printf "%.0f", $1/1000}'`"kb" || die "awk failed"
${BEDTOOLS} intersect -wa -a $PREFIX.windows.$THRESHOLD.bed -b asm.ends.bed > $PREFIX.windows.$THRESHOLD.$ENDS.ends.bed || die "bedtools intersect failed"

test -e ${BED} || die "telomere_analysis.sh failed"

echo -e "FILE\t${FILE}"
echo -e -n "total telomere regions at the end of contigs:\t"
cat $BED | wc -l
echo ""
echo ""
cat $BED | cut -f 1  | sort | uniq -c | awk 'BEGIN{t1=0;t2=0;t3=0}{if($1==1){t1+=1}else if($1==2){t2+=1} else {t3+=1}} END{print "telo in one end:\t"t1"\ntelo in two ends:\t"t2"\ntelo more than 2 (must be 0):\t"t3"\n"}' || die "awk failed"





