#!/bin/bash

set -o pipefail

die(){
    echo $1
    exit 1
}

[ $# -ne 1 ] && die "Usage: $0 <file>"

FILE=$1

test -f $FILE || die "File $FILE not found"
PREFIX=$(basename $FILE .fa)
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

cornetto telomere --patterns $FILE | awk '{print $1"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' - > $PREFIX.telomere
cornetto sdust $FILE > $PREFIX.sdust

test -f ${PREFIX}.fa.fai || samtools faidx ${PREFIX}.fa || die "samtools faidx on ${PREFIX}.fa failed"
awk '{print $1"\t"$2}' ${PREFIX}.fa.fai > ${PREFIX}.lens  || die "awk failed"

cornetto telomere --windows $PREFIX.telomere 99.9 0.1 > $PREFIX.windows
cornetto telomere --breaks $PREFIX.lens $PREFIX.sdust $PREFIX.telomere > $PREFIX.breaks
cornetto telomere --windows $PREFIX.telomere 99.9 $THRESHOLD > $PREFIX.windows.$THRESHOLD

echo "Merge telomere motifs in 100bp"
cat $PREFIX.windows.$THRESHOLD | awk '{print $2"\t"$(NF-2)"\t"$(NF-1)}' | sed 's/>//g' | bedtools merge -d 100  > $PREFIX.windows.$THRESHOLD.bed
echo

echo "Find those at end of scaffolds, within < $ENDS"
cat $PREFIX.lens | awk -v ends=$ENDS '{if ($2>(ends*2)) {print $1"\t0\t"ends"\n"$1"\t"($NF-ends)"\t"$NF} else {print $1"\t0\t"$NF}}' > asm.ends.bed

ENDS=`echo $ENDS | awk '{printf "%.0f", $1/1000}'`"kb"
bedtools intersect -wa -a $PREFIX.windows.$THRESHOLD.bed -b asm.ends.bed > $PREFIX.windows.$THRESHOLD.$ENDS.ends.bed

test -e ${BED} || die "telomere_analysis.sh failed"

echo -e "FILE\t${FILE}"
echo -e -n "total telomere regions at the end of contigs:\t"
cat $BED | wc -l
echo ""
echo ""
cat $BED | cut -f 1  | sort | uniq -c | awk 'BEGIN{t1=0;t2=0;t3=0}{if($1==1){t1+=1}else if($1==2){t2+=1} else {t3+=1}} END{print "telo in one end:\t"t1"\ntelo in two ends:\t"t2"\ntelo more than 2 (must be 0):\t"t3"\n"}'





