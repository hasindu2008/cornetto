#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 1 ] || die "1 argument required, $# provided. Usage: create_cornetto4.sh <assembly.fa>"

# dont run in parallel in same directory

FASTA=$1
ASSBED=${FASTA}.bed
BGTOTAL=${FASTA%.fasta}.cov-total.bg
BGMQ20=${FASTA%.fasta}.cov-mq20.bg
LOWQ=${FASTA%.fasta}.bp.p_ctg.lowQ.bed


test -f ${FASTA} || die "File ${FASTA} not found"
test -f ${BGTOTAL} || die "File ${BGTOTAL} not found"
test -f ${BGMQ20} || die "File ${BGMQ20} not found"
test -f ${LOWQ} || die "File ${LOWQ} not found"

samtools faidx ${FASTA} || die "samtools faidx failed"
awk '{print $1"\t0\t"$2}' ${FASTA}.fai > ${ASSBED} || die "awk failed"

#1# print all interesting windows with:
# low coverage: [<0.4x] genome average
# high coverage: [>2.5x] genome average
# low mappability: [mean MQ20 cov for window is < 0.4 x mean coverage for the window]

/home/hasindu/hasindu2008.git/cornetto/cornetto funbits -H 2.5 -L 0.4 -Q 0.4 ${BGTOTAL} -q ${BGMQ20} | awk '{if ($4!=".") print $1"\t"$2"\t"$3}' > 1_tmp.bed || die "cornetto failed"

#2# merge these interesting windows - overlapping or adjacent (within 1000bp)
cat 1_tmp.bed | sort -k1,1 -k2,2n | bedtools merge -d 1000 > 2_tmp.bed || die "bedtools merge failed"

#3# remove any merged intervals from (2) that are shorter than <30kb
awk '($3-$2)>=30000' 2_tmp.bed > 3_tmp.bed || die "awk failed"

#4# get regions longer than 8kb which were labelled by hifiasm as "low quality" (not sure what the definition of this is exactly)
cat ${LOWQ} | awk '($3-$2)>=8000' | cut -f 1-3 > lowQ_tmp.bed || die "awk failed"

#5# combine the funbits from (3) and (4) and extend them by +40kb in either direction
cat 3_tmp.bed lowQ_tmp.bed | sort -k1,1 -k2,2n | awk '{if($2>40000){print $1"\t"$2-40000"\t"$3+40000} else {print $0}}' > funbits.bed || die "awk failed"

#6# add new windows to the file from (5), which are 200kb intervals at the left edge and right edge of the contig
awk '{if(($3-$2)>200000) {print $1"\t0\t200000\n"$1"\t"$3-200000"\t"$3}}' ${ASSBED} >> funbits.bed || die "awk failed"

#7# merge any overlapping or adjacent (within 200kb) intervals from (6)
bedtools sort -i funbits.bed | bedtools merge -d 200000 > funbits_merged.bed || die "bedtools merge failed"

#8# subtract merged windows from (7) from the whole genome assembly
bedtools subtract -a ${ASSBED} -b funbits_merged.bed > boringbits_tmp.bed || die "bedtools subtract failed"

#9# subtract any contigs shorter than 800Mbase
awk '{if(($3-$2)<800000) print $0}' ${ASSBED} > short.bed || die "awk failed"
bedtools subtract -a boringbits_tmp.bed -b short.bed > boringbits.bed || die "bedtools subtract failed"

#10# if 'boring bits' are <50% of a single contig/scaffold, remove all boring bits on the whole scaffold. Use the below horrible inefficient code snippet for now
## i.e. if the contig is more than 50% interesting, capture the whole thing
INPUT=boringbits.bed
cut -f 1 ${INPUT}  | uniq > boring_ctg.tmp || die "cut failed"
while read p;
do
ctg_len=$(grep "$p" ${ASSBED} | cut -f 3)
ctg_boring=$(grep "$p" ${INPUT} | awk '{sum+=$3-$2}END{ print sum}')
fac=$(echo "$ctg_boring*100/$ctg_len" | bc)
if [ "$fac" -gt "50" ];then
    grep "$p" ${INPUT}
fi
done < boring_ctg.tmp > ${FASTA%.fasta}.boringbits.bed || die "while loop failed"

#11# print the size of the boring_bits panel as a % of human genome size
echo -n -e "${FASTA}\t"
cat ${FASTA%.fasta}.boringbits.bed | awk '{sum+=($3-$2)}END{print sum/3100000000*100}'




