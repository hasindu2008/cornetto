#!/bin/bash

set -o pipefail

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 1 ] || die "1 argument required, $# provided. Usage: create_cornetto.sh <assembly.fa>"

test -z ${CORNETTO} && CORNETTO=cornetto
${CORNETTO} --version || die "cornetto executable not found! Either put cornetto under path or set CORNETTO variable, e.g.,export CORNETTO=/path/to/cornetto"
test -z ${BEDTOOLS} && BEDTOOLS=bedtools
$BEDTOOLS --version > /dev/null 2>&1 || die "bedtools not found!. Either put bedtools under path or set BEDTOOLS variable, e.g.,export BEDTOOLS=/path/to/bedtools"

FASTA=$1
BGTOTAL=${FASTA%.fasta}.cov-total.bg
BGMQ20=${FASTA%.fasta}.cov-mq20.bg
LOWQ=${FASTA%.fasta}.bp.p_ctg.lowQ.bed

test -f ${FASTA} || die "File ${FASTA} not found"
test -f ${BGTOTAL} || die "File ${BGTOTAL} not found"
test -f ${BGMQ20} || die "File ${BGMQ20} not found"
test -f ${LOWQ} || die "File ${LOWQ} not found"

BASENAME=$(basename ${FASTA})
TMPOUT=tmp_create_cornetto
test -d ${TMPOUT} && die "Directory ${TMPOUT} already exists. Please remove it before running this script or change to a different working directory"
mkdir ${TMPOUT} || die "mkdir failed"

ASSBED=${TMPOUT}/${BASENAME}.bed
${CORNETTO} fa2bed ${FASTA} > ${ASSBED} || die "fa2bed failed"
# test -f ${FASTA}.fai || samtools faidx ${FASTA} || die "samtools faidx on ${FASTA} failed"
# awk '{print $1"\t0\t"$2}' ${FASTA}.fai > ${ASSBED} || die "awk failed"

#1# print all interesting windows with:
# low coverage: [<0.4x] genome average
# high coverage: [>2.5x] genome average
# low mappability: [mean MQ20 cov for window is < 0.4 x mean coverage for the window]
${CORNETTO} noboringbits -H 2.5 -L 0.4 -Q 0.4 ${BGTOTAL} -q ${BGMQ20} | awk '{if ($4!=".") print $1"\t"$2"\t"$3}' > ${TMPOUT}/1_tmp.bed || die "cornetto failed"

#2# merge these interesting windows - overlapping or adjacent (within 1000bp)
cat ${TMPOUT}/1_tmp.bed | sort -k1,1 -k2,2n | ${BEDTOOLS} merge -d 1000 > ${TMPOUT}/2_tmp.bed || die "bedtools merge failed"

#3# remove any merged intervals from (2) that are shorter than <30kb
awk '($3-$2)>=30000' ${TMPOUT}/2_tmp.bed > ${TMPOUT}/3_tmp.bed || die "awk failed"

#4# get regions longer than 8kb which were labelled by hifiasm as "low quality" (not sure what the definition of this is exactly)
cat ${LOWQ} | awk '($3-$2)>=8000' | cut -f 1-3 > ${TMPOUT}/lowQ_tmp.bed || die "awk failed"

#5# combine the funbits from (3) and (4) and extend them by +40kb in either direction
cat ${TMPOUT}/3_tmp.bed ${TMPOUT}/lowQ_tmp.bed | sort -k1,1 -k2,2n | awk '{if($2>40000){print $1"\t"$2-40000"\t"$3+40000} else {print $0}}' > ${TMPOUT}/funbits.bed || die "awk failed"

#6# add new windows to the file from (5), which are 200kb intervals at the left edge and right edge of the contig
awk '{if(($3-$2)>200000) {print $1"\t0\t200000\n"$1"\t"$3-200000"\t"$3}}' ${ASSBED} >> ${TMPOUT}/funbits.bed || die "awk failed"

#7# merge any overlapping or adjacent (within 200kb) intervals from (6)
${BEDTOOLS} sort -i ${TMPOUT}/funbits.bed | ${BEDTOOLS} merge -d 200000 > ${TMPOUT}/funbits_merged.bed || die "bedtools merge failed"

#8# subtract merged windows from (7) from the whole genome assembly
${BEDTOOLS} subtract -a ${ASSBED} -b ${TMPOUT}/funbits_merged.bed > ${TMPOUT}/boringbits_tmp.bed || die "bedtools subtract failed"

#9# subtract any contigs shorter than 800Mbase
awk '{if(($3-$2)<800000) print $0}' ${ASSBED} > ${TMPOUT}/short.bed || die "awk failed"
${BEDTOOLS} subtract -a ${TMPOUT}/boringbits_tmp.bed -b ${TMPOUT}/short.bed > ${TMPOUT}/boringbits.bed || die "bedtools subtract failed"

#10# if 'boring bits' are <50% of a single contig/scaffold, remove all boring bits on the whole scaffold.   Then create readfish targets
${CORNETTO} bigenough ${ASSBED} ${TMPOUT}/boringbits.bed -r ${BASENAME%.fasta}.boringbits.txt > ${BASENAME%.fasta}.boringbits.bed || die "cornetto bigenough failed"

# old crappy code to do the same thing. Remove after testing
#10# if 'boring bits' are <50% of a single contig/scaffold, remove all boring bits on the whole scaffold. Use the below horrible inefficient code snippet for now
# INPUT=${TMPOUT}/boringbits.bed
# cut -f 1 ${INPUT}  | uniq > ${TMPOUT}/boring_ctg.tmp || die "cut failed"
# while read p;
# do
# ctg_len=$(grep "$p" ${ASSBED} | cut -f 3)
# ctg_boring=$(grep "$p" ${INPUT} | awk '{sum+=$3-$2}END{ print sum}')
# fac=$(echo "$ctg_boring*100/$ctg_len" | bc)
# if [ "$fac" -gt "50" ];then
#     grep "$p" ${INPUT}
# fi
# done < ${TMPOUT}/boring_ctg.tmp > ${BASENAME%.fasta}.boringbits.bed || die "while loop failed"

#11# print the size of the boring_bits panel as a % of human genome size
# echo -n -e "${BASENAME%.fasta}\t"
# cat ${BASENAME%.fasta}.boringbits.bed | awk '{sum+=($3-$2)}END{print sum/3100000000*100}'

#12# create readfish targets
# cat ${BASENAME%.fasta}.boringbits.bed  | awk '{print $1","$2","$3",+"}' > ${TMPOUT}/plus_tmp
# cat ${BASENAME%.fasta}.boringbits.bed  | awk '{print $1","$2","$3",-"}' > ${TMPOUT}/minus_tmp
# cat ${TMPOUT}/plus_tmp ${TMPOUT}/minus_tmp | sort > ${BASENAME%.fasta}.boringbits.txt


