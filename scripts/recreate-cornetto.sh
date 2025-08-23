#!/bin/bash

set -o pipefail

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 1 ] || die "1 argument required, $# provided. Usage: recreate_cornetto.sh <assembly.fa>"

FASTA=$1
test -f ${FASTA} || die "Assembly FASTA not found"

test -z ${CORNETTO} && CORNETTO=cornetto
${CORNETTO} --version || die "cornetto executable not found! Either put cornetto under path or set CORNETTO variable, e.g.,export CORNETTO=/path/to/cornetto"
test -z ${BEDTOOLS} && BEDTOOLS=bedtools
$BEDTOOLS --version > /dev/null 2>&1 || die "bedtools not found!. Either put bedtools under path or set BEDTOOLS variable, e.g.,export BEDTOOLS=/path/to/bedtools"

PREFIX=$(basename $FASTA .fa)
PREFIX=$(basename $FASTA .fasta)

TMPOUT=tmp_recreate_cornetto
test -d ${TMPOUT} && die "Directory ${TMPOUT} already exists. Please remove it before running this script or change to a different working directory"
mkdir ${TMPOUT} || die "mkdir failed"

## generate CHROMBED file
${CORNETTO} fa2bed ${FASTA} | sort -k3,3nr > ${TMPOUT}/${PREFIX}.chroms.bed || die "fa2bed failed"
# test -f ${FASTA}.fai || samtools faidx ${FASTA} || die "samtools faidx failed"
# awk '{print $1"\t0\t"$2}' ${FASTA}.fai | sort -k3,3nr > ${TMPOUT}/${PREFIX}.chroms.bed || die "awk failed"

#1# get regions longer than 7.5kb which were labelled by hifiasm as "low quality" (not sure what the definition of this is exactly)
cat ${PREFIX}.bp.p_ctg.lowQ.bed | awk '($3-$2)>=7500' | cut -f 1-3 > ${TMPOUT}/lowQ_tmp.bed || die "awk failed"

#2# extend the lowQ regions from (1) by +50kb in either direction
cat ${TMPOUT}/lowQ_tmp.bed | sort -k1,1 -k2,2n | awk '{if($2>50000){print $1"\t"$2-40000"\t"$3+50000} else {print $0}}' > ${TMPOUT}/funbits.bed || die "awk failed"

#3# add new windows to the file from (5), which are 200kb intervals at the left edge and right edge of the contig
awk '{if(($3-$2)>200000) {print $1"\t0\t200000\n"$1"\t"$3-200000"\t"$3}}' ${TMPOUT}/${PREFIX}.chroms.bed >> ${TMPOUT}/funbits.bed || die "awk failed"

#4# merge any overlapping or adjacent (within 200kb) intervals from (3)
${BEDTOOLS} sort -i ${TMPOUT}/funbits.bed | ${BEDTOOLS} merge -d 200000 > ${TMPOUT}/funbits_merged.bed || die "bedtools merge failed"

#5# subtract merged windows from (4) from the whole genome assembly
${BEDTOOLS} subtract -a ${TMPOUT}/${PREFIX}.chroms.bed -b ${TMPOUT}/funbits_merged.bed > ${TMPOUT}/boringbits_tmp.bed || die "bedtools subtract failed"

#6# subtract any contigs shorter than 1Mbase
awk '{if(($3-$2)<1000000) print $0}' ${TMPOUT}/${PREFIX}.chroms.bed > ${TMPOUT}/short.bed || die "awk failed"
${BEDTOOLS} subtract -a ${TMPOUT}/boringbits_tmp.bed -b ${TMPOUT}/short.bed > ${TMPOUT}/boringbits.bed || die "bedtools subtract failed"

#7# if 'boring bits' are <50% of a single contig/scaffold, remove all boring bits on the whole scaffold.  Then create readfish targets
${CORNETTO} bigenough ${TMPOUT}/${PREFIX}.chroms.bed ${TMPOUT}/boringbits.bed -r ${PREFIX}.boringbits.txt > ${PREFIX}.boringbits.bed || die "cornetto bigenough failed"

# old crappy code to do the same thing. Remove after testing
#7# if 'boring bits' are <50% of a single contig/scaffold, remove all boring bits on the whole scaffold. Use the below horrible inefficient code snippet for now
# INPUT=${TMPOUT}/boringbits.bed
# cut -f 1 ${INPUT}  | uniq > ${TMPOUT}/boring_ctg.tmp
# while read p;
# do
# 	ctg_len=$(grep "$p" ${TMPOUT}/${PREFIX}.chroms.bed | cut -f 3)
# 	ctg_boring=$(grep "$p" ${INPUT} | awk '{sum+=$3-$2}END{ print sum}')
# 	fac=$(echo "$ctg_boring*100/$ctg_len" | bc)
# 	if [ "$fac" -gt "50" ];then
# 		    grep "$p" ${INPUT}
# 	fi
# done < ${TMPOUT}/boring_ctg.tmp > ${PREFIX}.boringbits.bed

#8# print the size of the boring_bits panel as a % of human genome size
# cat ${PREFIX}.boringbits.bed | awk '{sum+=($3-$2)}END{print sum/3100000000*100}'

#9# create readfish targets
# cat ${PREFIX}.boringbits.bed | awk '{print $1","$2","$3",+"}' > ${TMPOUT}/plus_tmp
# cat ${PREFIX}.boringbits.bed | awk '{print $1","$2","$3",-"}' > ${TMPOUT}/minus_tmp
# cat ${TMPOUT}/plus_tmp ${TMPOUT}/minus_tmp | sort > ${PREFIX}.boringbits.txt
