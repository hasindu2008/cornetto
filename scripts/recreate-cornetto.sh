#!/bin/bash

set -o pipefail

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 1 ] || die "1 argument required, $# provided. Usage: recreate_cornetto.sh <assembly.fa>"

FASTA=$1
test -f ${FASTA} || die "Assembly FASTA not found"

BASENAME=$(basename ${FASTA})
PREFIX=${FASTA%.fasta}

TMPOUT=tmp_recreate_cornetto
test -d ${TMPOUT} && die "Directory ${TMPOUT} already exists. Please remove it before running this script or change to a different working directory"
mkdir ${TMPOUT} || die "mkdir failed"

## generate CHROMBED file
test -f ${FASTA}.fai || samtools faidx ${FASTA} || die "samtools faidx failed"
awk '{print $1"\t0\t"$2}' ${FASTA}.fai | sort -k3,3nr > ${TMPOUT}/${PREFIX}.chroms.bed || die "awk failed"

#1# get regions longer than 7.5kb which were labelled by hifiasm as "low quality" (not sure what the definition of this is exactly)
cat ${PREFIX}.bp.p_ctg.lowQ.bed | awk '($3-$2)>=7500' | cut -f 1-3 > ${TMPOUT}/lowQ_tmp.bed

#2# extend the lowQ regions from (1) by +50kb in either direction
cat ${TMPOUT}/lowQ_tmp.bed | sort -k1,1 -k2,2n | awk '{if($2>50000){print $1"\t"$2-40000"\t"$3+50000} else {print $0}}' > ${TMPOUT}/funbits.bed

#3# add new windows to the file from (5), which are 200kb intervals at the left edge and right edge of the contig
awk '{if(($3-$2)>200000) {print $1"\t0\t200000\n"$1"\t"$3-200000"\t"$3}}' ${TMPOUT}/${PREFIX}.chroms.bed >> ${TMPOUT}/funbits.bed

#4# merge any overlapping or adjacent (within 200kb) intervals from (3)
bedtools sort -i ${TMPOUT}/funbits.bed | bedtools merge -d 200000 > ${TMPOUT}/funbits_merged.bed

#5# subtract merged windows from (4) from the whole genome assembly
bedtools subtract -a ${TMPOUT}/${PREFIX}.chroms.bed -b ${TMPOUT}/funbits_merged.bed > ${TMPOUT}/boringbits_tmp.bed

#6# subtract any contigs shorter than 1Mbase
awk '{if(($3-$2)<1000000) print $0}' ${TMPOUT}/${PREFIX}.chroms.bed > ${TMPOUT}/short.bed
bedtools subtract -a ${TMPOUT}/boringbits_tmp.bed -b ${TMPOUT}/short.bed > ${TMPOUT}/boringbits.bed

#7# if 'boring bits' are <50% of a single contig/scaffold, remove all boring bits on the whole scaffold. Use the below horrible inefficient code snippet for now
## i.e. if the contig is more than 50% interesting, capture the whole thing
INPUT=${TMPOUT}/boringbits.bed
cut -f 1 ${INPUT}  | uniq > ${TMPOUT}/boring_ctg.tmp
while read p;
do
	ctg_len=$(grep "$p" ${TMPOUT}/${PREFIX}.chroms.bed | cut -f 3)
	ctg_boring=$(grep "$p" ${INPUT} | awk '{sum+=$3-$2}END{ print sum}')
	fac=$(echo "$ctg_boring*100/$ctg_len" | bc)
	if [ "$fac" -gt "50" ];then
		    grep "$p" ${INPUT}
	    fi
done < ${TMPOUT}/boring_ctg.tmp > ${PREFIX}.boringbits.bed

#8# print the size of the boring_bits panel as a % of human genome size
cat ${PREFIX}.boringbits.bed | awk '{sum+=($3-$2)}END{print sum/3100000000*100}'

#9# create readfish targets
cat ${PREFIX}.boringbits.bed | awk '{print $1","$2","$3",+"}' > ${TMPOUT}/plus_tmp
cat ${PREFIX}.boringbits.bed | awk '{print $1","$2","$3",-"}' > ${TMPOUT}/minus_tmp
cat ${TMPOUT}/plus_tmp ${TMPOUT}/minus_tmp | sort > ${PREFIX}.boringbits.txt
