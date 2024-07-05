#!/bin/bash
#PBS -P ox63
#PBS -N getstat
#PBS -q normal
#PBS -l ncpus=4
#PBS -l mem=16GB
#PBS -l walltime=48:00:00
#PBS -l wd
#PBS -l storage=gdata/if89+scratch/wv19+gdata/wv19+gdata/ox63

###################################################################

###################################################################

# Make sure to change:
# 1. wv19 and ox63 to your own projects

# to run:
# qsub -v FASTQ=/path/to/duplex_reads.fastq,BAM=/path/to/duplex_reads.hg002v1.0.1_pat.bam ./fishembly.pbs.sh

###################################################################

CHR_LIST="chr3_PATERNAL chr6_PATERNAL chr11_PATERNAL chr12_PATERNAL  chr18_PATERNAL chr20_PATERNAL"
OUTPUT="duplex_reads_good.fastq"

usage() {
	echo "Usage: qsub -v FASTQ=/path/to/duplex_reads.fastq,BAM=/path/to/duplex_reads.hg002v1.0.1_pat.bam ./fishembly.pbs.sh" >&2
	echo
	exit 1
}


#asm
[ -z "${FASTQ}" ] && usage
#ref
[ -z "${BAM}" ] && usage

module load samtools/1.19

###################################################################

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

samtools --version || die "Could not find samtools"

test -e ${FASTQ} || die "FASTQ file not found"
test -e ${BAM} || die "BAM file not found"

test -e ${FASTQ}.fai || samtools faidx ${FASTQ} || die "samtools faidx failed"
test -e ${BAM}.bai || samtools index ${BAM} || die "samtools index failed"

test -e t2t_rid.txt && rm t2t_rid.txt

for chr in ${CHR_LIST}; do
	samtools view ${BAM} ${chr} | cut -f 1 >> t2t_rid.txt || die "samtools view failed"
done

sort -u t2t_rid.txt > t2t_rid_uniq.txt || die "sort failed"
cut -f 1 ${FASTQ}.fai > all_reads.txt || die "cut failed"
grep -v -F -f  t2t_rid_uniq.txt all_reads.txt > good.txt || die "grep failed"
samtools fqidx -r good.txt ${FASTQ}  > ${OUTPUT} || die "samtools fqidx failed"
