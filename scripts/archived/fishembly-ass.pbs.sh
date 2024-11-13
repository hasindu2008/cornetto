#!/bin/bash
#PBS -P ox63
#PBS -N fishsembly
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=128GB
#PBS -l walltime=6:00:00
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
	echo "Usage: qsub -v ASM=/path/to/asm.fasta,FASTQ=/path/to/duplex_reads.fastq ./fishembly.pbs.sh" >&2
	echo
	exit 1
}


#asm
[ -z "${ASM}" ] && usage
#fastq
[ -z "${FASTQ}" ] && usage

module load samtools/1.19
module load minimap2/2.24

###################################################################

set -x

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

minimap2 --version || die "Could not find minimap2"
samtools --version || die "Could not find samtools"

test -e ${FASTQ} || die "FASTQ file not found"
test -e ${ASM} || die "FASTA file not found"
test -e ${FASTQ}.fai || samtools fqidx ${FASTQ} || die "samtools failed to index FASTQ"

ASM_DIR=$(dirname ${ASM})
ASM_BASE=$(basename ${ASM} .fasta)
FASTQ_BASE=$(basename ${FASTQ} .fastq)

test -e ${ASM_DIR}/${ASM_BASE}/${ASM_BASE}.windows.0.4.50kb.ends.bed || die "ASM windows not found. please run getstat.sh first"

cat ${ASM_DIR}/${ASM_BASE}/${ASM_BASE}.windows.0.4.50kb.ends.bed | awk '{if($2==0) { print $1"\t0\t100000" } else if ($3-100000>100000 ) { print $1"\t"$3-100000"\t"$3 } }' > telocorners.txt || die "Could not extract telomere corners"


minimap2 -ax map-ont --secondary=no -t20 ${ASM} ${FASTQ} | samtools sort - -o ${FASTQ_BASE}.bam && samtools index ${FASTQ_BASE}.bam || die "Could not align reads"
samtools view ${FASTQ_BASE}.bam -L telocorners.txt | cut -f 1 >  t2t_rid.txt || die "Could not extract read ids"
sort -u t2t_rid.txt > t2t_rid_uniq.txt

cut -f 1 ${FASTQ}.fai > all_reads.txt
grep -v -F -f  t2t_rid_uniq.txt all_reads.txt > good.txt
samtools fqidx -r good.txt ${FASTQ}  > ${FASTQ_BASE}.clean.fastq

