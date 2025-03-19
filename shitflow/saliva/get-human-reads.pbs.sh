#!/bin/bash
#PBS -P ox63
#PBS -q normal
#PBS -N centrifuge
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/ox63+scratch/ox63+scratch/if89+gdata/if89+gdata/te53
#PBS -l mem=190GB
#PBS -l ncpus=48
#PBS -l wd

set -o pipefail

THREADS=48

usage() {
	echo "Usage: qsub -v FASTQ=/path/to/reads.fastq ./script.pbs.sh" >&2
	echo
	exit 1
}

BASENAME=$(basename ${FASTQ})
PREFIX=${BASENAME%%.*}
FQ=${PREFIX}_saliva.fastq

module load samtools/1.19

## run centrifuge on FASTQ
export CENTRIFUGE_HOME=/g/data/te53/ontsv/Figures/Review/microbial_contamination/centrifuge
INDEX_DIR=/g/data/te53/ontsv/Figures/Review/microbial_contamination/database/
REPORT=${PREFIX}.centrifuge_fastq_report.tsv
CLASSIFICATION=${PREFIX}.centrifuge_fastq_classification.tsv

$CENTRIFUGE_HOME/centrifuge -p ${THREADS} -q -x ${INDEX_DIR}/p_compressed+h+v -U ${FASTQ} -S ${CLASSIFICATION} --report-file ${REPORT}

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

test -e ${FQ} || cp ${FASTQ} ${FQ} || die "copying file failed"
test -e ${FQ}.fai || samtools faidx ${FQ} || die "samtools faidx failed"
cat ${CLASSIFICATION} | awk '$3!=9606' | cut -f 1 | sort -u > ${PREFIX}_nonhuman_reads.txt
cut -f 1 ${FQ}.fai > ${PREFIX}_all_reads.txt || die "cut failed"
grep -v -F -f ${PREFIX}_nonhuman_reads.txt ${PREFIX}_all_reads.txt > ${PREFIX}_human_reads.txt || die "grep failed"
samtools fqidx -r ${PREFIX}_human_reads.txt ${FQ} > ${PREFIX}_human_reads.fastq || die "samtools fqidx failed"
echo "non human removal success"

