#!/bin/bash
#PBS -P ox63
#PBS -q normal
#PBS -N centrifuge
#PBS -l walltime=48:00:00
#PBS -l storage=gdata/ox63+scratch/ox63+scratch/if89+gdata/if89+gdata/te53
#PBS -l mem=190GB
#PBS -l ncpus=48
#PBS -l wd

THREADS=48

usage() {
	echo "Usage: qsub -v ASM_PREFIX=F_0_RGBX240267_saliva,FASTQ_PREFIX=/path/to/reads ./xcript.pbs.sh" >&2
	echo
	exit 1
}
set -o pipefail

module load samtools/1.19
ASM=${ASM_PREFIX}
FASTQ=${FASTQ_PREFIX}

CHROMBED=${ASM}.chroms.bed
CHROMSIZES=${ASM}.chromsizes.tsv

## run centrifuge on assembly
export CENTRIFUGE_HOME=/g/data/te53/ontsv/Figures/Review/microbial_contamination/centrifuge
INDEX_DIR=/g/data/te53/ontsv/Figures/Review/microbial_contamination/database/
REPORT=${ASM}.centrifuge_report.tsv
CLASSIFICATION=${ASM}.centrifuge_classification.tsv

$CENTRIFUGE_HOME/centrifuge -p ${THREADS} -f -x ${INDEX_DIR}/p_compressed+h+v -U ${ASM}.fasta -S ${CLASSIFICATION} --report-file ${REPORT} || die "Failed to run centrifuge"

export PATH=$PATH:/g/data/te53/ontsv/sv_parsing/scripts/
FLATTEN=/g/data/te53/ontsv/sv_parsing/scripts/flattenFasta.pl
FETCH=/g/data/te53/ontsv/sv_parsing/scripts/fetchSubset.pl

## centrifuge files
FA_REPORT=${ASM}.centrifuge_report.tsv
FA_CLASS=${ASM}.centrifuge_classification.tsv
FQ_REPORT=${FASTQ}.centrifuge_fastq_report.tsv
FQ_CLASS=${FASTQ}.centrifuge_fastq_classification.tsv

## identify nonhuman species with minimum 100 reads
cat ${FQ_REPORT} | sed 's/ /-/g' | sort -k5,5nr | awk '$2 != 9606' | awk '$5 >= 100' | cut -f 2 | sort -u | awk '$1 != "taxID"' > ${ASM}_nonhuman_species_high_count.txt || die "Failed to identify nonhuman species"

## find the contigs that correspond to nonhuman species those nonhuman species
${FETCH} ${FA_CLASS} ${ASM}_nonhuman_species_high_count.txt 3 1 | cut -f 1 | sort -u > ${ASM}_nonhuman_species_high_count_contig_ids.txt || die "Failed to fetch nonhuman contigs"

## fetch the contig sequences for those species and create a FASTA file just those contigs
${FLATTEN} -tab ${ASM}.fasta > ${ASM}_tmp.tab || die "Failed to flatten ${ASM}.fasta"
${FETCH} ${ASM}_tmp.tab ${ASM}_nonhuman_species_high_count_contig_ids.txt 1 1 > ${ASM}_nonhuman_contig_seqs.tab || die "Failed to fetch nonhuman contigs"
${FLATTEN} -fa ${ASM}_nonhuman_contig_seqs.tab > ${ASM}.nonhuman_contigs.fasta || die "Failed to create ${ASM}.nonhuman_contigs.fasta"

## create a BED file covering full contigs for each nonhuman contig to be excluded during readfish
cat nonhuman_contig_seqs.tab | awk '{print $1"\t0\t"length($2)}' > ${ASM}.nonhuman_contigs.bed || die "Failed to create ${ASM}.nonhuman_contigs.bed"

