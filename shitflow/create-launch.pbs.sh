#!/bin/bash
#PBS -P ox63
#PBS -q normal
#PBS -N generate
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/ox63+scratch/ox63+scratch/if89+gdata/if89
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l wd

set -o pipefail

usage() {
	echo "Usage: qsub -v FQ=/g/data/ox63/cornetto/data/gtg_internal/HG002/RGBX240039_HG002.hifi.fastq.gz,ASM=RGBX240039_HG002.hifiasm-hifi_1x ./create-launch.pbs.sh" >&2
	echo
	exit 1
}

die () {
	echo >&2 "$@"
	exit 1
}

[ -z "${FQ}" ] && usage
[ -z "${ASM}" ] && usage

test -f ${ASM}.fasta || die "Assembly FASTA not found"
test -f ${FQ} || die "Input FASTQ not found"

export MODULEPATH=$MODULEPATH:/g/data/if89/apps/modulefiles/
module load minimap2/2.24 || die "loading minimap2/2.24 module failed"
module load samtools/1.19 || die "loading samtools/1.19 module failed"
module load kentutils/0.0 || die "loading kentutils/0.0 module failed"
module load bedtools/2.28.0 || die "loading bedtools/2.28.0 module failed"

export CORNETTO_BIN=/g/data/ox63/hasindu/cornetto/cornetto/cornetto
export SF_SCRIPT_DIR=/g/data/ox63/hasindu/cornetto/cornetto/shitflow

minimap2 --version || die "Could not find minimap2"
samtools --version || die "Could not find samtools"
bedtools --version || die "Could not find bedtools"
which bedGraphToBigWig  || die "Could not find bedGraphToBigWig"
${CORNETTO_BIN} --version || die "Could not find cornetto"

THREADS=24
BAM=${ASM}.realigned.bam
CHROMBED=${ASM}.chroms.bed
CHROMSIZES=${ASM}.chromsizes.tsv
TOTCOV=${ASM}.cov-total
MQ20COV=${ASM}.cov-mq20

## index the assembly FASTA
samtools faidx ${ASM}.fasta || die "samtools faidx failed"
minimap2 -x map-ont -d ${ASM}.fasta.idx ${ASM}.fasta || die "minimap2 failed"

## generate CHROMBED and CHROMSIZES files
awk '{print $1"\t0\t"$2}' ${FASTA}.fai | sort -k3,3nr > ${CHROMBED} || die "awk failed"
cat ${CHROMBED} | awk '{print $1"\t"$3}' > ${CHROMSIZES} || die "awk failed"

## align starting hifi FASTQ reads back to the assembly they generated
minimap2 -t ${THREADS} --secondary=no --MD -ax map-ont ${ASM}.fasta ${FQ} > tmp.sam || die "minimap2 failed"
samtools view -@ ${THREADS} -b tmp.sam | samtools sort -@ ${THREADS} > ${BAM} || die "samtools sort failed"
samtools index ${BAM} || die "samtools index failed"

## get per base coverage tracks for total alignemts (MQ>=0) and unique alignemnts (MQ>=20)
samtools depth -@ ${THREADS} -b ${CHROMBED} -aa ${BAM} | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > ${TOTCOV}.bg || die "samtools depth failed"
samtools depth -@ ${THREADS} -Q 20 -b ${CHROMBED} -aa ${BAM} | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > ${MQ20COV}.bg || die "samtools depth failed"

## convert to bigwig format
bedGraphToBigWig ${TOTCOV}.bg ${CHROMSIZES} ${TOTCOV}.bw || die "bedGraphToBigWig failed"
bedGraphToBigWig ${MQ20COV}.bg ${CHROMSIZES} ${MQ20COV}.bw || die "bedGraphToBigWig failed"

qsub -v ASM=${ASM} ${SF_SCRIPT_DIR}/create-core.pbs.sh || die "create.pbs.sh failed"
