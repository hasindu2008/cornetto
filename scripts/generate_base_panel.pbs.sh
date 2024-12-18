#!/bin/bash
#PBS -P ox63
#PBS -q normal
#PBS -N generate
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/ox63+scratch/ox63+scratch/if89+gdata/if89
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l wd

usage() {
	echo "Usage: qsub -v FQ=/g/data/ox63/cornetto/data/gtg_internal/HG002/RGBX240039_HG002.hifi.fastq.gz,ASM=RGBX240039_HG002.hifiasm-hifi_1x ./generate_base_panel.pbs.sh" >&2
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
module load minimap2/2.24
module load samtools/1.19
module load kentutils/0.0
module load bedtools/2.28.0
export PATH=$PATH:/g/data/ox63/ira/scripts/
FLATTEN=/g/data/ox63/ira/scripts/flattenFasta.pl
export CORNETTO_BIN=/g/data/ox63/hasindu/cornetto/cornetto/cornetto
export SCRIPT_DIR=/g/data/ox63/hasindu/cornetto/cornetto/scripts

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
${FLATTEN} -tab ${ASM}.fasta | awk '{print $1"\t0\t"length($2)}' | sort -k3,3nr > ${CHROMBED} || die "flattenFasta.pl failed"
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

${SCRIPT_DIR}/create_cornetto4.sh ${ASM}.fasta|| die "create_cornetto4.sh failed"
mkdir diploid || die "mkdir failed"
cd diploid || die "cd failed"
${SCRIPT_DIR}/create_hapnetto.sh ${ASM} || die "create_hapnetto.sh failed"

echo "all done. f ya."