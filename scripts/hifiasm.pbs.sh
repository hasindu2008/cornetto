#!/bin/bash
#PBS -P ox63
#PBS -q express
#PBS -N hifiasm
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/ox63+scratch/ox63+scratch/if89+gdata/if89+gdata/te53
#PBS -l mem=180GB
#PBS -l jobfs=200GB
#PBS -l ncpus=24
#PBS -l wd
#PBS -M i.deveson@garvan.org.au

usage() {
	echo "Usage: qsub -v BASE_FASTQ=/path/to/RGBX240039_HG002.hifi.fastq.gz,OUT_PREFIX=hg002-cornetto-A_1(optional),DUPLEX1=A_1_QGXHXX240275(optional),DUPLEX2=.... ./getstat.pbs.sh" >&2
	echo
	exit 1
}

#output prefix
[ -z "${OUT_PREFIX}" ] && OUT_PREFIX=hg002-cornetto-A_1
#pacbio fastq
[ -z "${BASE_FASTQ}" ] && usage

## inputs
HIFI_0=${BASE_FASTQ}

ONT_DATADIR=/g/data/ox63/hasindu/cornetto/autocall/
DUP1=${ONT_DATADIR}/${DUPLEX1}/${DUPLEX1}.duplex_reads.fastq
DUP2=${ONT_DATADIR}/${DUPLEX2}/${DUPLEX2}.duplex_reads.fastq
DUP3=${ONT_DATADIR}/${DUPLEX3}/${DUPLEX3}.duplex_reads.fastq
DUP4=${ONT_DATADIR}/${DUPLEX4}/${DUPLEX4}.duplex_reads.fastq

FASTQ_LIST=${HIFI_0}
test -z "${DUP1}" || FASTQ_LIST="${FASTQ_LIST} ${DUP1}"
test -z "${DUP2}" || FASTQ_LIST="${FASTQ_LIST} ${DUP2}"
test -z "${DUP3}" || FASTQ_LIST="${FASTQ_LIST} ${DUP3}"
test -z "${DUP4}" || FASTQ_LIST="${FASTQ_LIST} ${DUP4}"


## outputs
ASM=${PREFIX}
CHROMBED=${ASM}.chroms.bed
CHROMSIZES=${ASM}.chromsizes.tsv

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

## load software
export MODULEPATH=$MODULEPATH:/g/data/if89/apps/modulefiles/
module load hifiasm/0.19.8
module load minimap2/2.24
module load samtools/1.19
module load kentutils/0.0
module load quast/5.1.0rc1

GFATOOLS=/g/data/ox63/ira/adaptive_assembly/gfatools/gfatools
FLATTEN=/g/data/te53/ontsv/sv_parsing/scripts/flattenFasta.pl
REFERENCE=/g/data/ox63/cornetto/data/reference/hg002v1.0.1_pat.fa
GETSTAT_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/scripts/getstat.pbs.sh
THREADS=${PBS_NCPUS}

#########################################

## generate assembly with hifiasm
/usr/bin/time -v hifiasm -t ${THREADS} --hg-size 3g -o ${ASM} ${FASTQ_LIST}

## convert assembly graph to FASTA format
${GFATOOLS} gfa2fa ${ASM}.bp.p_ctg.gfa > ${ASM}.fasta
${GFATOOLS} gfa2fa ${ASM}.bp.hap1.p_ctg.gfa > ${ASM}.hap1.fasta
${GFATOOLS} gfa2fa ${ASM}.bp.hap2.p_ctg.gfa > ${ASM}.hap2.fasta

## index the assembly FASTA
samtools faidx ${ASM}.fasta
minimap2 -d ${ASM}.fasta.idx ${ASM}.fasta

## generate CHROMBED and CHROMSIZES files
${FLATTEN} -tab ${ASM}.fasta | awk '{print $1"\t0\t"length($2)}' | sort -k3,3nr > ${CHROMBED}
cat ${CHROMBED} | awk '{print $1"\t"$3}' > ${CHROMSIZES}

## run quast to evaluate assemblies
ASMPATH=`realpath ${ASM}.fasta`
quast.py -t ${THREADS} -o ${ASM}.quast_out -l ${ASM} --large ${ASMPATH}

## run hasindu chromosomes stats script
qsub -v REF=${REFERENCE},ASM=${ASMPATH} ${GETSTAT_SCRIPT}


