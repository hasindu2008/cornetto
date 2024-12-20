#!/bin/bash
#PBS -P ox63
#PBS -q normal
#PBS -N hifiasm
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/ox63+scratch/ox63+scratch/if89+gdata/if89+gdata/te53
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l wd

usage() {
	echo "Usage: qsub -v BASE_FASTQ=/path/to/RGBX240039_HG002.hifi.fastq.gz,FISH_PREV=A_1_QGXHXX240275:A_2_QGXHXX240279(optional),FISH_NOW=A_3_QGXHXX240283(optional),OUT_PREFIX=hg002-cornetto-A_3(optional) ./hifiasm.pbs.sh" >&2
	echo
	exit 1
}

#output prefix
[ -z "${OUT_PREFIX}" ] && OUT_PREFIX=hg002-cornetto
#pacbio fastq
[ -z "${BASE_FASTQ}" ] && usage

## inputs
HIFI_0=${BASE_FASTQ}
FASTQ_LIST=${HIFI_0}

ONT_DATADIR=/g/data/ox63/hasindu/cornetto/shitflow/
if [ -n "${FISH_PREV}" ]; then
	LIST=$(echo "$FISH_PREV" | tr ':' ' ' )
	for DUPLEX in ${LIST}; do
		DUP=${ONT_DATADIR}/${DUPLEX}/${DUPLEX}.duplex_reads.fastq
		FASTQ_LIST="${FASTQ_LIST} ${DUP}"
	done
fi

if [ -n "${FISH_NOW}" ]; then
	DUP=${ONT_DATADIR}/${FISH_NOW}/${FISH_NOW}.duplex_reads.fastq
	FASTQ_LIST="${FASTQ_LIST} ${DUP}"
fi

## outputs
ASM=${OUT_PREFIX}
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
module load minimap2/2.24
module load samtools/1.19
module load kentutils/0.0
module load quast/5.1.0rc1

GFATOOLS=/g/data/ox63/ira/adaptive_assembly/gfatools/gfatools
FLATTEN=/g/data/te53/ontsv/sv_parsing/scripts/flattenFasta.pl
REFERENCE=/g/data/ox63/cornetto/data/reference/hg002v1.0.1_pat.fa
GETSTAT_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/shitflow/getstat.pbs.sh
GENERATE_PANEL_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/shitflow/generate_panel.pbs.sh
HIFIASM=/g/data/ox63/install/hifiasm-0.22.0/hifiasm
THREADS=${PBS_NCPUS}

#########################################

echo "Running hifiasm with ${THREADS} threads, outprefix ${ASM}, and input fastq list ${FASTQ_LIST}" > hifiasm.log

## generate assembly with hifiasm
/usr/bin/time -v ${HIFIASM} --ont -t ${THREADS} --hg-size 3g -o ${ASM} ${FASTQ_LIST}
echo "hifiasm completed" >> hifiasm.log

## convert assembly graph to FASTA format
${GFATOOLS} gfa2fa ${ASM}.bp.p_ctg.gfa > ${ASM}.fasta
${GFATOOLS} gfa2fa ${ASM}.bp.hap1.p_ctg.gfa > ${ASM}.hap1.fasta
${GFATOOLS} gfa2fa ${ASM}.bp.hap2.p_ctg.gfa > ${ASM}.hap2.fasta
echo "gfa2fa completed" >> hifiasm.log

## index the assembly FASTA
samtools faidx ${ASM}.fasta
minimap2 -d ${ASM}.fasta.idx ${ASM}.fasta
echo "indexing completed" >> hifiasm.log

## generate CHROMBED and CHROMSIZES files
${FLATTEN} -tab ${ASM}.fasta | awk '{print $1"\t0\t"length($2)}' | sort -k3,3nr > ${CHROMBED}
cat ${CHROMBED} | awk '{print $1"\t"$3}' > ${CHROMSIZES}

## run quast to evaluate assemblies
ASMPATH=`realpath ${ASM}.fasta`
quast.py -t ${THREADS} -o ${ASM}.quast_out -l ${ASM} --large ${ASMPATH}
echo "quast completed" >> hifiasm.log

## run hasindu chromosomes stats script
qsub -v REF=${REFERENCE},ASM=${ASMPATH} ${GETSTAT_SCRIPT}
echo "getstat.pbs.sh submitted" >> hifiasm.log

## run generate panel script
qsub -v FISH_NOW=${FISH_NOW},PREFIX=${OUT_PREFIX} ${GENERATE_PANEL_SCRIPT}
echo "generate_panel.pbs.sh submitted" >> hifiasm.log