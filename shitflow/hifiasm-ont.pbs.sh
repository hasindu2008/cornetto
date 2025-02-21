#!/bin/bash
#PBS -P ox63
#PBS -q hugemem
#PBS -N hifiasm
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/ox63+scratch/ox63+scratch/if89+gdata/if89+gdata/te53
#PBS -l mem=1470GB
#PBS -l ncpus=48
#PBS -l wd

usage() {
	echo "Usage: qsub -v BASE_FASTQ=/path/to/A0_XX.fastq,FISH_PREV=A_1_QGXHXX240275:A_2_QGXHXX240279(optional),FISH_NOW=A_3_QGXHXX240283(optional),OUT_PREFIX=hg002-cornetto-A_3(optional) ./hifiasm.pbs.sh" >&2
	echo
	exit 1
}

#output prefix
[ -z "${OUT_PREFIX}" ] && OUT_PREFIX=hg002-cornetto
#ont base fastq
[ -z "${BASE_FASTQ}" ] && usage

## inputs
ONT_0=${BASE_FASTQ}
FASTQ_LIST=${ONT_0}

ONT_DATADIR=/g/data/ox63/hasindu/cornetto/shitflow/
if [ -n "${FISH_PREV}" ]; then
	LIST=$(echo "$FISH_PREV" | tr ':' ' ' )
	for DUPLEX in ${LIST}; do
		DUP=${ONT_DATADIR}/${DUPLEX}/${DUPLEX}.fastq
		FASTQ_LIST="${FASTQ_LIST} ${DUP}"
	done
fi

if [ -n "${FISH_NOW}" ]; then
	DUP=${ONT_DATADIR}/${FISH_NOW}/${FISH_NOW}.fastq
	FASTQ_LIST="${FASTQ_LIST} ${DUP}"
fi

## outputs
ASM=${OUT_PREFIX}
#CHROMBED=${ASM}.chroms.bed
#CHROMSIZES=${ASM}.chromsizes.tsv

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

## load software
export MODULEPATH=$MODULEPATH:/g/data/if89/apps/modulefiles/
module load minimap2/2.24 || die "loading minimap2/2.24 module failed"
module load samtools/1.19 || die "loading samtools/1.19 module failed"
module load quast/5.1.0rc1 || die "loading quast/5.1.0rc1 module failed"

#FLATTEN=/g/data/te53/ontsv/sv_parsing/scripts/flattenFasta.pl
REFERENCE=/g/data/ox63/cornetto/data/reference/hg002v1.0.1_pat.fa

GETSTAT_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/shitflow/getstat.pbs.sh
CREATE_PANEL_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/shitflow/create-launch.pbs.sh
RECREATE_PANEL_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/shitflow/recreate.pbs.sh
QUAST_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/shitflow/quast.pbs.sh

HIFIASM=/g/data/ox63/install/hifiasm-0.22.0/hifiasm
GFATOOLS=/g/data/ox63/ira/adaptive_assembly/gfatools/gfatools

THREADS=${PBS_NCPUS}

#########################################

echo "Running hifiasm with ${THREADS} threads, outprefix ${ASM}, and input fastq list ${FASTQ_LIST}" > hifiasm.log

## generate assembly with hifiasm
/usr/bin/time -v ${HIFIASM} --ont -t ${THREADS} --hg-size 3g -o ${ASM} ${FASTQ_LIST} || die "hifiasm failed"
echo "hifiasm completed" >> hifiasm.log

## convert assembly graph to FASTA format
${GFATOOLS} gfa2fa ${ASM}.bp.p_ctg.gfa > ${ASM}.fasta || die "gfa2fa failed"
${GFATOOLS} gfa2fa ${ASM}.bp.hap1.p_ctg.gfa > ${ASM}.hap1.fasta || die "gfa2fa failed"
${GFATOOLS} gfa2fa ${ASM}.bp.hap2.p_ctg.gfa > ${ASM}.hap2.fasta || die "gfa2fa failed"
echo "gfa2fa completed" >> hifiasm.log

# ## index the assembly FASTA
# samtools faidx ${ASM}.fasta
# minimap2 -d ${ASM}.fasta.idx ${ASM}.fasta
# echo "indexing completed" >> hifiasm.log

# ## generate CHROMBED and CHROMSIZES files
# ${FLATTEN} -tab ${ASM}.fasta | awk '{print $1"\t0\t"length($2)}' | sort -k3,3nr > ${CHROMBED}
# cat ${CHROMBED} | awk '{print $1"\t"$3}' > ${CHROMSIZES}

## run quast to evaluate assemblies
ASMPATH=`realpath ${ASM}.fasta`
qsub -v ASM=${ASMPATH},OUT_DIR=${ASM}.quast_out ${QUAST_SCRIPT} || die "quast submission failed"
echo "quast.pbs.sh submitted for primary" >> hifiasm.log

HAP1_PATH=`realpath ${ASM}.hap1.fasta`
qsub -v ASM=${HAP1_PATH},OUT_DIR=${ASM}.hap1.quast_out ${QUAST_SCRIPT} || die "quast submission failed"
echo "quast.pbs.sh submitted for hap 1" >> hifiasm.log

HAP2_PATH=`realpath ${ASM}.hap2.fasta`
qsub -v ASM=${HAP2_PATH},OUT_DIR=${ASM}.hap2.quast_out ${QUAST_SCRIPT} || die "quast submission failed"
echo "quast.pbs.sh submitted for hap 2" >> hifiasm.log

## run hasindu chromosomes stats script
qsub -v REF=${REFERENCE},ASM=${ASMPATH} ${GETSTAT_SCRIPT} || die "getstat submission failed"
echo "getstat.pbs.sh submitted for primary" >> hifiasm.log

qsub -v REF=${REFERENCE},ASM=${HAP1_PATH} ${GETSTAT_SCRIPT} || die "getstat submission failed"
echo "getstat.pbs.sh submitted for hap 1" >> hifiasm.log

qsub -v REF=${REFERENCE},ASM=${HAP2_PATH} ${GETSTAT_SCRIPT} || die "getstat submission failed"
echo "getstat.pbs.sh submitted for hap 2" >> hifiasm.log

## run generate panel script
if [ -z "${FISH_NOW}" ]; then
	qsub -v FQ=${BASE_FASTQ},ASM=${ASM} ${CREATE_PANEL_SCRIPT} || die "create-launch submission failed"
	echo "create-launch.pbs.sh submitted" >> hifiasm.log
else
	qsub -v ASM=${ASM} ${RECREATE_PANEL_SCRIPT} || die "recreate submission failed"
	echo "recreate.pbs.sh submitted" >> hifiasm.log
fi
