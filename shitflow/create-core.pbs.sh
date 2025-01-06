#!/bin/bash
#PBS -P ox63
#PBS -q normal
#PBS -N create-cornetto
#PBS -l walltime=3:00:00
#PBS -l storage=gdata/ox63+scratch/ox63+scratch/if89+gdata/if89
#PBS -l mem=64GB
#PBS -l ncpus=16
#PBS -l wd

usage() {
	echo "Usage: qsub -v ASM=RGBX240039_HG002.hifiasm-hifi_1x ./create.pbs.sh" >&2
	echo
	exit 1
}

die () {
	echo >&2 "$@"
	exit 1
}

[ -z "${ASM}" ] && usage

export MODULEPATH=$MODULEPATH:/g/data/if89/apps/modulefiles/
module load samtools/1.19
module load bedtools/2.28.0
module load minimap2/2.24

export CORNETTO_BIN=/g/data/ox63/hasindu/cornetto/cornetto/cornetto
export SCRIPT_DIR=/g/data/ox63/hasindu/cornetto/cornetto/scripts

minimap2 --version || die "Could not find minimap2"
samtools --version || die "Could not find samtools"
bedtools --version || die "Could not find bedtools"
${CORNETTO_BIN} --version || die "Could not find cornetto"

${SCRIPT_DIR}/create-cornetto.sh ${ASM}.fasta|| die "create-cornetto.sh failed"

${SCRIPT_DIR}/create-hapnetto.sh ${ASM} || die "create-hapnetto.sh failed"

echo "all done. f ya."