#!/bin/bash
#PBS -P ox63
#PBS -N getstat
#PBS -q normal
#PBS -l ncpus=16
#PBS -l mem=64GB
#PBS -l walltime=48:00:00
#PBS -l wd
#PBS -l storage=gdata/if89+scratch/wv19+gdata/wv19+gdata/ox63

###################################################################

###################################################################

# Make sure to change:
# 1. wv19 and ox63 to your own projects

# to run:
# qsub -v REF=/path/to/ref.fa,ASM=/path/to/asm.fa ./getstat.pbs.sh

###################################################################

usage() {
	echo "Usage: qsub -v REF=/path/to/ref.fa,ASM=/path/to/asm.fa ./getstat.pbs.sh" >&2
	echo
	exit 1
}


#asm
[ -z "${ASM}" ] && usage
#ref
[ -z "${REF}" ] && usage

module load minimap2/2.24
module load samtools/1.12

###################################################################

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

minimap2 --version || die "Could not find minimap2"
samtools --version || die "Could not find samtools"

export MINIDOT=/g/data/ox63/install/miniasm/minidot
export PATH=$PATH:/g/data/ox63/install/datamash
export TELO_SCRIPT_PATH=/g/data/ox63/install/vgp-pipeline/telomere/telomere_analysis.sh
export SCRIPT_DIR=/g/data/ox63/hasindu2008.git/cornetto/scripts

${SCRIPT_DIR}/minidotplot.sh ${REF} ${ASM}  || die "Failed to run minidotplot"

${SCRIPT_DIR}/telostas.sh ${ASM}  || die "Failed to run telostas"

${SCRIPT_DIR}/asmstats.sh ${ASM} > ${ASM}.asmstats.txt || die "Failed to run asmstats"
