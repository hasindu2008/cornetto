#!/bin/bash
#PBS -P ox63
#PBS -N getstat
#PBS -q normal
#PBS -l ncpus=16
#PBS -l mem=64GB
#PBS -l walltime=12:00:00
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
module load bedtools/2.28.0

###################################################################

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

minimap2 --version || die "Could not find minimap2"

export CORNETTO=/g/data/ox63/hasindu/cornetto/cornetto/cornetto
export SCRIPT_DIR=/g/data/ox63/hasindu/cornetto/cornetto/scripts

test -z $REF && die "Reference file not provided"
test -z $ASM && die "Assembly file not provided"

FILENAME=$(basename ${ASM})
test -e ${FILENAME} || cp ${ASM} ${FILENAME}
test -e ${FILENAME} ||  die "Assembly file not available in the current directory"

${SCRIPT_DIR}/minidotplot.sh ${REF} ${FILENAME}  || die "Failed to run minidotplot"

${SCRIPT_DIR}/telostats.sh ${FILENAME} > ${FILENAME}.telostats.txt || die "Failed to run telostats"

${SCRIPT_DIR}/asmstats.sh ${FILENAME} > ${FILENAME}.asmstats.txt || die "Failed to run asmstats"
