#!/bin/bash
#PBS -P ox63
#PBS -q normal
#PBS -N quast
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/ox63+scratch/ox63+scratch/if89+gdata/if89
#PBS -l mem=192GB
#PBS -l ncpus=48
#PBS -l wd

###################################################################

###################################################################

# Make sure to change:
# 1. ox63 to your own project

# to run:
# qsub -v ASM=/path/to/asm.fasta,OUT_DIR=/path/to/dir ./quast.pbs.sh

###################################################################


usage() {
	echo "Usage: qsub -v ASM=/path/to/asm.fasta,OUT_DIR=/path/to/dir ./quast.pbs.sh" >&2
	echo
	exit 1
}

[ -z "${ASM}" ] && usage
[ -z "${OUT_DIR}" ] && usage

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

export MODULEPATH=$MODULEPATH:/g/data/if89/apps/modulefiles/
module load quast/5.2.0
THREADS=${PBS_NCPUS}

#########################################

test -d ${OUT_DIR} && die "Output directory ${OUT_DIR} already exists. Please delete it first or give an alternate location. Exiting."

## run quast to evaluate assemblies
quast.py -t ${THREADS} -o ${OUT_DIR} -l ${ASM} --large ${ASM} -r /g/data/ox63/hasindu/cornetto/ref/hg002v1.0.1.fasta.gz || die "quast failed"





