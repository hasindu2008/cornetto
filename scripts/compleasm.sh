#!/bin/bash
#PBS -P ox63
#PBS -q normal
#PBS -N compleasm
#PBS -l walltime=12:00:00
#PBS -l storage=gdata/ox63+scratch/ox63+scratch/if89+gdata/if89+gdata/te53
#PBS -l mem=64GB
#PBS -l ncpus=16
#PBS -l wd

usage() {
	echo "Usage: qsub -v ASM=/path/to/asm.fa,OUT_DIR=busco/ ./compleasm.pbs.sh" >&2
	echo
	exit 1
}

#input assembly
[ -z "${ASM}" ] && usage

[ -z "${OUT_DIR}" ] && usage

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

## load software
module load python3/3.10.4

source /g/data/ox63/install/compleasm-0.2.6/venv/bin/activate
export PATH=/g/data/ox63/install/compleasm-0.2.6/:$PATH

THREADS=${PBS_NCPUS}

test -e ${ASM} || die "Assembly file not found: ${ASM}"
test -e ${OUT_DIR} && die "Output directory already there. Delete that first: ${OUT_DIR}"

#########################################

compleasm run -a ${ASM} -o ${OUT_DIR} -t 8 -l primates -L /g/data/ox63/install/compleasm-0.2.6/mb_downloads