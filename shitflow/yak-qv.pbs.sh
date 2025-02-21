#!/bin/bash
#PBS -P ox63
#PBS -N yak
#PBS -q normal
#PBS -l ncpus=16
#PBS -l mem=64GB
#PBS -l walltime=48:00:00
#PBS -l wd
#PBS -l storage=gdata/if89+scratch/wv19+gdata/wv19+gdata/ox63

usage() {
	echo "Usage: qsub -v REF=/path/to/ref.fasta,ASM=/path/to/asm.fa ./yak.pbs.sh" >&2
	echo
	exit 1
}

die () {
    echo "$1" >&2
    echo
    exit 1
}

YAK=/g/data/ox63/install/yak/yak
THREADS=${PBS_NCPUS}

#input assembly
[ -z "${ASM}" ] && usage
[ -z "${REF}" ] && usage

${YAK} version || die "yak not found"

test -e ${ASM} || die "Assembly file not found: ${ASM}"
test -e ${REF} || die "Reference file not found: ${REF}"

test -e ${REF}.yak || ${YAK} count -K1.5g -t ${THREADS} ${REF} -o ${REF}.yak || die "yak count failed"

/usr/bin/time -v ${YAK} qv ${REF}.yak ${ASM} -t ${THREADS} > ${ASM}.yak.txt || die "yak qv failed"