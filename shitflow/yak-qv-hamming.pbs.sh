#!/bin/bash
#PBS -P ox63
#PBS -N yak
#PBS -q normal
#PBS -l ncpus=48
#PBS -l mem=192GB
#PBS -l walltime=48:00:00
#PBS -l wd
#PBS -l storage=gdata/if89+scratch/wv19+gdata/wv19+gdata/ox63

usage() {
	echo "Usage: qsub -v ASM=/path/to/asm.fa,REF=/path/to/ref.yak(optional),MAT=/path/to/mat.yak(optional),PAT=/path/to/pat.yak(optional) ./yak.pbs.sh" >&2
	echo "       optional ones if not provided defaults to hg002"
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

# creating yak index ${YAK} count -K1.5g -t ${THREADS} ${REF} -o ${REF}.yak
#input assembly
[ -z "${ASM}" ] && usage
[ -z "${REF}" ] && REF=/g/data/ox63/hasindu/cornetto/ref/hg002v1.0.1.fasta.gz.yak
[ -z "${MAT}" ] && MAT=/g/data/ox63/hasindu/cornetto/ref/hg002_parents/mat.HG004.yak
[ -z "${PAT}" ] && PAT=/g/data/ox63/hasindu/cornetto/ref/hg002_parents/pat.HG003.yak

${YAK} version || die "yak not found"

test -e ${ASM} || die "Assembly file not found: ${ASM}"
test -e ${REF} || die "Reference yak not found: ${REF}"
test -e ${MAT} || die "Maternal yak not found: ${MAT}"
test -e ${PAT} || die "Paternal yak not found: ${PAT}"


/usr/bin/time -v ${YAK} qv ${REF} ${ASM} -t ${THREADS} > ${ASM}.yak.txt || die "yak qv failed"

/usr/bin/time -v ${YAK} trioeval ${PAT} ${MAT} ${ASM} -t ${THREADS} >> ${ASM}.yak.txt || die "yak hamming failed"