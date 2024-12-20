#!/bin/bash
#PBS -P ox63
#PBS -N generate_panel
#PBS -q normal
#PBS -l ncpus=4
#PBS -l mem=16GB
#PBS -l walltime=48:00:00
#PBS -l wd
#PBS -l storage=gdata/if89+scratch/wv19+gdata/wv19+gdata/ox63

usage() {
	echo "Usage: qsub -v FISH_NOW=A_3_QGXHXX240283,PREFIX=hg002-cornetto-A_3 ./generate_panel.pbs.sh" >&2
	echo
	exit 1
}

[ -z "${PREFIX}" ] && usage
[ -z "${FISH_NOW}" ] && usage

module load bedtools/2.28.0
module load minimap2/2.24
module load samtools/1.19

ONT_DATADIR=/g/data/ox63/hasindu/cornetto/shitflow/
export SCRIPT_DIR=/g/data/ox63/hasindu/cornetto/cornetto/scripts

cd ${ONT_DATADIR}/${FISH_NOW}

minimap2 -x map-ont -d ${ASM}.fasta.idx ${ASM}.fasta || die "minimap2 failed"
${SCRIPT_DIR}/recreate-cornetto.sh ${ASM} || die "create-cornetto.sh failed"

