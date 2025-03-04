
#PBS -P ox63
#PBS -q normal
#PBS -N index
#PBS -l walltime=4:00:00
#PBS -l storage=gdata/ox63+scratch/ox63+scratch/if89+gdata/if89+gdata/te53
#PBS -l mem=36GB
#PBS -l ncpus=8
#PBS -l wd

THREADS=${PBS_NCPUS}
module load minimap2/2.24

set -o pipefail

die () {
    echo >&2 "$@"
    exit 1
}

export PATH=$PATH:/g/data/te53/ontsv/sv_parsing/scripts/

usage() {
	echo "Usage: qsub -v HUMAN=human-only/,NONHUMAN=all-crap/centrifuge/,ASM=saliva-Q_0 ./xcript.pbs.sh" >&2
	echo
	exit 1
}

HUMAN_ASM=${HUMAN}/${ASM}.fasta
HUMAN_BORING_BED=${HUMAN}/${ASM}_dip.boringbits.bed

NONHUMAN_CTG_FA=${NONHUMAN}/${ASM}.nonhuman_contigs.fasta
NONHUMAN_CTG_BED=${NONHUMAN}/${ASM}.nonhuman_contigs.bed

COMBINED_FA=${ASM}.plus_nonhuman_ctg.fasta
COMBINED_BORING_BED=${ASM}_dip.boringbits.plus_nonhuman_ctg.bed
COMBINED_BORING_TXT=${ASM}_dip.boringbits.plus_nonhuman_ctg.txt

flattenFasta.pl -tab ${NONHUMAN_CTG_FA} | awk '{print $1"_nonhuman\t"$2}' > tmp1.tab || die "Failed to flatten ${NONHUMAN_CTG_FA}"
flattenFasta.pl -fa tmp1.tab > ${ASM}_nonhuman_ctg_renamed.fasta || die "Failed to create ${ASM}_nonhuman_ctg_renamed

cat ${NONHUMAN_CTG_BED} | awk '{print $1"_nonhuman\t"$2"\t"$3}' > ${ASM}_nonhuman_ctg_renamed.bed || die "Failed to create ${ASM}_nonhuman_ctg_renamed.bed"

cat ${HUMAN_ASM} ${ASM}_nonhuman_ctg_renamed.fasta > ${COMBINED_FA} || die "Failed to create ${COMBINED_FA}"
cat ${HUMAN_BORING_BED} ${ASM}_nonhuman_ctg_renamed.bed > ${COMBINED_BORING_BED} || die "Failed to create ${COMBINED_BORING_BED}"

cat ${COMBINED_BORING_BED} | awk '{print $1","$2","$3",+"}' > plus_tmp1 || die "Failed to create plus_tmp1"
cat ${COMBINED_BORING_BED} | awk '{print $1","$2","$3",-"}' > minus_tmp1 || die "Failed to create minus_tmp1"
cat plus_tmp1 minus_tmp1 | sort > ${COMBINED_BORING_TXT} || die "Failed to create ${COMBINED_BORING_TXT}"

minimap2 -x map-ont -d ${COMBINED_FA}.idx ${COMBINED_FA} || die "Failed to index ${COMBINED_FA}"