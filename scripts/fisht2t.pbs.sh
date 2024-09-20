#!/bin/bash
#PBS -P ox63
#PBS -N FISHT2T
#PBS -q normal
#PBS -l ncpus=32
#PBS -l mem=128GB
#PBS -l walltime=6:00:00
#PBS -l wd
#PBS -l storage=gdata/if89+scratch/wv19+gdata/wv19+gdata/ox63

###################################################################

###################################################################

# Make sure to change:
# 1. wv19 and ox63 to your own projects

# to run:
# qsub -v REF=/path/to/ref.fa,ASM=/path/to/asm.fa ./getstat.pbs.sh

###################################################################

#ASM_LIST="A_1_QGXHXX240275:A_2_QGXHXX240279:A_3_QGXHXX240283:A_4_QGXHXX240293:A_5_QGXHXX240298:A_6_QGXHXX240304:A_7_QGXHXX240308:A_8_QGXHXX240317:A_9_QGXHXX240325"
ASM_WORK_DIR=/g/data/ox63/hasindu/cornetto/autocall

usage() {
	echo "Usage: qsub -v ASM_LIST=A_1_QGXHXX240275:A_2_QGXHXX240279a.REF=/path/to/ref.fa ./fisht2t.pbs.sh" >&2
	echo
	exit 1
}


#asm
[ -z "${ASM_LIST}" ] && usage
#ref
[ -z "${REF}" ] && REF=/g/data/ox63/cornetto/data/reference/hg002v1.0.1_pat.fa
#ASM_NAME_PREFIX
[ -z "${ASM_NAME_PREFIX}" ] && ASM_NAME_PREFIX=hg002-cornetto-
#min contig len
[ -z "${MIN_CONTIG_LEN}" ] && MIN_CONTIG_LEN=40000000

module load minimap2/2.24
module load samtools/1.12

###################################################################

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}


REFERENCE=${REF}
ASM_DIR_LIST=$(echo $ASM_LIST | tr ':' ' ')
GETSTAT_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/scripts/getstat.pbs.sh
QUAST_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/scripts/quast.pbs.sh
THREADS=${PBS_NCPUS}


COUNT=$(echo $ASM_DIR_LIST | wc -w)
echo "Total number of assemblies: $COUNT"
[ $COUNT -lt 2 ] && die "Need at least two assemblies to work on"

MINIMAP2=minimap2
SAMTOOLS=samtools

CHECK_SHIT(){

	${MINIMAP2} --version || die "Could not find minimap2"
	${SAMTOOLS} --version || die "Could not find samtools"

	test -e ${GETSTAT_SCRIPT} || die "${GETSTAT_SCRIPT} not found"
	test -e ${QUAST_SCRIPT} || die "${QUAST_SCRIPT} not found"
	test -e ${REFERENCE} || die "${REFERENCE} not found"

	for ASM_DIR_NAME in $ASM_DIR_LIST
	do
		ASM_DIR=${ASM_WORK_DIR}/${ASM_DIR_NAME}
		LETTER=$(echo $ASM_DIR_NAME | cut -d'_' -f1)
		NUM=$(echo $ASM_DIR_NAME | cut -d'_' -f2)
		ASM_NAME=${ASM_NAME_PREFIX}${LETTER}_${NUM}
		ASM=$ASM_DIR/${ASM_NAME}.fasta
		TELO_END=$ASM_DIR/${ASM_NAME}/${ASM_NAME}.windows.0.4.50kb.ends.bed


		test -d $ASM_DIR || die "Assembly directory not found: $ASM_DIR"
		test -e $ASM || die "Assembly file not found: $ASM"
		test -e $TELO_END || die "Telo ends file not found: $TELO_END"

	done

}

CHECK_SHIT

EXTRACT_CONTIG(){

	NAME=$1
	samtools faidx ${ASM} -r ${ASM_DIR_NAME}.${NAME}.txt > ${ASM_DIR_NAME}.${NAME}.fasta || die "Could not extract the long t2t"
	awk -v prefix=${LETTER}_${NUM}_${NAME} '{if($1 ~ /^>/) print ">"prefix"_"substr($1,2); else print $0}' ${ASM_DIR_NAME}.${NAME}.fasta > ${ASM_DIR_NAME}.${NAME}.renamed.fasta || die "Could not rename the long t2t" || die "Could not rename the long t2t"
	grep "^>" ${ASM_DIR_NAME}.${NAME}.renamed.fasta | tr -d '>' > ${ASM_DIR_NAME}.${NAME}.renamed.txt || die "Could not get the names of renamed contigs"

}

GET_NEWFOUND_LIST(){

	NAME=$1
	rm -f ${ASM_DIR_NAME}.${NAME}.renamed.newfound.txt
	${MINIMAP2} -t ${THREADS} -K4G --eqx -cx asm5 ${OUTPUT_NAME}.t2t.fasta ${ASM_DIR_NAME}.${NAME}.renamed.fasta > ${ASM_DIR_NAME}.${NAME}.paf || die "minimap2 failed"

	while read p;
	do
		if grep -q $p ${ASM_DIR_NAME}.${NAME}.paf
		then
			grep $p ${ASM_DIR_NAME}.${NAME}.paf | awk 'BEGIN{sum=0} {sum+=($4-$3)} END{if(sum/$2<0.5){print $1}}'  >> ${ASM_DIR_NAME}.${NAME}.renamed.newfound.txt
		else
			echo "            $p is not in PAF, adding to newfound: ${ASM_DIR_NAME}.${NAME}.renamed.newfound.txt"
			echo $p >> ${ASM_DIR_NAME}.${NAME}.renamed.newfound.txt
		fi
	done <  ${ASM_DIR_NAME}.${NAME}.renamed.txt

}

CURR=1;
T2T_FOUND=0;
for ASM_DIR_NAME in $ASM_DIR_LIST
do
	ASM_DIR=${ASM_WORK_DIR}/${ASM_DIR_NAME}
	LETTER=$(echo $ASM_DIR_NAME | cut -d'_' -f1)
	NUM=$(echo $ASM_DIR_NAME | cut -d'_' -f2)
	ASM_NAME=${ASM_NAME_PREFIX}${LETTER}_${NUM}
	ASM=$ASM_DIR/${ASM_NAME}.fasta
	TELO_END=$ASM_DIR/${ASM_NAME}/${ASM_NAME}.windows.0.4.50kb.ends.bed
	OUTPUT_NAME=${ASM_NAME_PREFIX}${LETTER}

	echo "Doing ${ASM_DIR_NAME}"

	N_T2T=$(cat $TELO_END | cut -f 1 | sort | uniq -c | sort -k1,1 -n -r  | awk '{if($1==2) print $2}' | wc -l)

	if [ ${N_T2T} -gt 0 ]
	then
		echo "    ${ASM_DIR_NAME} has ${N_T2T} T2T"

		# get the t2t > 40M
		cat $TELO_END | cut -f 1 | sort | uniq -c | sort -k1,1 -n -r  | awk '{if($1==2) print $2}' > ${ASM_DIR_NAME}.t2t.tmp.txt || die "Could not get the all t2t"
		test -e ${ASM}.fai || samtools faidx ${$ASM} || die "Could not create the faidx"
		cat ${ASM}.fai  | awk -v LEN=${MIN_CONTIG_LEN} '{if($2>LEN) print $1}' > ${ASM_DIR_NAME}.longcontig.txt || die "grabbing the longcontig failed"
		grep ${ASM_DIR_NAME}.t2t.tmp.txt -F -f ${ASM_DIR_NAME}.longcontig.txt >  ${ASM_DIR_NAME}.t2t.txt || die "grabbing the t2t failed"

		# extract the long T2T to a file with the contif names appended with LETTER_NUM
		EXTRACT_CONTIG t2t

		# if there are any t2t, do the following

		if [ ${T2T_FOUND} -eq 0 ] # if no telo as been found so far, create the base t2t
		then
			echo "        Creating the base t2t: ${OUTPUT_NAME}.t2t.fasta"
			cp ${ASM_DIR_NAME}.t2t.renamed.fasta ${OUTPUT_NAME}.t2t.fasta || die "Could not copy the base t2t"
		else  # otherwise map the t2t to the base t2t; append any new t2t (which map <50% to any in base t2t) to the base t2t
			echo "        Finding newfound t2t"

			GET_NEWFOUND_LIST t2t

			if [ -s ${ASM_DIR_NAME}.t2t.renamed.newfound.txt ]
			then
				echo "            Found newfound t2t. See ${ASM_DIR_NAME}.t2t.renamed.newfound.txt"
				samtools faidx ${ASM_DIR_NAME}.t2t.renamed.fasta || die "Could not index the  fasta"
				samtools faidx ${ASM_DIR_NAME}.t2t.renamed.fasta -r ${ASM_DIR_NAME}.t2t.renamed.newfound.txt >> ${OUTPUT_NAME}.t2t.fasta || die "Could not append the newfound t2t"
			else
				echo "            No newfound t2t"
			fi

		fi

		T2T_FOUND=1

		# otherwise, do nothing

	else
		echo "    ${ASM_DIR_NAME} has no T2T"
	fi

	# if the last one, append any none t2t contigs to the base asm, which will be the final reference
	if [ $CURR -eq $COUNT ]
	then
		echo "    Last one. Appending none T2T contigs to the base asm"
		cut -f1 ${ASM}.fai  | grep -v -F -f ${ASM_DIR_NAME}.t2t.txt >  ${ASM_DIR_NAME}.nont2t.txt || die "grabbing the non t2t failed"

		EXTRACT_CONTIG nont2t
		GET_NEWFOUND_LIST nont2t

		samtools faidx ${ASM_DIR_NAME}.nont2t.renamed.fasta || die "Could not index the  fasta"
		samtools faidx ${ASM_DIR_NAME}.nont2t.renamed.fasta -r ${ASM_DIR_NAME}.nont2t.renamed.newfound.txt > ${OUTPUT_NAME}.nont2t.fasta || die "Could not append the newfound t2t"

		cat ${OUTPUT_NAME}.t2t.fasta ${OUTPUT_NAME}.nont2t.fasta > ${OUTPUT_NAME}.fasta || die "Could not append the non t2t to the t2t"
	fi


	CURR=$((CURR+1))
done

echo "Running the quast and asm stats"
## run quast to evaluate assemblies
ASMPATH=`realpath ${OUTPUT_NAME}.fasta`
qsub -v ASM=${ASMPATH},OUT_DIR=${ASMPATH}.quast_out ${QUAST_SCRIPT}
qsub -v REF=${REFERENCE},ASM=${ASMPATH} ${GETSTAT_SCRIPT}

ASMPATH=`realpath ${OUTPUT_NAME}.t2t.fasta`
qsub -v ASM=${ASMPATH},OUT_DIR=${ASMPATH}.quast_out ${QUAST_SCRIPT}
qsub -v REF=${REFERENCE},ASM=${ASMPATH} ${GETSTAT_SCRIPT}
