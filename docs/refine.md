# Further Refinements

For primary assemblies, you may use the following algorithm to further refine your assemblies.

Suppose your base assembly is called `asm-0.fasta` and the cornetto iterations are named `asm-1.fasta`, `asm-2.fasta`, `asm-3.fasta`, ..., `asm-n.fasta`.

1. Extracts contigs from `asm-1.fasta` that are longer than the expected minimum chromosome size and have telomeres in both ends.

2.

```
CURR=1;
T2T_FOUND=0;

for asm-1.fasta asm-2.fasta asm-3.fasta ... asm-n.fasta
do

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
```