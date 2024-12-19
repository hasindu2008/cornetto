#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 1 ] || die "1 argument required, $# provided. Usage: create_hapnetto.sh <assembly prefix>"

ASSNAME=$1
FASTA=../${ASSNAME}.fasta
ASSBED=${FASTA}.bed

test -f ${FASTA} || die "File ${FASTA} not found. Did you run create cornetto4 and are inside a directory inside that?"
test -f ${ASSBED} || die "File ${ASSBED} not found. Did you run create cornetto4 and are inside a directory inside that?"
test -f ../3_tmp.bed || die "File ../3_tmp.bed not found. Did you run create cornetto4 and are inside a directory inside that?"
test -f ../lowQ_tmp.bed || die "File ../lowQ_tmp.bed not found. Did you run create cornetto4 and are inside a directory inside that?"

test -f ../${ASSNAME}.hap1.fasta || die "File ../${ASSNAME}.hap1.fasta not found."
test -f ../${ASSNAME}.hap2.fasta || die "File ../${ASSNAME}.hap2.fasta not found."

#1# align the hapX assemblies to the primary assembly
minimap2 -t16 --eqx -cx asm5 ${FASTA} ../${ASSNAME}.hap1.fasta > ${ASSNAME}_hap1_to_asm.paf || die "minimap2 failed"
minimap2 -t16 --eqx -cx asm5 ${FASTA} ../${ASSNAME}.hap2.fasta > ${ASSNAME}_hap2_to_asm.paf || die "minimap2 failed"

GET_HAP_X_FUN () {
    HAP=$1

    # get the necessary columns from the paf file and sort (not assential, but for manual inspection if needed)
    cut -f 1-10 ${ASSNAME}_${HAP}_to_asm.paf | sort -k7,7nr -nk8,8 > ${HAP}.txt || die "cut failed"

    # go through every contig in hapX assembly separately while merging them if on same target contig within 1M bp
    rm -f ${HAP}_tmp.bed
    cut -f 1 ${HAP}.txt  | sort -u | while read ctg
    do
        grep $ctg ${HAP}.txt | awk '{print $6"\t"$8"\t"$9}' | bedtools sort | bedtools merge -d 1000000  >> ${HAP}_tmp.bed || die "awk failed"
    done

    # fun1: get the gaps on the primary assembly, that are not covered by hapX contigs
    bedtools subtract -a ${ASSBED} -b ${HAP}_tmp.bed > ${HAP}_tmp2.bed || die "bedtools subtract failed"

    # fun2:get the corners of the hapX contigs on the primary assembly (500 bp flank)
    awk '{if($2>=500){print $1"\t"$2-500"\t"$2+500} if($3>=500){print $1"\t"$3-500"\t"$3+500}}' ${HAP}_tmp.bed >> ${HAP}_tmp2.bed || die "awk failed"

    # merge the fun2 and func2 bit
    bedtools sort -i ${HAP}_tmp2.bed | bedtools merge > ${HAP}_funbits.bed || die "bedtools merge failed"
}

GET_HAP_X_FUN hap1
GET_HAP_X_FUN hap2

cat hap1_funbits.bed hap2_funbits.bed | bedtools sort | bedtools merge > hap1_hap2_funbits.bed || die "bedtools merge failed"

# now do the rest from cornetto4
#5# combine the funbits from (3) and (4) and extend them by +40kb in either direction
cat ../3_tmp.bed ../lowQ_tmp.bed hap1_hap2_funbits.bed | sort -k1,1 -k2,2n | awk '{if($2>40000){print $1"\t"$2-40000"\t"$3+40000} else {print $0}}' > funbits.bed || die "awk failed"

#6# add new windows to the file from (5), which are 200kb intervals at the left edge and right edge of the contig
awk '{if(($3-$2)>200000) {print $1"\t0\t200000\n"$1"\t"$3-200000"\t"$3}}' ${ASSBED} >> funbits.bed || die "awk failed"

#7# merge any overlapping or adjacent (within 200kb) intervals from (6)
bedtools sort -i funbits.bed | bedtools merge -d 200000 > funbits_merged.bed || die "bedtools merge failed"

#8# subtract merged windows from (7) from the whole genome assembly
bedtools subtract -a ${ASSBED} -b funbits_merged.bed > boringbits_tmp.bed || die "bedtools subtract failed"

#9# subtract any contigs shorter than 800Mbase
awk '{if(($3-$2)<800000) print $0}' ${ASSBED} > short.bed || die "awk failed"
bedtools subtract -a boringbits_tmp.bed -b short.bed > boringbits.bed || die "bedtools subtract failed"

#10# if 'boring bits' are <50% of a single contig/scaffold, remove all boring bits on the whole scaffold. Use the below horrible inefficient code snippet for now
## i.e. if the contig is more than 50% interesting, capture the whole thing
INPUT=boringbits.bed
cut -f 1 ${INPUT}  | uniq > boring_ctg.tmp || die "cut failed"
while read p;
do
ctg_len=$(grep "$p" ${ASSBED} | cut -f 3)
ctg_boring=$(grep "$p" ${INPUT} | awk '{sum+=$3-$2}END{ print sum}')
fac=$(echo "$ctg_boring*100/$ctg_len" | bc)
if [ "$fac" -gt "50" ];then
    grep "$p" ${INPUT}
fi
done < boring_ctg.tmp > ${ASSNAME}_dip.boringbits.bed || die "while loop failed"

#11# print the size of the boring_bits panel as a % of human genome size
echo -n -e "${ASSNAME}\t"
cat ${ASSNAME}_dip.boringbits.bed | awk '{sum+=($3-$2)}END{print sum/3100000000*100}'