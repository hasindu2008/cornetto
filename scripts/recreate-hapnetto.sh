#!/bin/bash

set -o pipefail

die () {
    echo >&2 "$@"
    exit 1
}

[ "$#" -eq 1 ] || die "1 argument required, $# provided. Usage: create_hapnetto.sh <assembly prefix>"

test -z ${CORNETTO} && CORNETTO=cornetto
${CORNETTO} --version || die "cornetto executable not found! Either put cornetto under path or set CORNETTO variable, e.g.,export CORNETTO=/path/to/cornetto"
test -z ${BEDTOOLS} && BEDTOOLS=bedtools
$BEDTOOLS --version > /dev/null 2>&1 || die "bedtools not found!. Either put bedtools under path or set BEDTOOLS variable, e.g.,export BEDTOOLS=/path/to/bedtools"
test -z ${MINIMAP2} && MINIMAP2=minimap2
$MINIMAP2 --version > /dev/null 2>&1 || die "minimap2 not found!. Either put minimap2 under path or set MINIMAP2 variable, e.g.,export MINIMAP2=/path/to/minimap2"

TMPOUT=tmp_recreate_hapnetto
test -d ${TMPOUT} && die "Directory ${TMPOUT} already exists. Please remove it before running this script or change to a different working directory"
mkdir ${TMPOUT} || die "mkdir failed"

TMPOUT_PREV=tmp_recreate_cornetto
test -d ${TMPOUT_PREV} || die "Directory ${TMPOUT_PREV} not found. Did you run recreate cornetto under the current directory?"

ASSNAME=$1
FASTA=${ASSNAME}.fasta
ASSBED=${TMPOUT_PREV}/${ASSNAME}.chroms.bed
test -f ${FASTA} || die "File ${FASTA} not found. Did you run recreate cornetto and are inside a directory inside that?"
test -f ${ASSBED} || die "File ${ASSBED} not found. Did you run create cornetto4 and are inside a directory inside that?"
test -f ${TMPOUT_PREV}/lowQ_tmp.bed || die "File ${TMPOUT_PREV}/lowQ_tmp.bed not found. Did you run create cornetto4 and are inside a directory inside that?"

test -f ${ASSNAME}.hap1.fasta || die "File ${ASSNAME}.hap1.fasta not found."
test -f ${ASSNAME}.hap2.fasta || die "File ${ASSNAME}.hap2.fasta not found."

#1# align the hapX assemblies to the primary assembly
${MINIMAP2} -t16 --eqx -cx asm5 ${FASTA} ${ASSNAME}.hap1.fasta > ${TMPOUT}/${ASSNAME}_hap1_to_asm.paf || ${MINIMAP2} -t16 --eqx -x asm5 ${FASTA} ${ASSNAME}.hap1.fasta > ${TMPOUT}/${ASSNAME}_hap1_to_asm.paf || die "minimap2 failed"
${MINIMAP2} -t16 --eqx -cx asm5 ${FASTA} ${ASSNAME}.hap2.fasta > ${TMPOUT}/${ASSNAME}_hap2_to_asm.paf || ${MINIMAP2} -t16 --eqx -x asm5 ${FASTA} ${ASSNAME}.hap2.fasta > ${TMPOUT}/${ASSNAME}_hap2_to_asm.paf || die "minimap2 failed"

GET_HAP_X_FUN () {
    HAP=$1

    # get the necessary columns from the paf file and sort (not essential, but for manual inspection if needed)
    cut -f 1-10 ${TMPOUT}/${ASSNAME}_${HAP}_to_asm.paf | sort -k7,7nr -nk8,8 > ${TMPOUT}/${HAP}.txt || die "cut failed"

    # go through every contig in hapX assembly separately while merging them if on same target contig within 1M bp
    # TODO would benefit from a good implementation
    rm -f ${TMPOUT}/${HAP}_tmp.bed
    cut -f 1 ${TMPOUT}/${HAP}.txt  | sort -u | while read ctg
    do
        cat ${TMPOUT}/${HAP}.txt | awk -v ctg=$ctg '{if($1==ctg){print $6"\t"$8"\t"$9}}' | ${BEDTOOLS} sort | ${BEDTOOLS} merge -d 1000000  >> ${TMPOUT}/${HAP}_tmp.bed || die "awk failed"
    done

    # fun1: get the gaps on the primary assembly, that are not covered by hapX contigs
    ${BEDTOOLS} subtract -a ${ASSBED} -b ${TMPOUT}/${HAP}_tmp.bed > ${TMPOUT}/${HAP}_tmp2.bed || die "bedtools subtract failed"

    # fun2:get the corners of the hapX contigs on the primary assembly (500 bp flank)
    awk '{if($2>=500){print $1"\t"$2-500"\t"$2+500} if($3>=500){print $1"\t"$3-500"\t"$3+500}}' ${TMPOUT}/${HAP}_tmp.bed >> ${TMPOUT}/${HAP}_tmp2.bed || die "awk failed"

    # merge the fun1 and fun2 bit
    ${BEDTOOLS} sort -i ${TMPOUT}/${HAP}_tmp2.bed | ${BEDTOOLS} merge > ${TMPOUT}/${HAP}_funbits.bed || die "bedtools merge failed"
}

GET_HAP_X_FUN hap1
GET_HAP_X_FUN hap2

cat ${TMPOUT}/hap1_funbits.bed ${TMPOUT}/hap2_funbits.bed | ${BEDTOOLS} sort | ${BEDTOOLS} merge > ${TMPOUT}/hap1_hap2_funbits.bed || die "bedtools merge failed"

# now do the rest from cornetto4
#5# combine the funbits from (3) and (4) and extend them by +40kb in either direction
cat ${TMPOUT_PREV}/lowQ_tmp.bed  ${TMPOUT}/hap1_hap2_funbits.bed | sort -k1,1 -k2,2n | awk '{if($2>40000){print $1"\t"$2-40000"\t"$3+40000} else {print $0}}' >  ${TMPOUT}/funbits.bed || die "awk failed"

#6# add new windows to the file from (5), which are 200kb intervals at the left edge and right edge of the contig
awk '{if(($3-$2)>200000) {print $1"\t0\t200000\n"$1"\t"$3-200000"\t"$3}}' ${ASSBED} >>  ${TMPOUT}/funbits.bed || die "awk failed"

#7# merge any overlapping or adjacent (within 200kb) intervals from (6)
${BEDTOOLS} sort -i  ${TMPOUT}/funbits.bed | ${BEDTOOLS} merge -d 200000 >  ${TMPOUT}/funbits_merged.bed || die "bedtools merge failed"

#8# subtract merged windows from (7) from the whole genome assembly
${BEDTOOLS} subtract -a ${ASSBED} -b  ${TMPOUT}/funbits_merged.bed >  ${TMPOUT}/boringbits_tmp.bed || die "bedtools subtract failed"

#9# subtract any contigs shorter than 800Mbase
awk '{if(($3-$2)<800000) print $0}' ${ASSBED} >  ${TMPOUT}/short.bed || die "awk failed"
${BEDTOOLS} subtract -a  ${TMPOUT}/boringbits_tmp.bed -b  ${TMPOUT}/short.bed >  ${TMPOUT}/boringbits.bed || die "bedtools subtract failed"

#10# if 'boring bits' are <50% of a single contig/scaffold, remove all boring bits on the whole scaffold. Then create readfish targets
${CORNETTO} bigenough ${ASSBED} ${TMPOUT}/boringbits.bed -r ${ASSNAME}_dip.boringbits.txt > ${ASSNAME}_dip.boringbits.bed || die "cornetto bigenough failed"

# old crappy code to do the same thing. Remove after testing
#10# if 'boring bits' are <50% of a single contig/scaffold, remove all boring bits on the whole scaffold. Use the below horrible inefficient code snippet for now
# INPUT=${TMPOUT}/boringbits.bed
# cut -f 1 ${INPUT}  | uniq >  ${TMPOUT}/boring_ctg.tmp || die "cut failed"
# while read p;
# do
# ctg_len=$(grep "$p" ${ASSBED} | cut -f 3)
# ctg_boring=$(grep "$p" ${INPUT} | awk '{sum+=$3-$2}END{ print sum}')
# fac=$(echo "$ctg_boring*100/$ctg_len" | bc)
# if [ "$fac" -gt "50" ];then
#     grep "$p" ${INPUT}
# fi
# done <  ${TMPOUT}/boring_ctg.tmp >  ${ASSNAME}_dip.boringbits.bed || die "while loop failed"

#11# print the size of the boring_bits panel as a % of human genome size
# echo -n -e "${ASSNAME}\t"
# cat ${ASSNAME}_dip.boringbits.bed | awk '{sum+=($3-$2)}END{print sum/3100000000*100}'

#12# create readfish targets
# cat ${ASSNAME}_dip.boringbits.bed  | awk '{print $1","$2","$3",+"}' > ${TMPOUT}/plus_tmp
# cat ${ASSNAME}_dip.boringbits.bed  | awk '{print $1","$2","$3",-"}' > ${TMPOUT}/minus_tmp
# cat ${TMPOUT}/plus_tmp ${TMPOUT}/minus_tmp | sort > ${ASSNAME}_dip.boringbits.txt

