#!/bin/bash


LIST="1 2 3 4 5 6 7 8 9"
ASM=hg002-cornetto-A
LINEAGE=primates


for EACH in ${LIST}
do
    cd ${EACH}
    qsub -v ASM=${ASM}.fasta,LINEAGE=${LINEAGE},OUT_DIR=${ASM}.busco_out  ~/cornetto-hasindu/cornetto/shitflow/compleasm.pbs.sh
    cd ..
done



