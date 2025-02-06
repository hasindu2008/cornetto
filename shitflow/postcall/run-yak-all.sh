#!/bin/bash


LIST="1 2 3 4 5 6 7 8 9"
REF=/g/data/ox63/hasindu/cornetto/ref/hg002v1.0.1.fasta.gz
ASM=hg002-cornetto-A

for EACH in ${LIST}
do
    cd ${EACH}
    qsub -v REF=${REF},ASM=${ASM}.fasta ~/cornetto-hasindu/cornetto/shitflow/yak.pbs.sh
    cd ..
done



