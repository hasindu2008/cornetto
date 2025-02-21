#!/bin/bash

ASM=hg002

qsub -v REF=/g/data/ox63/hasindu/cornetto/ref/hg002v1.0.1.fasta.gz,ASM=${ASM}.fasta ~/cornetto-hasindu/cornetto/shitflow/yak-qv.pbs.sh

cat ${ASM}.hap1.fasta ${ASM}.hap2.fasta > ${ASM}.hap1+hap2.fasta
qsub -v ASM=${ASM}.hap1+hap2.fasta ~/cornetto-hasindu/cornetto/shitflow/yak-qv-hamming.pbs.sh

