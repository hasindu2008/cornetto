#!/bin/bash

set -x

ASM=hg002
LINEAGE=primates
qsub -v ASM=${ASM}.fasta,LINEAGE=${LINEAGE},OUT_DIR=${ASM}.busco_out ~/cornetto-hasindu/cornetto/shitflow/compleasm.pbs.sh
qsub -v ASM=${ASM}.hap1.fasta,LINEAGE=${LINEAGE},OUT_DIR=${ASM}.hap1.busco_out ~/cornetto-hasindu/cornetto/shitflow/compleasm.pbs.sh
qsub -v ASM=${ASM}.hap2.fasta,LINEAGE=${LINEAGE},OUT_DIR=${ASM}.hap2.busco_out ~/cornetto-hasindu/cornetto/shitflow/compleasm.pbs.sh


