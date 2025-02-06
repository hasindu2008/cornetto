#!/bin/bash

HAP1=assembly_verkko_v1.3_trio.maternal.fasta
HAP2=assembly_verkko_v1.3_trio.paternal.fasta

set -x

qsub -v ASM=${HAP1},OUT_DIR=${HAP1%%.fasta}.quast_out /g/data/ox63/hasindu/cornetto/cornetto/shitflow/quast.pbs.sh
qsub -v ASM=${HAP2},OUT_DIR=${HAP2%%.fasta}.quast_out /g/data/ox63/hasindu/cornetto/cornetto/shitflow/quast.pbs.sh

qsub -v REF=~/cornetto-hasindu/ref/hg002v1.0.1_pat.fa,ASM=${HAP1} /g/data/ox63/hasindu/cornetto/cornetto/shitflow/getstat.pbs.sh
qsub -v REF=~/cornetto-hasindu/ref/hg002v1.0.1_pat.fa,ASM=${HAP2} /g/data/ox63/hasindu/cornetto/cornetto/shitflow/getstat.pbs.sh

LINEAGE=primates
qsub -v ASM=${HAP1},LINEAGE=${LINEAGE},OUT_DIR=${HAP1%%.fasta}.busco_out /g/data/ox63/hasindu/cornetto/cornetto/shitflow/compleasm.pbs.sh
qsub -v ASM=${HAP2},LINEAGE=${LINEAGE},OUT_DIR=${HAP2%%.fasta}.busco_out /g/data/ox63/hasindu/cornetto/cornetto/shitflow/compleasm.pbs.sh

cat ${HAP1} ${HAP2} >  ${HAP1%%.fasta}"+"${HAP2%%.fasta}.fasta
qsub -v ASM=${HAP1%%.fasta}"+"${HAP2%%.fasta}.fasta /g/data/ox63/hasindu/cornetto/cornetto/shitflow/yak-qv-hamming.pbs.sh
