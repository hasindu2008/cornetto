#!/bin/bash

REF=/g/data/ox63/hasindu/cornetto/ref/hg002v1.0.1.fasta.gz

while read p; do
    DIR=$(dirname $p)
    FILE=$(basename $p)
    cd $DIR
    qsub -v REF=${REF},ASM=$FILE ~/cornetto-hasindu/cornetto/shitflow/yak-qv.k21.pbs.sh
done