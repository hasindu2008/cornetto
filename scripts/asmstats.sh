#!/bin/bash

die() {
    echo "$1"
    exit 1
}

#DEBUG
# cd /home/hasindu/scratch/cornetto/ass/1x_hifi
# ./asmstats.sh RGBX240039_HG002.hifiasm.primary_asm.fasta

[ $# -ne 1 ] && die "Usage: $0 <FASTA>"

FASTA=$1
PREFIX=$(basename $FASTA .fa)
PREFIX=$(basename $PREFIX .fasta)
FILE=$PREFIX.fasta.tmp.paf

test -z ${CORNETTO} && CORNETTO=cornetto
${CORNETTO} --version > /dev/null || die "cornetto executable not found! Either put cornetto under path or set CORNETTO variable, e.g.,export CORNETTO=/path/to/cornetto"

test -e $FILE  || die "File $FILE does not exist. Did you run minidotplot.sh?"
test -e ${PREFIX}.fasta.report.tsv || die "File ${PREFIX}.fasta.report.tsv does not exist. Did you run minidotplot.sh?"
test -e ${PREFIX}/${PREFIX}.windows.0.4.50kb.ends.bed || die "File ${PREFIX}/${PREFIX}.windows.0.4.50kb.ends.bed does not exist. Did you run telostats.sh?"

echo "${FASTA}"
echo ""
${CORNETTO} asmstats ${FILE} ${PREFIX}/${PREFIX}.windows.0.4.50kb.ends.bed -r ${PREFIX}.fasta.report.tsv -H -T
