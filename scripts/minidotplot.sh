#!/bin/bash

set -o pipefail

die() {
    echo "$1"
    exit 1
}

[ $# -ne 2 ] && die "Usage: $0 <reference> <myassembly>"

REF=$1
ASM=$2

[ ! -f $REF ] && die "Reference $REF not found"
[ ! -f $ASM ] && die "Assembly $ASM not found"

test -z ${CORNETTO} && CORNETTO=cornetto
${CORNETTO} --version > /dev/null 2>&1 || die "cornetto executable not found! Either put cornetto under path or set CORNETTO variable, e.g.,export CORNETTO=/path/to/cornetto"
test -z ${MINIMAP2} && MINIMAP2=minimap2
$MINIMAP2 --version > /dev/null 2>&1 || die "minimap2 not found!. Either put minimap2 under path or set MINIMAP2 variable, e.g.,export MINIMAP2=/path/to/minimap2"


PREFIX=$(basename $ASM .fa)
PREFIX=$(basename $ASM .fasta)
TEMPDIR=tmp_${PREFIX}_minidot

mkdir -p $TEMPDIR || die "mkdir $TEMPDIR failed"

${MINIMAP2} -t16 --eqx -cx asm5 -I8G $REF $ASM > ${PREFIX}.paf || die "minimap2 failed"

${CORNETTO} fixasm ${ASM} ${PREFIX}.paf --report ${PREFIX}.report.tsv -w $TEMPDIR/${PREFIX}.fix.paf > $TEMPDIR/${PREFIX}.fix.fasta || die "cornetto failed"

${CORNETTO} minidot $TEMPDIR/${PREFIX}.fix.paf -f 2  > ${PREFIX}.eps || die "minidot failed"

echo "yey, all done for minidotplot"
