#!/bin/bash

die(){
    echo $1
    exit 1
}

[ $# -ne 1 ] && die "Usage: $0 <file>"

FILE=$1

test -z ${TELO_SCRIPT_PATH} && TELO_SCRIPT_PATH=/install/vgp-pipeline/telomere/telomere_analysis.sh

test -f $FILE || die "File $FILE not found"
PREFIX=$(basename $FILE .fasta)
BED=${PREFIX}/${PREFIX}.windows.0.4.50kb.ends.bed

${TELO_SCRIPT_PATH} ${PREFIX} 0.4 50000 ${FILE}
test -e ${BED} || die "telomere_analysis.sh failed"

echo -e "FILE\t${FILE}"
echo -e -n "total telomere regions at the end of contigs:\t"
cat $BED | wc -l
echo ""
echo ""
cat $BED | cut -f 1  | sort | uniq -c | awk 'BEGIN{t1=0;t2=0;t3=0}{if($1==1){t1+=1}else if($1==2){t2+=1} else {t3+=1}} END{print "telo in one end:\t"t1"\ntelo in two ends:\t"t2"\ntelo more than 2 (must be 0):\t"t3"\n"}'





