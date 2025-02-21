#!/bin/bash
SCRIPT_DIR=$(pwd)
set -o pipefail
VGP_PIPELINE=/path/to/vgp-assembly/pipeline
die(){
    echo $1
    exit 1
}

[ $# -ne 1 ] && die "Usage: $0 <file>"

FILE=$1

test -f $FILE || die "File $FILE not found"
PREFIX=$(basename $FILE .fa)
BED=${PREFIX}/${PREFIX}.windows.0.4.50kb.ends.bed

#================================================================================
#module load bedtools
#module load java/12.0.1
cpus=$SLURM_CPUS_PER_TASK # todo

mkdir -p $PREFIX
cd $PREFIX

if [[ -z $cpus ]]; then
    cpus=1
fi

threshold=0.4
ends=50000

echo "genome: $PREFIX"
echo "threshold: $threshold"
echo "ends: $ends"
echo "asm: $FILE"

if [[ -z $FILE ]]; then
    echo "Failed to find fasta file."
    exit -1
fi

#ln -s $FILE 2> /dev/null

#echo "Received $FILE"

#================================================================================================

cornetto telomere --patterns $FILE | awk '{print $1"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' - > $PREFIX.telomere
sdust $FILE > $PREFIX.sdust
test -f ${PREFIX}.fa.fai || samtools faidx ${PREFIX}.fa || die "samtools faidx on ${PREFIX}.fa failed"
awk '{print $1"\t"$2}' ${PREFIX}.fa.fai > ${PREFIX}.lens  || die "awk failed"

# lowering threshold to 0.10 (10%) from the initial 0.40 (40%)
java -cp $VGP_PIPELINE/telomere/telomere.jar FindTelomereWindows $PREFIX.telomere 99.9 0.1 > $PREFIX.windows.vgp
cornetto telomere --windows $PREFIX.telomere 99.9 0.1 > $PREFIX.windows.cornetto

java -cp $VGP_PIPELINE/telomere/telomere.jar FindTelomereBreaks $PREFIX.lens $PREFIX.sdust $PREFIX.telomere > $PREFIX.breaks.vgp
cornetto telomere --breaks $PREFIX.lens $PREFIX.sdust $PREFIX.telomere > $PREFIX.breaks.cornetto
java -cp $VGP_PIPELINE/telomere/telomere.jar FindTelomereWindows $PREFIX.telomere 99.9 $threshold > $PREFIX.windows.$threshold
