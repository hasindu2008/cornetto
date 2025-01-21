#!/bin/bash

usage (){
    echo "Usage: $0 <D_1_PGXXXX240596:D_2_PGXXXX250005:D_3_PGXXXX250015> <saliva_cornetto>"
    exit 1
}

die() {
    echo "$1"
    exit 1
}

test $# -eq 2 || usage

LIST=$1
PREFIX=$2

FISHT2T_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/shitflow/fisht2t.pbs.sh

LIST_SPACE=$(echo $LIST | tr ':' ' ')

NUM=1
for EACH in ${LIST_SPACE}
do
    CURR_LIST=$(echo $LIST | cut -d ':' -f 1-$NUM)
    echo $CURR_LIST
    mkdir ${NUM} || die "Could not create directory ${NUM}"
    cd ${NUM} || die "Could not cd to ${NUM}"
    qsub -v ASM_LIST=${CURR_LIST},ASM_NAME_PREFIX=${PREFIX}- ${FISHT2T_SCRIPT} || die "Could not submit job"
    cd .. || die "Could not cd back"
    NUM=$((NUM+1))
done



