#!/bin/bash

LIST=H_1_QGXHXX240451:H_2_QGXHXX240452:H_3_QGXHXX240458:H_4_QGXHXX240462:H_5_QGXHXX240474:H_6_QGXHXX240476

LIST_SPACE=$(echo $LIST | tr ':' ' ')

NUM=1
for EACH in ${LIST_SPACE}
do
    CURR_LIST=$(echo $LIST | cut -d ':' -f 1-$NUM)
    echo $CURR_LIST
    mkdir ${NUM}
    cd ${NUM}
    qsub -v ASM_LIST=${CURR_LIST},ASM_NAME_PREFIX=saliva_cornetto- ~/cornetto-hasindu/cornetto/shitflow/fisht2t.pbs.sh
    cd ..
    NUM=$((NUM+1))
done



