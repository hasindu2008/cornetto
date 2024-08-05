
# protocol for basecall 

## On fridge: run inside a screen session

PREFIX=A_1
SAMPLE=QGXHXX240275

cd /data3/cornetto
slow5tools merge /data/${SAMPLE}/*/*/slow5/ -o ${PREFIX}_${SAMPLE}.blow5
slow5tools stats ${PREFIX}_${SAMPLE}.blow5

## on brenner

cd /directflow/KCCGGenometechTemp/projects/iradev/operation_cornetto/temp_hasindu
mkdir ${PREFIX}_${SAMPLE} && cd ${PREFIX}_${SAMPLE}
qsub ./fridge_dload.sge.sh /data3/cornetto/${PREFIX}_${SAMPLE}.blow5 ./
qsub ../../split_dorado_duplex.sge.sh ${PREFIX}_${SAMPLE}
