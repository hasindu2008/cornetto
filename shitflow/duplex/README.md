
# shitflow for Duplex reads

NOTE: As some restructuring was done, the documentation below might be inconsistent. The documentation left there if would become necessary later. Duplex is not there now anyways.

To create the base assembly and panel you may use `hifiasm.pbs.sh` instead of `hifiasm-ont.pbs.sh` similar to [here](../README.md).

To directly launch an iteration from the fridge just call `shitflow/duplex-shitflow.sh`:
```
# first iteration
shitflow/duplex-shitflow.sh -b /g/data/ox63/cornetto/cichlid/GSU_1.fastq.gz -o cichlid-cornetto-C_1 C_1 QGXHXX240408
# second iteration
shitflow/duplex-shitflow.sh -b /g/data/ox63/cornetto/cichlid/GSU_1.fastq.gz -p C_1_QGXHXX240408 -o cichlid-cornetto-C_2 C_2 QGXHXX240418
# third iteration
shitflow/duplex-shitflow.sh -b /g/data/ox63/cornetto/cichlid/GSU_1.fastq.gz -p C_1_QGXHXX240408:C_2_QGXHXX240418 -o cichlid-cornetto-C_3 C_3 QGXHXX240421
# foruth iteration ...
shitflow/duplex-shitflow.sh -b /g/data/ox63/cornetto/cichlid/GSU_1.fastq.gz -p C_1_QGXHXX240408:C_2_QGXHXX240418:C_3_QGXHXX240421  -o cichlid-cornetto-C_4 C_4 QGXHXX240440
```

## Manual Steps

```
PREFIX=A_1
SAMPLE=QGXHXX240275
```

### On fridge: run inside a screen session

```
cd /data3/cornetto
slow5tools merge /data/${SAMPLE}/*/*/slow5/ -o ${PREFIX}_${SAMPLE}.blow5
slow5tools stats ${PREFIX}_${SAMPLE}.blow5
```
### on brenner
```
cd /directflow/KCCGGenometechTemp/projects/iradev/operation_cornetto/temp_hasindu
mkdir ${PREFIX}_${SAMPLE} && cd ${PREFIX}_${SAMPLE}
qsub fridge_dload.sge.sh /data3/cornetto/${PREFIX}_${SAMPLE}.blow5 ./
qsub split_dorado_duplex.sge.sh ${PREFIX}_${SAMPLE}
```

### on gadi

You may use `hifiasm.pbs.sh` instead of `hifiasm-ont.pbs.sh` similar to [here](../README.md).