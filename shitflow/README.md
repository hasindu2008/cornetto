# shitflow (shell-based internode transfer flow)

## Base assembly and creating the panel

On gadi:

```
qsub -v BASE_FASTQ=/path/to/A0_XX.fastq,OUT_PREFIX=hg002-cornetto-A_1 shitflow/hifiasm-ont.pbs.sh
```

##  Adaptive assembly and recreating the panel

```
# first one
qsub -v BASE_FASTQ=/path/to/A0_XX.fastq,FISH_PREV=A_1_QGXHXX240275:A_2_QGXHXX240279,FISH_NOW=A_3_QGXHXX240283,OUT_PREFIX=hg002-cornetto-A_3 ./hifiasm-ont.pbs.sh

# onwards
qsub -v BASE_FASTQ=/path/to/A0_XX.fastq,FISH_PREV=A_1_QGXHXX240275:A_2_QGXHXX240279,FISH_NOW=A_3_QGXHXX240283,OUT_PREFIX=hg002-cornetto-A_3 ./hifiasm-ont.pbs.sh
```

## pacbio base, followed by ont-duplex

On fridge:
```
shitflow/duplex-shitflow.sh -b /g/data/ox63/cornetto/cichlid/GSU_1.fastq.gz -o cichlid-cornetto-C_1 C_1 QGXHXX240408
shitflow/duplex-shitflow.sh -b /g/data/ox63/cornetto/cichlid/GSU_1.fastq.gz -p C_1_QGXHXX240408 -o cichlid-cornetto-C_2 C_2 QGXHXX240418
shitflow/duplex-shitflow.sh -b /g/data/ox63/cornetto/cichlid/GSU_1.fastq.gz -p C_1_QGXHXX240408:C_2_QGXHXX240418 -o cichlid-cornetto-C_3 C_3 QGXHXX240421
shitflow/duplex-shitflow.sh -b /g/data/ox63/cornetto/cichlid/GSU_1.fastq.gz -p C_1_QGXHXX240408:C_2_QGXHXX240418:C_3_QGXHXX240421  -o cichlid-cornetto-C_4 C_4 QGXHXX240440
```

On gadi use, `hifiasm.pbs.sh`, similar to above.