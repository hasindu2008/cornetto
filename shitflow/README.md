# shitflow (shell-based internode transfer flow)

## Base assembly and creating the panel

On Gadi launch the following:
```
qsub -v BASE_FASTQ=/path/to/A0_XX.fastq,OUT_PREFIX=hg002-cornetto-A_1 shitflow/hifiasm-ont.pbs.sh
```

##  Adaptive assembly and recreating the panel

On Gadi if data is already copied:
```
# first one
qsub -v BASE_FASTQ=/path/to/A0_XX.fastq,FISH_NOW=A_1_QGXHXX240275,OUT_PREFIX=hg002-cornetto-A_3 ./hifiasm-ont.pbs.sh

# onwards
qsub -v BASE_FASTQ=/path/to/A0_XX.fastq,FISH_PREV=A_1_QGXHXX240275:A_2_QGXHXX240279,FISH_NOW=A_3_QGXHXX240283,OUT_PREFIX=hg002-cornetto-A_3 ./hifiasm-ont.pbs.sh
```

## pacbio base, followed by ont-duplex

See [here](duplex/README.md)



