# Shitflow (**Sh**ell-based **i**nternode **t**ransfer **f**low)

This directory contains the scripts we used semi-automate the execution and analysis on a combination of in-house computer systems and the Australia NCI Gadi supercomputer nodes. These scripts are for our own use with various hardcoded paths, so you better keep away. Just documenting here in case one is curious what the heck this.


## Base assembly and creating the first panel

On Gadi launch the following:
```
qsub -v BASE_FASTQ=/path/to/D_0_PGXXSX240470.fastq,OUT_PREFIX=hg002-cornetto-D_0 shitflow/hifiasm-ont.pbs.sh
```

If saliva, On gadi:
```
# get human-only reads
qsub -v FASTQ=/path/to/Q_O_PGXXXX250135_saliva.fastq shitflow/saliva/get-human-reads.pbs.sh
# assembly with human-only reads
qsub -v BASE_FASTQ=/path/to/Q_O_PGXXXX250135_saliva_human_reads.fastq,OUT_PREFIX=saliva-Q_0 shitflow/hifiasm-ont.pbs.sh

# assembly with all reads
qsub -v BASE_FASTQ=/path/to/Q_O_PGXXXX250135_saliva.fastq,OUT_PREFIX=saliva-Q_0 shitflow/hifiasm-ont.pbs.sh
# get all non-human contigs
qsub -v ASM_PREFIX=saliva-Q_0,FASTQ_PREFIX=~/cornetto-hasindu/data/saliva/Q_O_PGXXXX250135_pass_sup_500_saliva shitflow/saliva/get-nonhuman-contigs.pbs.sh

# combine
qsub -v HUMAN=human-only/,NONHUMAN=all-crap/centrifuge/,ASM=saliva-Q_0 ./create-combined-panel.pbs.sh
```

Then copy the panel and in the index to the fridge and launch fishing.

##  Adaptive assembly and recreating the panel

On fridge:

```
# first iteration
shitflow/simplex-shitflow.sh -b /gadi/path/to/D_0_PGXXSX240470.fastq -o hg002-cornetto-D_1 D_1 PGXXXX240596

# second iteration
shitflow/simplex-shitflow.sh -b /gadi/path/to/D_0_PGXXSX240470.fastq -p D_1_PGXXXX240596 -o hg002-cornetto-D_2 D_2 PGXXXX250005

# third iteration
shitflow/simplex-shitflow.sh -b /gadi/path/to/D_0_PGXXSX240470.fastq -p D_1_PGXXXX240596:D_2_PGXXXX250005 -o hg002-cornetto-D_3 D_3 PGXXXX250015
```

If saliva, following additional step needed on Gadi:
```
qsub -v ASM=saliva-cornetto-Q_1,NONHUMAN_PREFIX=/g/data/ox63/hasindu/cornetto/base_assemblies/saliva_Q_0/all-crap/centrifuge/saliva-Q_0 ~/cornetto-hasindu/cornetto/shitflow/saliva/recreate-combined-panel.pbs.sh
```

Some intermediate steps if need to intervene

On Gadi if FASTQ data is already copied to the shitflow directory `/g/data/ox63/hasindu/cornetto/shitflow`:
```
# first one
qsub -v BASE_FASTQ=/path/to/D_0_PGXXSX240470.fastq,FISH_NOW=D_1_PGXXXX240596,OUT_PREFIX=hg002-cornetto-D_1 ./hifiasm-ont.pbs.sh

# second one
qsub -v BASE_FASTQ=/path/to/D_0_PGXXSX240470.fastq,FISH_PREV=D_1_PGXXXX240596,FISH_NOW=D_2_PGXXXX250005,OUT_PREFIX=hg002-cornetto-D_2 ./hifiasm-ont.pbs.sh

# onwards
qsub -v BASE_FASTQ=/path/to/D_0_PGXXSX240470.fastq,FISH_PREV=D_1_PGXXXX240596:D_2_PGXXXX250005,FISH_NOW=D_3_QGXHXX240283,OUT_PREFIX=hg002-cornetto-D_3 ./hifiasm-ont.pbs.sh
```


On gta100 if the blow5 is already copied to the shitflow directory `/data/hasindu/shitflow/`

```
# example for second iteration
shitflow/simplex/basecall-gta100.sh D2_PGXXXX250005 BASE_FASTQ=/g/data/ox63/hasindu/cornetto/data/D_0_PGXXSX240470_pass_sup_500.fastq,FISH_PREV=D_1_PGXXXX240596,FISH_NOW=D2_PGXXXX250005,OUT_PREFIX=hg002-cornetto-D_2
```



## pacbio base, followed by ont-duplex

See [here](duplex/README.md)



