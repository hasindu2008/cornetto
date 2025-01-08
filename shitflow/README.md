# shitflow (shell-based internode transfer flow)

## Base assembly and creating the first panel

On Gadi launch the following:
```
qsub -v BASE_FASTQ=/path/to/D_0_PGXXSX240470.fastq,OUT_PREFIX=hg002-cornetto-D_1 shitflow/hifiasm-ont.pbs.sh
```

Then copy the panel and in the index to the fridge and launch fishing.

##  Adaptive assembly and recreating the panel

On fridge:

```
# first iteration
shitflow/simplex-shitflow.sh -b /gadi/path/to/D_0_PGXXSX240470.fastq -o hg002-cornetto-D_1 D1 PGXXXX240596

# second iteration
shitflow/simplex-shitflow.sh -b /gadi/path/to/D_0_PGXXSX240470.fastq -p D_1_PGXXXX240596 -o hg002-cornetto-D_2 D2 PGXXXX250005
```

Some intermediate steps if need to intervene

On Gadi if FASTQ data is already copied to the shitflow directory `/g/data/ox63/hasindu/cornetto/shitflow`:
```
# first one
qsub -v BASE_FASTQ=/path/to/D_0_PGXXSX240470.fastq,FISH_NOW=D_1_PGXXXX240596,OUT_PREFIX=hg002-cornetto-D_1 ./hifiasm-ont.pbs.sh

# second one
qsub -v BASE_FASTQ=/path/to/D_0_PGXXSX240470.fastq,FISH_PREV=D_1_PGXXXX240596,FISH_NOW=D2_PGXXXX250005,OUT_PREFIX=hg002-cornetto-D_2 ./hifiasm-ont.pbs.sh

# onwards
qsub -v BASE_FASTQ=/path/to/D_0_PGXXSX240470.fastq,FISH_PREV=D_1_PGXXXX240596:D2_PGXXXX250005,FISH_NOW=D_3_QGXHXX240283,OUT_PREFIX=hg002-cornetto-D_3 ./hifiasm-ont.pbs.sh
```


On gta100 if the blow5 is already copied to the shitflow directory `/data/hasindu/shitflow/`

```
# example for second iteration
shitflow/simplex/basecall-gta100.sh D2_PGXXXX250005 BASE_FASTQ=/g/data/ox63/hasindu/cornetto/data/D_0_PGXXSX240470_pass_sup_500.fastq,FISH_PREV=D_1_PGXXXX240596,FISH_NOW=D2_PGXXXX250005,OUT_PREFIX=hg002-cornetto-D_2
```

## pacbio base, followed by ont-duplex

See [here](duplex/README.md)



