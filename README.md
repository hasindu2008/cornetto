# cornetto

Cornetto is a method for adaptive genome assembly using nanopore sequencing. This repository documents the cornetto bioinformatics protocol and contains the source code for cornetto (a collection of shell scripts and a C programme).

## Prerequisites

* This documentation assumes that you are well-versed in nanopore bioinformatics and genome assembly. Thus, not every tiny step is explained in details.
* This documentation also assumes that you have an ONT nanopore sequencer connected to a Linux host installed with MinKNOW software.
* For adaptive sampling, we will use the open-source [readfish](https://github.com/LooseLab/readfish) software, so get it installed on the host running MinKNOW. ONT MInKNOW's inbuilt adaptive sampling also should work, but we sticked to the open-source readfish.
* Get the cornetto C programme compiled as instructed in the section [below](#compiling-the-cornetto-c-programme).
* For assembling the genomes, you will need [hifiasm](https://github.com/chhylp123/hifiasm). For ONT-only assemblies, make sure you have a newer version that supports ONT data (`--ont` option).
* For the creating the cornetto readfish panels, we need the following software
    - [gfatools](https://github.com/lh3/gfatools)
    - [minimap2](https://github.com/lh3/minimap2/)
    - [samtools](https://www.htslib.org/download/)
    - [bedtools](https://github.com/arq5x/bedtools2)
- For evaluating assemblies, you may use methods that you may wish. Our method requires the following software:
    - [quast](https://quast.sourceforge.net)
    - [compleasm](https://github.com/huangnengCSU/compleasm)
    - [yak](https://github.com/lh3/yak)

### Compiling the cornetto C programme

Building the cornetto C programme requires a compiler that supports C99 standard (with X/Open 7 POSIX 2008 extensions), which is widely available. To build:

```bash
git clone https://github.com/hasindu2008/cornetto
cd cornetto
make
```

## Creating a base assembly and initial cornetto panel

### Step 1: generating the base assembly

Use hifiasm to generate the base assembly. Depending on the type of your input data for the initial assembly, pick the example command below. Make sure to change the parameters such as the number of threads and the genome size as necessary.

```bash
# If based on pacbio hifi data
hifiasm -t 48 --hg-size 3g -o asm-0 reads-0.fastq
# If based on ONT simplex data (R10.4.1 flowcell, LSK114 kit, super accuracy basecalls)
hifiasm --ont -t 48 --hg-size 3g -o asm-0 reads-0.fastq
```

Now convert the assemblies from gfa format to fasta:

```bash
# primary assembly
gfatools gfa2fa asm-0.bp.p_ctg.gfa > asm-0.fasta

# haplotype assemblies
# only needed for diploid assemblies with ONT simplex data
gfatools gfa2fa asm-0.bp.hap1.p_ctg.gfa > asm-0.hap1.fasta
gfatools gfa2fa asm-0.bp.hap2.p_ctg.gfa > asm-0.hap2.fasta
```

You may want to evaluate the quality of the assembly to see if it is sane. For this, please refer to the section on [evaluation](#evaluating-assemblies).

### Step 2: generating the cornetto panel

Align the starting input FASTQ reads back to the primary assembly we generated:
```bash
# if pacbio-based base assembly
minimap2 -t 24 --secondary=no --MD -ax map-hifi asm-0.fasta reads-0.fastq -o asm-0.realigned.sam
# if ONT simplex-based base assembly
minimap2 -t 24 --secondary=no --MD -ax map-ont asm-0.fasta reads-0.fastq

# sort and index
samtools sort -@ $24 asm-0.realigned.sam -o asm-0.realigned.bam
samtools index -o asm-0.realigned.bam
```

Get the per base coverage information for total alignments (mapq>=0) and unique alignments (mapq>=20)
```bash
samtools faidx asm-0.fasta
awk '{print $1"\t0\t"$2}' asm-0.fasta.fai | sort -k3,3nr > asm-0.chroms.bed
samtools depth -@ 24 -b asm-0.chroms.bed -aa asm-0.realigned.bam  | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > asm-0.cov-total.bg
samtools depth -@ 24 -Q 20 -b asm-0.chroms.bed -aa asm-0.realigned.bam  | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > asm-0.cov-mq20.bg
```

Now to create the initial cornetto panel, you can launch the script at [scripts/create-cornetto.sh](scripts/create-cornetto.sh):

```bash
scripts/create-cornetto.sh asm-0.fasta
```

See comments inside [scripts/create-cornetto.sh](scripts/create-cornetto.sh) to understand what the script is doing.
Running this script will generate two files `asm-0.boringbits.bed` and `asm-0.boringbits.txt`.


### Step 3: generating the haplotype-aware cornetto panel

This step is only required for diploid assemblies using ONT simplex data (ONT simplex for both the base assembly and cornetto iterations). Note that you should have already run step 2 above, before running this step. For primary assembly using pacbio data for the base assembly, followed by ONT duplex cornetto iterations, this step can be skipped.

Run the script at [scripts/create-hapnetto.sh](scripts/create-hapnetto.sh):

```
scripts/create-hapnetto.sh asm-0
```

See comments inside [scripts/create-hapnetto.sh](scripts/create-hapnetto.sh) to understand what the script is doing. The final outputs we want are the two files `asm-0_dip.boringbits.bed ` and `asm-0_dip.boringbits.txt`.

### Step 4: configuring readfish

Now create the minimap2 index for the primary assembly `asm-0.fasta` to be used with readfish:
```bash
minimap2 -x map-ont -d asm-0.fasta.idx asm-0.fasta
```

Create a readfish toml file named `asm-1.boringbits.toml` as per the example below. The example below assumes you are using the diploid assembly panel `asm-0_dip.boringbits.txt` for readfish. For using the primary assembly panel (for pacbio-based base assembly + ONT duplex cornetto), change it to `asm-0.boringbits.txt`.

```
[caller_settings]
config_name = "dna_r10.4.1_e8.2_400bps_5khz_fast_prom"
host = "ipc:///tmp/.guppy"
port = 5555
align_ref = "/path/to/asm-0.fasta.idx"

[conditions]
reference = "/path/to/asm-0.fasta.idx"

[conditions.0]
name = "asm-0_dip.boringbits"
control = false
min_chunks = 0
max_chunks = 16
targets = "asm-0_dip.boringbits.txt"
single_on = "unblock"
multi_on = "unblock"
single_off = "stop_receiving"
multi_off = "stop_receiving"
no_seq = "proceed"
no_map = "proceed"
```

## Running a cornetto iteration

### Step 1: Run adaptive sampling using readfish

Do the sequencing and run readfish with the `asm-1.boringbits.toml` we created above. Example command is given below. You will have to change parameters as appropriate.

```bash
readfish targets --device ${DEVICE_ID} --experiment-name asm-1 --toml asm-1.boringbits.toml --port 9502 --cache-size 3000 --batch-size 3000 --channels 1 3000 --log-file my.log
```

### Step 2: Basecall your data

If you are using pacbio base + ONT duplex cornetto, basecall your data with the duplex super accuracy basecalling model. Then extract only the duplex reads in to a file called `reads-1.fastq`.

If you are using ONT simplex for the base and cornetto, basecall your data with the simplex super accuracy basecalling model. Then extract reads which are longer than a threshold (30 kbases) into a file called `reads-1.fastq`.

### Step 3: Assemble

Now launch hifiasm with the base FASTQ and the fastq from the cornetto iteration. Commands are very similar to what we used before when generating the base assembly:

```bash
# if pacbio base + ONT duplex cornetto
hifiasm -t 48 --hg-size 3g -o asm-1 reads-0.fastq reads-1.fastq
# if ONT simplex for the base and cornetto
hifiasm --ont -t 48 --hg-size 3g -o asm-1 reads-0.fastq reads-1.fastq
```

Now convert the assemblies from gfa format to fasta, using commands similar to what we used for base assembly. Let us say out primary assembly is `asm-1.fasta` and the haplotype assemblies are `asm-1.hap1.fasta` and `asm-1.hap2.fasta`.

You may want to evaluate the quality of the assembly to see if it has improved. For this, please refer to the section on [evaluation](#evaluating-assemblies).

### Step 4: Create the cornetto panel for the new cornetto iteration

Now to create the cornetto panel for the next iteration, you can launch the script at [scripts/recreate-cornetto.sh](scripts/recreate-cornetto.sh):

```bash
scripts/recreate-cornetto.sh asm-1.fasta
```

See comments inside [scripts/recreate-cornetto.sh](scripts/create-cornetto.sh) to understand what the script is doing.
Running this script will generate two files `asm-1.boringbits.bed` and `asm-1.boringbits.txt`.


### Step 5: Creating the diploid cornetto panel for the new cornetto iteration

This step is only required for diploid assemblies using ONT simplex data (ONT simplex for both the base assembly and cornetto iterations). Note that you should have already run step 4 above, before running this step. For primary assembly using pacbio data for the base assembly, followed by ONT duplex cornetto iterations, this step can be skipped.

Run the script at [scripts/recreate-hapnetto.sh](scripts/recreate-hapnetto.sh):

```
scripts/create-hapnetto.sh asm-1
```

See comments inside [scripts/recreate-hapnetto.sh](scripts/recreate-hapnetto.sh) to understand what the script is doing. The final outputs we want are the two files `asm-1_dip.boringbits.bed ` and `asm-1_dip.boringbits.txt`.

### Step 6: configuring readfish

Just as before when configuring readfish for the base assembly [above](#step-4-configuring-readfish), now create the minimap2 index for the primary assembly `asm-1.fasta` to be used with readfish. Then similarly create a readfish configuration for the next iteration (`asm-2.boringbits.toml`) with `asm-1_dip.boringbits.txt` or `asm-1_dip.boringbits.txt` created with step 4/5 above.

Now repeat from step 1 [above](#running-a-cornetto-iteration) to start another cornetto iteration (asm-2, asm-3, and so on). For assembling at each iteration, use the base FASTQ file and the FASTQ files from all previous iterations.


## Evaluating assemblies

[todo]

```
scripts/minidotplot.sh ref.fasta assembly.fasta
scripts/telostat.sh assembly.fasta
scripts/asmstats.sh assembly.fasta
```


## Post analysis

[todo]

## t2t-aware iterative assembly

Launch `scripts/fisht2t.pbs.sh`.
For cumulative assemblies, you may use `scripts/postcall_all/run_fisht2t_all.sh`.


## Usage of C programme

[todo]

See [C programme commands and options](docs/command.md).


## shitflow (shell-based internode transfer flow)

The [shitflow](shitflow/README.md) directory in the repository contains the scripts we used semi-automate the execution and analysis on a combination of in-house computer systems and the Australia NCI Gadi supercomputer nodes. These scripts are for our own use with various hardcoded paths, so you better keep away. Just documenting here in case one is curious what the heck this is.


## Notes

- Our scripts and the C programme is not tested on non Linux platforms, so might need some adjustments.

## Acknowledgement

- minidot programme in src/minidot is from https://github.com/lh3/miniasm/ under the MIT license
