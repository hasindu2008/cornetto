# cornetto

Cornetto is a method for adaptive genome assembly using nanopore sequencing. This repository documents the cornetto bioinformatics protocol and contains the source code for cornetto (a collection of shell scripts and a C programme).

## Prerequisites

* This documentation assumes that you are well-versed in nanopore bioinformatics and genome assembly. Thus, not every tiny step is explained in details.
* This documentation also assumes that you have an ONT nanopore sequencer connected to a Linux host installed with MinKNOW software.
* For adaptive sampling, we will use the open-source [readfish](https://github.com/LooseLab/readfish) software, so get it installed on the host running MinKNOW. ONT MInKNOW's inbuilt adaptive sampling also should work, but we sticked to the open-source readfish.
* Get the cornetto C programme compiled as instructed in the section [below](#compiling-the-cornetto-c-programme).
* For assembling the genomes, you will need [hifiasm](https://github.com/chhylp123/hifiasm). For ONT-only assemblies, make sure you have a newer version that supports ONT data.
* For the creating the cornetto readfish panels, we need the following software
    - gfatools
    - minimap2
    - samtools
    - bedtools
- For evaluating assemblies, you may use methods that you may wish. our method requires the following software:
    - quast
    - compleasm
    - yak

### Compiling the cornetto C programme

Building the cornetto C programme requires a compiler that supports C99 standard (with X/Open 7 POSIX 2008 extensions), which is widely available. To build:

```bash
git clone https://github.com/hasindu2008/cornetto
cd cornetto
make
```

## Creating a base assembly and initial panel

### Step 1: generating the base assembly

Use hifiasm to generate the base assembly. Depending on the type of your input data for the initial assembly, pick the example command below. Make sure to change the parameters as necessary.

```bash
# If based on pacbio hifi data
hifiasm -t 48 --hg-size 3g -o asm-0 run-0.fastq
# If based on ONT simplex data (R10.4.1 flowcell, LSK114 kit, super accuracy basecalls)
hifiasm --ont -t 48 --hg-size 3g -o asm-0 run-0.fastq
```

Now convert the assemblies from gfa format to fasta:

```bash
# primary assembly
gfatools gfa2fa asm-0.bp.p_ctg.gfa > asm-0.fasta

# haplotype assemblies
# only needed for haplotype-aware adaptive sampling with ONT simplex data for the base assembly
gfatools gfa2fa asm-0.bp.hap1.p_ctg.gfa > asm-0.hap1.fasta
gfatools gfa2fa asm-0.bp.hap2.p_ctg.gfa > asm-0.hap2.fasta
```

You may want to evaluate the quality of the assembly to see if it is sane. For this please refer to the section on [evaluation](#evaluating-assemblies).

### Step 2: generating the cornetto panel

Align starting input FASTQ reads back to the primary assembly they generated:
```bash
# if pacbio base
minimap2 -t 24 --secondary=no --MD -ax map-hifi asm-0.fasta run-0.fastq -o asm-0.realigned.sam
# if ONT base
minimap2 -t 24 --secondary=no --MD -ax map-ont asm-0.fasta run-0.fastq

# sort and index
samtools sort -@ $24 asm-0.realigned.sam -o asm-0.realigned.bam
samtools index -o asm-0.realigned.bam
```

Get per base coverage information for total alignments (MQ>=0) and unique alignments (MQ>=20)
```bash
samtools faidx asm-0.fasta
awk '{print $1"\t0\t"$2}' asm-0.fasta.fai | sort -k3,3nr > asm-0.chroms.bed
samtools depth -@ 24 -b asm-0.chroms.bed -aa asm-0.realigned.bam  | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > asm-0.cov-total.bg
samtools depth -@ 24 -Q 20 -b asm-0.chroms.bed -aa asm-0.realigned.bam  | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > asm-0.cov-mq20.bg
```

Now to create the initial adaptive sampling panel, you can launch the script at [scripts/create-cornetto.sh](scripts/create-cornetto.sh):

```bash
scripts/create-cornetto.sh asm-0.fasta
```

See [Creating a new cornetto panel](docs/create.md) or the comments inside the script to know what the script is doing.
This script will generate two files `asm-0.boringbits.bed` and `asm-0.boringbits.txt`.


### Step 3: generating the haplotype-aware cornetto panel

This step is only required for haplotype-aware adaptive assembly using an ONT base assembly. Even for such, You should have still run step 2 above.

Run the script at [scripts/create-hapnetto.sh](scripts/create-hapnetto.sh):

```
scripts/create-hapnetto.sh asm-0
```
The final outputs we want are the two files `asm-0_dip.boringbits.bed ` and `asm-0_dip.boringbits.txt`.

### Step 4: configuring readfish

Now create the minimap2 index for the primary assembly `asm-0.fasta` to be used with readfish:
```bash
minimap2 -x map-ont -d asm-0.fasta.idx asm-0.fasta
```

Create a readfish toml file named `asm-1_dip.boringbits.toml` as per the example below. The example below assumes you are using the `asm-0_dip.boringbits.txt` for readfish, if you are focused on improving only the diploid assembly. If you are only working on primary assembly for pacbio base, change it so that `asm-0.boringbits.txt` is used.

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

Step 1: Run adaptive sampling using readfish

Do the sequencing and run readfish with the `asm-1_dip.boringbits.toml` we created above. Example command is given below. You will have to change parameters as appropriate.

```bash
readfish targets --device ${DEVICE_ID} --experiment-name asm-1_dip.boringbits --toml asm-1_dip.boringbits.toml --port 9502 --cache-size 3000 --batch-size 3000 --channels 1 3000 --log-file my.log
```

Step 2:  Basecall your data

Step 3: Assemble

Step 4: Create the next cornetto panel


## Evaluating assemblies

```
scripts/minidotplot.sh ref.fasta assembly.fasta
scripts/telostat.sh assembly.fasta
scripts/asmstats.sh assembly.fasta
```

See [Evaluation](docs/eval.md) for more details.

## Post analysis

## t2t-aware iterative assembly

Launch `scripts/fisht2t.pbs.sh`.
For cumulative assemblies, you may use `scripts/postcall_all/run_fisht2t_all.sh`.

The older methods can be found under [archived](archived.md).

## Usage of C programme

See [C programme commands and options](docs/command.md).


## shitflow (shell-based internode transfer flow)

see [here](shitflow/README.md).

## Notes

- Our scripts and the C programme is not tested on non Linux platforms, so might need some adjustments.

## Acknowledgement

- minidot programme in src/minidot is from https://github.com/lh3/miniasm/ under the MIT license
