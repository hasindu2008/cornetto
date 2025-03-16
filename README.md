# cornetto

Cornetto is a method for adaptive genome assembly using nanopore sequencing. This repository documents the cornetto bioinformatics protocol and contains the source code for cornetto (a collection of shell scripts and a C programme).

## Prerequisites

* This documentation assumes that you are well-versed in nanopore bioinformatics and genome assembly. Thus, not every tiny step is explained in details.

* This documentation also assumes that you have an ONT nanopore sequencer connected to a Linux host installed with MinKNOW software.

* For adaptive sampling, we will use the open-source [readfish](https://github.com/LooseLab/readfish) software, so get it installed on the host running MinKNOW. ONT MInKNOW's inbuilt adaptive sampling also should work, but we sticked to readfish.

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

If your base assembly is based on pacbio hifi data, use the following example hifiasm command to assemble the genome. Change the parameters as necessary.

```bash
hifiasm -t 48 --hg-size 3g -o asm hifi-run-0.fastq
```

If your base assembly is based on ONT simplex data (R10.4.1 flowcell, LSK114 kit), make sure you have basecalled them with super accuracy basecalling models. The hifiasm command example is like:

```bash
hifiasm --ont -t 48 --hg-size 3g -o asm ont-run-0.fastq
```

Convert the assemblies (primary assembly and the two haplotypes) from gfa format to fasta:
```bash
gfatools gfa2fa asm.bp.p_ctg.gfa > asm.fasta
gfatools gfa2fa asm.bp.hap1.p_ctg.gfa > asm.hap1.fasta
gfatools gfa2fa asm.bp.hap2.p_ctg.gfa > asm.hap2.fasta
```

Align starting FASTQ reads back to the assembly they generated:
```bash
minimap2 -t 24 --secondary=no --MD -ax map-hifi/map-ont asm.fasta hifi-run-0.fastq/ont-run-0.fastq | samtools sort -@ $24 -o asm.realigned.bam
samtools index -o asm.realigned.bam
```

Get per base coverage tracks for total alignments (MQ>=0) and unique alignments (MQ>=20)
```bash
samtools faidx asm.fasta
awk '{print $1"\t0\t"$2}' asm.fasta.fai | sort -k3,3nr > asm.chroms.bed
samtools depth -@ 24 -b asm.chroms.bed -aa asm.realigned.bam  | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > asm.cov-total.bg
samtools depth -@ 24 -Q 20 -b asm.chroms.bed -aa asm.realigned.bam  | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > asm.cov-mq20.bg
```

Now to create the initial adaptive sampling panel, you can launch the script at [scripts/create-cornetto.sh](scripts/create-cornetto.sh):

```bash
scripts/create-cornetto.sh asm.fasta
```

See [Creating a new cornetto panel](docs/create.md) or the comments inside the script to know what the script is doing.

Running the create-cornetto script will create many temporary files, but the final outputs we want are the two files `asm.boringbits.bed` and `asm.boringbits.txt`. The `asm.boringbits.txt` file will be used for readfish, if you are focused on improving only the primary assembly.

If you want to focus on improving the diploid assembly, run the script at [scripts/create-hapnetto.sh](scripts/create-hapnetto.sh):

```
scripts/create-hapnetto.sh asm
```
The final outputs we want are the two files `asm_dip.boringbits.bed ` and `asm_dip.boringbits.txt`. The `asm_dip.boringbits.txt` file will be used for readfish, if you are focused on improving only the diploid assembly.


Now create the minimap2 index for `asm.fasta` to be used with readfish:
```bash
minimap2 -x map-ont -d asm.fasta.idx asm.fasta
```


Create the readfish toml file:

```

```

## Running a cornetto iteration

Run readfish




Create the readfish toml file:

```

```

Now repeat ^^

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
