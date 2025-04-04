# Cornetto

Cornetto is a method for adaptive genome assembly using nanopore sequencing from Oxford Nanopore Technologies (ONT). This repository documents the Cornetto bioinformatics protocol and contains the source code for Cornetto (a programme written in C and a collection of shell scripts).

## Table of Contents

- [Prerequisites](#prerequisites)
- [Creating the base assembly and initial Cornetto panel](#creating-the-base-assembly-and-initial-cornetto-panel)
- [Running a Cornetto iteration](#running-a-cornetto-iteration)
- [Evaluating assemblies](#evaluating-assemblies)
- [Additional refinements](#additional-refinements)
- [Usage of C programme](#usage-of-c-programme)
- [Notes](#notes)
- [Acknowledgement](#acknowledgement)

## Prerequisites

* This documentation assumes that you are well-versed in nanopore bioinformatics and genome assembly. Thus, not every tiny step is explained in detail.
* It is expected that you have an ONT nanopore sequencer connected to a Linux host installed with MinKNOW software.
* It is expected that you have read the Cornetto manuscript thoroughly, as information there has not been repeated here.
* For adaptive sampling, we use the open-source [readfish](https://github.com/LooseLab/readfish) software, so get it installed on the host running MinKNOW. ONT MinKNOW's inbuilt adaptive sampling also should work, but we relied on the open-source readfish.
* Get the Cornetto C programme compiled as instructed in the section [below](#compiling-the-cornetto-c-programme).
* For assembling the genomes, you will need [hifiasm](https://github.com/chhylp123/hifiasm). For ONT-only assemblies, make sure you have a newer version (>= 0.22.0) that supports ONT data (`--ont` option).
* For creating the Cornetto readfish panels, we need following additional software.
    - [gfatools](https://github.com/lh3/gfatools)
    - [minimap2](https://github.com/lh3/minimap2/)
    - [samtools](https://www.htslib.org/download/)
    - [bedtools](https://github.com/arq5x/bedtools2)
    - [seqkit](https://bioinf.shenwei.me/seqkit) (only for ONT-only simplex)
    - [centrifuge](https://ccb.jhu.edu/software/centrifuge) (only for saliva samples)
- For evaluating assemblies, you may use methods of your choice. Our suggested method requires following additional software:
    - [quast](https://quast.sourceforge.net)
    - [compleasm](https://github.com/huangnengCSU/compleasm)
    - [yak](https://github.com/lh3/yak)

### Compiling the Cornetto C programme

Building the Cornetto C programme requires a compiler that supports C99 standard (with X/Open 7 POSIX 2008 extensions), which is widely available. To build:

```bash
sudo apt-get install zlib1g-dev   #install zlib development libraries
git clone https://github.com/hasindu2008/cornetto
cd cornetto
make
```
The commands to zlib development libraries on some popular distributions :
```bash
On Debian/Ubuntu : sudo apt-get install zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install zlib-devel
On OS X : brew install zlib
```

## Creating the base assembly and initial Cornetto panel

### Step 1: Generating the base assembly

Use hifiasm to generate the base assembly. Depending on the type of your input data for the initial assembly, pick the example command below. Make sure to change the parameters such as the number of threads and the genome size as necessary. If your sample is saliva, make sure you first remove non-human reads by following steps [here](docs/saliva.md#removing-non-human-reads-from-saliva-samples)

```bash
# If based on PacBio HiFi data
hifiasm -t 48 --hg-size 3g -o asm-0 reads-0.fastq
# If based on ONT simplex data (R10.4.1 flowcell, LSK114 kit, super accuracy basecalls)
hifiasm --ont -t 48 --hg-size 3g -o asm-0 reads-0.fastq
```

For ONT simplex data, the command we used for basecalling was:
```bash
slow5-dorado basecaller -x cuda:all dna_r10.4.1_e8.2_400bps_sup@v5.0.0 reads-0.blow5 --emit-fastq --min-qscore 10  > reads-0.fastq
```

Once the assembly is done, convert the assemblies from GFA format to FASTA format:

```bash
# primary assembly
gfatools gfa2fa asm-0.bp.p_ctg.gfa > asm-0.fasta

# haplotype assemblies
# only needed for diploid assemblies with ONT simplex data
gfatools gfa2fa asm-0.bp.hap1.p_ctg.gfa > asm-0.hap1.fasta
gfatools gfa2fa asm-0.bp.hap2.p_ctg.gfa > asm-0.hap2.fasta
```

You may want to evaluate the quality of the assembly to see if it is sane. For this, please refer to the section on [evaluation](#evaluating-assemblies).


### Step 2: Generating the Cornetto panel

Align the starting input FASTQ reads back to the primary assembly we generated:
```bash
# if PacBio-based base assembly
minimap2 -t 24 --secondary=no --MD -ax map-hifi asm-0.fasta reads-0.fastq -o asm-0.realigned.sam
# if ONT simplex-based base assembly
minimap2 -t 24 --secondary=no --MD -ax map-ont asm-0.fasta reads-0.fastq -o asm-0.realigned.sam

# sort and index
samtools sort -@ 24 asm-0.realigned.sam -o asm-0.realigned.bam
samtools index asm-0.realigned.bam
```

Get the per base coverage information for total alignments (mapq>=0) and unique alignments (mapq>=20)
```bash
cornetto fa2bed asm-0.fasta | sort -k3,3nr > asm-0.chroms.bed
samtools depth -@ 24 -b asm-0.chroms.bed -aa asm-0.realigned.bam  | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > asm-0.cov-total.bg
samtools depth -@ 24 -Q 20 -b asm-0.chroms.bed -aa asm-0.realigned.bam  | awk '{print $1"\t"$2-1"\t"$2"\t"$3}' > asm-0.cov-mq20.bg
```

Now to create the initial Cornetto panel, you can launch the script at [scripts/create-cornetto.sh](scripts/create-cornetto.sh):

```bash
scripts/create-cornetto.sh asm-0.fasta
```

See comments inside [scripts/create-cornetto.sh](scripts/create-cornetto.sh) to understand what the script does.
Running this script will generate two files `asm-0.boringbits.bed` and `asm-0.boringbits.txt`.


### Step 3: Generating the haplotype-aware Cornetto panel

This step is only required for diploid assemblies using ONT simplex data (ONT simplex for both the base assembly and Cornetto iterations). Note that you should have already run step 2 above, before running this step. For primary assembly using PacBio data for the base assembly, followed by ONT duplex Cornetto iterations, this step can be skipped.

Run the script at [scripts/create-hapnetto.sh](scripts/create-hapnetto.sh):

```
scripts/create-hapnetto.sh asm-0
```

See comments inside [scripts/create-hapnetto.sh](scripts/create-hapnetto.sh) to understand what the script does. The final outputs we want are the two files `asm-0_dip.boringbits.bed` and `asm-0_dip.boringbits.txt`.

### Step 4: Only for human saliva samples

You need to append any non-human contigs to the panel if we are using a human saliva sample. First, follow the additional instructions [here](docs/saliva.md#get-the-non-human-contigs) to generate two files `asm-all-0.nonhuman_contigs.fasta` and `asm-all-0.nonhuman_contigs.bed`.

Then, append `asm-all-0.nonhuman_contigs.fasta` to primary assembly `asm-0.fasta`. Append the `asm-all-0.nonhuman_contigs.bed` to the `asm-0.boringbits.bed` or `asm-0_dip.boringbits.bed` based on what you are after. Generate the `asm-0.boringbits.txt` or `asm-0_dip.boringbits.txt`based on the bed file.


### Step 5: Configuring readfish

Now create the minimap2 index for the primary assembly `asm-0.fasta` to be used with readfish:
```bash
minimap2 -x map-ont asm-0.fasta -d asm-0.fasta.idx
```

Create a readfish toml file named `asm-1.boringbits.toml` as per the example below. The example below assumes you are using the diploid assembly panel `asm-0_dip.boringbits.txt` for readfish. For using the primary assembly panel (for PacBio-based base assembly + ONT duplex Cornetto), change it to `asm-0.boringbits.txt`.

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

## Running a Cornetto iteration

### Step 1: Running adaptive sampling using readfish

Do the sequencing and run readfish with the `asm-1.boringbits.toml` we created above. Example command is given below. You will have to change parameters as appropriate.

```bash
readfish targets --device ${DEVICE_ID} --experiment-name asm-1 --toml asm-1.boringbits.toml --port 9502 --cache-size 3000 --batch-size 3000 --channels 1 3000 --log-file my.log
```

### Step 2: Basecalling your data

If you are using PacBio base + ONT duplex Cornetto, basecall your data with the duplex super accuracy basecalling model. Then extract only the duplex reads in to a file called `reads-1.fastq`. The command we used for duplex basecalling was:
```bash
slow5-dorado duplex dna_r10.4.1_e8.2_400bps_sup@v4.2.0 reads-1.blow5  > reads-1.bam
```

If you are using ONT simplex for the base and Cornetto, basecall your data with the simplex super accuracy basecalling model. Then extract reads which are longer than a threshold (30 kbases) into a file called `reads-1.fastq`. The commands we used for simplex basecalling was:
```bash
slow5-dorado basecaller -x cuda:all dna_r10.4.1_e8.2_400bps_sup@v5.0.0 reads-1.blow5 --emit-fastq --min-qscore 10  > reads-1_all.fastq
seqkit seq -m 30000 reads-1_all.fastq -o reads-1.fastq
```

### Step 3: Assembling

Now launch hifiasm with the base FASTQ and the FASTQ from the Cornetto iteration. Commands are very similar to what we used before when generating the base assembly:

```bash
# if PacBio base + ONT duplex Cornetto
hifiasm -t 48 --hg-size 3g -o asm-1 reads-0.fastq reads-1.fastq
# if ONT simplex for the base and Cornetto
hifiasm --ont -t 48 --hg-size 3g -o asm-1 reads-0.fastq reads-1.fastq
```

Now convert the assemblies from GFA format to FASTA, using commands similar to what we used for base assembly. Let us say out primary assembly is `asm-1.fasta` and the haplotype assemblies are `asm-1.hap1.fasta` and `asm-1.hap2.fasta`.

You may want to evaluate the quality of the assembly to see if it has improved. For this, please refer to the section on [evaluation](#evaluating-assemblies).

### Step 4: Create the Cornetto panel for the new Cornetto iteration

Now to create the Cornetto panel for the next iteration, you can launch the script at [scripts/recreate-cornetto.sh](scripts/recreate-cornetto.sh):

```bash
scripts/recreate-cornetto.sh asm-1.fasta
```

See comments inside [scripts/recreate-cornetto.sh](scripts/create-cornetto.sh) to understand what the script is doing.
Running this script will generate two files `asm-1.boringbits.bed` and `asm-1.boringbits.txt`.


### Step 5: Creating the diploid Cornetto panel for the new Cornetto iteration

This step is only required for diploid assemblies using ONT simplex data (ONT simplex for both the base assembly and Cornetto iterations). Note that you should have already run step 4 above, before running this step. For primary assembly using PacBio data for the base assembly, followed by ONT duplex Cornetto iterations, this step can be skipped.

Run the script at [scripts/recreate-hapnetto.sh](scripts/recreate-hapnetto.sh):

```
scripts/recreate-hapnetto.sh asm-1
```

See comments inside [scripts/recreate-hapnetto.sh](scripts/recreate-hapnetto.sh) to understand what the script is doing. The final outputs we want are the two files `asm-1_dip.boringbits.bed`  and `asm-1_dip.boringbits.txt`.

### Step 6: Only for human saliva samples

Append the non-human contigs we generated in [step 4 of creating the base assembly](#step-4-only-for-human-saliva-samples), into the `asm-1.fasta`, `asm-1.boringbits.bed`/`asm-1_dip.boringbits.bed` and `asm-1.boringbits.txt`/`asm-1_dip.boringbits.txt`

### Step 7: Configuring readfish

Just as before when configuring readfish for the base assembly [above](#step-5-configuring-readfish), now create the minimap2 index for the primary assembly `asm-1.fasta` to be used with readfish. Then similarly create a readfish configuration for the next iteration (`asm-2.boringbits.toml`) with `asm-1_dip.boringbits.txt` or `asm-1_dip.boringbits.txt` created with step 4/5 above.

Now repeat from step 1 [above](#running-a-cornetto-iteration) to start another Cornetto iteration (asm-2, asm-3, and so on). For assembly at each iteration, use the base FASTQ file along with the FASTQ files from all previous iterations.


## Evaluating assemblies

In this section we document the methods you may use for evaluating the assemblies and some scripts are provided where relevant.

### Evaluating BUSCO completeness

We used [compleasm](https://github.com/huangnengCSU/compleasm/) for getting the BUSCO scores. The command we used was:
```
compleasm run -a asm.fasta -o asm.busco_out -t 8 -l ${LINEAGE} -L /path/to/mb_downloads
# LINEAGE: primates for human, actinopterygii_odb10 for cichlid, tetrapoda_odb10 for birds and turtles
```

### Evaluating using Quast

You can use quast as follows:
```
quast.py -t 24 -o asm.quast_out -l asm.fasta --large asm.fasta
```

### Additional evaluations for HG002

If the sample is HG002, there is high quality reference data that we can use to evaluate our assemblies.

First download the Q100 HG002 assembly and extract the two haplotypes:

```bash
# download and index
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.0.1.fasta.gz
samtools faidx hg002v1.0.1.fasta.gz

# paternal haplotype
grep "PATERNAL\|chrEBV\|chrM\|chrX\|chrY" hg002v1.0.1.fasta.gz.fai | cut -f 1 > paternal.txt
samtools faidx hg002v1.0.1.fasta.gz -r paternal.txt -o hg002v1.0.1_pat.fasta

# maternal haplotype
grep "MATERNAL\|chrEBV\|chrM\|chrX\|chrY" hg002v1.0.1.fasta.gz.fai | cut -f 1 > maternal.txt
samtools faidx hg002v1.0.1.fasta.gz -r maternal.txt -o hg002v1.0.1_mat.fasta
```

To generate the dotplot use the `scripts/minidotplot.sh`. This script requires minimap2, samtools.

```bash
scripts/minidotplot.sh hg002v1.0.1_pat.fasta asm.fasta
```

To get the telomere statistics use the `scripts/telostats.sh`. This script uses cornetto subcommands that implements  functionality of [teleomere analysis scripts from the VGP project](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere).


```bash
scripts/telostats.sh asm.fasta
```

To get per-chromosome statistics use the `scripts/asmstats.sh`. Make sure you have already run `scripts/minidotplot.sh` and `scripts/telostats.sh` before running this script. This is because the files generated in those steps are reused by this script. We plan to make the functionality of this bash script into the Cornetto programme in future.

```bash
scripts/asmstats.sh asm.fasta
```

To calculate the QV for the assembly, we can use yak.

```bash
# create the yak index using the HG002 Q100 assembly, with a k-mer size of 21
yak count -k 21 -K1.5g -t 16 hg002v1.0.1.fasta.gz -o hg002v1.0.1.k21.yak
# to get the QV of primary ASM
yak qv hg002v1.0.1.k21.yak asm.fasta -t 24

# to get the QV of the combined haplotypes
cat asm.hap1.fasta asm.hap2.fasta > asm.hap1+hap2.fasta
yak qv hg002v1.0.1.k21.yak asm.hap1+hap2.fasta -t 24
```

To calculate the Hamming and switch error, we can use yak too. But first we need parental yak indexes for HG002. We downloaded them from the [human-pangenomics project](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/6040D518-FE32-4CEB-B55C-504A05E4D662--HG002_PARENTAL_YAKS/HG002_PARENTS_FULL/) mentioned [here](https://github.com/chhylp123/hifiasm/issues/262#issuecomment-1081303007). After doing so, you can do:

``` bash
yak trioeval pat.HG003.yak mat.HG004.yak asm.hap1+hap2.fasta -t 16
```

## Additional refinements

If you are focused on primary assemblies, you may use following the approach documented [here](docs/refine.md) for further refinements.

## Usage of C programme

Our Cornetto C programme contains a number of subtools that are used by the above explained scripts. If you want to use those subtools in your scripts, see the [manual page](docs/command.md).

## Notes

- Our scripts and the C programme is not tested on non-Linux platforms, so might need some adjustments.
- A very useful article that explains various assembly-related terms: [concepts-in-phased-assemblies](https://lh3.github.io/2021/04/17/concepts-in-phased-assemblies)

## Acknowledgement

- cornetto uses klib [https://github.com/attractivechaos/klib] which is under the MIT license
- `minidot` from [miniasm](https://github.com/lh3/miniasm) package (MIT license) has been integrated as `cornetto minidot` at [src/minidot](src/minidot).
- `sdust` programme from https://github.com/lh3/sdust is integrated as `cornetto sdust` at [src/sdust](src/sdust).
- `telofind`, `telobreaks` and `telowin` subcommands are implemented based on [VGP telomere scripts](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere)





