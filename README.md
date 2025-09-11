# Cornetto

Cornetto is a method for iterative genome assembly using nanopore sequencing from Oxford Nanopore Technologies (ONT). This repository documents the Cornetto bioinformatics protocol and the Cornetto toolkit (a programme written in C and a collection of shell scripts). Cornetto toolkit also features some useful commands for evaluating assemblies generated from any other method.

**Cornetto is under development and there could be interface changes and changes to default parameters. Do not hesitate to open an issue if you found a bug, something is not clear or for any feature requests.**

Documentation: https://hasindu2008.github.io/cornetto <br>

[![CI](https://github.com/hasindu2008/cornetto/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/hasindu2008/cornetto/actions/workflows/c-cpp.yml)

## Publication and data

- [Preprint link (bioxiv)](https://doi.org/10.1101/2025.03.31.646505)
- [Talk video link (London Calling 2025)](https://youtu.be/ci0OoM6VbsA)
- [Raw data and assemblies for the preprint](docs/data.md)

## Cornetto bioinformatics protocol

### Detailed documentation

See [here](docs/protocol.md).

### Summary

1. Generate a base assembly: Sequence on a single flowcell as normal and use hifiasm to assemble.
2. Generate a cornetto panel for adaptive sampling: Identify the 'boringbits' from the assembly and create a adaptive sampling target panel for rejecting reads.
3. Run adaptive sampling using readfish with the panel created in step 2 for around 24 hours.
4. Create the next iteration of assembly: Use collective data from step 1 and step 3 to create an assembly with hifiasm.
5. Repeat from step 2.

## Cornetto toolkit

### Detailed documentation

See [here](docs/toolkit.md).


### Summary

Using helper scripts:

```bash
scripts/telostats.sh asm.fasta # prints the telomere counts
scripts/minidotplot.sh ref.fasta asm.fasta # creates reference vs assembly dotplot. output: assembly.eps
scripts/asmstats.sh asm.fasta  # prints chromosome-wise assembly to reference report
```

Using individual commands:

```bash
# creating a dot plot
## 1. align assembly to reference
minimap2 --eqx -cx asm5 ref.fasta asm.fasta > asm.paf
## 2. fix the +/- directions to match the reference
cornetto fixasm asm.fasta asm.paf -r asm.report.tsv -w asm.fix.paf > asm.fix.fasta
## 3. dot plot
cornetto minidot asm.fix.paf -f 2 > asm.eps

# per-chromosome evaluation. Note: asm.windows.0.4.50kb.ends.bed is from `scripts/telostats.sh`
cornetto asmstats asm.paf asm.windows.0.4.50kb.ends.bed -r asm.report.tsv -s ref.fasta

# miscellaneous commands
cornetto fa2bed asm.fasta > asm.bed  # create a bed file with assembly contig lengths
cornetto seq -m 10000 reads.fastq > long.fastq # extract reads >=10kb
```

## Acknowledgement

- cornetto uses [klib](https://github.com/attractivechaos/klib) (MIT licensed).
- `cornetto minidot` is simply `minidot` from [miniasm](https://github.com/lh3/miniasm) (MIT licensed).
- `cornetto sdust` is simply [sdust](https://github.com/lh3/sdust).
- `cornetto telofind`, `telobreaks` & `telowin` are implemented based on [VGP](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere) telomere scripts (BSD licensed).


## Fun fact

Regarding the name – the method is named by Ira Deveson after the ‘Cornetto’ ice-cream cone. There was a [famous TV advertisement in Australia](https://www.youtube.com/watch?v=NInHrq02LmI&ab_channel=terryschopshop) in 2000s that promoted the Cornetto as having ‘no boring bits’. In our cornetto workflow, we reject the 'boring bits', the bits which are too easy for the assemblers to resolve.

![cornetto](cover.png)

