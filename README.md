# Cornetto

Cornetto is a method for adaptive genome assembly using nanopore sequencing from Oxford Nanopore Technologies (ONT). This repository documents the Cornetto bioinformatics protocol and contains the source code for Cornetto (a programme written in C and a collection of shell scripts).

Documentation: https://hasindu2008.github.io/cornetto <br>

[![CI](https://github.com/hasindu2008/cornetto/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/hasindu2008/cornetto/actions/workflows/c-cpp.yml)

## Publication and data

- [Preprint link](https://doi.org/10.1101/2025.03.31.646505)
- [Talk video link](https://youtu.be/ci0OoM6VbsA)
- [Datasets](docs/data.md)

## Cornetto bioinformatics protocol

### Summary

1. Generate a base assembly
2. Generating a cornetto panel for adaptive sampling
3. Run adaptive sampling
4. Create the next iteration of assembly
5. Repeat from step 2

### Detailed documentation

See [here](docs/protocol.md).

## Cornetto toolkit

### Summary

Using helper scripts:

```bash
scripts/telostats.sh asm.fasta # prints the telomere counts
scripts/minidotplot.sh ref.fasta asm.fasta # creates reference vs assembly dotplot in assembly.eps
scripts/asmstats.sh asm.fasta  # prints chromosome-wise assembly to reference report
```

Using individual commands:

```bash
# creating a dot plot
## 1. align assembly to reference
minimap2 -t16 --eqx -cx asm5 ref.fasta asm.fasta > asm.paf
## 2. fix the +/- directions to match the reference
cornetto fixasm asm.fasta asm.paf --report asm.report.tsv -w asm.fix.paf > asm.fix.fasta
## 3. dot plot
cornetto minidot asm.fix.paf -f 2 > asm.eps
```

### Detailed documentation

See [here](docs/toolkit.md).

## Acknowledgement

- cornetto uses [klib](https://github.com/attractivechaos/klib) which is under the MIT license.
- `minidot` from [miniasm](https://github.com/lh3/miniasm) package (MIT license) has been integrated as `cornetto minidot` at [src/minidot](src/minidot).
- `sdust` programme from https://github.com/lh3/sdust is integrated as `cornetto sdust` at [src/sdust](src/sdust).
- `telofind`, `telobreaks` and `telowin` subcommands are implemented based on [VGP telomere scripts](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere) license under BSD.





