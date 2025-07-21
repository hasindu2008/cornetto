# Cornetto

Cornetto is a method for adaptive genome assembly using nanopore sequencing from Oxford Nanopore Technologies (ONT). This repository documents the Cornetto bioinformatics protocol and contains the source code for Cornetto (a programme written in C and a collection of shell scripts).

Documentation: https://hasindu2008.github.io/cornetto <br>

[![CI](https://github.com/hasindu2008/cornetto/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/hasindu2008/cornetto/actions/workflows/c-cpp.yml)

## Publication and data

- [Preprint link](https://doi.org/10.1101/2025.03.31.646505)
- [Talk video link](https://youtu.be/ci0OoM6VbsA)
- [Datasets](docs/data.md)

## Cornetto Bioinformatics protocol

See [here](docs/protocol.md)

## Cornetto Toolkit

See [here](docs/toolkit.md)

## Acknowledgement

- cornetto uses [klib](https://github.com/attractivechaos/klib) which is under the MIT license.
- `minidot` from [miniasm](https://github.com/lh3/miniasm) package (MIT license) has been integrated as `cornetto minidot` at [src/minidot](src/minidot).
- `sdust` programme from https://github.com/lh3/sdust is integrated as `cornetto sdust` at [src/sdust](src/sdust).
- `telofind`, `telobreaks` and `telowin` subcommands are implemented based on [VGP telomere scripts](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere) license under BSD.





