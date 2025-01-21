# cornetto

Cornetto is a collection of shell scripts and a C programme for adaptive assembly.

## Compiling the C programme

Building the cornetto C programme requires a compiler that supports C99 standard (with X/Open 7 POSIX 2008 extensions), which is widely available. To build:

```
git clone https://github.com/hasindu2008/cornetto
cd cornetto
make
```

## Usage of C programme

See [C programme commands and options](docs/command.md).

## Tools required for shell scripts

- minimap2
- samtools
- bedtools

## Shell scripts

### Creating a new cornetto panel

```
scripts/create-cornetto base_assembly.fasta
```

See [Creating a new cornetto panel](docs/create.md) for more details.


## t2t-aware iterative assembly

Launch `scripts/fisht2t.pbs.sh`.
For cumulative assemblies, you may use `scripts/postcall_all/run_fisht2t_all.sh`.

The older methods can be found under [archived](archived.md)

### Evaluating

```
scripts/minidotplot.sh ref.fasta assembly.fasta
scripts/telostat.sh assembly.fasta
scripts/assembly.fasta
```

See [Evaluation](docs/eval.md) for more details.

## shitflow (shell-based internode transfer flow)

see [here](shitflow/README.md).

## Acknowledgement

- minidot programme in src/minidot is from https://github.com/lh3/miniasm/ under the MIT license
