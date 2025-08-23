# Cornetto Toolkit

## Compiling the Cornetto C programme

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

## Using helper scripts

### dotplot

To generate the dotplot use the `scripts/minidotplot.sh`. This script requires minimap2, samtools.

```bash
scripts/minidotplot.sh hg002v1.0.1_pat.fasta asm.fasta
```

Output:
- asm.tmp.paf: minimap2 alignment
- asm.report.tsv: which chromosome is the best match per each contig, and if needs to reverse complement
- asm.fix.tmp.paf: renamed and reverse complemented (when necessary)
- asm.tmp.renamed.fasta: renamed and reverse complemented (when necessary)
- *asm.eps*: the dot plot


### Telomere stats

To get the telomere statistics use the `scripts/telostats.sh`. This script uses cornetto subcommands that implements  functionality of teleomere analysis scripts from the [VGP project](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere).

```bash
scripts/telostats.sh asm.fasta
```

Output:
- asm.windows.0.4.50kb.ends.bed: telomeres at ends of contigs
- stdout: counts of number of contigs with 1 telemore at the end, 2 telomores at the end

### Assembly stats

To get per-chromosome statistics use the `scripts/asmstats.sh`. Make sure you have already run `scripts/minidotplot.sh` and `scripts/telostats.sh` before running this script. This is because the files generated in those steps are reused by this script.

```bash
scripts/asmstats.sh asm.fasta
```

Output is explained [here](asmstats.md)


## Using individual commands

TBD


## Usage of C programme

Our Cornetto C programme contains a number of subtools that are used by the above explained scripts. If you want to use those subtools in your scripts, see the [manual page](command.md).


