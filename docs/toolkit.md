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

### Telomere stats

To get the telomere statistics of your assembly (vertebrates with `TTAGGG` telomere sequence) use the script `scripts/telostats.sh`. This script uses cornetto subcommands that implements  functionality of telomere analysis scripts from the [VGP project](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere). If your assembly is `asm.fasta`:

```bash
scripts/telostats.sh asm.fasta
```

Output:
- `asm.windows.0.4.50kb.ends.bed`: telomeres at the ends of contigs
- `stdout`: counts of number of contigs with 1 telomere at the end, 2 telomeres at the end, more than 2 telomeres at the end

Example:
If you run this on the [HG002 Q100](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.0.1.fasta.gz) genome, you should see 46 contigs with 2 telomeres at the end.
```bash
scripts/telostats.sh hg002v1.0.1.fasta

$ contigs with 2 telo:       46
```

### dotplot

To generate a dotplot of the mapping between your assembly and a reference genome, use the `scripts/minidotplot.sh`. This script requires minimap2, samtools. If your assembly is `ref.fasta` and the reference is `asm.fasta`:

```bash
scripts/minidotplot.sh ref.fasta asm.fasta
```

Output:
- `asm.paf`: minimap2 alignment
- `asm.report.tsv`: for each assembly contig: which chromosome in the reference is the best match; if it had to be reverse complemented to match the reference; new name
- `asm.eps`: the dot plot in eps format

Examples:

```bash
# dot plot of a hifiasm primary assembly against the chm13
awk '/^S/{print ">"$2;print $3}' asm.p_ctg.gfa > asm.fasta
scripts/minidotplot.sh chm13.fa asm.fasta

# dot plot of a hifiasm primary assembly against a haploid assembly of HG002 Q100
samtools faidx hg002v1.0.1.fasta
grep "PATERNAL\|chrEBV\|chrM\|chrX\|chrY" hg002v1.0.1.fasta.gz.fai | cut -f 1 > paternal.txt
samtools faidx hg002v1.0.1.fasta.gz -r paternal.txt -o hg002v1.0.1_pat.fasta
scripts/minidotplot.sh hg002v1.0.1_pat.fasta asm.fasta

# dot plot of hifiasm hap1+hap2 against HG002 Q100 diploid reference
awk '/^S/{print ">"$2;print $3}' asm.hap1.p_ctg.gfa > asm.hap1.fasta
awk '/^S/{print ">"$2;print $3}' asm.hap2.p_ctg.gfa > asm.hap2.fasta
cat asm.hap1.fasta asm.hap2.fasta > asm.hap1+hap2.fasta
scripts/minidotplot.sh hg002v1.0.1.fasta asm.hap1+hap2.fasta

# dot plot of a hifiasm hap1 against the chm13
scripts/minidotplot.sh hg002v1.0.1.fasta asm.hap1.fasta
```


### Assembly stats

To get per-chromosome statistics use the script `scripts/asmstats.sh`. Make sure you have already run `scripts/telostats.sh` and `scripts/minidotplot.sh` before running this script. This is because the files generated in those steps are reused by this script.

```bash
scripts/asmstats.sh asm.fasta
```

Output which is printed to the `stdout` is explained [here](asmstats.md).

Examples:
```bash
# for a primary assembly called asm.fasta against chm13.fa
scripts/telostats.sh asm.fasta
scripts/minidotplot.sh chm13.fa asm.fasta
scripts/asmstats.sh asm.fasta

# for a primary assembly called asm.fasta against a haploid hg002 we created before
scripts/telostats.sh asm.fasta
scripts/minidotplot.sh hg002v1.0.1_pat.fasta asm.fasta
scripts/asmstats.sh asm.fasta

# for a diploid assembly called asm.hap1+hap2.fasta against hg002v1.0.1.fasta diploid reference
scripts/telostats.sh asm.hap1+hap2.fasta
scripts/minidotplot.sh hg002v1.0.1.fasta asm.hap1+hap2.fasta
scripts/asmstats.sh asm.hap1+hap2.fasta

# for haplotype 1 assembly called asm.hap1.fasta against the chm13
scripts/telostats.sh asm.hap1.fasta
scripts/minidotplot.sh chm13.fa asm.hap1.fasta
scripts/asmstats.sh asm.hap1.fasta
```


## Using individual commands

What the aforementioned scripts perform, can also be done manually through individual commands.

### Generate a dot plot step by step

1. First map your assembly to the reference using minimap2.
```bash
minimap2 --eqx -cx asm5 ref.fasta asm.fasta > asm.paf
```
`ref.fasta` and `asm.fasta` are what we detailed under the `minidotplot.sh` section above. You may want to change the pre-set to asm10 and tune other minimap2 options.

2. fix the +/- directions in the paf file we generated in step 1 to match the direction on the reference
```bash
cornetto fixasm asm.fasta asm.paf -r asm.report.tsv -w asm.fix.paf > asm.fix.fasta
```
Here, the inputs are `asm.fasta` and `asm.paf`. `asm.report.tsv`, `asm.fix.paf` and `asm.fix.fasta` are outputs:

- `asm.fix.paf` will contain the `asm.paf` with +/- directions fixed to match those in the reference.
- `asm.fix.fasta` will contain the contigs in `asm.fasta` but with the directions corrected to match the reference. The contigs are also renamed here based on the chromosome for which the majority of a contig mapped to.
- `asm.report.tsv`will contain a report like the example below that tells for each assembly contig: the chromosome which the majority of the contig mapped; if the direction was swapped; and, the new name for the contig we assigned.
```
ptg000001l	chr3	+	chr3_0
ptg000002l	chr13	+	chr13_0
ptg000003l	chr15	+	chr15_0
ptg000004l	chr5	-	chr5_0
```

3. dot plot
```
cornetto minidot asm.fix.paf -f 2 > asm.eps
```
You may convert the eps file to pdf to view.


### per-chromosome evaluation

```bash
cornetto asmstats asm.paf asm.windows.0.4.50kb.ends.bed -r asm.report.tsv
```
Here, `asm.paf` and `asm.report.tsv` is what we generated in the section "Generate a dot plot step by step" above. ` asm.windows.0.4.50kb.ends.bed` is an output from `scripts/telostats.sh`.

Output which is printed to the `stdout` is explained [here](asmstats.md).

### Usage of C programme

Our Cornetto C programme contains a number of subtools in addition to the ones explained above. For more details of each and every command, see the [manual page](command.md).






