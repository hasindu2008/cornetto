# Usage of the C programme

## create panels

### noboringbits

This programme loads the whole depth file to memory, thus would need tens of gigabytes of RAM. It is not memory-optimized because the assembly process already requires several hundred gigabytes of RAM. Therefore, the user is expected to have access to a computer with a large amount of RAM.

**Options:**

* `-q FILE`:       depth file with high mapq read coverage
* `-w INT`:        window size [default: 2500]
* `-i INT`:        window increment [default: 50]
* `-L FLOAT`:      low coverage threshold factor [default: 0.4]
* `-H FLOAT`:      high coverage threshold factor [default: 2.5]
* `-Q FLOAT`:      mapq low coverage threshold factor [default: 0.4]
* `-m INT`:        minimum contig length [default: 1000000]
* `-e INT`:        edge length to ignore [default: 100000]
* `-h`:            help
* `--verbose INT`: verbosity level [default: 4]
* `--version`:     print version


Cornetto noboringbits prints coordinate windows that meet any of the following:

1. contigs < 1Mbase in size
2. 100kbase edge regions at each
3. Windows that meet any of the following criteria:
   - windows with low coverage: < [0.4]x genome average
   - windows with high coverage: > [2.5]x genome average
   - windows with low mappability: mean MAPQ 20 coverage for window is < [0.4]x mean coverage for the window


**Example usage:**
```
./cornetto noboringbits test/cov-total.bg -q test/cov-mq20.bg > noboringbits.txt
```


### bigenough

```bash
cornetto bigenough [options] <assembly.bed> <boring.bed>
```

For each contig, if the total length of the regions listed in <boring.bed> covers more than T% of the contig's total length in <assembly.bed>, include all of that contig’s regions from <boring.bed> in the output.

Options:
* `-r FILE`:  also output in readfish format to FILE
* `-T INT`:   percentage threshold to consider as sufficient boring bits on a contig [default: 50]

---

## dotplot related

### fixasm

This programme processes a FASTA file and a PAF alignment file to fix the direction of contigs based on the total base length being more positive or negative. It outputs the direction fixed FASTA to `stdout` and logs missing sequences to `stderr`.

**Input:**

* `<assembly.fa>`: Input FASTA file containing the assembly to be fixed for contig direction.
* `<asm_to_ref.paf>`: Input PAF file containing alignments of the assembly to a reference.

**Output:**

* fixed FASTA is written to `stdout`.

**Options:**

* `-m FILE`:   write missing contig names to FILE
* `-r FILE`:   write report to FILE
* `-w FILE`:   write fixed PAF to FILE


**Algorithm:**

1. Parse the PAF file to calculate the total positive and negative alignment lengths for each contig.
2. Reverse complement contigs with a higher negative alignment length.
3. Write the fixed contigs to `stdout`.
4. Log sequences missing from the PAF file if requested.
5. Write the fixed PAF file if requested.

**Example usage:**

```bash
./cornetto fixasm assembly.fa asm2ref.paf -m missing_sequences.txt -r report.tsv -w fixed.paf > fixed_contigs.fasta
```


### minidot

This subcommand generates a dot plot from a PAF file. From https://github.com/lh3/miniasm.

**Options:**

* `-m INT`:        minimum match length [default: 100]
* `-i FLOAT`:      minimum identity [default: 0.1]
* `-s INT`:        minimum span [default: 1000]
* `-w INT`:        image width [default: 600]
* `-f INT`:        font size [default: 11]
* `-L`:            do not print labels
* `-D`:            do not align hits to the diagonal

**Example usage:**

```bash
cornetto minidot -m 500 -i 0.9 -s 2000 -w 800 input.paf > output.eps
```

---

## Telemere evaluation


### telowin

This subcommand analyses telomere windows in a genome assembly.

**Options:**

* `<input_file>`:  Input file containing telomere regions.
* `<identity>`:    Identity percentage (e.g., 99.9).
* `<threshold>`:   Threshold for telomere detection.

**Example usage:**

```bash
cornetto telowin input.telomere 99.9 0.4 > output.windows
```


### telobreaks

This subcommand identifies telomere breaks in a genome assembly.

**Inputs:**

* `<lens_file>`:   File containing contig lengths.
* `<sdust_file>`:  File containing low-complexity regions.
* `<telomere_file>`: File containing telomere regions.

**Example usage:**

```bash
cornetto telobreaks assembly.lens assembly.sdust assembly.telomere > output.breaks
```


### telofind

This subcommand identifies telomere sequences in a FASTA file.

**Options:**

* `<input.fasta>`: Input FASTA file.
* `[sequence]`:    Optional sequence to search for (default: `TTAGGG`).

**Example usage:**

```bash
cornetto telofind input.fasta > output.telomere
```

### sdust

This subcommand identifies low-complexity regions in a FASTA file using the symmetric DUST algorithm. From https://github.com/lh3/sdust.

**Options:**

* `-w INT`:        window size [default: 64]
* `-t INT`:        threshold [default: 20]

**Example usage:**

```bash
cornetto sdust -w 64 -t 20 input.fa > output.sdust
```

## Misc

### fa2bed

create a bed file with assembly contig lengths

**Example usage:**
```
cornetto asmbed <assembly.fasta>
```

### seq

extract reads equal to larger than a threshold from a fastq

**Options:**
- `-m INT`:        min length [30000]

**Example usage:**
```
cornetto seq <reads.fastq>
```