## Usage of the C programme

### noboringbits

**This programme loads the whole depth file to memory, thus would need tens of gigabytes of RAM. It is not memory-optimized because the assembly process already requires several hundred gigabytes of RAM. Therefore, the user is expected to have access to a computer with a large amount of RAM.**

Options:

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


Example usage:
```
./cornetto noboringbits test/cov-total.bg -q test/cov-mq20.bg > noboringbits.txt
```

---

### fixasm

**This programme processes a FASTA file and a PAF alignment file to fix the direction of contigs based on the total base length being more positive or negative. It outputs the corrected FASTA to `stdout` and logs missing sequences to `stderr`.**

Options:

* `<incorrect_assembly.fa>`: Input FASTA file containing the assembly to be corrected.
* `<a.paf>`: Input PAF file containing alignments of the assembly to a reference.

**Output:**

- Corrected FASTA is written to `stdout`.
- Missing sequences are logged to `stderr`.

**Algorithm:**

1. Parse the PAF file to calculate the total positive and negative alignment lengths for each contig.
2. Reverse complement contigs with a higher negative alignment length.
3. Write the corrected contigs to `stdout`.
4. Log sequences missing from the PAF file to `stderr`.

**Example usage:**

```bash
./cornetto fixasm incorrect_assembly.fa a.paf > corrected_contigs.fasta 2> missing_sequences.log
```

**Output example:**

- `stdout` (corrected FASTA):
  ```
  >contig1
  ATCGTACGATCG
  >contig2
  CGTACGATCGTA
  ```

- `stderr` (missing sequences log):
  ```
  contig3
  contig4
  ```

### bigenough

```bash
cornetto bigenough [options] <assembly.bed> <boring.bed>
```

For each contig, if the total length of the regions listed in <boring.bed> covers more than T% of the contig's total length in <assembly.bed>, include all of that contig’s regions from <boring.bed> in the output.

Options:
* `-r FILE`:  also output in readfish format to FILE
* `-T INT`:   percentage threshold to consider as sufficient boring bits on a contig [default: 50]