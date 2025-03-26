# Additional Refinements

For primary assemblies, you may use the following algorithm to further refine your assemblies.

Suppose your base assembly is called `asm-0.fasta` and the cornetto iterations are named `asm-1.fasta`, `asm-2.fasta`, `asm-3.fasta`, ..., `asm-n.fasta`.

1. Starting from `asm-1.fasta`, go through `asm-2.fasta`, `asm-3.fasta`, ..., `asm-n.fasta` until any contigs are found longer than the expected minimum chromosome size and have telomeres in both ends (let us call them near-complete chromosomes). Suppose we found such in `asm-k.fasta`. Now extract such contigs from `asm-k.fasta` into a file called `asm.fasta`.

2. Now starting from `asm-(k+1).fasta`, iterate till `asm-n.fasta` (including `asm-n.fasta`). At each iteration, map any near-complete chromosomes in the assembly to `asm-k.fasta`. Append any newly found near-complete chromosomes to `asm.fasta` - those contigs which map <50% of their length to a contig into `asm.fasta`.

3. If we are at the last iteration (`asm-n.fasta`), map any other contigs (that are not considered near-complete chromosomes) to `asm.fasta`. Append any contigs which map <50% of their length to a contig into `asm.fasta`.

4. At the end, `asm.fasta` is the final curated primary assembly.
