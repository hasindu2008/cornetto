# asmstats output

Output from the `asmstats` subtool is explained on this page. First of all pay attention to the following definition.

Definition: A contig is considered to 'belong' to a reference chromosome if classified so by the minidotplot.sh script. minidotplot.sh uses `cornetto fixasm` subtool, which will classify so when >50% of a contig's length is aligned to the chromosome.

Output of `asmstats` contains four tables as explained below.

## 1. Contigs containing one or more telomeres at the ends:

Using belonging contigs with at least one telomere at the ends, a table with 4 columns below will created, with one row per each reference chromosome.
- chr: chromosome
- T2T?: `y` for two telomeres at the ends, `n` if only one telomere. Comma separated 'y'/'n' chacaters per each assembly contig that belongs to the chromosome
- NTelo: total number of telomeres in contigs that belong to the chromosome
- Telocontiglen: comma separated contig lengths (same order as T2T?)


## 2. Contigs whose majority is mapped to the corresponding chromosome

Using belonging contigs (irrespective of telomeres present or not), a table with following columns will be created, with one row per each reference chromosome.
- chr: chromosome
- Ncontigsofsize>=KMbasealignedtochr: Number of contigs of size >= K Mbase aligned to the reference chromosome
- %ofchrsequencecoveredbycontigsofsize>=KMbase: Percentage of the reference chromosome covered by contigs of size >= K Mbase

K values reported are 0Mbase, 0.1Mbase, 1Mbase, 5Mbase  and 10Mbase

## 3. LX of Contigs whose majority is mapped to the corresponding chromosome

Using belonging contigs (irrespective of telomeres present or not), a table with following columns will be created, with one row per each chromosome on the reference.
- chr: chromosome
- L50: Min N major-mapping contigs that cover >= 50% of chromosome [n=X]
- L90: Min N major-mapping contigs that cover >= 90% of chromosome [n=X]
- L95: Min N major-mapping contigs that cover >= 95% of chromosome [n=X]
- L99: Min N major-mapping contigs that cover >= 99% of chromosome [n=X]
- CumCovN5: Cumulative % covered by n=5 largest contigs [n%,n%,n%,n%,n%]

## 4. Contigs whose majority is mapped to another chromosome. Higher values in this table means potential misassemblies.

Using non-belonging contigs (irrespective of telomeres present or not), a table with following columns will be created, with one row per each chromosome on the reference.
- chr: chromosome
- Ncontigsofsize>=KMbasealignedtochr: Number of contigs of size >= K Mbase aligned to the chr
- %ofchrsequencecoveredbycontigsofsize>=KMbase: Pecentage of the chromosome sequence covered by contigs of size >= K Mbase
