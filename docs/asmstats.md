# asmstats output

Output will be as follows

Definition: A contig will taken to belong to a chromosome if classified so by minidotplot script above (cornetto fixasm, 50% of its length is actually aligned to the chromosome).

1. Contigs containing one or more telomeres at the end:

Using belonging contigs, a table with 4 columns below will created, with one row per each chromosome on the reference.
    - chr: chromosome
    - T2T?: `y` for two telomeres at the ends, `n otherwise. Comma separated chr per each contig belonging to the chromosome
    - NTelo: total number of telomeres in contigs belonging to the chromosome
    - Telocontiglen: comma separated contig lengths


2. Contigs whose majority is mapped to the corresponding chromosome

Using belonging contigs, a table with following columns will be created, with one row per each chromosome on the reference.
    - chr: chromosome
    - Ncontigsofsize>=KMbasealignedtochr: Number of contigs of size >= K Mbase aligned to the chr
    - %ofchrsequencecoveredbycontigsofsize>=KMbase: Pecentage of the chromosome sequence covered by contigs of size >= K Mbase

K values are 0Mbase, 0.1Mbase, 1Mbase, 5Mbase  and 10Mbase

3. LX of Contigs whose majority is mapped to the corresponding chromosome

Using belonging contigs, a table with following columns will be created, with one row per each chromosome on the reference.
- chr: chromosome
- L50: Min N major-mapping contigs that cover >= 50% of chromosome [n=X]
- L90: Min N major-mapping contigs that cover >= 90% of chromosome [n=X]
- L95: Min N major-mapping contigs that cover >= 95% of chromosome [n=X]
- L99: Min N major-mapping contigs that cover >= 99% of chromosome [n=X]
- CumCovN5: Cumulative % covered by n=5 largest contigs [n%,n%,n%,n%,n%]

4. Contigs whose majority is mapped to another chromosome

Using non-belonging contigs, a table with following columns will be created, with one row per each chromosome on the reference.
    - chr: chromosome
    - Ncontigsofsize>=KMbasealignedtochr: Number of contigs of size >= K Mbase aligned to the chr
    - %ofchrsequencecoveredbycontigsofsize>=KMbase: Pecentage of the chromosome sequence covered by contigs of size >= K Mbase
