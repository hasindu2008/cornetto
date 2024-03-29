# cornetto

## Building

```
sudo apt-get install zlib1g-dev   #install zlib development libraries
scripts/install-hts.sh  # download and compile the htslib
make
```

The commands to zlib __development libraries__ on some popular distributions :
```sh
On Debian/Ubuntu : sudo apt-get install zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install zlib-devel
On OS X : brew install zlib
```

## Usage

```
Usage: cornetto boringbits/funbits cov-total.bg -q cov-mq20.bg

basic options:
   -q FILE                    depth file with high mapq read coverage
   -w INT                     window size [2500]
   -i INT                     window increment [50]
   -l FLOAT                   low coverage threshold factor [0.6]
   -h FLOAT                   high coverage threshold factor [1.6]
   -L FLOAT                   mapq low coverage threshold factor [0.6]
   -m INT                     minimum contig length [1000000]
   -e INT                     edge length to ignore [100000]
   -h                         help
   --verbose INT              verbosity level [4]
   --version                  print version
```

If running `boringbits`:

It prints windows that meet 1, 2 and 2
1. contigs must be > [1MBase] in size
2. excluding [100kb] edge regions at each
3. Does not fall into any of the below category
   - windows with low coverage: [<0.6x] genome average
   - windows with high coverage: [>1.6x] genome average
   - windows with low mappability: [mean MQ20 cov for window is < 0.6 x mean coverage for the window]


To merge the windows that are overlapping or adjacent (within 500bp), you may use bedtools as follows:
```
./cornetto boringbits test/cov-total.bg -q test/cov-mq20.bg  | cut -f 1,2,3 | bedtools sort | bedtools merge -d 500 > boringbits.bed
```

If running `funbits`:


It prints windows that meet any of the following
1. contigs < [1MBase] in size
2. [100kb] edge regions at each
3. Windows that meet any of the following criteria
   - windows with low coverage: [<0.6x] genome average
   - windows with high coverage: [>1.6x] genome average
   - windows with low mappability: [mean MQ20 cov for window is < 0.6 x mean coverage for the window]
```
./cornetto funbits test/cov-total.bg -q test/cov-mq20.bg > fun.txt
```

Create a full genome assembly bed:
```
samtools faidx assembly.fa
awk '{print $1"\t0\t"$2}' assembly.fa.faidx > assembly.bed
```

1. print all interesting windows with :
   - low coverage: [<0.6x] genome average
   - high coverage: [>1.6x] genome average
   - low mappability: [mean MQ20 cov for window is < 0.6 x mean coverage for the window]

```
./cornetto funbits test/cov-total.bg -q test/cov-mq20.bg | awk '{if ($4!=".") print $1"\t"$2"\t"$3}' > 1.bed
```
2. merge these interesting windows - overlapping or adjacent (within 500bp)
```
 bedtools sort 1.bed | bedtools merge -d 500 > 2.bed
```
3. remove any merged intervals from (2) that are shorter than <10kb
```
awk '($3-$2)>=10000'  2.bed > 3.bed
```
4. extend merged intervals from (3) by +10kb in either direction
```
awk '{print $1"\t"$2-10000"\t"$3+10000}' 3.bed > funbits.bed
```

5. add new windows to the file from (4), which are 100kb intervals at the left edge and right edge of the contig

```
awk '{if(($3-$2)>100000) {print $1"\t0\t100000\n"$1"\t"$3-100000"\t"$3}}'  assembly.bed >> funbits.bed
```

6. merge any overlapping or adjacent (within 10kb) intervals from (5)
```
bedtools sort -i funbits.bed | bedtools merge -d 10000 > funbits_merged.bed
```

8. subtract merged windows from (6) from the whole genome assembly
```
bedtools subtract -a assembly.bed -b funbits_merged.bed > boringbits_tmp.bed
```

9. subtract any contigs shorter than 1Mbase
```
awk '{if(($3-$2)<1000000) print $0}' short.bed
bedtools subtract -a boringbits_tmp.bed -b funbits.bed> boringbits.bed
```


## Notes

Only tested briefly on a tiny dataset. Could be full of bugs yet. Lazily loads the whole bloody depth file to memory, thus memory inefficient. Later, let us use the htslib to get depth from the BAM.
