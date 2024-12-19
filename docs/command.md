## Usage of the C programme

**Lazily loads the whole bloody depth file to memory, thus memory inefficient. Later, let us use the htslib to get depth from the BAM.**

```
Usage: cornetto noboringbits cov-total.bg -q cov-mq20.bg

basic options:
   -q FILE                    depth file with high mapq read coverage
   -w INT                     window size [2500]
   -i INT                     window increment [50]
   -L FLOAT                   low coverage threshold factor [0.6]
   -H FLOAT                   high coverage threshold factor [1.6]
   -Q FLOAT                   mapq low coverage threshold factor [0.6]
   -m INT                     minimum contig length [1000000]
   -e INT                     edge length to ignore [100000]
   -h                         help
   --verbose INT              verbosity level [4]
   --version                  print version
```

It prints windows that meet any of the following
1. contigs < [1MBase] in size
2. [100kb] edge regions at each
3. Windows that meet any of the following criteria
   - windows with low coverage: [<0.6x] genome average
   - windows with high coverage: [>1.6x] genome average
   - windows with low mappability: [mean MQ20 cov for window is < 0.6 x mean coverage for the window]
```
./cornetto noboringbits test/cov-total.bg -q test/cov-mq20.bg > fun.txt
```


## deprecated

Instead of `noboringbits`, there is a deprecated `boringbits` which:

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

