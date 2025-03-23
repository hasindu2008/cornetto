## Usage of the C programme

### noboringbits

**This programme loads the whole depth file to memory, thus would need tens of gigabytes of RAM. Not memory optimised because the assembly process anyway needs a couple of hundred gigabytes of RAM and thus the user would have access to a computer with large amount of RAM.
**

Options:

* `-q FILE`:       depth file with high mapq read coverage
* `-w INT`:        window size [default: 2500]
* `-i INT`:        window increment [default: 50]
* `-L FLOAT`:      low coverage threshold factor [default: 0.6]
* `-H FLOAT`:      high coverage threshold factor [default: 1.6]
* `-Q FLOAT`:      mapq low coverage threshold factor [default: 0.6]
* `-m INT`:        minimum contig length [default: 1000000]
* `-e INT`:        edge length to ignore [default: 100000]
* `-h`:            help
* `--verbose INT`: verbosity level [4]
* `--version`:     print version


Cornetto noboringbits prints coordinate windows that meet any of the following
1. contigs < 1MBase in size
2. 100kbase edge regions at each
3. Windows that meet any of the following criteria
   - windows with low coverage: <0.6x genome average
   - windows with high coverage: >1.6x genome average
   - windows with low mappability: mean MQ20 cov for window is < 0.6 x mean coverage for the window


Example usage:
```
./cornetto noboringbits test/cov-total.bg -q test/cov-mq20.bg > noboringbits.txt
```



