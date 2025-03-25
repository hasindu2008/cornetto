## Usage of the C programme

### noboringbits

**This programme loads the whole depth file to memory, thus would need tens of gigabytes of RAM. It is not memory-optimized because the assembly process already requires several hundred gigabytes of RAM. Therefore, the user is expected to have access to a computer with a large amount of RAM.
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
* `--verbose INT`: verbosity level [default: 4]
* `--version`:     print version


Cornetto noboringbits prints coordinate windows that meet any of the following:
1. contigs < 1Mbase in size
2. 100kbase edge regions at each
3. Windows that meet any of the following criteria:
   - windows with low coverage: <0.6x genome average
   - windows with high coverage: >1.6x genome average
   - windows with low mappability: mean MAPQ 20 coverage for window is < 0.6 x mean coverage for the window


Example usage:
```
./cornetto noboringbits test/cov-total.bg -q test/cov-mq20.bg > noboringbits.txt
```



