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
Usage: cornetto boringbits cov-total.bg -q cov-mq20.bg

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

Currently simply prints windows that meet 1, 2 and 2
1. contigs must be > [1MBase] in size
2. excluding [100kb] edge regions at each
3. Does not fall into any of the below category
   - windows with low coverage: [<0.6x] genome average
   - windows with high coverage: [>1.6x] genome average
   - windows with low mappability: [mean MQ20 cov for window is < 0.6 x mean coverage for the window]

Only tested briefly on a tiny dataset. Could be full of bugs yet. Lazily loads the whole bloody depth file to memory, thus memory inefficient. Later, let us use the htslib to get depth from the BAM.

To merge the regions you may use bedtools as follows:
```
./cornetto boringbits test/cov-total.bg -q test/cov-mq20.bg  | cut -f 1,2,3 > all_regions.bed
bedtools sort -i all_regions.bed | bedtools merge -d 500


#### Can we please try this new merging strategy?

## merge any overlapping or adjacent (within 500bp) instances of the sliding window, then:
## remove any merged intervals that are < 10kb, then:
## extend merged intervals by 10kb in either direction, in order to capture the surrounding 'boring' flanking regions
bedtools sort -i all_regions.bed | bedtools merge -d 500 | awk '($3-$2)>=10000' | awk '{print $1"\t"$2-10000"\t"$3+10000}' > merged_windows1.bed

## then perform second merging, to merge any overlapping or adjacent (within 10kb) windows
bedtools sort -i merged_windows1.bed | bedtools merge -d 10000 > merged_windows2.bed

```
