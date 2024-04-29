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
   -L FLOAT                   low coverage threshold factor [0.6]
   -H FLOAT                   high coverage threshold factor [1.6]
   -Q FLOAT                   mapq low coverage threshold factor [0.6]
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

## Method to get eventual boring bits

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
 bedtools sort -i 1.bed | bedtools merge -d 500 > 2.bed
```
3. remove any merged intervals from (2) that are shorter than <10kb
```
awk '($3-$2)>=20000'  2.bed > 3.bed
```
4. extend merged intervals from (3) by +10kb in either direction
```
awk '{if($2>10000){print $1"\t"$2-10000"\t"$3+10000} else {print $0}}' 3.bed > funbits.bed
```

5. add new windows to the file from (4), which are 100kb intervals at the left edge and right edge of the contig

```
awk '{if(($3-$2)>100000) {print $1"\t0\t100000\n"$1"\t"$3-100000"\t"$3}}'  assembly.bed >> funbits.bed
```

6. merge any overlapping or adjacent (within 10kb) intervals from (5)
```
bedtools sort -i funbits.bed | bedtools merge -d 10000 > funbits_merged.bed
```

7. subtract merged windows from (6) from the whole genome assembly
```
bedtools subtract -a assembly.bed -b funbits_merged.bed > boringbits_tmp.bed
```

8. subtract any contigs shorter than 1Mbase
```
awk '{if(($3-$2)<1000000) print $0}'  assembly.bed  > short.bed
bedtools subtract -a boringbits_tmp.bed -b short.bed > boringbits.bed
```

9. if 'boring bits' are <50% of a single contig/scaffold, remove all boring bits on the whole scaffold. Use the below horrible inefficient code snippet for now

```
INPUT=boringbits_10k.bed
cut -f 1 ${INPUT}  | uniq > boring_ctg.tmp
while read p; 
do
ctg_len=$(grep "$p" assembly.bed | cut -f 3)
ctg_boring=$(grep "$p" ${INPUT} | awk '{sum+=$3-$2}END{ print sum}') 
fac=$(echo "$ctg_boring*100/$ctg_len" | bc)
if [ "$fac" -gt "50" ];then
	grep "$p" ${INPUT}
fi
done < boring_ctg.tmp > boringbits_10k_final.bed
```

## Notes

 Lazily loads the whole bloody depth file to memory, thus memory inefficient. Later, let us use the htslib to get depth from the BAM.

## Evaluation methods

A useful article https://lh3.github.io/2021/04/17/concepts-in-phased-assemblies

### Dot plots

1. Grab the HG2 

```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.0.1.fasta.gz
```

2. Minimap

```
minimap2 -t16 --eqx -cx asm5 hg002v1.0.1.fasta.gz RGBX240039_HG002.hifiasm.primary_asm.fasta > a.paf
```

3. fix the directions

```
cut -f 1 a.paf  | sort -u > ctg.list

while read p;
do 
  grep $p  a.paf | awk 'BEGIN{sump=0;sumn=0} {if($5=="-"){sumn+=($9-$8)}else{sump+=($9-$8)}} END{if(sump>sumn){print $1"\t+"}else{print $1"\t-"}}'
done < ctg.list > ctg_dir.txt

grep "+" ctg_dir.txt | cut -f 1 > ctg_plus.txt
grep "-" ctg_dir.txt | cut -f 1 > ctg_mins.txt

samtools faidx  RGBX240039_HG002.hifiasm.primary_asm.fasta -r ctg_plus.txt > RGBX240039_HG002_fixed.hifiasm.primary_asm.fasta
samtools faidx  RGBX240039_HG002.hifiasm.primary_asm.fasta -r ctg_mins.txt -i >> RGBX240039_HG002_fixed.hifiasm.primary_asm.fasta
```

4. Map again

```
minimap2 -t16 --eqx -cx asm5 hg002v1.0.1.fasta.gz RGBX240039_HG002_fixed.hifiasm.primary_asm.fasta > b.paf
```

5. Plottttt
```
/install/miniasm/minidot b.paf -f 4  > a.eps
```

### assembly QV

Using yak:

1. Create a k-mer count profile using the HG2 assemblt (we should use illumina reads instead - the gapless assembly paper does so?)
```
/install/yak/yak count -K1.5g -t32 -o hg002v1.0.1.yak hg002v1.0.1.fasta.gz
```
2. Get the QV
```
/install/yak/yak qv hg002v1.0.1.yak RGBX240039_HG002.hifiasm.primary_asm.fasta > qv.txt
```

Using Merqury:

1. Get Meryl dbs from Illumina WGS and hapmers are available (wtf is a hapmer?)

```
mkdir merydbs
wget https://obj.umiacs.umd.edu/marbl_publications/merqury/HG002/HG002.k21.meryl.tar.gz
wget https://obj.umiacs.umd.edu/marbl_publications/merqury/HG002/HG002.mat.hapmers.meryl.tar.gz
wget https://obj.umiacs.umd.edu/marbl_publications/merqury/HG002/HG002.pat.hapmers.meryl.tar.gz

tar xf HG002.k21.meryl.tar.gz
tar xf HG002.mat.hapmers.meryl.tar.gz
tar xf HG002.pat.hapmers.meryl.tar.gz

/install/merqury/merqury.sh HG002.k21.meryl HG002.mat.hapmers.meryl pat.hapmers.meryl RGBX240039_HG002.hifiasm.primary_asm.fasta RGBX240039_HG002_meryl

```

### Hamming errors

Need trio fastq it seems?

https://hifiasm.readthedocs.io/en/latest/trio-assembly.html

### Gene stuff

```
source /install/compleasm-0.2.6/venv3/bin/activate
compleasm run -a RGBX240039_HG002.hifiasm.primary_asm.fasta -o output_dir -t 8 -l primates
```

### T2T counts

```
# /install/vgp-pipeline/telomere/find_telomere.sh RGBX240039_HG002.hifiasm.primary_asm.fasta
/install/vgp-pipeline/telomere/find_telomere.sh 0.4 50000 RGBX240039_HG002.hifiasm.primary_asm.fasta telo
```
