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

## Cornetto 2 Pannel

1. print all interesting windows with :
   - low coverage: [<0.5x] genome average
   - high coverage: [>2.2x] genome average
   - low mappability: [mean MQ20 cov for window is < 0.5 x mean coverage for the window]

```
./cornetto funbits -H 2.2 -L 0.5 -Q 0.5 test/cov-total.bg -q test/cov-mq20.bg | awk '{if ($4!=".") print $1"\t"$2"\t"$3}' > 1.bed
```
2. merge these interesting windows - overlapping or adjacent (within 500bp)
```
 bedtools sort -i 1.bed | bedtools merge -d 500 > 2.bed
```
3. remove any merged intervals from (2) that are shorter than <20kb
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

6. merge any overlapping or adjacent (within 50kb) intervals from (5)
```
bedtools sort -i funbits.bed | bedtools merge -d 50000 > funbits_merged.bed
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
INPUT=boringbits.bed
cut -f 1 ${INPUT}  | uniq > boring_ctg.tmp
while read p; 
do
ctg_len=$(grep "$p" assembly.bed | cut -f 3)
ctg_boring=$(grep "$p" ${INPUT} | awk '{sum+=$3-$2}END{ print sum}') 
fac=$(echo "$ctg_boring*100/$ctg_len" | bc)
if [ "$fac" -gt "50" ];then
	grep "$p" ${INPUT}
fi
done < boring_ctg.tmp > boringbits_cornetto2.bed
```

### cornetto 3 panel

```
#1# get contigs longer than 1Mbase and make BED file covering whole contig excluding 200kb at each end
cat assembly.bed | awk '$3 >= 1000000' | awk '{print $1"\t"200000"\t"$3-200000}' > long_contigs.200kb-excluded.bed

#2# get regions longer than 5kb which were labelled by hifiasm as "low quality" (not sure what the definition of this is exactly)
cat RGBX240039_HG002.hifiasm.primary_asm.bp.p_ctg.lowQ.bed | awk '($3-$2)>=5000' > lowQ-5kbplus.bed

#3# exclude each low-Q region by 50kb on either side
cat lowQ-5kbplus.bed | awk '{print $1"\t"$2-50000"\t"$3+50000}' | awk -vOFS='\t' '{for(i=1;i<=NF;i++)if($i<0||$i=="-nan")$i=0}1' > lowQ-5kbplus_extended.bed

#4# subtract extended lowQ regions from the BED file above ins step #2#
bedtools subtract -a long_contigs.200kb-excluded.bed -b lowQ-5kbplus_extended.bed > boringbits_cornetto3.bed

```

### cornetto 4 panel

```
scripts/create_cornetto4.sh
```

## Notes

 Lazily loads the whole bloody depth file to memory, thus memory inefficient. Later, let us use the htslib to get depth from the BAM.


## Evaluation methods

A useful article https://lh3.github.io/2021/04/17/concepts-in-phased-assemblies


### Dot plots

1. Grab the HG2 

```
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.0.1.fasta.gz
samtools faidx hg002v1.0.1.fasta.gz
grep "PATERNAL\|chrEBV\|chrM\|chrX\|chrY" hg002v1.0.1.fasta.gz.fai | cut -f 1 > paternal.txt
grep "MATERNAL\|chrEBV\|chrM\|chrX\|chrY" hg002v1.0.1.fasta.gz.fai | cut -f 1 > maternal.txt
samtools faidx hg002v1.0.1.fasta.gz -r paternal.txt -o hg002v1.0.1_pat.fa
samtools faidx hg002v1.0.1.fasta.gz -r maternal.txt -o hg002v1.0.1_mat.fa
```

2. Minimap

```
minimap2 -t16 --eqx -cx asm5 hg002v1.0.1_pat.fa RGBX240039_HG002.hifiasm.primary_asm.fasta > a.paf
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

A shit script that I wrote quickly (ultra-inefficient) is in scripts/minidotplot.sh which you can use as `minidotplot.sh hg002v1.0.1_pat.fa RGBX240039_HG002.hifiasm.primary_asm.fasta` for instance.


### assembly QV

Get the trio data to generate the yak databases.

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

WARNING: Meryl seem to take ages to run

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
cat RGBX240039_HG002.hifiasm.primary_asm.hap1.fasta RGBX240039_HG002.hifiasm.primary_asm.hap2.fasta > RGBX240039_HG002.hifiasm.primary_asm.fasta
/install/vgp-pipeline/telomere/telomere_analysis.sh RGBX240039_HG002 0.4 50000 RGBX240039_HG002.hifiasm.primary_asm.fasta
cat telomere/RGBX240039_HG002.hifiasm.primary_asm.windows.0.4.50kb.ends.bed | cut -f 1  | sort | uniq -c | awk 'BEGIN{t1=0;t2=0;t3=0}{if($1==1){t1+=1}else if($1==2){t2+=1} else {t3+=1}} END{print "telo in one end:\t"t1"\ntelo in two ends:\t"t2"\ntelo more than 2 (must be 0):\t"t3"\n"}'
```
