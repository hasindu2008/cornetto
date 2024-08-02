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
#1# print all interesting windows with:
# low coverage: [<0.5x] genome average
# high coverage: [>2.5x] genome average
# low mappability: [mean MQ20 cov for window is < 0.5 x mean coverage for the window]
./cornetto funbits -H 2.5 -L 0.5 -Q 0.5 test/cov-total.bg -q test/cov-mq20.bg | awk '{if ($4!=".") print $1"\t"$2"\t"$3}' > 1.bed

#2# merge these interesting windows - overlapping or adjacent (within 1000bp)
bedtools sort -i 1.bed | bedtools merge -d 1000 > 2.bed

#3# remove any merged intervals from (2) that are shorter than <25kb
awk '($3-$2)>=25000' 2.bed > 3.bed

#4# get regions longer than 10kb which were labelled by hifiasm as "low quality" (not sure what the definition of this is exactly)
cat RGBX240039_HG002.hifiasm.primary_asm.bp.p_ctg.lowQ.bed | awk '($3-$2)>=10000' | awk '{print $1"\t"$2"\t"$3}' > lowQ-10kbplus.bed

#5# combine the funbits from (3) and (4) and extend them by +50kb in either direction
cat 3.bed lowQ-10kbplus.bed | bedtools sort | awk '{if($2>50000){print $1"\t"$2-50000"\t"$3+50000} else {print $0}}' > funbits.bed

#6# add new windows to the file from (5), which are 200kb intervals at the left edge and right edge of the contig
awk '{if(($3-$2)>200000) {print $1"\t0\t200000\n"$1"\t"$3-200000"\t"$3}}' assembly.bed >> funbits.bed

#7# merge any overlapping or adjacent (within 50kb) intervals from (6)
bedtools sort -i funbits.bed | bedtools merge -d 50000 > funbits_merged.bed

#8# subtract merged windows from (7) from the whole genome assembly
bedtools subtract -a assembly.bed -b funbits_merged.bed > boringbits_tmp.bed

#9# subtract any contigs shorter than 1Mbase
awk '{if(($3-$2)<1000000) print $0}' assembly.bed > short.bed
bedtools subtract -a boringbits_tmp.bed -b short.bed > boringbits.bed

#10# if 'boring bits' are <50% of a single contig/scaffold, remove all boring bits on the whole scaffold. Use the below horrible inefficient code snippet for now
## i.e. if the contig is more than 50% interesting, capture the whole thing
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
done < boring_ctg.tmp > boringbits_cornetto4.bed
```

Yu may use the script under `scripts/create_cornetto4.sh`.


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
/install/vgp-pipeline/telomere/telomere_analysis.sh RGBX240039_HG002 0.4 50000 RGBX240039_HG002.hifiasm.primary_asm.fasta
cat telomere/RGBX240039_HG002.hifiasm.primary_asm.windows.0.4.50kb.ends.bed | cut -f 1  | sort | uniq -c | awk 'BEGIN{t1=0;t2=0;t3=0}{if($1==1){t1+=1}else if($1==2){t2+=1} else {t3+=1}} END{print "telo in one end:\t"t1"\ntelo in two ends:\t"t2"\ntelo more than 2 (must be 0):\t"t3"\n"}'
```

You may use the `scripts/telostats.sh`.

### Per chr stats

Launch `scripts/asmstats.sh` e.g. RGBX240039_HG002.hifiasm.primary_asm.fasta.

### PAF to depth

```
cat HG002_asm.hifiasm-cornetto5-2.fasta.fix.tmp.paf | awk '{print $6"\t"$8"\t"$9"\t"$1"\t"$12"\t"$5}' > HG002_asm.hifiasm-cornetto5-2.fasta.fix.tmp.paf.bed
```

### t2t-aware fishembly

Ref based:
```
scp gadi:/g/data/ox63/cornetto/HG002/evaluation/mapping_to_reference/PGXHXX240192_0.5.duplex_reads_out/PGXHXX240192_0.5.duplex_reads.hg002v1.0.1_pat.bam .
samtools index PGXHXX240192_0.5.duplex_reads.hg002v1.0.1_pat.bam
scp gadi:/g/data/ox63/cornetto/data/gtg_internal/HG002/PGXHXX240192_0.5.duplex_reads.fastq .
samtools fqidx PGXHXX240192_0.5.duplex_reads.fastq

samtools view PGXHXX240192_0.5.duplex_reads.hg002v1.0.1_pat.bam chr3_PATERNAL | cut -f 1 > t2t_rid.txt
samtools view PGXHXX240192_0.5.duplex_reads.hg002v1.0.1_pat.bam chr11_PATERNAL | cut -f 1 >> t2t_rid.txt
sort -u t2t_rid.txt > t2t_rid_uniq.txt
cut -f 1 PGXHXX240192_0.5.duplex_reads.fastq.fai > all_reads.txt
grep -v -F -f  t2t_rid_uniq.txt all_reads.txt > good.txt
samtools fqidx -r good.txt PGXHXX240192_0.5.duplex_reads.fastq  > PGXHXX240192_0.5.duplex_reads_good2.fastq
```

Ass based:
```
scp gadi:/g/data/ox63/cornetto/data/gtg_internal/HG002/PGXHXX240192_0.5.duplex_reads.fastq .
samtools fqidx PGXHXX240192_0.5.duplex_reads.fastq
scp gadi:/g/data/ox63/cornetto/HG002/cornetto_assemblies/HG002_asm.hifiasm-cornetto4-5/HG002_asm.hifiasm-cornetto4-5.fasta .
minimap2 -ax map-ont --secondary=no -t20 HG002_asm.hifiasm-cornetto4-5.fasta PGXHXX240192_0.5.duplex_reads.fastq | samtools sort - -o PGXHXX240192_0.5.duplex_reads.cornetto4-5.bam
cat HG002_asm.hifiasm-cornetto4-5/HG002_asm.hifiasm-cornetto4-5.windows.0.4.50kb.ends.bed  | cut -f 1  | sort | uniq -c | sort -nr -k1,1

samtools index PGXHXX240192_0.5.duplex_reads.cornetto4-5.bam
samtools view PGXHXX240192_0.5.duplex_reads.cornetto4-5.bam ptg000019l | cut -f 1 > t2t_rid.txt
samtools view PGXHXX240192_0.5.duplex_reads.cornetto4-5.bam ptg000002l | cut -f 1 >> t2t_rid.txt
sort -u t2t_rid.txt > t2t_rid_uniq.txt
cut -f 1 PGXHXX240192_0.5.duplex_reads.fastq.fai > all_reads.txt
grep -v -F -f  t2t_rid_uniq.txt all_reads.txt > good.txt
samtools fqidx -r good.txt PGXHXX240192_0.5.duplex_reads.fastq  > PGXHXX240192_0.5.duplex_reads_good.fastq
```

Ass based remove all telo reads:
```
HiFi FASTQ (A_0):
/g/data/ox63/cornetto/data/gtg_internal/HG002/RGBX240039_HG002.hifi.fastq.gz

HiFi base assembly (A_0):
/g/data/ox63/cornetto/HG002/base_assemblies/hifi-hifiasm/RGBX240039_HG002.hifiasm-hifi_1x/A_0-RGBX240039_HG002.fasta

Cornetto Duplex reads (A_1):
/g/data/ox63/cornetto/data/gtg_internal/HG002/A_1-QGXHXX240262.duplex_reads.fastq.gz

Cornetto duplex assembly (created without removing any telomere reads) (A_1):
/directflow/KCCGGenometechTemp/projects/iradev/operation_cornetto/A_1/hg002-cornetto-A_1.fasta

1. Run telomere finding script on HiFi assembly (A_0) to find contigs with a telomere at one corner (or both corners).
~/hasindu2008.git/cornetto/scripts/telostas.sh A_0-RGBX240039_HG002.fasta

2. If a contig has a corner with a telomere, create a 100kb window covering that corner.
cat A_0-RGBX240039_HG002/A_0-RGBX240039_HG002.windows.0.4.50kb.ends.bed | awk '{if($2==0) { print $1":0-100000" } else if ($3-100000>100000 ) { print $1":"$3-100000"-"$3 } }' > telocorners.txt

3. Align duplex FASTQ file (A_1) to HiFi assembly (A_0) and identify any read that is a unique (MapQ>=5), primary alignment within a telomere corner window; save IDs.
minimap2 -ax map-ont --secondary=no -t20 A_0-RGBX240039_HG002.fasta A_1-QGXHXX240262.duplex_reads.fastq.gz | samtools sort - -o A_1-QGXHXX240262.duplex_reads.A_0.bam && samtools index A_1-QGXHXX240262.duplex_reads.A_0.bam 
samtools view A_1-QGXHXX240262.duplex_reads.A_0.bam -r telocorners.txt | cut -f 1 >  t2t_rid.txt
sort -u t2t_rid.txt > t2t_rid_uniq.txt

4. Remove reads labelled in (3) from the Duplex FASTQ file (A_1) to create a new FASTQ file: A_1-clean.
cut -f 1 PGXHXX240192_0.5.duplex_reads.fastq.fai > all_reads.txt
grep -v -F -f  t2t_rid_uniq.txt all_reads.txt > good.txt
gunzip A_1-QGXHXX240262.duplex_reads.fastq.gz
samtools fqidx A_1-QGXHXX240262.duplex_reads.fastq
samtools fqidx -r good.txt A_1-QGXHXX240262.duplex_reads.fastq  > A_1-QGXHXX240262_clean.duplex_reads.fastq
```

