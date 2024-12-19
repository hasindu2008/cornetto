## Cornetto 2 Pannel

1. print all interesting windows with :
   - low coverage: [<0.5x] genome average
   - high coverage: [>2.2x] genome average
   - low mappability: [mean MQ20 cov for window is < 0.5 x mean coverage for the window]

```
./cornetto noboringbits -H 2.2 -L 0.5 -Q 0.5 test/cov-total.bg -q test/cov-mq20.bg | awk '{if ($4!=".") print $1"\t"$2"\t"$3}' > 1.bed
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
cat A_0-RGBX240039_HG002/A_0-RGBX240039_HG002.windows.0.4.50kb.ends.bed | awk '{if($2==0) { print $1"\t0\t100000" } else if ($3-100000>100000 ) { print $1"\t"$3-100000"\t"$3 } }' > telocorners.txt

3. Align duplex FASTQ file (A_1) to HiFi assembly (A_0) and identify any read that is a unique (MapQ>=5), primary alignment within a telomere corner window; save IDs.
minimap2 -ax map-ont --secondary=no -t20 A_0-RGBX240039_HG002.fasta A_1-QGXHXX240262.duplex_reads.fastq.gz | samtools sort - -o A_1-QGXHXX240262.duplex_reads.A_0.bam && samtools index A_1-QGXHXX240262.duplex_reads.A_0.bam
samtools view A_1-QGXHXX240262.duplex_reads.A_0.bam -L telocorners.txt | cut -f 1 >  t2t_rid.txt
sort -u t2t_rid.txt > t2t_rid_uniq.txt

4. Remove reads labelled in (3) from the Duplex FASTQ file (A_1) to create a new FASTQ file: A_1-clean.

gunzip A_1-QGXHXX240262.duplex_reads.fastq.gz
samtools fqidx A_1-QGXHXX240262.duplex_reads.fastq
cut -f 1 A_1-QGXHXX240262.duplex_reads.fastq.fai > all_reads.txt
grep -v -F -f  t2t_rid_uniq.txt all_reads.txt > good.txt
samtools fqidx -r good.txt A_1-QGXHXX240262.duplex_reads.fastq  > A_1-QGXHXX240262_clean.duplex_reads.fastq
```