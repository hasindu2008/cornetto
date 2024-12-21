## Method to get eventual boring bits

Create a full genome assembly bed:
```
samtools faidx assembly.fa
awk '{print $1"\t0\t"$2}' assembly.fa.faidx > assembly.bed
```

### cornetto panel

```
#1# print all interesting windows with:
# low coverage: [<0.5x] genome average
# high coverage: [>2.5x] genome average
# low mappability: [mean MQ20 cov for window is < 0.5 x mean coverage for the window]
./cornetto noboringbits -H 2.5 -L 0.5 -Q 0.5 test/cov-total.bg -q test/cov-mq20.bg | awk '{if ($4!=".") print $1"\t"$2"\t"$3}' > 1.bed

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

You may use the script under `scripts/create-cornetto.sh`.

For previous deprecated cornetto panels refer to [archived](archived.md).