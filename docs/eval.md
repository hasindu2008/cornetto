## Evaluation methods

A useful article https://lh3.github.io/2021/04/17/concepts-in-phased-assemblies.

On gadi you can use scripts/getstat.pbs.sh to get dotplots, telomere and per-chr stats. Use `scripts/compleasm.pbs.sh` for buscos.

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

Need trio fastq to create the yaks: https://hifiasm.readthedocs.io/en/latest/trio-assembly.html.
You can download the yak for trio HG002 as suggested in https://github.com/chhylp123/hifiasm/issues/262, at https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=submissions/6040D518-FE32-4CEB-B55C-504A05E4D662--HG002_PARENTAL_YAKS/HG002_PARENTS_FULL/.

```
yak trioeval -t16 pat.yak mat.yak assembly.fa
```

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
