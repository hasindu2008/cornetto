#!/bin/bash

die() {
    echo "$1"
    exit 1
}

#DEBUG
# cd /home/hasindu/scratch/cornetto/ass/1x_hifi
# ./asmstats.sh RGBX240039_HG002.hifiasm.primary_asm.fasta

[ $# -ne 1 ] && die "Usage: $0 <FASTA>"

datamash --version > /dev/null || die "Datamash not found"

FASTA=$1
PREFIX=$(basename $FASTA .fasta)
FILE=$PREFIX.fasta.fix.tmp.paf
test -e $FILE  || die "File $FILE does not exist. Did you run minidotplot.sh?"
test -e ${PREFIX}.report.tsv || die "File ${PREFIX}.report.tsv does not exist. Did you run minidotplot.sh?"
awk '{print $1"\t"$4}' ${PREFIX}.report.tsv > ${PREFIX}.chr.rename.txt || die "awk failed"
test -e ${PREFIX}/${PREFIX}.windows.0.4.50kb.ends.bed || die "File ${PREFIX}/${PREFIX}.windows.0.4.50kb.ends.bed does not exist. Did you run telostats.sh?"

awk '{print "s/"$1"/"$2"/g"}' ${PREFIX}.fasta.chr.rename.txt | sed -f - ${PREFIX}/${PREFIX}.windows.0.4.50kb.ends.bed | sort -k1,1 | cut -f1 | uniq -c | sort -k1,1 -r -n | awk '{print $2"\t"$1}' > ${PREFIX}/$PREFIX.fasta.fix.tmp.telostat.txt

echo "${FASTA}"
echo ""
echo -e "chr\tT2T?\tNTelo\tTelocontiglen"
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do
T2T=$(grep $chr"_" ${PREFIX}/$PREFIX.fasta.fix.tmp.telostat.txt | awk '{if($2==2){print "y"} else {print "n"}}' |  tr '\n' ',')
NTelo=$(grep $chr"_" ${PREFIX}/$PREFIX.fasta.fix.tmp.telostat.txt | datamash sum 2)
cat ${PREFIX}/$PREFIX.fasta.fix.tmp.telostat.txt  | grep $chr"_" | cut -f 1 > ${PREFIX}/$PREFIX.fasta.fix.tmp.telostat.$chr.txt || die "Error in extracting telomere contigs"
TELLEN=$(cat $FILE | cut -f 1,2 | grep -w -F -f ${PREFIX}/$PREFIX.fasta.fix.tmp.telostat.$chr.txt | uniq -c | cut -f 2 | tr '\n' ',')
echo -e $chr"\t"$T2T"\t"$NTelo"\t"$TELLEN
done

echo ""
echo ""
echo "Contigs whose majority is mapped to the corresponding chromosome"
echo -e "\tNcontigsofsize>=KMbasealignedtochr\t\t\t\t\t%ofchrsequencecoveredbycontigsofsize>=KMbase"
echo -e "chr\t0Mbase\t0.1Mbase\t1Mbase\t5Mbase\t10Mbase\t0Mbase\t0.1Mbase\t1Mbase\t5Mbase\t10Mbase"
for chr in chr1_PATERNAL chr2_PATERNAL chr3_PATERNAL chr4_PATERNAL chr5_PATERNAL chr6_PATERNAL chr7_PATERNAL chr8_PATERNAL chr9_PATERNAL chr10_PATERNAL chr11_PATERNAL chr12_PATERNAL chr13_PATERNAL chr14_PATERNAL chr15_PATERNAL chr16_PATERNAL chr17_PATERNAL chr18_PATERNAL chr19_PATERNAL chr20_PATERNAL chr21_PATERNAL chr22_PATERNAL chrX_MATERNAL chrY_PATERNAL
do
chr_prefix=$(echo $chr | awk -F'_' '{print $1}')
LEN=$(awk -v chr="$chr" '{if($6==chr) print $1"\t"$7}' $FILE | head -1 | cut -f 2 ) # length of the chromosome

# 1. get the contigs aligned to the chromosome with the alignment block length on the chr
# 2. sort the contigs by contig name
# 3. get the sum of the alignment block length for each contig
# 4. filter the contigs that are >50% of its length are actually aligned to the chromosome
# 5. get the counts and stats
awk -v chr="$chr" '{if($6==chr) print $1"\t"$9-$8}' $FILE | sort -k1,1 | datamash -g 1 sum 2 | grep "$chr_prefix""_" | awk -v chr=$chr_prefix -v len="$LEN" 'BEGIN{c0=0;s0=0;c01=0;s01=0;c1=0;s1=0;c5=0;s5=0;c10=0;c10=0;}{if($2>0){c0+=1; s0+=$2} if($2>=100000){c01+=1; s01+=$2} if($2>=1000000){c1+=1; s1+=$2} if($2>=5000000){c5+=1; s5+=$2} if($2>=10000000){c10+=1; s10+=$2}  }END{print chr"\t"c0"\t"c01"\t"c1"\t"c5"\t"c10"\t"s0/len*100"\t"s01/len*100"\t"s1/len*100"\t"s5/len*100"\t"s10/len*100}'
done


# L50: Min N major-mapping contigs that cover >= 50% of chromosome [n=X]
# L90: Min N major-mapping contigs that cover >= 90% of chromosome [n=X]
# L95: Min N major-mapping contigs that cover >= 95% of chromosome [n=X]
# L99: Min N major-mapping contigs that cover >= 99% of chromosome [n=X]
# CumCovN5: Cumulative % covered by n=5 largest contigs [n%,n%,n%,n%,n%]
echo ""
echo ""
echo "LX of Contigs whose majority is mapped to the corresponding chromosome"
echo -e "\tL50\tL90\tL95\tL99\tCumCovN5"
for chr in chr1_PATERNAL chr2_PATERNAL chr3_PATERNAL chr4_PATERNAL chr5_PATERNAL chr6_PATERNAL chr7_PATERNAL chr8_PATERNAL chr9_PATERNAL chr10_PATERNAL chr11_PATERNAL chr12_PATERNAL chr13_PATERNAL chr14_PATERNAL chr15_PATERNAL chr16_PATERNAL chr17_PATERNAL chr18_PATERNAL chr19_PATERNAL chr20_PATERNAL chr21_PATERNAL chr22_PATERNAL chrX_MATERNAL chrY_PATERNAL
do
chr_prefix=$(echo $chr | awk -F'_' '{print $1}')
LEN=$(awk -v chr="$chr" '{if($6==chr) print $1"\t"$7}' $FILE | head -1 | cut -f 2 ) # length of the chromosome

# 1. get the contigs aligned to the chromosome with the alignment block length on the chr
# 2. sort the contigs by contig name
# 3. get the sum of the alignment block length for each contig
# 4. filter the contigs that are >50% of its length are actually aligned to the chromosome
# 5. get the counts and stats
awk -v chr="$chr" '{if($6==chr) print $1"\t"$9-$8}' $FILE | sort -k1,1 | datamash -g 1 sum 2 | grep "$chr_prefix""_" | cut -f 2 | sort -r -n | awk -v chr=$chr_prefix -v len="$LEN" 'BEGIN{l50=0;l90=0;l95=0;l99=0;sum=0;CumCov1=0;CumCov2=0;CumCov3=0;CumCov4=0;CumCov5=0;} { sum+=$1; if(sum>=len*0.50 && l50==0){l50=NR} if(sum>=len*0.90 && l90==0){l90=NR} if(sum>=len*0.95 && l95==0){l95=NR}  if(sum>=len*0.99 && l99==0){l99=NR} if(NR<=1) {CumCov1+=$1} if(NR<=2) {CumCov2+=$1} if(NR<=3) {CumCov3+=$1} if(NR<=4) {CumCov4+=$1} if(NR<=5) {CumCov5+=$1} } END{print chr"\t"l50"\t"l90"\t"l95"\t"l99"\t"CumCov1/len*100","CumCov2/len*100","CumCov3/len*100","CumCov4/len*100","CumCov5/len*100}'
done


echo ""
echo ""
echo "Contigs whose majority is mapped to another chromosome"
# now do the same with grep -v
echo -e "\tNcontigsofsize>=KMbasealignedtochr\t\t\t\t\t%ofchrsequencecoveredbycontigsofsize>=KMbase				"
echo -e "chr\t0Mbase\t0.1Mbase\t1Mbase\t5Mbase\t10Mbase\t0Mbase\t0.1Mbase\t1Mbase\t5Mbase\t10Mbase"
for chr in chr1_PATERNAL chr2_PATERNAL chr3_PATERNAL chr4_PATERNAL chr5_PATERNAL chr6_PATERNAL chr7_PATERNAL chr8_PATERNAL chr9_PATERNAL chr10_PATERNAL chr11_PATERNAL chr12_PATERNAL chr13_PATERNAL chr14_PATERNAL chr15_PATERNAL chr16_PATERNAL chr17_PATERNAL chr18_PATERNAL chr19_PATERNAL chr20_PATERNAL chr21_PATERNAL chr22_PATERNAL chrX_MATERNAL chrY_PATERNAL
do
chr_prefix=$(echo $chr | awk -F'_' '{print $1}')
LEN=$(awk -v chr="$chr" '{if($6==chr) print $1"\t"$7}' $FILE | head -1 | cut -f 2 ) # length of the chromosome

# 1. get the contigs aligned to the chromosome with the alignment block length on the chr
# 2. sort the contigs by contig name
# 3. get the sum of the alignment block length for each contig
# 4. filter the contigs that are >50% of its length are actually aligned to the chromosome
# 5. get the counts and stats
awk -v chr="$chr" '{if($6==chr) print $1"\t"$9-$8}' $FILE | sort -k1,1 | datamash -g 1 sum 2 | grep -v "$chr_prefix""_" | awk -v chr=$chr_prefix -v len="$LEN" 'BEGIN{c0=0;s0=0;c01=0;s01=0;c1=0;s1=0;c5=0;s5=0;c10=0;c10=0;}{if($2>0){c0+=1; s0+=$2} if($2>=100000){c01+=1; s01+=$2} if($2>=1000000){c1+=1; s1+=$2} if($2>=5000000){c5+=1; s5+=$2} if($2>=10000000){c10+=1; s10+=$2}  }END{print chr"\t"c0"\t"c01"\t"c1"\t"c5"\t"c10"\t"s0/len*100"\t"s01/len*100"\t"s1/len*100"\t"s5/len*100"\t"s10/len*100}'
done
