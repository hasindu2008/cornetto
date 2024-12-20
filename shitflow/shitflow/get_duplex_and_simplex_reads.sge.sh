#!/bin/bash
#$ -cwd
#$ -N merge
#$ -S /bin/bash
#$ -b y
#$ -l mem_requested=24G
#$ -l h_vmem=24G
#$ -l tmp_requested=24G
#$ -pe smp 8
#$ -l h='!epsilon*'
#$ -V

THREADS=8

SAMPLE=${1}

BAM_LIST=${SAMPLE}_duplex_out/split_bam.list

find ${SAMPLE}_duplex_out/split_bam/ -name *.bam > ${BAM_LIST}

export PATH=$PATH:/home/iradev/scripts/
export MODULEPATH=$MODULEPATH:/share/ClusterShare/Modules/modulefiles/noarch:/share/ClusterShare/Modules/modulefiles/centos6.2_x86_64:/share/ClusterShare/Modules/modulefiles/contrib/:/usr/share/Modules/modulefiles:/etc/modulefiles
module load centos7.8/aletan/samtools/1.14
module load centos6.10/gi/pigz/2.3
module load centos7.8/qiadu/minimap2/2.22
module load centos6.10/hasgam/slow5tools-0.8.0
module load centos6.10/kseskv/seqtk/1.3

DUPLEX_FQ=${SAMPLE}.duplex_reads.fastq
SIMPLEX_FQ=${SAMPLE}.simplex_reads.fastq
SIMPLEX_FILT=${SAMPLE}.simplex-min10kb.fastq

rm -f ${DUPLEX_FQ}
rm -f ${SIMPLEX_FQ}

while read SEG
do
## get duplex reads from the dorado output
samtools view -H ${SEG} > ${SAMPLE}_tmp.sam
samtools view -@ ${THREADS} ${SEG} | awk 'length($1)==73' >> ${SAMPLE}_tmp.sam
samtools fastq ${SAMPLE}_tmp.sam >> ${DUPLEX_FQ}

## get list of original individual readIDs that were successfully incorporated into duplexes
samtools view ${SAMPLE}_tmp.sam | cut -f 1 | tr ";" "\n" > ${SAMPLE}_duplex_indv_ids.txt

## get simplex reads from the dorado output
samtools view -@ ${THREADS} ${SEG} | awk 'length($1)==36' > ${SAMPLE}_tmp.sam
samtools view -H ${SEG} > ${SAMPLE}_tmp2.sam
removeSubset.pl ${SAMPLE}_tmp.sam ${SAMPLE}_duplex_indv_ids.txt 1 1 >> ${SAMPLE}_tmp2.sam
samtools fastq ${SAMPLE}_tmp2.sam >> ${SIMPLEX_FQ}
done < ${BAM_LIST}

## exclude simplex reads shorter than 10kb
seqtk seq -L 10000 ${SIMPLEX_FQ} > ${SIMPLEX_FILT}

