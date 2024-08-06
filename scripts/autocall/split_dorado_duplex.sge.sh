#!/bin/bash
#$ -cwd
#$ -N SLOW5DORADO
#$ -S /bin/bash
#$ -b y
#$ -l mem_requested=8G
#$ -l h_vmem=8G
#$ -l nvgpu=4
#$ -l tmp_requested=150G
#$ -pe smp 32
#$ -V

die() {
    echo "$@" >&2
    exit 1
}

SAMPLE=${1}
BLOW5=/directflow/KCCGGenometechTemp/projects/iradev/operation_cornetto/autocall_hasindu/${SAMPLE}/${SAMPLE}.blow5

#don't edit stuff below here
OUTDIR="split_bam/"
SLOW5TOOLS=/share/ClusterShare/software/contrib/hasgam/slow5tools-v1.1.0/slow5tools

SLOW5_DORADO=/share/ClusterShare/software/contrib/hasgam/slow5-dorado-0.3.4/bin/slow5-dorado
MODEL_SUP=/share/ClusterShare/software/contrib/hasgam/slow5-dorado-0.3.4/models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0

mkdir ${SAMPLE}_duplex_out
cd ${SAMPLE}_duplex_out

mkdir split_lists || die "Could not create split directory"
mkdir $TMPDIR/split_blow5 || die "Could not create split blow5 directory"
mkdir split_bam || die "Could not create split bam directory"
mkdir split_blow5_failed || die "Could not create split blow5 directory"

BLOW5_LOCAL=$TMPDIR/reads.blow5

echo "Copying SLOW5 file"
cp $BLOW5 $BLOW5_LOCAL || die "Copying $BLOW5 to $BLOW5_LOCAL  failed"

echo "generating split lists"
${SLOW5TOOLS} skim -t40 $BLOW5_LOCAL | awk -v c="channel_number" 'NR==1{for (i=1; i<=NF; i++) if ($i==c){p=i; break};} {print $1"\t"$p}' | tail -n+2 | awk '{print $1 > "split_lists/"int($2/50)".txt"}'

NUM1=$(cat split_lists/*.txt | wc -l)
NUM2=$(${SLOW5TOOLS} skim --rid ${BLOW5_LOCAL} | wc -l)

if [ "$NUM1" -ne "$NUM2" ]; then
    echo "ERROR: Number of reads in split files does not match number of reads in original file"
    exit 1
fi

echo "going through each channel"
for i in split_lists/*.txt; do
    channel=$(basename $i .txt)
    echo "Extracting channel group: $channel"
    ${SLOW5TOOLS} get -t40 ${BLOW5_LOCAL} --list $i  -o $TMPDIR/split_blow5/${channel}.blow5 || die "Could not extract channel $channel"
    echo "basecalling channel group: $channel"
    /usr/bin/time -v ${SLOW5_DORADO} duplex ${MODEL_SUP} $TMPDIR/split_blow5/${channel}.blow5  > split_bam/sup_duplex_unmapped_${channel}.bam
    status=$?
    if [ $status -eq 0 ]
    then
        echo "basecalling succeeded for channel group $channel"
    else
        echo "ERROR: basecallings failed for channel group $channel"
        cp $TMPDIR/split_blow5/${channel}.blow5 split_blow5_failed/
        rm split_bam/sup_duplex_unmapped_${channel}.bam
    fi
done
echo "all done"


