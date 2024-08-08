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
BLOW5_DIR=/directflow/KCCGGenometechTemp/projects/iradev/operation_cornetto/autocall_hasindu/${SAMPLE}/${SAMPLE}_duplex_out/split_blow5_failed

SLOW5_DORADO=/share/ClusterShare/software/contrib/hasgam/slow5-dorado-0.3.4/bin/slow5-dorado
MODEL_SUP=/share/ClusterShare/software/contrib/hasgam/slow5-dorado-0.3.4/models/dna_r10.4.1_e8.2_400bps_sup@v4.2.0


ls -1 -A ${BLOW5_DIR} | while read file
do
    echo "Reprocessing failed channel: $file"
    channel=$(echo $file | cut -d'.' -f1)
    BLOW5=${BLOW5_DIR}/${file}
    BLOW5_LOCAL=$TMPDIR/reads.blow5
    echo "Copying SLOW5 file $BLOW5"
    cp $BLOW5 $BLOW5_LOCAL || die "Copying $BLOW5 to $BLOW5_LOCAL  failed"

    /usr/bin/time -v ${SLOW5_DORADO} duplex ${MODEL_SUP} ${BLOW5_LOCAL}  > sup_duplex_unmapped_${channel}.bam
    status=$?
    if [ $status -eq 0 ]
    then
        echo "basecalling succeeded for channel group $channel"
        mv sup_duplex_unmapped_${channel}.bam ${BLOW5_DIR}/../split_bam/
    else
        echo "ERROR: basecallings failed for channel group $channel"
        rm sup_duplex_unmapped_${channel}.bam
    fi

    rm $TMPDIR/reads.blow5

done






