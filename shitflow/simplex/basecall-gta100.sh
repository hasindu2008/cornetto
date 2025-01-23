#!/bin/bash

die() {
    echo "$1" 1>&2
    exit 1
}

usage () {
    echo "Usage: $0 <prefix_sample> [<GADI_PBS_ARGS>(optional)]"
    exit 1
}

if [ $# -lt 1 ] ; then
    usage
elif [ $# -gt 2 ]; then
    usage
fi

NAME=$1
if [ $# -eq 2 ]; then
    GADI_PBS_ARGS=$2
fi

GTA100_DATA=/data/hasindu/shitflow/
GADI_DATA=/g/data/ox63/hasindu/cornetto/shitflow
GADI_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/shitflow/hifiasm-ont.pbs.sh

SCRIPT_REALPATH=$(realpath "$0")
SCRIPT_PATH=$(dirname "$SCRIPT_REALPATH")

test -d $GTA100_DATA || die "$GTA100_DATA not found"
test -d $GTA100_DATA/${NAME} || die "$GTA100_DATA/${NAME} not found"
cd $GTA100_DATA/${NAME} || die "cd failed"

DEVICES=

GET_FREE_GPU () {

    nvidia-smi || die "nvidia-smi command failed"
    DEVICES_ALL=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd "," -)
    [ -z "$DEVICES_ALL" ] && die "No GPU found using nvidia-smi --query-gpu=index --format=csv,noheader"
    echo "All GPUs: $DEVICES_ALL"

    while [ -z "$DEVICES" ]; do
        # this needs some beautification and catch errors
        DEVICES=$(for i in {1..10}; do nvidia-smi --query-gpu=index,utilization.gpu,utilization.memory --format=csv,noheader,nounits; done | awk -F',' '{gpu[$1]+=$2; mem[$1]+=$3} END{for(i=0; i<length(gpu); i++){if(gpu[i]==0 && mem[i]==0) print i} }' | paste -sd ',')
        if [ -z "$DEVICES" ]; then
            echo "No free GPU. Waiting for 5 minutes"
            sleep 300
        fi
    done
    echo "Using GPUs: $DEVICES"
}

GET_FREE_GPU

export CUDA_VISIBLE_DEVICES=${DEVICES}

/usr/bin/time -v /install/slow5-dorado-0.8.3/bin/slow5-dorado basecaller -x cuda:all /install/dorado-0.8.3/models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0/ ${NAME}.blow5 --emit-fastq --min-qscore 10  > ${NAME}.basecalls.fastq || die "Failed to basecall"

/install/seqkit-v2.3.0/seqkit stats --threads 16 -o ${NAME}.seqkit_stats.tsv --all --tabular ${NAME}.basecalls.fastq || die "Failed to generate seqkit stats"
/install/seqkit-v2.3.0/seqkit seq -m 30000 ${NAME}.basecalls.fastq -o ${NAME}.fastq || die "Failed to filter reads"

ssh gadi "mkdir ${GADI_DATA}/${NAME}" || die "gadi ssh failed"
scp ${NAME}.fastq gadi-dm:${GADI_DATA}/${NAME}/ || die "copying to Gadi failed"

if [ -n "${GADI_PBS_ARGS}" ]; then
    GADI_COMMAND="cd ${GADI_DATA}/${NAME}/ && qsub -v ${GADI_PBS_ARGS} ${GADI_SCRIPT}"
    echo "Running on gadi: ${GADI_COMMAND}"
    ssh gadi "${GADI_COMMAND}" || die "gadi qsub failed"
    echo "Handed over work to gadi"
else
    echo "No GADI_PBS_ARGS provided. Not running on gadi"
fi






