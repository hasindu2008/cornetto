#!/bin/bash

die() {
    echo "$1"
    exit 1
}

## On fridge: run inside a screen session

#PREFIX=D_1
#SAMPLE=QGXHXX240275

BASE_FASTQ=
FISH_PREV=
OUT_PREFIX=

usage(){
    echo "Usage: $0 [-b /path/to/D_0.fastq -p D_1_QGXHXX240275:D_2_QGXHXX240279 -o hg002-cornetto-D_3] <D_3> <QGXHXX240283>"
    echo "Options:"
    echo "  -b  Base fastq file"
    echo "  -p  Previous fish data"
    echo "  -o  Output prefix"
}

while getopts "o:b:p:h" option; do
   case $option in
        b) BASE_FASTQ=$OPTARG;;
        p) FISH_PREV=$OPTARG;;
        o) OUT_PREFIX=$OPTARG;;
        h) usage; exit 0;;
        \?) echo "Invalid option: $OPTARG" >&2; exit 1;;
   esac
done
shift $(($OPTIND - 1))

test $# -eq 2 || { usage; exit 1; }
PREFIX=$1
SAMPLE=$2

FRIDGE_TMP=/data3/cornetto
GTA100_SCRIPT=/home/hasindu/hasindu2008.git/cornetto/shitflow/simplex/basecall-gta100.sh
GTA100_DATA=/data/hasindu/shitflow/
GADI_DATA=/g/data/ox63/hasindu/cornetto/shitflow
GADI_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/shitflow/hifiasm-ont.pbs.sh

NAME=${PREFIX}_${SAMPLE}
GADI_PBS_ARGS=
echo "Launching shitflow for ${NAME}"

checkshit() {
    test -d ${FRIDGE_TMP} || die "${FRIDGE_TMP} not found on fridge"
    ssh gta100 "test -x ${GTA100_SCRIPT}" || die "${GTA100_SCRIPT} not found on gta100"
    ssh gta100 "test -d ${GTA100_DATA}" || die "${GTA100_DATA} not found on gta100"
    ssh gta100 "test -d ${GTA100_DATA}/${NAME}" && die "${GTA100_DATA}/${NAME} already exists on gta100. Delete that shit first"
    ssh gadi "test -d ${GADI_DATA}" || die "${GADI_DATA} not found on gadi"
    ssh gadi "test -d ${GADI_DATA}/${NAME}" && die "${GADI_DATA}/${NAME} already exists on gadi. Delete that shit first"
    ssh gadi "test -x ${GADI_SCRIPT}" || die "${GADI_SCRIPT} not found on gadi"

    if [ -n "${BASE_FASTQ}" ]
    then
        ssh gadi "test -e ${BASE_FASTQ}" || die "${BASE_FASTQ} not found on gadi"
        GADI_PBS_ARGS="BASE_FASTQ="${BASE_FASTQ}
    fi
    if [ -n "${FISH_PREV}" ]; then
        LIST=$(echo "$FISH_PREV" | tr ':' ' ' )
        for SIMPLEX in ${LIST}; do
            DUP=${GADI_DATA}/${SIMPLEX}/${SIMPLEX}.fastq
            ssh gadi "test -e ${DUP}" || die "${DUP} not found on gadi"
        done
        GADI_PBS_ARGS=${GADI_PBS_ARGS}",FISH_PREV="${FISH_PREV}
    fi
}

checkshit

if [ -n "${GADI_PBS_ARGS}" ]; then
    GADI_PBS_ARGS=${GADI_PBS_ARGS}",FISH_NOW="${NAME}
    if [ -n "${OUT_PREFIX}" ]; then
        GADI_PBS_ARGS=${GADI_PBS_ARGS}",OUT_PREFIX="${OUT_PREFIX}
    fi
fi

cd ${FRIDGE_TMP} || die "Could not cd to ${FRIDGE_TMP}"
slow5tools merge /data/${SAMPLE}/*/*/slow5/ -o ${PREFIX}_${SAMPLE}.blow5 || die "Could not merge slow5 files"
slow5tools stats ${PREFIX}_${SAMPLE}.blow5 || die "Could not get stats"

ssh gta100 "mkdir ${GTA100_DATA}/${NAME}" || die "gta100 ssh failed"
scp ${NAME}.blow5 gta100:${GTA100_DATA}/${NAME}/ || die "copying to Gadi failed"

COMMAND="source /etc/profile; screen -S shitflow_${PREFIX}_${SAMPLE} -d -m -L ${GTA100_SCRIPT} ${NAME} ${GADI_PBS_ARGS}"
echo "$COMMAND"
ssh gta100 "$COMMAND"
echo ""
echo "Handed the work to the gta100"
