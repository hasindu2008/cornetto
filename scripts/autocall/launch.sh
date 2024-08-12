#!/bin/bash

die() {
    echo "$1"
    exit 1
}

## On fridge: run inside a screen session

#PREFIX=A_1
#SAMPLE=QGXHXX240275

BASE_FASTQ=
FISH_PREV=
OUT_PREFIX=

usage(){
    echo "Usage: $0 [-b /path/to/RGBX240039_HG002.hifi.fastq.gz -p A_1_QGXHXX240275:A_2_QGXHXX240279 -o hg002-cornetto-A_3] <A_3> <QGXHXX240283>"
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
BRENNER_SCRIPT=/home/hasgam/hasindu2008.git/cornetto/scripts/autocall/autocall.sh
BRENNER_DATA=/directflow/KCCGGenometechTemp/projects/iradev/operation_cornetto/autocall_hasindu/
GADI_DATA=/g/data/ox63/hasindu/cornetto/autocall
GADI_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/scripts/hifiasm.pbs.sh

NAME=${PREFIX}_${SAMPLE}
GADI_PBS_ARGS=
echo "Launching autocall for ${NAME}"

checkshit() {
    test -d ${FRIDGE_TMP} || die "${FRIDGE_TMP} not found on fridge"
    ssh brenner-fpga "test -x ${BRENNER_SCRIPT}" || die "${BRENNER_SCRIPT} not found on brenner-fpga"
    ssh brenner-fpga "test -d ${BRENNER_DATA}" || die "${BRENNER_DATA} not found on brenner-fpga"
    ssh brenner-fpga "test -d ${BRENNER_DATA}/${NAME}" && die "${BRENNER_DATA}/${NAME} already exists on brenner-fpga. Delete that shit first"
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
        for DUPLEX in ${LIST}; do
            DUP=${GADI_DATA}/${DUPLEX}/${DUPLEX}.duplex_reads.fastq
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

COMMAND="source /etc/profile; screen -S autocall_${PREFIX}_${SAMPLE} -d -m -L ${BRENNER_SCRIPT} ${NAME} ${GADI_PBS_ARGS}"
echo "$COMMAND"
ssh brenner-fpga "$COMMAND"
echo ""
echo "Handed the work to the brenner-fpga"
