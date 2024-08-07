#!/bin/bash

die() {
    echo "$1"
    exit 1
}

## On fridge: run inside a screen session

#PREFIX=A_1
#SAMPLE=QGXHXX240275

test $# -eq 3 || die "Usage: $0 <prefix:A_1> <sample:QGXHXX240275> <pacbiofastq>"

PREFIX=$1
SAMPLE=$2
BASEFASTQ=$3

FRIDGE_TMP=/data3/cornetto
BRENNER_SCRIPT=/home/hasgam/hasindu2008.git/cornetto/scripts/autocall/autocall.sh
BRENNER_DATA=/directflow/KCCGGenometechTemp/projects/iradev/operation_cornetto/autocall_hasindu/
GADI_DATA=/g/data/ox63/hasindu/cornetto/autocall
GADI_BASEDATA=/g/data/ox63/cornetto/data/gtg_internal/HG002/

NAME=${PREFIX}_${SAMPLE}

checkshit() {
    test -d ${FRIDGE_TMP} || die "${FRIDGE_TMP} not found on fridge"
    ssh brenner-fpga "test -x ${BRENNER_SCRIPT}" || die "${BRENNER_SCRIPT} not found on brenner-fpga"
    ssh brenner-fpga "test -d ${BRENNER_DATA}" || die "${BRENNER_DATA} not found on brenner-fpga"
    ssh brenner-fpga "test -d ${BRENNER_DATA}/${NAME}" && die "${BRENNER_DATA}/${NAME} already exists on brenner-fpga. Delete that shit first"
    ssh gadi "test -d ${GADI_DATA}" || die "${GADI_DATA} not found on gadi"
    ssh gadi "test -d ${GADI_DATA}/${NAME}" && die "${GADI_DATA}/${NAME} already exists on gadi. Delete that shit first"
    ssh gadi "test -e ${GADI_BASEDATA}/BASEFASTQ" || die "${GADI_BASEDATA}/BASEFASTQ not found on gadi"
}

checkshit

cd ${FRIDGE_TMP} || die "Could not cd to ${FRIDGE_TMP}"
slow5tools merge /data/${SAMPLE}/*/*/slow5/ -o ${PREFIX}_${SAMPLE}.blow5 || die "Could not merge slow5 files"
slow5tools stats ${PREFIX}_${SAMPLE}.blow5 || die "Could not get stats"

COMMAND="source /etc/profile; screen -S autocall_${PREFIX}_${SAMPLE} -d -m -L ${BRENNER_SCRIPT} ${NAME}"
echo "$COMMAND"
ssh brenner-fpga "$COMMAND"

echo "Handed the work to the brenner-fpga"
