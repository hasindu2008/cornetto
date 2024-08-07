#!/bin/bash

die() {
    echo "$1"
    exit 1
}

[ $# -ne 1 ] && die "Usage: $0 <prefix_sample>"

NAME=$1
BRENNER_DATA=/directflow/KCCGGenometechTemp/projects/iradev/operation_cornetto/autocall_hasindu/
FRIDGE_TMP=/data3/cornetto/
GADI_DATA=/g/data/ox63/hasindu/cornetto/autocall

SCRIPT_REALPATH=$(realpath "$0")
SCRIPT_PATH=$(dirname "$SCRIPT_REALPATH")

test -e $SCRIPT_PATH/fridge_dload.sge.sh || die "Could not find $SCRIPT_PATH/fridge_dload.sge.sh"
test -e $SCRIPT_PATH/split_dorado_duplex.sge.sh || die "Could not find $SCRIPT_PATH/split_dorado_duplex.sge.sh"
test -e $SCRIPT_PATH/get_duplex_and_simplex_reads.sge.sh || die "Could not find $SCRIPT_PATH/get_duplex_and_simplex_reads.sge.sh"

cd $BRENNER_DATA/ || die "Could not cd to $BRENNER_DATA"
mkdir ${NAME} || die "Could not create ${NAME}"
cd ${NAME} || die "Could not cd to ${NAME}"

qsub -sync y $SCRIPT_PATH/fridge_dload.sge.sh ${FRIDGE_TMP}/${NAME}.blow5 ./ || die "Could not download the data"
qsub -sync y $SCRIPT_PATH/split_dorado_duplex.sge.sh ${NAME} || die "Could not basecall the data"

#handle those failed basecall files if any?
# these qsub scripts are not propoerly error checked

qsub -sync y $SCRIPT_PATH/get_duplex_and_simplex_reads.sge.sh ${NAME} || die "Could not basecall the data"

ssh gadi "mkdir ${GADI_DATA}/${NAME}" || die "gadi ssh failed"
scp ${NAME}.duplex_reads.fastq gadi-dm:${GADI_DATA}/${NAME}/ || die "copying to Gadi failed"

echo "Done"
