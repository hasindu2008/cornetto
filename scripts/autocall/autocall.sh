#!/bin/bash

die() {
    echo "$1"
    exit 1
}

[ $# -ne 1 ] && die "Usage: $0 <prefix_sample>"

NAME=$1
DATA_PATH=/directflow/KCCGGenometechTemp/projects/iradev/operation_cornetto/autocall_hasindu/
FRIDGE_PATH=/data3/cornetto/
SCRIPT_PATH=$(dirname "$0")


test -e $SCRIPT_PATH/fridge_dload.sge.sh || die "Could not find $SCRIPT_PATH/fridge_dload.sge.sh"
test -e $SCRIPT_PATH/split_dorado_duplex.sge.sh || die "Could not find $SCRIPT_PATH/split_dorado_duplex.sge.sh"

cd $DATA_PATH/ || die "Could not cd to $DATA_PATH"
mkdir ${NAME} || die "Could not create ${NAME}"
cd ${NAME} || die "Could not cd to ${NAME}"

qsub -sync y $SCRIPT_PATH/fridge_dload.sge.sh ${FRIDGE_PATH}/${NAME}.blow5 ./ || die "Could not download the data"
qsub -sync y $SCRIPT_PATH/split_dorado_duplex.sge.sh ${NAME} || die "Could not basecall the data"

echo "Done"