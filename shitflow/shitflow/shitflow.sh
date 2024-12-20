#!/bin/bash

die() {
    echo "$1"
    exit 1
}

usage () {
    echo "Usage: $0 <prefix_sample> <GADI_PBS_ARGS>"
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

BRENNER_DATA=/directflow/KCCGGenometechTemp/projects/iradev/operation_cornetto/shitflow_hasindu/
FRIDGE_TMP=/data3/cornetto/
GADI_DATA=/g/data/ox63/hasindu/cornetto/shitflow
GADI_SCRIPT=/g/data/ox63/hasindu/cornetto/cornetto/shitflow/hifiasm.pbs.sh

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


if [ -z "$( ls -A ${NAME}_duplex_out/split_blow5_failed/ )" ]; then
   echo "All channels are done"
else
    echo "Some channels failed"
    qsub -sync y $SCRIPT_PATH/dorado_duplex_retry.sge.sh ${NAME} || echo "even the retry for the failed channels failed"
fi

qsub -sync y $SCRIPT_PATH/get_duplex_and_simplex_reads.sge.sh ${NAME} || die "Could not basecall the data"

ssh gadi "mkdir ${GADI_DATA}/${NAME}" || die "gadi ssh failed"
scp ${NAME}.duplex_reads.fastq gadi-dm:${GADI_DATA}/${NAME}/ || die "copying to Gadi failed"

if [ -n "${GADI_PBS_ARGS}" ]; then
    GADI_COMMAND="cd ${GADI_DATA}/${NAME}/ && qsub -v ${GADI_PBS_ARGS} ${GADI_SCRIPT}"
    echo "Running on gadi: ${GADI_COMMAND}"
    ssh gadi "${GADI_COMMAND}" || die "gadi qsub failed"
    echo "Handed over work to gadi"
else
    echo "No GADI_PBS_ARGS provided. Not running on gadi"
fi


