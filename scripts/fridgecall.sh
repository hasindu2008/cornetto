#!/bin/bash

# run inside a screen session

PREFIX=A_1
SAMPLE=QGXHXX240275

cd /data3/cornetto
slow5tools merge /data/${SAMPLE}/*/*/slow5/ -o ${PREFIX}_${SAMPLE}.blow5



