#!/bin/bash

die () {
    echo "$1" >&2
    echo
    exit 1
}

source /install/truvari/bin/activate || die "Failed to activate truvari environment"
truvari bench -b ../hg002q100/structural_split.vcf.gz -c structural_split.vcf.gz -f /genome/chm13v2.0.fa -o output_dir/ || die "truvari bench failed"

