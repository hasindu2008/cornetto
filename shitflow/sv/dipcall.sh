#!/bin/bash

set -o pipefail

DIPCALL=/install/dipcall.kit/run-dipcall
BGZIP=bgzip
TABIX=tabix
BCFTOOLS=bcftools

REF=/genome/chm13v2.0.fa

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <pat.fa> <mat.fa>"
    exit 1
fi

PAT=${1}
MAT=${2}

die() {
    echo "$1" >&2
    echo
    exit 1
}

test -x ${DIPCALL} || die "DIPCALL not found: ${DIPCALL}"
${BGZIP} --version || die "BGZIP not found: ${BGZIP}"
${TABIX} --version || die "TABIX not found: ${TABIX}"
${BCFTOOLS} --version || die "BCFTOOLS not found: ${BCFTOOLS}"

test -e ${REF} || die "Reference genome not found: ${REF}"
test -e ${PAT} || die "Parental assembly not found: ${PAT}"
test -e ${MAT} || die "Maternal assembly not found: ${MAT}"

test -e dip.mak && die "dip.mak already exists. Delete the content in the current directory first."

${DIPCALL} dip ${REF} ${PAT} ${MAT}  > dip.mak || die "Failed to create dip.mak"

make -j2 -f dip.mak || die "Failed to run make on dip.mak"

VCF=dip.dip.vcf.gz
test -e ${VCF} || die "VCF file not found: ${VCF}"

${BCFTOOLS} norm -m-any ${VCF} > split.vcf  || die "Failed to normalize VCF"
cat split.vcf | grep "^#" > structural_split.vcf || die "Failed to create header for structural.vcf"
cat split.vcf | grep -v "^#" | awk '{if(length($4)>50 || length($5)>50) print $0}' >> structural_split.vcf || die "Failed to filter structural variants"
${BGZIP} structural_split.vcf || die "Failed to compress structural_split.vcf"
${TABIX} structural_split.vcf.gz || die "Failed to index structural_split.vcf.gz"

echo "Structural variants have been filtered and indexed in structural_split.vcf.gz"

