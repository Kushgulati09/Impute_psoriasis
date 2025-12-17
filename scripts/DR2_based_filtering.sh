#!/usr/bin/env bash
set -euo pipefail

###########################################################################
# Script to filter imputed VCF based on Stratified Beagle DR2 threshold
# 
# Usage:
# bash scripts/DR2_based_filtering.sh \
#   results/imputed_chr22_biallelic.vcf.gz \
#   results/imputed_chr22_biallelic_dr2filt.vcf.gz \
#   0.8
###########################################################################
IN_VCF=$1
OUT_VCF=$2
THRESHOLD=${3:-0.8} # Default threshold is 0.8 if not provided

# Filter VCF based on DR2 threshold
echo "Filtering variants with DR2 >= ${THRESHOLD} from ${IN_VCF} to ${OUT_VCF}"
bcftools view -i "INFO/DR2>=${THRESHOLD}" "$IN_VCF" -Oz -o "$OUT_VCF"
tabix -p vcf "$OUT_VCF"

# Summary statistics
echo "Filtering complete."
echo "Variants retained:"
bcftools view -H "$OUT_VCF" | wc -l
