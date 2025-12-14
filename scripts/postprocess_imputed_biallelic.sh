#!/usr/bin/env bash
set -euo pipefail

################################################################################
# postprocess_imputed_biallelic.sh
#
# Post-process Beagle-imputed VCF to retain only biallelic variants.
#
# Usage:
#   bash scripts/postprocess_imputed_biallelic.sh \
#       results/imputed_chr22.vcf.gz \
#       results/imputed_chr22_biallelic.vcf.gz
# 
################################################################################
# Input arguments
IN_VCF=$1
OUT_VCF=$2

# Command
bcftools view -m2 -M2 "$IN_VCF" -Oz -o "$OUT_VCF"
tabix -p vcf "$OUT_VCF"

echo "Biallelic VCF written to $OUT_VCF"