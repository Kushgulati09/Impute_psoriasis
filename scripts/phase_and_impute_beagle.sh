#!/usr/bin/env bash
set -euo pipefail   #stopping on errors, undefined variables, or failed pipes

################################################################################
# phase_and_impute_beagle.sh
#
# Run Beagle to phase and impute the pseudo-array using the
# reference VCF created earlier (ref/ref_chr22.vcf.gz).
#
# Usage:
#   bash scripts/phase_and_impute_beagle.sh \
#       --target target/pseudoarray_maf5_final.vcf.gz \
#       --ref ref/ref_chr22.vcf.gz \
#       --out results/imputed_chr22 \
#       --beagle tools/beagle.jar
#
################################################################################

# Input arguments
TARGET_VCF="target/pseudoarray_maf5_final.vcf.gz"
REF_VCF="ref/ref_chr22.vcf.gz"
OUT_PREFIX="results/imputed_chr22"
BEAGLE_JAR="tools/beagle.jar"
JAVA_MEM="4g" # Java memory allocation for Beagle (change as available RAM allows)

# Function to display usage information
usage() {
  echo "Usage: $0 [--target <target_vcf>] [--ref <ref_vcf>] [--out <out_prefix>] [--beagle <beagle_jar>] [--mem <java_mem>]"
  exit 1
}

# Parsing command-line arguments
while [[ $# -gt 0 ]]; do # While there are arguments left to process
  case $1 in #define cases for each argument expected
    --target) TARGET_VCF="$2"; shift 2;; # Shift to next argument after processing
    --ref) REF_VCF="$2"; shift 2;;
    --out) OUT_PREFIX="$2"; shift 2;;
    --beagle) BEAGLE_JAR="$2"; shift 2;;
    --mem) JAVA_MEM="$2"; shift 2;;
    *) usage ;; # If unknown argument, show usage and exit
  esac
done

# Checking if input files exist
if [ ! -f "$TARGET_VCF" ]; then echo "Target VCF not found: $TARGET_VCF"; exit 1; fi
if [ ! -f "$REF_VCF" ]; then echo "Reference VCF not found: $REF_VCF"; exit 1; fi
if [ ! -f "$BEAGLE_JAR" ]; then echo "Beagle jar not found: $BEAGLE_JAR"; exit 1; fi

mkdir -p "$(dirname "$OUT_PREFIX")"

echo ">>> Running Beagle imputation"
echo "Target: $TARGET_VCF"
echo "Reference: $REF_VCF"
echo "Output prefix: $OUT_PREFIX"
echo "Beagle jar: $BEAGLE_JAR"

# Running Beagle:
# gt=target genotypes, 
# ref=reference, 
# out=prefix for all output files
# Beagle outputs VCF + .gp, with DS/GP fields (dosage/probabilities)
# change nthreads as availability of CPU cores
java -Xmx${JAVA_MEM} -jar "$BEAGLE_JAR" \
     gt="$TARGET_VCF" \
     ref="$REF_VCF" \
     out="$OUT_PREFIX" \
     nthreads=4

echo ">>> Beagle finished. Output (gzipped VCF): ${OUT_PREFIX}.vcf.gz"