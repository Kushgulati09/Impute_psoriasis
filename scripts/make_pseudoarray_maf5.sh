#!/usr/bin/env bash
set -euo pipefail   #stopping on errors, undefined variables, or failed pipes

################################################################################
# make_pseudoarray_maf5.sh
#
# Create a pseudo-array from a dense target VCF by selecting variants with
# MAF > 5%, applying QC, removing ambiguous SNPs, and exporting a final VCF.
#
# Usage:
#   bash scripts/make_pseudoarray_maf5.sh target/target_dense_chr22.vcf.gz
#
################################################################################
# Input argument check
if [ "$#" -ne 1 ]; then  # Expect exactly one argument 
    echo "Usage: $0 <target_dense_vcf.gz>"
    exit 1 # Exit with error code
fi

# Input file
TARGET_DENSE_VCF=$1

# Check if input file exists
if [ ! -f "$TARGET_DENSE_VCF" ]; then
    echo "Error: dense VCF not found at $TARGET_DENSE_VCF"
    exit 1
fi

echo ">>> Running pseudo-array creation for: $TARGET_DENSE_VCF"

# Output files
OUT_DIR="target"
mkdir -p "$OUT_DIR"

# Define output file paths
VCF_WITH_MAF="${OUT_DIR}/target_withMAF.vcf.gz"
VCF_PSEUDO="${OUT_DIR}/target_pseudoarray_maf5.vcf.gz"
PLINK_BASE="${OUT_DIR}/pseudoarray_maf5"
PLINK_QC="${OUT_DIR}/pseudoarray_maf5_qc"
VCF_QC="${OUT_DIR}/pseudoarray_maf5_qc.vcf.gz"
VCF_BIALLELIC="${OUT_DIR}/pseudoarray_maf5_qc_biallelic.vcf.gz"
AMBIG_LIST="${OUT_DIR}/ambiguous_snps_maf5.txt"
VCF_FINAL="${OUT_DIR}/pseudoarray_maf5_final.vcf.gz"

# bcftools +fill-tags: computes variant summary statistics from samples and writes them into the VCF INFO field.
# -t MAF: specifically adds MAF (minor allele frequency) for each site
# Output compressed (-Oz) and indexed with (tabix) for quicker query in downstream steps.
echo ">>> Step 1: Annotate MAF using bcftools +fill-tags"
bcftools +fill-tags "$TARGET_DENSE_VCF" -Oz -o "$VCF_WITH_MAF" -- -t MAF
tabix -p vcf "$VCF_WITH_MAF"

# Results: VCF containing only common variants — i.e. pseudo-array candidate set
# Arrays are enriched for common, well-tagging SNPs. Using high-MAF SNPs gives a realistic set for imputation tests.
echo ">>> Step 2: Select variants with INFO/MAF > 0.05"
bcftools view -i 'INFO/MAF>0.05' "$VCF_WITH_MAF" -Oz -o "$VCF_PSEUDO"
tabix -p vcf "$VCF_PSEUDO"

# How many variant records passed the MAF>5% filter
echo ">>> Variant count after MAF>5% filter:"
bcftools view -H "$VCF_PSEUDO" | wc -l

# Convert to PLINK binary format (.bed/ .bim/ .fam) for QC steps using --make-bed 
    # .bed:	stores Actual genotype data (compact binary matrix)
    # .bim:	stores Variant information (chrom, pos, alleles)
    # .fam:	stores Sample metadata information (IDs, sex, phenotype)
echo ">>> Step 3: Convert pseudo-array VCF to PLINK for QC"
plink --vcf "$VCF_PSEUDO" --make-bed --out "$PLINK_BASE"

echo ">>> Step 4: Basic QC (sample missingness, variant missingness, MAF>=0.01)"
plink --bfile "$PLINK_BASE" \
       --mind 0.05 --geno 0.05 --maf 0.01 \
       --make-bed --out "$PLINK_QC"
# Explanation of QC thresholds:
    # --mind 0.05 = removes samples with >5% missing genotypes to avoid unreliable individuals in phasing and imputation.
    # --geno 0.05 = removes SNPs missing in >5% of samples to ensure markers used for imputation have high call quality.
    # --maf 0.01 = ensures eliminating monomorphic / unstable SNPs after missingness filtering
    #              and drops extremely rare variants (<1% MAF) that provide little tagging power and often reduce imputation accuracy.

# Note: We first selected MAF>5% from the dense truth. After removing some samples (by --mind) and variants (by --geno), 
# allele counts can change and (--maf 0.01) is a safety net to remove sites that became almost monomorphic because of missingness 
# or small sample size. It doesn’t contradict the earlier MAF>0.05.

# Back to a compressed VCF and indexed it. This is the QC'd pseudo-array VCF.
echo ">>> Step 5: Convert QC'd PLINK files back to VCF"
plink --bfile "$PLINK_QC" --recode vcf --out "${OUT_DIR}/pseudoarray_maf5_qc"
bgzip "${OUT_DIR}/pseudoarray_maf5_qc.vcf"
tabix -p vcf "$VCF_QC"

echo ">>> Step 6: Keep only biallelic SNPs"
bcftools view -v snps -m 2 -M 2 "$VCF_QC" -Oz -o "$VCF_BIALLELIC"
tabix -p vcf "$VCF_BIALLELIC"
# -v snps : keep only SNPs (drop indels)
# -m 2 -M 2 : keep only sites with exactly 2 alleles (REF + 1 ALT)
# (Multiallelic sites complicate allele matching and dosage calculation)

echo ">>> Step 7: Identify ambiguous SNPs (A/T, T/A, C/G, G/C)"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "$VCF_BIALLELIC" | \
  awk '($4$5=="AT"||$4$5=="TA"||$4$5=="CG"||$4$5=="GC"){print $1"\t"$2}' > "$AMBIG_LIST"
# This command lists sites where the REF/ALT allele pair is A/T, T/A, C/G, or G/C — these are strand-ambiguous SNPs.
# and writes chromosome + position pairs into (AMBIG_LIST) for removal in the next step.

AMBIG_COUNT=$(wc -l < "$AMBIG_LIST")
echo ">>> Found $AMBIG_COUNT ambiguous SNPs"

if [ "$AMBIG_COUNT" -gt 0 ]; then
    echo ">>> Removing ambiguous SNPs"
    bcftools view -T ^"$AMBIG_LIST" "$VCF_BIALLELIC" -Oz -o "$VCF_FINAL"
    # (-T ^file) tells bcftools to exclude the positions in AMBIG_LIST file
else
    echo ">>> No ambiguous SNPs found, copying file"
    cp "$VCF_BIALLELIC" "$VCF_FINAL"
fi
tabix -p vcf "$VCF_FINAL"
# compressing

echo ">>> Final pseudo-array variant count:"
bcftools view -H "$VCF_FINAL" | wc -l

echo ">>> Pseudo-array creation complete."
echo ">>> Final output: $VCF_FINAL"


