#!/usr/bin/env python3

################################################
# Script to evaluate imputation performance

# For each variant that exists in both:
# - the true dense target VCF
# - the imputed VCF

# It computes: Empirical dosage R²
# i.e. how well the imputed dosages match the true dosages across individuals.

# Usage:
# 
# mkdir -p results/evaluation
# 
# python scripts/evaluate_imputation.py \
#   --truth target/target_dense_chr22.vcf.gz \
#   --imputed results/imputed_chr22_biallelic.vcf.gz \
#   --out results/evaluation/imputation_accuracy.tsv
################################################

import argparse # For parsing command-line arguments
import numpy as np
import pandas as pd
import cyvcf2 # For reading VCF files
from scipy.stats import pearsonr # For computing Pearson correlation coeff r, ( its sq = r² is standard metric for imputation)

def dosage_from_gt(gt):
    """
    Convert GT array (e.g., [0,1], [1,1]) to dosage (0,1,2).
    Missing -> np.nan
    """
    if gt is None or -1 in gt: # missing genotype
        return np.nan # marking them as NaN so we can mask them later
    return sum(gt) #giving dosage as sum of alleles


def main(args):
    truth_vcf = cyvcf2.VCF(args.truth) #loading both truth and imputed VCFs
    imp_vcf = cyvcf2.VCF(args.imputed)

    # Index imputed variants by (chr, pos, ref, alt)
    imp_map = {} #dict for quick lookup of imputed verson of variant
    for v in imp_vcf:
        key = (v.CHROM, v.POS, v.REF, v.ALT[0]) #better practise than using ID
        ds = v.format("DS") #get imputed dosages
        dr2 = v.INFO.get("DR2", np.nan) #Beagle’s internal imputation quality estimate per variant (not ground truth)
        if ds is not None:
            imp_map[key] = (ds.flatten(), dr2) #flattening for easy alignment with true dosages

    records = []

    for v in truth_vcf:
        key = (v.CHROM, v.POS, v.REF, v.ALT[0])
        if key not in imp_map:
            continue #skip if variant not in imputed VCF

        imp_ds, dr2 = imp_map[key] #tuple unpacking

        # true dosages from GT
        true_ds = np.array([dosage_from_gt(gt[:2]) for gt in v.genotypes]) #converts GT to dosage (0,1,2)
        # v.genotypes[:, :2] → take only allele indices, ignore phased info

        # align lengths & remove missing
        mask = ~np.isnan(true_ds) & ~np.isnan(imp_ds)
        if mask.sum() < 3: # 3 is the safe minimum to compute correlation 
            continue

        true_ds = true_ds[mask]
        imp_ds = imp_ds[mask]

        # skip monomorphic
        # If all samples have the same genotype → variance = 0, i.e. no info to compute correlation
        if np.std(true_ds) == 0 or np.std(imp_ds) == 0:
            continue

        r, _ = pearsonr(true_ds, imp_ds)
        r2 = r ** 2

        # Each row corresponds to one variant
        records.append({
            "chrom": v.CHROM,
            "pos": v.POS,
            "ref": v.REF,
            "alt": v.ALT[0],
            "n_samples": mask.sum(),
            "r2_empirical": r2,
            "dr2_beagle": dr2
        })

    df = pd.DataFrame(records)
    df.to_csv(args.out, sep="\t", index=False) #saving results to tsv
    print(f"Saved results to {args.out}")
    print(df.describe()) #printing summary statistics (mean r², distribution, range)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--truth", required=True, help="Dense truth VCF")
    parser.add_argument("--imputed", required=True, help="Imputed biallelic VCF")
    parser.add_argument("--out", required=True, help="Output TSV")
    args = parser.parse_args()
    main(args)

