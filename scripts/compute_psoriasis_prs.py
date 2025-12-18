#!/usr/bin/env python3

#################################################################
# Compute Polygenic Risk Scores (PRS) for Psoriasis
# 
# Usage:
# python scripts/compute_psoriasis_prs.py \
#   --vcf results/imputed_chr22_biallelic.vcf.gz \
#   --score data/PGS001312.txt \
#   --out results/psoriasis_prs.tsv
#################################################################

import argparse
import pandas as pd
import numpy as np
import cyvcf2

# defining strand-ambiguous SNP pairs
AMBIGUOUS = {("A","T"), ("T","A"), ("C","G"), ("G","C")}

def main(args):
    # Load PGS file (skip metadata lines starting with '#')
    pgs = pd.read_csv(
        args.score,
        sep="\t",
        comment="#" # Skip metadata lines starting with '#'
    )

    # Keep only needed columns
    pgs = pgs[[
        "chr_name", # chromosome name
        "chr_position", # position on chromosome
        "effect_allele", #the allele carrying the risk
        "other_allele", #the alternative allele
        "effect_weight" #how strongly the effect_allele influences risk
    ]]

    # Drop rows without valid genomic position (e.g. haplotypes)
    pgs = pgs.dropna(subset=["chr_name", "chr_position"])

    # Ensure positions are integers
    pgs["chr_position"] = pgs["chr_position"].astype(int)

    # Index by (chr, pos) | Building a look-up dictionary for fast access
    pgs_map = {
        (str(r.chr_name), int(r.chr_position)):
            (r.effect_allele, r.other_allele, float(r.effect_weight))
        for r in pgs.itertuples() # Iterate over rows
    }

    vcf = cyvcf2.VCF(args.vcf) # Load VCF file
    samples = vcf.samples # List of sample IDs in VCF
    prs = np.zeros(len(samples)) #one PRS per sample, initialized to zero, incremented below SNP-by-SNP

    used = 0 
    dropped_ambiguous = 0 
    dropped_mismatch = 0

    for v in vcf: # Each v is one imputed SNP
        key = (v.CHROM, v.POS) 
        if key not in pgs_map: # only variants that are imputed AND in the PGS file
            continue

        # Ensure biallelic again (just to be safe)
        if len(v.ALT) != 1:
            dropped_mismatch += 1
            continue
        
        # Get reference and alternate alleles
        ref = v.REF
        alt = v.ALT[0]

        # Drop ambiguous SNPs
        if (ref, alt) in AMBIGUOUS:
            dropped_ambiguous += 1
            continue

        effect_allele, other_allele, weight = pgs_map[key]

        # Important
        # Harmonize alleles using both effect and other allele
        if effect_allele == alt and other_allele == ref:
            signed_weight = weight
        elif effect_allele == ref and other_allele == alt:
            signed_weight = -weight
        else:
            dropped_mismatch += 1
            continue
        
        ds = v.format("DS") #dosage format
        if ds is None:
            continue

        prs += ds.flatten() * signed_weight # Update PRS for all samples
        used += 1
    
    # Output PRS results
    out = pd.DataFrame({
        "sample": samples,
        "PRS": prs
    })
    out.to_csv(args.out, sep="\t", index=False)

    print("PRS computation complete")
    print(f"Variants used: {used}")
    print(f"Dropped ambiguous SNPs: {dropped_ambiguous}")
    print(f"Dropped allele mismatches: {dropped_mismatch}")
    print(f"Output written to: {args.out}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser() #CLI (Command Line Interface) parser
    parser.add_argument("--vcf", required=True, help="Imputed biallelic VCF")
    parser.add_argument("--score", required=True, help="PGS score file (PGS001312)")
    parser.add_argument("--out", required=True, help="Output PRS file")
    args = parser.parse_args()
    main(args) # Run the main function with parsed arguments
