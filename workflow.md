## Download 1000G chr22 VCF
- Downloaded 1000G Phase 3 chr22 VCF (v5b) and index from official FTP.
- Verified v5b as the correct file based on FTP directory listing.
- Commands:
    - wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
    - wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbiwo

## Extract Sample list
- Created "all_sample.txt" from the downloaded genotype VCF file using:
    - bcftools query -l ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz > all_samples.txt

## Select 100 samples
- Created Python script "scripts/select_samples.py" to randomly select 100 samples (seed=42).
- Ran script to generate raw_data/chosen_samples.txt.
 
## 80 / 20 split 
- Deterministic split using commands:
    - head -n 80 chosen_samples.txt > ref_samples.txt
    - tail -n 20 chosen_samples.txt > target_samples.txt 

## Extract Ref. and Target Dense Truth VCFs
- Extract reference individuals (80)(these will act as the reference panel):
    - bcftools view -S ref_samples.txt -Oz -o ../ref/ref_chr22.vcf.gz ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
    - tabix -p vcf ../ref/ref_chr22.vcf.gz

- Extract target dense truth (full variants retained to later test accuracy of imputation):
    - bcftools view -S target_samples.txt -Oz -o ../target/target_dense_chr22.vcf.gz ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
    - tabix -p vcf ../target/target_dense_chr22.vcf.gz

## Construct Pseudo-array (MAF > 5%) and QC
- Created a bash script `scripts/make_pseudoarray_maf5.sh` which performs given tasks and maintains reproducibility.
- Created a reproducible pseudo-array dataset from the dense target chr22 VCF using an MAF > 5% selection strategy.  
- This simulates an array-like SNP panel without relying on external manifests and follows a transparent, reproducible pipeline.
- Steps Performed in bash script:
    1. Annotated the dense target VCF with allele frequency and MAF using `bcftools +fill-tags`.
    2. Selected variants with `INFO/MAF > 0.05` to produce an initial pseudo-array.
    3. Converted pseudo-array VCF to PLINK (PLINK version 1.9) and applied QC:
        - `--mind 0.05` (remove samples with >5% missingness)
        - `--geno 0.05` (remove SNPs missing in >5% samples)
        - `--maf 0.01` (remove rare / unstable SNPs after missingness filtering)
    4. Converted QC’ed PLINK set back to VCF for downstream phasing/imputation.
    5. Enforced biallelic SNPs only.
    6. Identified and removed ambiguous SNPs.
    7. Produced final pseudo-array VCF ready for phasing:
        - `target/pseudoarray_maf5_final.vcf.gz`


