# Imputation + Psoriasis PRS 

This repository implements a minimal, reproducible DNA genotype imputation and polygenic risk score (PRS) workflow using public human genotype data. Starting from dense DNA variant data (VCF), the pipeline simulates array-like genotypes, performs quality control, phasing, and imputation, and computes a disease-specific PRS (psoriasis) using published PGS Catalogue weights. The project is designed as an educational and benchmarking exercise to demonstrate core population-genetics and PRS concepts on a small dataset. This repository is a work in progress and is actively being updated. For a step-by-step explanation of each stage and the reasoning behind it, please follow the `workflow.md` file.

Repository to reproduce imputation benchmark and PRS workflow (chr22 subset).
Contents:
- `data/` : small datasets, DO NOT commit large VCFs here (use external storage)
- `raw_data/` : original downloads (ignored)
- `ref/`, `target/` : reference and target VCFs derived from 1000G
- `scripts/` : pipeline scripts (bash / python)
- `results/` : plots and evaluation outputs
- `docs/` : notes and methods
- `tools/` : tool dependencies (java .jar file) 

### Environment setup

All analyses were performed in a Conda environment documented in `environment.yml`.

To recreate the environment:
```bash
conda env create -f environment.yml
conda activate Impute
```