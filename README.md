# Imputation + Psoriasis PRS 

This repository implements a minimal, **reproducible DNA genotype imputation and polygenic risk score (PRS) workflow** using public human genotype data. Starting from dense DNA variant data (VCF), the pipeline simulates array-like genotypes, performs quality control, phasing, and imputation, and computes a disease-specific PRS (psoriasis) using published PGS Catalogue weights. 

---

The project is designed as an educational and benchmarking exercise to demonstrate core population-genetics and PRS concepts on a small dataset, restricted to **chromosome 22**. 

---

For a step-by-step explanation of each stage and the reasoning behind it, please follow the `workflow.md` file.

---

## Repository structure:
```text
├── environment.yml # Conda environment definition
├── README.md # Project overview (this file)
├── workflow.md # Detailed pipeline documentation
├── raw_data/ # Downloaded raw VCFs (not tracked)
├── ref/ # Reference panel data, VCFs (not tracked) 
├── target/ # Target data and pseudo-array VCFs
├── scripts/ # Pipeline and analysis scripts
├── results/ # Imputation outputs, evaluation, PRS
├── data/ # External score files (not tracked)
└── tools/ # External tools (e.g. Beagle JAR)
```

### Environment setup

All analyses were performed in a Conda environment documented in `environment.yml`.

To recreate the environment:
```bash
conda env create -f environment.yml
conda activate Impute
```