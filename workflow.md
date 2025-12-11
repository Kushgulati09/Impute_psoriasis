- Downloaded 1000G Phase 3 chr22 VCF (v5b) and index from official FTP.
- Verified v5b as the correct file based on FTP directory listing.
- Commands:
- wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
- wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbiwo

- Created "all_sample.txt" from the downloaded genotype VCF file using:
- bcftools query -l ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz > all_samples.txt
- wc -l ./all_samples.txt  (output: 2504 (as expected))

- Created Python script scripts/select_samples.py to randomly select 100 samples (seed=42).
- Ran script to generate raw_data/chosen_samples.txt.