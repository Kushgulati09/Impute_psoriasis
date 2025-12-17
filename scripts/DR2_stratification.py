import pandas as pd
import numpy as np

df = pd.read_csv("results/evaluation/imputation_accuracy.tsv", sep="\t")

# defining bins
bins = [0, 0.2, 0.4, 0.6, 0.8, 1.0]
df["dr2_bin"] = pd.cut(df["dr2_beagle"], bins=bins)

summary = (
    df.groupby("dr2_bin")
      .agg(
          n_variants=("r2_empirical", "count"),
          mean_r2=("r2_empirical", "mean"),
          median_r2=("r2_empirical", "median")
      )
)

summary.to_csv("results/evaluation/DR2_bin_imputation_accuracy.tsv", sep="\t", index=True)
print(summary)

