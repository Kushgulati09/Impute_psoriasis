import pandas as pd
import matplotlib.pyplot as plt

# Load the imputation accuracy results
df = pd.read_csv("results/evaluation/imputation_accuracy.tsv", sep="\t")

# Basic statistics
print("Number of variants evaluated:", len(df))
print("\nEmpirical R² summary:")
print(df["r2_empirical"].describe()) #printing summary statistics (mean r², distribution, range)

print("\nBeagle DR2 summary:")
print(df["dr2_beagle"].describe()) #printing summary statistics (mean DR2, distribution, range)

print("\nCorrelation between empirical R² and Beagle DR2:")
print(df[["r2_empirical", "dr2_beagle"]].corr()) #correlation matrix

#  Same information saved to a text file for reproducibility
with open("results/evaluation/imputation_accuracy_summary.txt", "w") as f:
    f.write(f"Number of variants evaluated: {len(df)}\n\n")
    f.write("Empirical R² summary:\n")
    f.write(df["r2_empirical"].describe().to_string()) 
    f.write("\n\nBeagle DR2 summary:\n")
    f.write(df["dr2_beagle"].describe().to_string())
    f.write("\n\nCorrelation matrix:\n")
    f.write(df[["r2_empirical", "dr2_beagle"]].corr().to_string())

# Scatter plot of Beagle DR2 vs empirical R²
plt.figure()
plt.scatter(df["dr2_beagle"], df["r2_empirical"], s=5, alpha=0.3) #s=5 for point size, alpha=0.3 for transparency
plt.xlabel("Beagle DR2")
plt.ylabel("Empirical dosage R²")
plt.title("Imputation accuracy: DR2 vs empirical R²")
plt.savefig("results/evaluation/imputation_accuracy_scatter.png", dpi=300) #saving as well
plt.show()
