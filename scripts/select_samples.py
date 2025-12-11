import random
all_samples_file = 'raw_data/all_samples.txt'
chosen_samples_file = 'raw_data/chosen_samples.txt'
with open(all_samples_file) as f:
    samples = [s.strip() for s in f if s.strip()]
random.seed(42)                                           # for reproducibility in sample selection
sel = random.sample(samples, 100)
with open(chosen_samples_file,'w') as out:
    out.write("\n".join(sel)+"\n")
print("Selected sample count:", len(sel))
print("Wrote:", chosen_samples_file)