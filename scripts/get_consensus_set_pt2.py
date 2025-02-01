import os

# Directory containing the files
gene_set="nde"
# gene_set="queen"
in_dir = f"/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/{gene_set}_gene_seqs/"
out_dir = f"/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/MSA_{gene_set}/"

# List to hold the keys from each file
all_keys = []

# Loop through each file in the directory
for filename in os.listdir(in_dir):
    if filename.endswith("Mphar_prot_matches.txt"):
        filepath = os.path.join(in_dir, filename)
        with open(filepath, "r") as file:
            keys = {line.strip().split()[0] for line in file}
            all_keys.append(keys)

# Find the intersection of all key sets
consensus_keys = set.intersection(*all_keys)

# Output the consensus keys
with open(f"{out_dir}consensus_Mphar_prot_set.txt", "w") as outfile:
    for key in consensus_keys:
        outfile.write(f"{key}\n")