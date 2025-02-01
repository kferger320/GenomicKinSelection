aimport sys, os, re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
from collections import defaultdict


MSA_dir = sys.argv[1] #MSA or MSA_queen
WD=f"/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/{MSA_dir}/proteins_all_taxa/aligned/"

protein_master_file=f"{WD}final_prot_set.txt"
protein_master_list = pd.read_csv(protein_master_file, header=None)[0].to_list()

##Step 5. Concatenate all (filled) alignments together to create a single alignment file
# List to store sequences from each file
sequences = defaultdict(list)
records = []

# Read sequences from each input file
for file in protein_master_list:
	ID=file.split(".nuc.final.fas")[0]
	with open(f"{WD}nuc_aligned/{ID}.nuc.final.consensus.fas", "r") as handle:
		for i, record in enumerate(SeqIO.parse(handle, "fasta")):
			seq_acc = re.search(r'GC[^_]*_\d+', record.id).group(0)
			sequences[f"{seq_acc}_{i}"].append(str(record.seq))

# for key in sequences.keys():
# 	print(key)
# print(len(sequences.keys()))
for acc, seq_list in sequences.items():
	#combine all sequences into one string, write a record for each seq
	records.append(SeqRecord(Seq(''.join(seq_list)), id=acc, description=''))

# Write the combined sequences to a single output file
with open(f"{WD}nuc_aligned/all_seqs.nuc.final.consensus.fas", "w") as output_handle:
	SeqIO.write(records, output_handle, "fasta")

