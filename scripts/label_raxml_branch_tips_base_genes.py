from ete3 import Tree
import sys, re, os, glob
import pandas as pd

group=sys.argv[1] #social_poly, polygyne, etc
gene_set=sys.argv[2] #worker, nde, queen
prots=sys.argv[3] ##eg. queen_gene_set_1.txt

print(group)
print(gene_set)
print(prots)

tree_dir = f"/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/RAxML_sbatch_dir/{gene_set}_genes/"
WD="/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/"

# ref_list = pd.read_csv(f"{WD}{group}_GC_refs.txt", sep="\t", header=None)[0].tolist()
# monogyne_list = pd.read_csv(f"{WD}monogyne_GC_refs.txt", sep="\t", header=None)[0].tolist()

##rerun with bin modifications 12/29/24
ref_list = pd.read_csv(f"{WD}{group}_GC_refs_redo.txt", sep="\t", header=None)[0].tolist()
monogyne_list = pd.read_csv(f"{WD}monogyne_GC_refs_redo.txt", sep="\t", header=None)[0].tolist()
prot_list = pd.read_csv(prots, sep="\t", header=None)[0].tolist()


def extract_id(input_str):
	# Define the regular expression pattern to match the ID
	pattern = r'lcl[0-9A-Za-z_.-]+cds_([0-9A-Za-z_.-]+)_([0-9A-Za-z_.-]+)_([0-9A-Za-z_.-]+)_([0-9A-Za-z_.-]+)'

	# Use re.search to find the pattern in the input string
	match = re.search(pattern, input_str)

	# If a match is found, return the captured group
	if match:
		return f"{match.group(1)}.{match.group(2)}"
	else:
		return None


#loop through each protein in input list
prot_labs = []
for prot in prot_list:
	# tree = Tree(f"{tree_dir}RAxML_bestTree.{prot}")
	prot_base = os.path.basename(prot)

	# pattern = f"{tree_dir}{prot_base}.*.nuc.final.consensus.nodummy.fas"
	pattern = f"{tree_dir}RAxML_bestTree.{prot_base}"
	matching_files = glob.glob(pattern)

	if len(matching_files > 0): #continue if RAxML tree file found with corresponding protein ID
		tree = Tree(matching_files[0])
		prot_lab = prot_base.split(".")[0]

		##triples since one path per bin (social poly, etc)
		if prot_lab in prot_labs:
			continue

		#prot lab not yet in list
		prot_labs.append(prot_lab)
		#check each leaf to see if that gene is from a reference species in our specified group
		for leaf in tree:

			if 'GC' in leaf.name:
				leaf_id = leaf.name.split('_')[0] + "_" + leaf.name.split('_')[1]
			
				# print(f"leaf id: {leaf_id}")

				#label genes found in test group as 'Test'
				# if leaf_id in gene_list:
				if leaf_id in ref_list:
					# Ensure the name is properly formatted with the {Test} outside the quotations
					leaf.name = f"'{leaf.name}'{{Test}}"
				
				##added, non-chemo 
				if leaf_id in monogyne_list:
					leaf.name = f"'{leaf.name}'{{Reference}}"
			

		##Write our final labeled tree 
		tree.write(format=1, outfile=f"{tree_dir}labeled/RAxML_bestTree.{prot_lab}-{gene_set}-{group}")
		# tree.write(format=1, outfile=f"{tree_dir}labeled/redo-12-29-24/RAxML_bestTree.{prot_lab}-{gene_set}-{group}")
	
	else: #corresponding tree file not found, skip
		print(f"Tree file not found for {pattern}; skipping")
		continue