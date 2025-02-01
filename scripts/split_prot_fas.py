import os, sys
import pandas as pd

gene_set = sys.argv[1] #worker, queen
genes_per_group = int(sys.argv[2])

##toggle if generating groups of conv probs gene subset 
conv_probs=True

if gene_set == "worker":
    MSA_DIR = "MSA"
else:
    MSA_DIR = f"MSA_{gene_set}"

base_dir = f"/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/"
WD = f"/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/{MSA_DIR}/proteins_all_taxa/aligned/nuc_aligned/"
tree_DIR=f"/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/RAxML_sbatch_dir/{gene_set}_genes/labeled/"


def write_group_files(files, genes_per_group, gene_set, species_set, conv_probs=True):

    if gene_set == "worker":
        MSA_DIR = "MSA"
    else:
        MSA_DIR = f"MSA_{gene_set}"

    WD = f"/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/{MSA_DIR}/proteins_all_taxa/aligned/nuc_aligned/"

    group = 1
    outline = []
    if conv_probs:
        out_prefix = f"{WD}{gene_set}_gene_gs_{genes_per_group}_{species_set}_conv_probs"
    else:
        out_prefix = f"{WD}{gene_set}_gene_gs_{genes_per_group}_{species_set}_no_conv_probs"

    for i, file in enumerate(files):
        if (i < genes_per_group * group) and (i >= genes_per_group * (group - 1)):
            outline.append(file)
            
        if i == group * genes_per_group:

            with open(f"{out_prefix}_set_{group}.txt", "w") as out:
                for line in outline:
                    out.write(line+"\n")

            group += 1
            outline = []
            #append last file
            outline.append(file)

    ##last file
    with open(f"{out_prefix}_set_{group}.txt", "w") as out:
        for line in outline:
            out.write(line+"\n")


def get_filenames(gene_list, gene_set):

    if gene_set == "worker":
        MSA_DIR = "MSA"
    else:
        MSA_DIR = f"MSA_{gene_set}"

    WD = f"/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/{MSA_DIR}/proteins_all_taxa/aligned/nuc_aligned/"

    files = []
    for root, dirs, filenames in os.walk(WD):
        # Check the depth of the current directory
        if root == WD:
            files.extend(os.path.join(root, file) for file in filenames)
            # Exit after processing the first level
            break

    gene_files = []
    for gene in gene_list:
        for file in files:
            if gene in file:
                gene_files.append(file)

    gene_files_final = [file for file in gene_files if "all_seqs" not in file and 
        "nuc.final.consensus.nodummy.fas" in file and 
        "reduced" not in file]
    
    return gene_files_final


def process_no_conv_probs(gene_set, species_set):
    conv_probs_genes = pd.read_csv(f"{base_dir}{gene_set}_{species_set}_conv_probs_genes.txt", sep="\t", header=None)[0].tolist()

    if gene_set == "worker":
        MSA_DIR = "MSA"
    else:
        MSA_DIR = f"MSA_{gene_set}"

    WD = f"/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/{MSA_DIR}/proteins_all_taxa/aligned/nuc_aligned/"

    files = []
    for root, dirs, filenames in os.walk(WD):
        # Check the depth of the current directory
        if root == WD:
            files.extend(os.path.join(root, file) for file in filenames)
            # Exit after processing the first level
            break

    gene_files_final = [file for file in files if "all_seqs" not in file and 
        "nuc.final.consensus.nodummy.fas" in file and 
        "reduced" not in file]    # print(len(files))
    
        ##added 
    no_conv_probs_files = []
    for file in gene_files_final: #check if gene experienced convergence probs
        found = False
        file_name = os.path.basename(file)
        query_gene = file_name.split('.')[0]

        if query_gene in conv_probs_genes:
            found=True

        if not found:
            no_conv_probs_files.append(file)
    
    return no_conv_probs_files


if conv_probs:

    ##separate processing for conv probs genes
    polygyne = pd.read_csv(f"{base_dir}{gene_set}_polygyne_conv_probs_genes.txt", sep="\t", header=None)[0].tolist()
    high_poly = pd.read_csv(f"{base_dir}{gene_set}_high_poly_conv_probs_genes.txt", sep="\t", header=None)[0].tolist()
    social_poly = pd.read_csv(f"{base_dir}{gene_set}_social_poly_conv_probs_genes.txt", sep="\t", header=None)[0].tolist()

    files = []

    for root, dirs, filenames in os.walk(WD):
        # Check the depth of the current directory
        if root == WD:
            files.extend(os.path.join(root, file) for file in filenames)
            # Exit after processing the first level
            break

    poly_files = get_filenames(polygyne, gene_set=gene_set)
    high_poly_files = get_filenames(high_poly, gene_set=gene_set)
    social_poly_files = get_filenames(social_poly, gene_set=gene_set)
    

    # ##generate files
    write_group_files(files=poly_files, conv_probs=True, genes_per_group=genes_per_group, gene_set=gene_set, species_set="polygyne")
    write_group_files(files=high_poly_files, conv_probs=True, genes_per_group=genes_per_group, gene_set=gene_set, species_set="high_poly")
    write_group_files(files=social_poly_files, conv_probs=True, genes_per_group=genes_per_group, gene_set=gene_set, species_set="social_poly")


else: # generate files for all genes instead

    poly_no_conv_probs_files = process_no_conv_probs(gene_set, "polygyne")
    high_poly_no_conv_probs_files = process_no_conv_probs(gene_set, "high_poly")
    social_poly_no_conv_probs_files = process_no_conv_probs(gene_set, "social_poly")

    # print(len(no_conv_probs_files))

    ##generate files
    # write_group_files(files, conv_probs=False, genes_per_group=genes_per_group, gene_set=gene_set, species_set=NULL)

    #added
    write_group_files(files=poly_no_conv_probs_files, conv_probs=False, genes_per_group=genes_per_group, gene_set=gene_set, species_set="polygyne")
    write_group_files(files=high_poly_no_conv_probs_files, conv_probs=False, genes_per_group=genes_per_group, gene_set=gene_set, species_set="high_poly")
    write_group_files(files=social_poly_no_conv_probs_files, conv_probs=False, genes_per_group=genes_per_group, gene_set=gene_set, species_set="social_poly")