#!/bin/bash
#SBATCH --job-name=sort_seqs
#SBATCH --account=fc_tsutsuifca
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --open-mode=append
#SBATCH --requeue

gene_set=nde
MSA_DIR=MSA_${gene_set}
BASE_DIR=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/
WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/${MSA_DIR}/proteins_all_taxa/aligned/
results=final_prot_set.txt

##Step 4. Sort all final prot files by ID to ensure identical order for later concatenation 
while file=; IFS=$' \t\r\n' read -r file || [[ $file ]]; do
	ID="${file%.nuc.final.fas}"
	seqkit sort -n ${WD}nuc_aligned/${ID}.nuc.final.consensus.fas > ${WD}nuc_aligned/${ID}.tmp
	mv ${WD}nuc_aligned/${ID}.tmp ${WD}nuc_aligned/${ID}.nuc.final.consensus.fas

done<${WD}${results}