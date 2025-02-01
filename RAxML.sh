#!/bin/bash
#SBATCH --job-name=RAxML
#SBATCH --account=fc_tsutsuifca
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --exclusive
#SBATCH --open-mode=append
#SBATCH -o raxml-log/slurm-%j.out
#SBATCH --requeue

MSA_DIR=MSA_${gene_set}
# MSA_DIR=MSA
WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/${MSA_DIR}/proteins_all_taxa/aligned/nuc_aligned/
base_DIR=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/RAxML_sbatch_dir/
gene_set=$1

##Step 1. Generate trees for each gene in input file
while file=; IFS=$' \t\r\n' read -r file || [[ $file ]]; do

	echo ${file}
	rm ${base_DIR}*${file}
	# raxmlHPC-PTHREADS -m GTRCAT -n ${gene_set}_genes_all_seqs.nuc.final.nodup.consensus -N 100 -p 14130 -s ${WD}all_seqs.nuc.final.consensus.fas -T 56 -x 633 -f a
	raxmlHPC-PTHREADS -m GTRCAT -n ${file} -N 100 -p 14130 -s ${WD}${file} -T 20 -x 633 -f a
	echo "Job completed successfully"

done<${gene_set}