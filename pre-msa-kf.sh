#!/bin/bash
#SBATCH --job-name=pre-msa
#SBATCH --account=fc_tsutsuifca
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --open-mode=append
#SBATCH -o pre-msa-kf-log/slurm-%j.out
#SBATCH --requeue

##input 
WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/
ref=$1
ref_basename_nopath="${ref##*/}"
ref_basename="${ref_basename_nopath%.*}"
gene_set=nde #worker or queen
genes=${WD}${gene_set}_gene_seqs/${ref_basename}_longest_iso_cds*.${gene_set}_genes.fna

##Step 1. translate input nt seqs from worker_genes, concatenate and output .faa (translated seqs) and .fna (corresponding nt seqs, divisible by 3, frame taken into account)

python ${WD}scripts/pre-msa-kf.py ${genes}
##outputs: ${gene_set}_gene_seqs/pre-msa/{basename}.nuc_prot.fna and {basename}.nuc_prot.faa


