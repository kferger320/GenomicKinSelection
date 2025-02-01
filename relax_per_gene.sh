#!/bin/bash
#SBATCH --job-name=hyphy-analyses
#SBATCH --account=fc_tsutsuifca
#SBATCH --partition=savio2_htc
#SBATCH --cpus-per-task=1 
#SBATCH --ntasks=1
#SBATCH --qos=savio_normal
#SBATCH --time=72:00:00
#SBATCH --open-mode=append
#SBATCH --requeue

module load gcc openmpi

bin=$1 #worker, queen, nde
test_group=$2 #high_poly_unicolonial, polygyne, social_poly
gene_set=$3

if [ "$bin" = "worker" ]; then
    MSA_DIR="MSA"
else
    MSA_DIR="MSA_${bin}"
fi

base_DIR=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/
WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/${MSA_DIR}/proteins_all_taxa/aligned/nuc_aligned/
tree_DIR=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/RAxML_sbatch_dir/${bin}_genes/labeled/
# tree_DIR=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/RAxML_sbatch_dir/${bin}_genes/labeled/redo-12-29-24/
output_DIR=${WD}relax-srv/
# output_DIR=${WD}relax-srv/redo-12-29-24/

# ##Step 1. run relax with 1CPU for each gene; parallelized with ${gene_set}_gene_set_${i}.txt (~30 genes/file, looping through)
while file=; IFS=$' \t\r\n' read -r file || [[ $file ]]; do

	base=$(basename "$file" | cut -d'.' -f1-2 |cut -d'.' -f1) ##expects full path to be included with filename  
	found=$(find ${output_DIR} -type f -name "RELAX-*-${base}-${bin}-${test_group}")

	if [[ -z "$found" ]]; then
	
		##if running the General descriptive model, change --models from 'Minimal' to 'All'
		hyphy CPU=1 relax --alignment ${file} \
			--models Minimal \
			--tree ${tree_DIR}RAxML_bestTree.${base}-${bin}-${test_group} \
			--test Test \
			--reference Reference \
			--starting-points 5 \
			--grid-size 500 \
			--srv Yes \
			--output ${output_DIR}RELAX-${SLURM_JOB_ID}-${base}-${bin}-${test_group}
			# --output ${output_DIR}RELAX-${SLURM_JOB_ID}-${base}-${bin}-${test_group}-gen-model-rerun

	fi

done<${gene_set}

echo "Job completed successfully"