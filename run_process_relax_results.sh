#!/bin/bash
#SBATCH --job-name=process_relax
#SBATCH --account=fc_tsutsuifca
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --open-mode=append
#SBATCH --requeue


script_DIR=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/scripts/

for bin in social_poly polygyne high_poly; do
    for gene_set in worker queen nde; do
        python ${script_DIR}process_relax_results.py ${gene_set} ${bin}
    done
done