scripts_DIR=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/scripts/
base_DIR=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/
conv_probs=$1

#####WORKER#####
gene_set=worker
MSA_DIR=MSA
gene_set_DIR=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/${MSA_DIR}/proteins_all_taxa/aligned/nuc_aligned/
out_DIR=${gene_set_DIR}relax-srv/
# out_DIR=${gene_set_DIR}relax-srv/redo-12-29-24/

##polygyne
species_bin=polygyne
if [[ ${conv_probs} == "True" ]]; then

	for i in {1..17}; do ##change this range based on number of input group files
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_conv_probs_set_${i}.txt
	done
else 
	for i in {1..13}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-no-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_no_conv_probs_set_${i}.txt
	done

fi


# ##high_poly
species_bin=high_poly
if [[ ${conv_probs} == "True" ]]; then

	for i in {1..18}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_conv_probs_set_${i}.txt
	done
else
	for i in {1..11}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-no-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_no_conv_probs_set_${i}.txt
	done

fi


##social_poly
species_bin=social_poly
if [[ ${conv_probs} == "True" ]]; then

	for i in {1..18}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_conv_probs_set_${i}.txt
	done
else
	for i in {1..12}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-no-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_no_conv_probs_set_${i}.txt
	done

fi


#####QUEEN#####
gene_set=queen
MSA_DIR=MSA_${gene_set}
gene_set_DIR=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/${MSA_DIR}/proteins_all_taxa/aligned/nuc_aligned/
out_DIR=${gene_set_DIR}relax-srv/
# out_DIR=${gene_set_DIR}relax-srv/redo-12-29-24/

##polygyne
species_bin=polygyne

if [[ ${conv_probs} == "True" ]]; then

	for i in {1..22}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_conv_probs_set_${i}.txt
	done
else
	for i in {1..24}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-no-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_no_conv_probs_set_${i}.txt
	done

fi


##high_poly
species_bin=high_poly
if [[ ${conv_probs} == "True" ]]; then

	for i in {1..28}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_conv_probs_set_${i}.txt
	done
else
	for i in {1..19}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-no-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_no_conv_probs_set_${i}.txt
	done

fi


##social_poly
species_bin=social_poly
if [[ ${conv_probs} == "True" ]]; then

	for i in {1..24}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_conv_probs_set_${i}.txt
	done
else
	for i in {1..23}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-no-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_no_conv_probs_set_${i}.txt
	done

fi


#####NDE#####
gene_set=nde
MSA_DIR=MSA_${gene_set}
gene_set_DIR=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/${MSA_DIR}/proteins_all_taxa/aligned/nuc_aligned/
out_DIR=${gene_set_DIR}relax-srv/
# out_DIR=${gene_set_DIR}relax-srv/redo-12-29-24/

# ##polygyne
species_bin=polygyne

if [[ ${conv_probs} == "True" ]]; then

	for i in {1..15}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_conv_probs_set_${i}.txt
	done
else
	for i in {1..11}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-no-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_no_conv_probs_set_${i}.txt
	done

fi 


##high_poly
species_bin=high_poly

if [[ ${conv_probs} == "True" ]]; then

	for i in {1..16}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_conv_probs_set_${i}.txt
	done
else
	for i in {1..10}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-no-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_no_conv_probs_set_${i}.txt
	done

fi


#social_poly
species_bin=social_poly

if [[ ${conv_probs} == "True" ]]; then

	for i in {1..14}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_conv_probs_set_${i}.txt
	done
else
	for i in {1..12}; do
		sbatch -o ${out_DIR}log/relax-${gene_set}-${species_bin}-gs${i}-no-conv-probs.log ${base_DIR}relax_per_gene.sh ${gene_set} ${species_bin} ${gene_set_DIR}${gene_set}_gene_gs_10_${species_bin}_no_conv_probs_set_${i}.txt
	done

fi
