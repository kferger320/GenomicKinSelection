gene_set=queen #queen, worker, nde
test_bin=social_poly #high_poly, polygyne, social_poly
MSA_DIR=MSA_${gene_set}
group_files=16 #this varies depeding on how many prot files per gene set there are/files per group
group_DIR=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/${MSA_DIR}/proteins_all_taxa/aligned/nuc_aligned/
WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/RAxML_sbatch_dir/${gene_set}_genes/


for ((i=1; i<=group_files; i++)); do
	##running RAxML command, 1 per gene set file
	sbatch ${WD}RAxML.sh ${group_DIR}${gene_set}_gene_set_${i}.txt
	
	##running branch tip labelling, 1 per gene set file
	python ${WD}scripts/label_raxml_branch_tips_base_genes.py ${test_bin} ${gene_set} ${group_DIR}${gene_set}_gene_set_${i}.txt

done
