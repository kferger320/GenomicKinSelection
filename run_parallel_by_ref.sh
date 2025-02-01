WD="/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/ref_genomes/unannotated/"
# WD="/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/ref_genomes/"

for file in ${WD}*genomic.fna; do 
	# echo ${file} >> all_refs.txt
	# ref_basename_nopath="${file##*/}"
	# ref_basename="${ref_basename_nopath%.*}"
	# mkdir /global/scratch/users/kailey_ferger/genomic_sig_kin_selection/ref_genomes/${ref_basename}
	# echo ${file}
	# sbatch augustus.sh ${file}
	# sbatch fetch_worker_genes.sh ${file}
	# sbatch blast_worker_genes.sh ${file}
	# sbatch blast_queen_genes.sh ${file}
	# sbatch extract_cds.sh ${file} unannotated
	# sbatch extract_cds.sh ${file} annotated
	# sbatch pre-msa-kf.sh ${file}
	sbatch get_consensus_gene_set.sh ${file} unannotated
	# sbatch get_consensus_gene_set.sh ${file} annotated

	# ref_basename_nopath="${file##*/}"
	# ref_basename="${ref_basename_nopath%.*}"
	# acc=$(echo "$ref_basename" | awk -F'_' '{print $1"_"$2}')
	# echo ${acc} >> org_unann_ref_accessions.tsv
	# echo ${acc} >> org_ann_ref_accessions.tsv

	# sbatch bedtools_get_anno_worker_genes_gff.sh ${file}
done