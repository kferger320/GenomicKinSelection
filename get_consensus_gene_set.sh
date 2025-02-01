#!/bin/bash
#SBATCH --job-name=get_con_genes
#SBATCH --account=fc_tsutsuifca
#SBATCH --partition=savio4_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=16:00:00
#SBATCH --open-mode=append
#SBATCH -o get_consensus_gene_set-log/slurm-%j.out
#SBATCH --requeue

WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/
ref=$1
annotated=$2 # "annotated" or "unannotated"
ref_basename_nopath="${ref##*/}"
ref_basename="${ref_basename_nopath%.*}"

module load bedtools

MSA_DIR=MSA_nde #or MSA
gene_set=nde #queen or worker
results=consensus_Mphar_prot_set.txt


##Step 1. Get matched gene names with corresponding M/pharaonis protein names, output matches to new file
# blast_file=${WD}blast_queen/${ref_basename}_longest_iso.uniq.blast 
blast_file=${WD}blast_nde/${ref_basename}_longest_iso.uniq.blast 

if [ "$annotated" = "annotated" ]; then
	gff_file=${WD}ref_genomes/${ref_basename}.gff
else
	gff_file=${WD}ref_genomes/${ref_basename}/${ref_basename}.fna.aug.gff3
fi

python ${WD}scripts/get_consensus_set.py ${blast_file} ${gff_file} ${annotated} ${gene_set}


##Step 3. Extract all nt sequences per ref and per consensus M.phar protein
rm ${WD}prot_not_mapped_to_gene.txt

while ID=; IFS=$' \t\r\n' read -r ID || [[ $ID ]]; do

	###remove old versions of prot_ID files since we're appending in the loop
	# stripped_name=$(echo $ref_basename | sed 's/genomic//')
	# echo $stripped_name
	if [ -f ${WD}${MSA_DIR}/${ref_basename}_${ID}.tmp ]; then
    	rm -f ${WD}${MSA_DIR}/${ref_basename}_${ID}.tmp
	fi

	if [ -f ${WD}${MSA_DIR}/${ref_basename}_${ID}.faa.tmp ]; then
    	rm -f ${WD}${MSA_DIR}/${ref_basename}_${ID}.faa.tmp
	fi

	##find associated gene for each M.phar prot ID from files generated in Step 1.
	matches=$(grep ${ID} ${WD}${gene_set}_gene_seqs/${ref_basename}_Mphar_prot_matches.txt | awk -F"\t" '{print $2}')

	if [[ -z "${matches// }" ]]; then
		echo ${ID} >> ${WD}${MSA_DIR}/prot_not_mapped_to_gene.txt
		continue
	fi
	# echo ${ID} ${res}

	# #remove tmp bedfile
	# rm ${WD}${MSA_DIR}/${ref_basename}.bed

	##Extract matching nt gene seqs to gene file in MSA/
	for res in ${matches}; do
		# search_string='${res}_'
		seqkit grep -r -n -p "${res} " ${WD}${gene_set}_gene_seqs/pre-msa/${ref_basename}_longest_iso_cds*.${gene_set}_genes.nuc_prot.fna | \
			seqkit replace -p "(.*)" -r '$1 | '"$ref_basename" >> ${WD}${MSA_DIR}/${ref_basename}_${ID}.tmp
	done

	##same thing for protein seqs
	for res in ${matches}; do
		seqkit grep -r -n -p "${res} " ${WD}${gene_set}_gene_seqs/pre-msa/${ref_basename}_longest_iso_cds*.${gene_set}_genes.nuc_prot.faa | \
			seqkit replace -p "(.*)" -r '$1 | '"$ref_basename" >> ${WD}${MSA_DIR}/${ref_basename}_${ID}.faa.tmp
	done

	echo ${WD}${MSA_DIR}/${ref_basename}_${ID}.tmp

	##Step 4. Extract the longest isoform for each Mphar protein and each reference- outputs have only 1 seq per protein per ref
	python ${WD}scripts/longest.py -i ${WD}${MSA_DIR}/${ref_basename}_${ID}.tmp -n 1 > ${WD}${ref_basename}_${ID}.long_iso.tmp
	python ${WD}scripts/longest.py -i ${WD}${MSA_DIR}/${ref_basename}_${ID}.faa.tmp -n 1 > ${WD}${ref_basename}_${ID}.long_iso.faa.tmp

	# #rename
	mv ${WD}${ref_basename}_${ID}.long_iso.tmp ${WD}${MSA_DIR}/${ref_basename}_${ID}.tmp
	mv ${WD}${ref_basename}_${ID}.long_iso.faa.tmp ${WD}${MSA_DIR}/${ref_basename}_${ID}.faa.tmp

done<${WD}${MSA_DIR}/${results}




#######RUN THESE STEPS ONLY AFTER ALL JOBS FOR ALL REFERENCES HAVE COMPLETED FOR THE COMMANDS ABOVE#####
#Run only once; not per-reference

# ##Step 5. Concatenate all nt and aa seqs for each protein to create 1 file per M.phar protein

# # ##Remove old versions of files
# while ID=; IS=$' \t\r\n' read -r ID || [[ $ID ]]; do
# 	# for file in ${WD}${MSA_DIR}/*${ID}.tmp; do
# 		# echo ${file}

# 	#remove old versions since appending in a loop
# 	if [ -f ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.fna ]; then
# 		rm -f ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.fna
# 	fi

	# done

# done<${WD}${MSA_DIR}/${results}

# while ID=; IFS=$' \t\r\n' read -r ID || [[ $ID ]]; do
# 	# for file in ${WD}${MSA_DIR}/*${ID}.faa.tmp; do
# 		# echo ${file}

# 	#remove old versions since appending in a loop
# 	if [ -f ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.faa ]; then
# 		rm -f ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.faa
# 	fi

# 	# done

# done<${WD}${MSA_DIR}/${results}



# ##concatenate all nt tmp files to one final file per protein ()
# while ID=; IFS=$' \t\r\n' read -r ID || [[ $ID ]]; do

# 	for file in ${WD}${MSA_DIR}/*${ID}.tmp; do
# 		# echo ${file}
# 		cat ${file} >> ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.fna

# 	done

# 	##remove duplicates found in files (cannot figure out how they got in there?)
# 	seqkit rmdup -s ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.fna > ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.tmp.fna
# 	mv ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.tmp.fna ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.fna


# done<${WD}${MSA_DIR}/${results}

# ##concatenate all protein (aa) tmp files to one final file per protein ()
# while ID=; IFS=$' \t\r\n' read -r ID || [[ $ID ]]; do

# 	for file in ${WD}${MSA_DIR}/*${ID}.faa.tmp; do
# 		cat ${file} >> ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.faa

# 	done

# 	##remove duplicates found in files (cannot figure out how they got in there?)
# 	seqkit rmdup -s ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.faa > ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.tmp.faa
# 	mv ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.tmp.faa ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.faa

# done<${WD}${MSA_DIR}/${results}