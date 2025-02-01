#!/bin/bash
#SBATCH --job-name=extract-exons
#SBATCH --account=fc_tsutsuifca
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --open-mode=append
#SBATCH -o extract-cds-log/slurm-%j.out
#SBATCH --requeue

WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/
ref=$1
annotated=$2 ##"annotated" or "unannotated"
ref_basename_nopath="${ref##*/}"
ref_basename="${ref_basename_nopath%.*}"
gene_set=nde ##worker, queen, nde

# gff=${WD}ref_genomes/${ref_basename}.gff
blast=${WD}blast_${gene_set}/${ref_basename}_longest_iso.uniq.blast #worker

module load bedtools

##extract longest isoforms for each gene; then extract only the cds corresponding to that isoform
if [ "$annotated" = "annotated" ]; then

	gff=${WD}ref_genomes/${ref_basename}.gff

	##Step 1. Extract the longest isoforms for each gene from each gff file
	ID_tags=$(python ${WD}scripts/get_the_longest_transcripts.py ${gff} | \
		awk -F"\t" '{extracted=substr($9, index($9, "ID=rna-")+3, index($9, ";")-index($9, "ID=rna-")-3); gsub("^rna-", "", extracted); print extracted}') #extract longest isoform tag
	

	# #appending in a loop
	rm ${WD}ref_genomes/${ref_basename}_${gene_set}_longest_iso_cds.bed

	##Step 2. Extract the CDS's from each longest isoform
	while IFS= read -r tag; do

		##use awk to filter by mRNA rows with matching tags; extract gene name and save to var
		m_lines=$(grep -F "${tag};" ${gff} | sed 's/ /\\t/g')
		# echo "${m_lines}\n" >> ${WD}m_lines.txt
		gene=$(awk -F"\t" -v p="mRNA" 'BEGIN{OFS="\t";} $3 == p { match($9, /Parent=gene-([^;]+)/, m); if(m[1]) print m[1] }' <<< "$m_lines")
		# echo ${gene}

		##output sepearate bedfile containing whole gene ranges for blast intersect command later
		awk -F"\t" -v p="mRNA" -v gene="$gene" 'BEGIN{OFS="\t";} $3 == p { print $1, $4, $5, gene" "$4"-"$5"-"$8, $6, $7 }' <<< "$m_lines" >> ${WD}ref_genomes/${ref_basename}_longest_iso_transcript.bed
		
		#extract cds sequences
		awk -F"\t" -v p="CDS" -v gene="$gene" 'BEGIN{OFS="\t";} $3 == p { print $1, $4, $5, gene" "$4"-"$5"-"$8, $6, $7 }' <<< "$m_lines" >> ${WD}ref_genomes/${ref_basename}_longest_iso_cds.tmp

	done <<< "$ID_tags"
	

	# ##Step 3. Filter extracted transcripts to those only occurring in (blasted) genes, extract corresponding CDS seqs
	bedtools intersect -a ${WD}ref_genomes/${ref_basename}_longest_iso_transcript.bed \
		-b ${blast} -wa > ${WD}ref_genomes/${ref_basename}_${gene_set}_gene_lines.bed #| awk -F "\t" 'BEGIN{OFS="\t";} {print $7, $8, $9, $10, $11, $12}'
	
	##Create file for unique gene matches to extract CDS lines later
	genes=$(awk -F"\t" '{print $4}' ${WD}ref_genomes/${ref_basename}_${gene_set}_gene_lines.bed | awk -F" " '{print $1}' |sort | uniq)
	# echo ${gene_entry}

	##appending in a loop
	rm ${WD}ref_genomes/${ref_basename}_${gene_set}_longest_iso_cds.bed

	##extract corresponding CDS seqs to create final CDS bedfile
	for gene in ${genes}; do
		# echo ${gene}
		grep ${gene} ${WD}ref_genomes/${ref_basename}_longest_iso_cds.tmp >> ${WD}ref_genomes/${ref_basename}_${gene_set}_longest_iso_cds.bed
	done

	#remove temporary file
	rm ${WD}ref_genomes/${ref_basename}_longest_iso_cds.tmp

	
	##sort by gene name to ensure proper order for concatenation
	sort -k4 ${WD}ref_genomes/${ref_basename}_${gene_set}_longest_iso_cds.bed > ${ref_basename}.tmp
	mv ${ref_basename}.tmp ${WD}ref_genomes/${ref_basename}_${gene_set}_longest_iso_cds.bed

	##Step 4. Extract nt sequences of matching CDS's
	bedtools getfasta -fi ${ref} \
		-bed ${WD}ref_genomes/${ref_basename}_${gene_set}_longest_iso_cds.bed \
		-s -name -fo ${WD}${gene_set}_gene_seqs/${ref_basename}_longest_iso_cds.${gene_set}_genes.fna
	


else #unannotated files have different format

	gff=${WD}ref_genomes/${ref_basename}/${ref_basename}.fna.aug.gff3

	##Step 1. Extract the longest isoforms for each gene from each gff file
	ID_tags=$(python ${WD}scripts/get_the_longest_transcripts.py ${gff} | \
		awk -F"\t" '{extracted=substr($9, index($9, "ID=")+3, index($9, ";")-index($9, "ID=")-3); print extracted}') #extract longest isoform tag

	#appending in a loop
	rm ${WD}ref_genomes/unannotated/${ref_basename}_longest_iso_transcript.bed
	rm ${WD}ref_genomes/unannotated/${ref_basename}_longest_iso_cds.tmp

	##Step 2. Extract the transcripts and CDS's from each longest isoform
	while IFS= read -r tag; do

		##use awk to filter by transcript rows with matching tags; extract gene name and save to var
		m_lines=$(grep -F "${tag}" ${gff} | sed 's/ /\\t/g')
		# echo "${m_lines}\n"

		gene=$(awk -F"\t" -v p="transcript" 'BEGIN{OFS="\t";} $3 == p { match($9, /Parent=([^;]+)/, m); if(m[1]) print m[1] }' <<< "$m_lines")

		#output sepearate bedfile containing whole gene ranges for blast intersect command later
		awk -F"\t" -v p="transcript" -v gene="$tag" 'BEGIN{OFS="\t";} $3 == p { print $1, $4, $5, gene" "$4"-"$5"-"$8, $6, $7 }' <<< "$m_lines" >> ${WD}ref_genomes/unannotated/${ref_basename}_longest_iso_transcript.bed

		#extract cds sequences
		awk -F"\t" -v p="CDS" -v gene="$tag" 'BEGIN{OFS="\t";} $3 == p { print $1, $4, $5, gene" "$4"-"$5"-"$8, $6, $7 }' <<< "$m_lines" >> ${WD}ref_genomes/unannotated/${ref_basename}_longest_iso_cds.tmp


	done <<< "$ID_tags"

	# ##Step 3. Filter extracted transcripts to those only occurring in (blasted) genes, extract corresponding CDS seqs
	bedtools intersect -a ${WD}ref_genomes/unannotated/${ref_basename}_longest_iso_transcript.bed \
		-b ${blast} -wa > ${WD}ref_genomes/unannotated/${ref_basename}_${gene_set}_gene_lines.bed #| awk -F "\t" 'BEGIN{OFS="\t";} {print $7, $8, $9, $10, $11, $12}'
	
	##Create file for unique gene matches to extract CDS lines later
	genes=$(awk -F"\t" '{print $4}' ${WD}ref_genomes/unannotated/${ref_basename}_${gene_set}_gene_lines.bed | awk -F" " '{print $1}' |sort | uniq)
	# echo ${gene_entry}

	##appending in a loop
	rm ${WD}ref_genomes/unannotated/${ref_basename}_${gene_set}_longest_iso_cds.bed

	##extract corresponding CDS seqs to create final CDS bedfile
	for gene in ${genes}; do
		# echo ${gene}
		grep ${gene} ${WD}ref_genomes/unannotated/${ref_basename}_longest_iso_cds.tmp >> ${WD}ref_genomes/unannotated/${ref_basename}_${gene_set}_longest_iso_cds.bed
	done

	#remove temporary file
	rm ${WD}ref_genomes/unannotated/${ref_basename}_longest_iso_cds.tmp

	
	##sort by gene name to ensure proper order for concatenation
	sort -k4 ${WD}ref_genomes/unannotated/${ref_basename}_${gene_set}_longest_iso_cds.bed > ${ref_basename}.tmp
	mv ${ref_basename}.tmp ${WD}ref_genomes/unannotated/${ref_basename}_${gene_set}_longest_iso_cds.bed

	##Step 4. Extract nt sequences of matching CDS's
	bedtools getfasta -fi ${ref} \
		-bed ${WD}ref_genomes/unannotated/${ref_basename}_${gene_set}_longest_iso_cds.bed \
		-s -name -fo ${WD}${gene_set}_gene_seqs/${ref_basename}_longest_iso_cds.aug3.${gene_set}_genes.fna

fi