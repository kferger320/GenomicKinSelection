#!/bin/bash
#SBATCH --job-name=fetch-genes
#SBATCH --account=fc_tsutsuifca
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --open-mode=append
#SBATCH --requeue

WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/
gene_group=$1 #queen, worker, nde


##Step 1. Extract protein IDs/cleaning steps using gene2accession.tsv 

##remove old files appended in a loop
rm ${WD}Warneretal2017_${gene_group}_gene_IDs.strip.txt
rm ${WD}ordered_Warneretal2017_${gene_group}_protein_IDs.txt

##Creating new geneID file with 'LOC' stripped (numbers only)
while read ID; do
	stripped_ID=$(echo "$ID" | sed 's/LOC//')
	printf "${stripped_ID}\n" >> ${WD}Warneretal2017_${gene_group}_gene_IDs.strip.txt
done<${WD}Warneretal2017_${gene_group}_gene_IDs.txt

sort -t '	' -k1,1 ${WD}Warneretal2017_${gene_group}_gene_IDs.strip.txt | uniq > ${WD}ordered_Warneretal2017_${gene_group}_gene_IDs.strip.txt

#Finding the associated proteinID from each geneID in gene2accession.tsv
while read p; do
	grep "$p" ${WD}gene2accession.tsv | cut -d$'\t' -f2 >> ${WD}ordered_Warneretal2017_${gene_group}_protein_IDs.txt
done<${WD}ordered_Warneretal2017_${gene_group}_gene_IDs.strip.txt

##cleaning file (removing '-' characters from presumably no match found)
grep -v '-' ${WD}ordered_Warneretal2017_${gene_group}_protein_IDs.txt > ${WD}ordered_Warneretal2017_${gene_group}_protein_IDs.cleaned.txt


## Step 2. fetch protein sequences using Entrez tools and deposit into singular .faa file: using webtool:
## http://www.ncbi.nlm.nih.gov/sites/batchentrez with ordered_Warneretal2017_${gene_group}_protein_IDs.cleaned.txt as input


##Optional step: fetch ${gene_group} gene sequences
while read ID; do
	echo ${ID}
	stripped_ID=$(echo "$ID" | sed 's/LOC//')
	esummary -db gene -id ${stripped_ID} \
	| xtract -pattern DocumentSummary -if GenomicInfoType -element Id \
		-block GenomicInfoType -element ChrAccVer ChrStart ChrStop \
	| while read -r gene_id chr_acc chr_start chr_stop ; do 
		efetch -db nuccore -id $chr_acc -chr_start $chr_start -chr_stop $chr_stop -format fasta ; 
	done >> ${WD}Mpharaonis_${gene_group}_genes.fasta
done<${WD}Warneretal2017_${gene_group}_gene_IDs.txt


##Step 3. Extract longest isoforms from protein fasta file
cd-hit -i ${WD}Mpharaonis_${gene_group}_proteins.fasta -o ${WD}Mpharaonis_${gene_group}_proteins_longest_iso.fasta -c 0.9