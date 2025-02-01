#!/bin/bash
#SBATCH --job-name=blast
#SBATCH --account=fc_tsutsuifca
#SBATCH --partition=savio2_htc
#SBATCH --qos=savio_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --open-mode=append
#SBATCH -o blast-log/slurm-%j.out
#SBATCH --requeue

WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/
ref=$1
ref_basename_nopath="${ref##*/}"
ref_basename="${ref_basename_nopath%.*}"
gene_group=nde #worker, queen, nde

##pre-step: make blastdb
makeblastdb -in ${ref} -dbtype nucl

##Step 1. Perform blast command against queen protein fasta using tblastn

#outputting in accordance with BED6 format: chr, start, stop, name, score, strand
tblastn -db ${ref} -evalue 1e-5 -max_hsps 1 -num_alignments 1 \
	-outfmt "6 sseqid sstart send qseqid evalue sstrand" \
	-query ${WD}Mpharaonis_${gene_group}_proteins_longest_iso.fasta \
	-out ${WD}blast_${gene_group}/${ref_basename}.blast


##Step 2. remove query isoforms and fix formatting issues 
# create unique file to eliminate query isoforms
awk '!a[$1 $2 $3]++ { print ;}' ${WD}blast_${gene_group}/${ref_basename}.blast > ${WD}blast_${gene_group}/${ref_basename}_longest_iso.uniq.blast

##fix formatting issues (change 'plus' into '+', 'minus' into '-', and switch 'start' and 'end' if on '-' strand (blast reverses these for - strand))
awk 'BEGIN {OFS="\t"} {if ($NF == "plus") {$NF = "+"} else {$NF = "-"; t=$2; $2=$3; $3=t}}1' ${WD}blast_${gene_group}/${ref_basename}_longest_iso.uniq.blast > ${ref_basename}.tmp
mv ${ref_basename}.tmp ${WD}blast_${gene_group}/${ref_basename}_longest_iso.uniq.blast


# #extract only the gene names
# awk -F"\t" '{print $2}' ${WD}blast/${ref_basename}_longest_iso.uniq.blast | sort | uniq > ${WD}blast/${ref_basename}.queen_match_geneIDs_longest_iso.uniq.txt
# awk -F"\t" '{print $1}' ${WD}blast/${ref_basename}_longest_iso.uniq.blast | sort | uniq > ${WD}blast/${ref_basename}.queen_match_geneIDs_longest_iso.uniq.txt

## extract queen gene sequences from full gene fasta
# seqtk subseq ${WD}ref_genomes/${ref_basename}.aug3.genes_full.fna ${WD}blast/${ref_basename}.queen_match_geneIDs_longest_iso.uniq.txt > ${WD}queen_gene_seqs/${ref_basename}.aug3.queen_genes.fna


