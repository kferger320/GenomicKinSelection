#!/bin/bash
#SBATCH --job-name=msa
#SBATCH --account=fc_tsutsuifca
#SBATCH --partition=savio3
#SBATCH --qos=savio_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH --open-mode=append
#SBATCH -o msa-log/slurm-%j.out
#SBATCH --requeue


WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/
hyphy_analyses_bin=/global/scratch/users/kailey_ferger/bin/hyphy-analyses/codon-msa/
gene_set=nde #queen, worker, nde
MSA_DIR=MSA_${gene_set} #or MSA
results=consensus_Mphar_prot_set.txt

source activate hyphy-redo

##appended in loop
rm ${WD}${MSA_DIR}/proteins_all_taxa/aligned/seqcounts.txt

while ID=; IFS=$' \t\r\n' read -r ID || [[ $ID ]]; do

	# #remove later, added bc run ended too early 
	# if ! [ -f ${WD}${MSA_DIR}/proteins_all_taxa/aligned/${ID}.nuc.final.fas ]; then

	##Step 1. Run pre-msa.bf on concatenated prot files to prep for msa
		hyphy ${hyphy_analyses_bin}pre-msa.bf --input ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.fna --E 0.0005

	##Step 2. Perform multiple sequence alignment on each protein output generated above 
		mafft --thread 32 ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.fna_protein.fas > ${WD}${MSA_DIR}/proteins_all_taxa/aligned/${ID}.fna_protein.mafft.fas
		# mafft ${WD}proteins_all_taxa/${ID}.faa > ${WD}proteins_all_taxa/aligned/${ID}.mafft.faa

	##Step 3. obtain a final nucleotide msa using frameshift corrected nucleotide sequences from step 1
		hyphy ${hyphy_analyses_bin}post-msa.bf ${WD}${MSA_DIR}/proteins_all_taxa/aligned/${ID}.fna_protein.mafft.fas \
			--nucleotide-sequences ${WD}${MSA_DIR}/proteins_all_taxa/${ID}.fna_nuc.fas \
			--output ${WD}${MSA_DIR}/proteins_all_taxa/aligned/${ID}.nuc.final.fas

##Step 4. count total number of seqs per file and output
	#appended in loop
	rm ${WD}${MSA_DIR}/proteins_all_taxa/aligned/${ID}.seqheaders.txt

	count=$(grep ">" ${WD}${MSA_DIR}/proteins_all_taxa/aligned/${ID}.nuc.final.fas | wc -l)
	##output seq counts
	printf "${ID}.nuc.final.fas\t${count}\n" >> ${WD}${MSA_DIR}/proteins_all_taxa/aligned/seqcounts.txt
	##output sequence headers
	grep ">" ${WD}${MSA_DIR}/proteins_all_taxa/aligned/${ID}.nuc.final.fas >> ${WD}${MSA_DIR}/proteins_all_taxa/aligned/${ID}.seqheaders.txt
	
	# fi

done<${WD}${MSA_DIR}/${results}