#!/bin/bash
#SBATCH --job-name=augustus
#SBATCH --account=co_rosalind
#SBATCH --partition=savio3
#SBATCH --qos=savio_lowprio
#SBATCH --time=36:00:00
#SBATCH --open-mode=append
#SBATCH -o augustus-log/slurm-%j.out
#SBATCH --requeue

module load augustus/3.3 
reference=$1
ref_basename_nopath="${reference##*/}"
ref_basename="${ref_basename_nopath%.*}"

WD="/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/"
mkdir ${WD}ref_genomes/${ref_basename}/

#temp param increases sensitivity- more exons sampled 
	# --codingseq=on \
	# --minmeanexonintronprob=0.1 \


# augustus --species=fly \
# 	--introns=on \
# 	--start=on \
# 	--stop=on \
# 	--exonnames=on \
# 	--sample=100 \
# 	--alternatives-from-sampling=true \
# 	--temperature=5 \
# 	--progress=true \
# 	${WD}ref_genomes/GCF_002006095.1_ASM200609v1_genomic.fna > \
# 	${WD}ref_genomes/GCF_002006095.1_ASM200609v1_predic_genes_aug.gff

##Step 1. Run AUGUSTUS with threading using the 'run_augustus_parallel' command
run_augustus_parallel -f ${reference} \
	-j 32 \
	-c 32 \
	-s ${ref_basename} \
	-p '--species=fly,--introns=on,--start=on,--stop=on,--exonnames=on,--sample=100,--alternatives-from-sampling=true,--temperature=5,--progress=true' \
	-o ${WD}ref_genomes/${ref_basename}/