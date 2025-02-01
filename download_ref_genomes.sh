#!/bin/bash
#SBATCH --job-name=download-refs
#SBATCH --account=co_rosalind
#SBATCH --partition=savio2_htc
#SBATCH --qos=rosalind_htc2_normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --open-mode=append
#SBATCH -o download-log/slurm-%j.out
#SBATCH --requeue

WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/ref_genomes/

##Step 1. download each reference genome using the 'datasets download' command (NCBI Datasets command-line tools (CLI)); unzip and move to single folder
while read acc; do 

	acc_basename_nopath="${acc##*/}"
	acc_basename="${acc_basename_nopath%.*}"
	acc_ID=$(echo ${acc_basename} |awk -F_ '{print $1"_"$2}')

	if [ ! -d ${WD}${acc_basename} ]; then
		mkdir -p ${WD}${acc_basename};
	fi

	#fetch
	#datasets download genome accession ${acc} --include genome --reference --filename ${WD}${acc}.zip
	# datasets download genome accession ${acc} --include protein,gff3 --reference --filename ${WD}${acc}/${acc}.zip
	datasets download genome accession ${acc_ID} --include genome,gff3 --reference --filename ${WD}${acc_basename}/${acc_basename}.zip


	#unzip into ref_genomes
	echo "N" | unzip -d ${WD}${acc_basename} -j ${WD}${acc_basename}/${acc_basename}.zip

	rm ${WD}${acc_basename}/${acc_basename}.zip
	rm ${WD}${acc_basename}/README.md

	#rename to include organism name
	# mv ${WD}${acc}/protein.faa ${WD}${acc}/${org}.faa
	mv ${WD}${acc_basename}/genomic.gff ${WD}${acc_basename}/${acc_basename}.gff


	#move genome and gff files back to WD
	mv ${WD}${acc_basename}/* ${WD}
	
	sleep 5

done</global/scratch/users/kailey_ferger/genomic_sig_kin_selection/GC_accessions.txt