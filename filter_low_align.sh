MSA_DIR=MSA_nde #or MSA
threshold=62

WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/${MSA_DIR}/proteins_all_taxa/aligned/

##Step 1. Find protein .nuc_final.fas files with <62 sequences aligned (coudl not align in pre-msa.bf step)
awk '$2 < 62 {print $1}' ${WD}seqcounts.txt > ${WD}low_align_seqs.txt
awk '$2 >= 62 {print $1}' ${WD}seqcounts.txt > ${WD}final_prot_set.txt

##Step 2. Move corresponding files to low_align folder, don't include them in rest of analysis 
while file=; IFS=$' \t\r\n' read -r file || [[ $file ]]; do
	mv ${WD}${file} ${WD}low_align/

done<${WD}low_align_seqs.txt