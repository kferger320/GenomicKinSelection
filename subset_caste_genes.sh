WD=/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/

##worker genes
awk -F"," -v p="Worker" 'tolower($11) == tolower(p)' ${WD}Warneretal2017_all_genes.csv > ${WD}Warneretal2017_worker_genes.csv
awk -F"," '{print $1}' ${WD}Warneretal2017_worker_genes.csv > ${WD}Warneretal2017_worker_gene_IDs.txt

##queen genes
awk -F"," -v p="Reproductive" 'tolower($11) == tolower(p)' ${WD}Warneretal2017_all_genes.csv > ${WD}Warneretal2017_queen_genes.csv
awk -F"," '{print $1}' ${WD}Warneretal2017_queen_genes.csv > ${WD}Warneretal2017_queen_gene_IDs.txt

##NDE genes:
#'NDE' in all life stages, BSnIPRE.class == 'neutral', and PS1 == 'cellular'
## && $26 == PS1_value
##        PS1_value = "cellular"

awk -F"," '
    BEGIN {
        # Define the specific values you want to filter for each column
        L2_value = "NDE"
        L3_value = "NDE"
        L4_value = "NDE"
        L5_value = "NDE"
        Head_value = "NDE"
        Gaster_value = "NDE"
        Larval_value = "NDE"
        Overall_value = "NDE"
        BSnIPRE_class_value = "neut"
    }

    # Skip the header line
    NR == 1 {
        print $0
        next
    }

    # Filter rows based on the specified values
    $4 == L2_value && $5 == L3_value && $6 == L4_value && $7 == L5_value && \
    $8 == Head_value && $9 == Gaster_value && $10 == Larval_value && $11 == Overall_value && \
    $18 == BSnIPRE_class_value {
        print $0
    }
' ${WD}Warneretal2017_all_genes.csv > ${WD}Warneretal2017_nde_genes.csv
awk -F"," '{print $1}' ${WD}Warneretal2017_nde_genes.csv > ${WD}Warneretal2017_nde_gene_IDs.txt