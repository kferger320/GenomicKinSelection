#!/usr/bin/env python
# coding: utf-8

import json, os, sys
import pandas as pd
from collections import defaultdict
import numpy as np
from scipy.stats import false_discovery_control

gene_set = sys.argv[1] #"worker"
gene_bin = sys.argv[2] #"high_poly_unicolonial"
#optional
gen_desc=False #default process hypothesis tests only; change to "True" to process general descriptive model results

if gene_set == "worker":
    MSA_dir = "MSA"
    
else:
    MSA_dir = f"MSA_{gene_set}"


base_dir="/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/"
gene_dir = f"{base_dir}{MSA_dir}/proteins_all_taxa/aligned/nuc_aligned/"    
WD = f"{gene_dir}/relax-srv/"
# WD = f"{gene_dir}relax-srv/redo-12-29-24/"

all_genes = pd.read_csv(f"{gene_dir}{gene_set}_genes_all.txt", sep="\t", header=None)[0].tolist()
check_genes = {}
for gene in all_genes:
    check_genes[gene] = False #initialize

# ref_list = pd.read_csv(f"{base_dir}{gene_bin}_GC_refs.txt", sep="\t", header=None)[0].tolist()
ref_list = pd.read_csv(f"{base_dir}{gene_bin}_GC_refs_redo.txt", sep="\t", header=None)[0].tolist()
# outlier_list = pd.read_csv(f"{base_dir}outlier_GC_refs.txt", sep="\t", header=None)[0].tolist()
outlier_list = pd.read_csv(f"{base_dir}outlier_GC_refs_redo.txt", sep="\t", header=None)[0].tolist()
# monogyne_list = pd.read_csv(f"{base_dir}monogyne_GC_refs.txt", sep="\t", header=None)[0].tolist()
monogyne_list = pd.read_csv(f"{base_dir}monogyne_GC_refs_redo.txt", sep="\t", header=None)[0].tolist()


def get_nested_value(data, key_path):
    keys = key_path.split('.')
    value = data
    try:
        for key in keys:
            value = value[key]
        return value
    except KeyError:
        return None


def is_sig(df):
    df['significant'] = np.where(df['adj-p-value'] < 0.05, "yes", "no")
    return df


def main():

    final_df = defaultdict(list)
    ##general descriptive model
    final_df_gen_desc = defaultdict(list)

    for filename in os.listdir(WD):
        # Construct full file path
        file_path = os.path.join(WD, filename)
        # print(file_path)
        file_base = os.path.basename(file_path)
        # print(file_base)
        
        # Check if it's a file (not a directory)
        if os.path.isfile(file_path):
            # Check if the file size is 0 and file is correct bin
            if (os.path.getsize(file_path) > 0) and (gene_bin in file_base):
                # (gene_bin in file_base) and \

            # if (os.path.getsize(file_path) > 0) and \
            #     # ("gen-model-rerun" in file_base):
            #     ("conv-probs" in file_base): ##conv-probs runs
                # print(file_path)
                
                ##run scraping
                with open(file_path, 'r') as file:
                    data = json.load(file)

                    ##record gene/protein
                    gene = file_base.split("-")[2]

                    if not check_genes[gene]: #don't store duplicates
                        final_df["gene"].append(gene)
                        final_df_gen_desc["gene"].append(gene)
                        check_genes[gene] = True
                    else:
                        print(f"duplicate gene: {gene}")
                        continue

                    ##first check for convergence problems
                    conv_path = "analysis.settings"
                    conv_values = get_nested_value(data, conv_path)

                    if len(conv_values) > 0:
    #                     print(conv_values)
                        final_df["convergence_probs"].append("yes")   
                        final_df_gen_desc["convergence_probs"].append("yes")
                        # print(f"{gene}\tconv probs")   
                    else:
                        final_df["convergence_probs"].append("no")
                        final_df_gen_desc["convergence_probs"].append("no")  
                        # print(f"{gene}\tno conv probs")

                    ##dN/dS fits
                    dnds_path = "fits.MG94xREV with separate rates for branch sets.Rate Distributions"
                    dnds_value = get_nested_value(data, dnds_path)
                    final_df["dNdS_ref"].append(dnds_value['non-synonymous/synonymous rate ratio for *Reference*'][0][0])
                    final_df["dNdS_test"].append(dnds_value['non-synonymous/synonymous rate ratio for *Test*'][0][0])
                    final_df["dNdS_unclass"].append(dnds_value['non-synonymous/synonymous rate ratio for *Unclassified*'][0][0])

                    ##hypothesis test
                    ref_omegas = {} #omega_1: [0.1, 1, 0.1]
                    test_omegas = {}

                    hyp_path = "fits.RELAX alternative.Rate Distributions"
                    hyp_values = get_nested_value(data, hyp_path)

                    for i in range(3): #3 omega classes

                        ##reference values
                        ref_omega_path = f"Reference.{i}.omega"
                        ref_omega_values = get_nested_value(hyp_values, ref_omega_path)
                        ref_omegas[f"omega_{i+1}"] = ref_omega_values
                        final_df[f"ref_omega_{i+1}"].append(ref_omega_values)

                        ##test values
                        test_omega_path = f"Test.{i}.omega"
                        test_omega_values = get_nested_value(hyp_values, test_omega_path)
                        test_omegas[f"omega_{i+1}"] = test_omega_values
                        final_df[f"test_omega_{i+1}"].append(test_omega_values)


                    ##test results
                    test_results = {}
                    test_res_path = "test results"
                    test_res_res = get_nested_value(data, test_res_path)

                    test_results["LRT"] = test_res_res["LRT"]
                    test_results["p-value"] = test_res_res["p-value"]
                    test_results["K"] = test_res_res["relaxation or intensification parameter"]

                    final_df["LRT"].append(test_results["LRT"])
                    final_df["p-value"].append(test_results["p-value"])
                    final_df["K"].append(test_results["K"])
                    
                    
                    if gen_desc: #run only if general descriptive model was run for current set of genes
                        ##General descriptive model 
                        #initialize final_df_gen_desc to ensure consistency with labels
                        
                        for lab in ref_list:
                            final_df_gen_desc[f"test_{lab}"].append(np.nan)

                        for lab in monogyne_list:
                            final_df_gen_desc[f"ref_{lab}"].append(np.nan)
                        
                        ##get labels for all other possible branches (other bins, outliers)
                        ##outlier branches
                        for lab in outlier_list:
                            final_df_gen_desc[f"unclass_{lab}"].append(np.nan)

                        if gene_bin != "polygyne":
                            # polygyne_list = pd.read_csv(f"{base_dir}polygyne_GC_refs.txt", sep="\t", header=None)[0].tolist()
                            polygyne_list = pd.read_csv(f"{base_dir}polygyne_GC_refs_redo.txt", sep="\t", header=None)[0].tolist()
                            for lab in polygyne_list:
                                final_df_gen_desc[f"unclass_{lab}"].append(np.nan)

                        if gene_bin != "social_poly":
                            # social_poly_list = pd.read_csv(f"{base_dir}social_poly_GC_refs.txt", sep="\t", header=None)[0].tolist()
                            social_poly_list = pd.read_csv(f"{base_dir}social_poly_GC_refs_redo.txt", sep="\t", header=None)[0].tolist()
                            for lab in social_poly_list:
                                final_df_gen_desc[f"unclass_{lab}"].append(np.nan)

                        if gene_bin != "high_poly_unicolonial":
                            # high_poly_unicolonial_list = pd.read_csv(f"{base_dir}high_poly_unicolonial_GC_refs.txt", sep="\t", header=None)[0].tolist()
                            high_poly_unicolonial_list = pd.read_csv(f"{base_dir}high_poly_GC_refs_redo.txt", sep="\t", header=None)[0].tolist()
                            for lab in high_poly_unicolonial_list:
                                final_df_gen_desc[f"unclass_{lab}"].append(np.nan)
                            
                        # print(final_df_gen_desc)
                        gen_desc_path = "branch attributes.0"
                        target_value = get_nested_value(data, gen_desc_path)

                        ##get k values per branch for test and ref branches
                        gen_desc_kvals_test = {}
                        gen_desc_kvals_ref = {}

                        for key, value in target_value.items():
                            GC_ref = "_".join(key.split("_")[0:2])

                            #if branch is test branch
                            if GC_ref in ref_list:
                                gen_desc_kvals_test[GC_ref] = value['k (general descriptive)']
                                #replace np.nan with real value if found, leave if not found
                                final_df_gen_desc[f"test_{GC_ref}"][-1] = gen_desc_kvals_test[GC_ref]

                            elif GC_ref in monogyne_list:
                                gen_desc_kvals_ref[GC_ref] = value['k (general descriptive)']
                                final_df_gen_desc[f"ref_{GC_ref}"][-1] = gen_desc_kvals_ref[GC_ref]
                            
                            elif "Node" not in GC_ref: #unclassified branches
                                # print(GC_ref)
                                gen_desc_kvals_ref[GC_ref] = value['k (general descriptive)']
                                final_df_gen_desc[f"unclass_{GC_ref}"][-1] = gen_desc_kvals_ref[GC_ref]
                            else:
                                continue
            else:
                continue

    ##write basic stats output
    df = pd.DataFrame.from_dict(final_df, orient="columns")

    ##Benjamini Hochberg p-value adjustment 
    df["adj-p-value"] = false_discovery_control(df["p-value"])

    df = is_sig(df)
    df.to_csv(f"{base_dir}RELAX-srv_{gene_set}_{gene_bin}_results_allgenes.txt", sep="\t", index=None)
    # df.to_csv(f"{base_dir}RELAX-srv_{gene_set}_{gene_bin}_results_redo-12-29-24_allgenes.txt", sep="\t", index=None)

    if gen_desc:
        #write general descriptive model output
        gen_desc_df = pd.DataFrame.from_dict(final_df_gen_desc, orient="columns")
        conv_col = gen_desc_df.pop("convergence_probs")
        gen_desc_df.insert(len(gen_desc_df.columns), 'convergence_probs', conv_col)
        gen_desc_df.to_csv(f"{base_dir}RELAX-srv_{gene_set}_{gene_bin}-rerun-gen_desc-results-allgenes.txt", sep="\t", index=None)


if __name__=="__main__":
    main()