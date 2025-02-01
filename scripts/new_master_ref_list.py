import sys, os, re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
from collections import defaultdict

MSA_dir = sys.argv[1] #MSA or MSA_queen or MSA_nde
WD=f"/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/{MSA_dir}/proteins_all_taxa/aligned/"

protein_master_file=f"{WD}final_prot_set.txt"
protein_master_list = pd.read_csv(protein_master_file, header=None)[0].to_list()

ref_master_dict = defaultdict(list) #ref_master_dict[GCA_12345] = [max_occurrences, seqlen]

def get_sorting_key(record):
	#sort records in fasta file to make sure they all align
	id_str = record.id
	# start_index = id_str.find('GC')
	# if start_index != -1:
	#     sorting_key = id_str[start_index:start_index+11]
	# else:
	#     sorting_key = 'zzz'  # Assign a high value for records without 'GC*'
	sorting_key = re.search(r'GC[^_]*_\d+', record.id).group(0)
	return sorting_key

def rearrange_string(input_string):
    pattern = r"(.*)_(GC[A|F])_(.*)"
    match = re.match(pattern, input_string)

    if match:
        return f"{match.group(2)}_{match.group(3)}_{match.group(1)}"
    return None

##keep track of all seq headers for dummy seqs later on
seqID_tracking = defaultdict(list)

for file in protein_master_list:
	ID=file.split(".nuc.final.fas")[0]
	##tracking of accession frequencies within-file
	within_file_ref_acc = {}

	for i, record in enumerate(sorted(SeqIO.parse(f"{WD}{file}", "fasta"), key=get_sorting_key)):
		# seqID_tracking[ID].append(record.id)

		#get ref accession. eg. GCA_001045655
		ref_acc = re.search(r'GC[^_]*_\d+', record.id).group(0)

		#tally GC[A/F]_... within each file
		if ref_acc not in within_file_ref_acc.keys():
			within_file_ref_acc[ref_acc] = 1
		else: #already found, replicate sequence
			within_file_ref_acc[ref_acc] += 1
			##append any duplicate names for later use 
			seqID_tracking[ref_acc].append(record.id)
			print(f"{ID}\tref_acc: {ref_acc}\trecord_id: {record.id}")


	for key, value in within_file_ref_acc.items():

		if key in ref_master_dict.keys(): #if GC[A/F]_... already in master list
			if value > ref_master_dict[key]: #if current ref has more gene duplicates than previously recorded
				ref_master_dict[key] = value #update new master value to reflect number of duplicates in current ref
		
		else: #GC[A/F]_... not yet recorded in master list
			ref_master_dict[key] = value
	

final_fastas = defaultdict(list)
##loop through query proteins again to add dummy seqs and change IDs
for file in protein_master_list:
	ID=file.split(".nuc.final.fas")[0]
	##tracking of accession frequencies within-file
	within_file_ref_acc = {}
	unique_GCs = [] #only adding series of potential duplicates once 

	for i, record in enumerate(sorted(SeqIO.parse(f"{WD}{file}", "fasta"), key=get_sorting_key)):
		seqlen = len(record.seq) #aligned sequence length, the same for all seqs in file

			#get ref accession. eg. GCA_001045655
		ref_acc = re.search(r'GC[^_]*_\d+', record.id).group(0)

		#tally GC[A/F]_... within each file
		if ref_acc not in within_file_ref_acc.keys():
			within_file_ref_acc[ref_acc] = 1
		else: #already found, replicate sequence
			within_file_ref_acc[ref_acc] += 1
		
		##append current record to final_fastas 
		#change name
		new_id = rearrange_string(record.id)
		final_fastas[ID].append(SeqRecord(record.seq, id=new_id, description=''))
	
	for key, value in ref_master_dict.items():
		#compare number of duplicates between current ref and master list
		if key in within_file_ref_acc.keys():
			if within_file_ref_acc[key] < value:
				# print(f"{key}: {within_file_ref_acc[key]}; master_ref_value:{value}")
				##how many short are we?
				short = value - within_file_ref_acc[key]
				possible_headers = [i for i in seqID_tracking[key] if key in i]

				for i in range(short):
					#Add dummy seqs to fill out duplicate seqs
					dummy_seq = ''.join(['-' for i in range(seqlen)])
					#add actual ref headers of missing duplicates from some master list, so that sorting will be consistent at end?
					dummy_id = possible_headers[i]
					final_dummy_id = rearrange_string(dummy_id)

					# dummy_id = f"{key}_dummy{i}" 
					final_fastas[ID].append(SeqRecord(Seq(dummy_seq), id=final_dummy_id, description=''))

		else: #reference accession not found in current file at all; insert dummy seq
			# # possible_headers = [i for i in seqID_tracking[key] if key in i]

			for i in range(ref_master_dict[key]):
				#loop through all possible master seqs
				#Add dummy seqs to fill out duplicate seqs
				dummy_seq = ''.join(['-' for i in range(seqlen)])
				#add actual ref headers of missing duplicates from some master list, so that sorting will be consistent at end?
				dummy_id = f"{key}_dummyseq{i}"
				final_fastas[ID].append(SeqRecord(Seq(dummy_seq), id=dummy_id, description=''))


##write the final files 

for ID, seq in final_fastas.items():

	output_file = f"{WD}nuc_aligned/{ID}.nuc.final.consensus.fas"
	SeqIO.write(seq, output_file, "fasta")
