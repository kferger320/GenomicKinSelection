##adapted from https://www.biostars.org/p/55851/

import sys, os, re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# WD="/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/worker_gene_seqs/"


input_file = sys.argv[1]  #{WD}GCF_001594115.1_Tsep1.0_genomic_longest_iso_cds.anno.worker_genes.fna
WD = f"{os.path.dirname(input_file)}/" # eg. /global/scratch/users/kailey_ferger/genomic_sig_kin_selection/worker_gene_seqs/
basename = os.path.splitext(os.path.basename(input_file))[0]
output_nt_file = f"{WD}pre-msa/{basename}.nuc_prot.fna"
output_aa_file = f"{WD}pre-msa/{basename}.nuc_prot.faa"

gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
      'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

basepairs = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

def translate_frameshifted( sequence ):
    translate = ''.join([gencode.get(sequence[3*i:3*i+3],'X') for i in range(len(sequence)//3)])
    translated_nt = [sequence[3*i:3*i+3] for i in range(len(sequence)//3)]
    return translate, ''.join(translated_nt)

def reverse_complement( sequence ):
    reversed_sequence = (sequence[::-1])
    rc = ''.join([basepairs.get(reversed_sequence[i], 'X') for i in range(len(sequence))])
    return rc


gene_name=""
gene_desc=""
modified_prot_records = []
modified_nt_records = []

concat_protein_records = []
concat_nt_records = []


for record in SeqIO.parse(input_file, "fasta"):
    frame = int(record.description.split('-')[2][0])
#     print(frame)
    
    if (gene_name != record.id) & (len(modified_prot_records) > 0):
        # print(f"Gene name {gene_name} != record.id {record.id}")
        #indicates switching to a new gene
        all_prot_seqs = ''.join([str(i.seq) for i in modified_prot_records])
        all_nt_seqs = ''.join([str(i.seq) for i in modified_nt_records])
        
        ##outputting concatenated information from previous gene
        strand_pattern = r'\([+-]\)'
        strand = re.search(strand_pattern, gene_desc).group()
        concat_protein_records.append(SeqRecord(Seq(all_prot_seqs), id=gene_name, description=strand))
        concat_nt_records.append(SeqRecord(Seq(all_nt_seqs), id=gene_name, description=strand))
        
#         print(all_prot_seqs)
#         print(all_nt_seqs)
        
        #initialize new storage of sequences/output new sequence
        modified_prot_records = []
        modified_nt_records = []
        #process all same seqs together and output as one
        gene_name = record.id
        gene_desc = record.description
    
    if len(modified_prot_records) == 0:
        #initialize
        gene_name = record.id
        gene_desc = record.description
    
    
    #translate the input nt string at the appropriate frame start pos
    protein_sequence, nt = translate_frameshifted(str(record.seq[frame:]))
#     print(protein_sequence)

    # Create a new SeqRecord with the modified sequence and the same header
    new_prot_record = SeqRecord(Seq(protein_sequence), id=record.id, description=record.description)
    new_nt_record = SeqRecord(Seq(nt), id=record.id, description=record.description)

    modified_prot_records.append(new_prot_record)
    modified_nt_records.append(new_nt_record)

#dealing with last entry
all_prot_seqs = ''.join([str(i.seq) for i in modified_prot_records])
all_nt_seqs = ''.join([str(i.seq) for i in modified_nt_records])

strand_pattern = r'\([+-]\)'
strand = re.search(strand_pattern, record.description).group()
concat_protein_records.append(SeqRecord(Seq(all_prot_seqs), id=record.id, description=strand))
concat_nt_records.append(SeqRecord(Seq(all_nt_seqs), id=record.id, description=strand))

    
#     print(modified_prot_records)
#     print(modified_nt_records)

# # Write the modified sequences to an output FASTA file
# #aa and nt output fastas
SeqIO.write(concat_nt_records[-1], "/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/test.fasta", "fasta")
# SeqIO.write(concat_protein_records, output_aa_file, "fasta")
# SeqIO.write(concat_nt_records, output_nt_file, "fasta")