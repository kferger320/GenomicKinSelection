import sys, os

blast_file = sys.argv[1]
gff_file = sys.argv[2]
annotated = sys.argv[3] #annotated or unannotated (differnet file formats)
gene_set = sys.argv[4]

out_dir = f"/global/scratch/users/kailey_ferger/genomic_sig_kin_selection/{gene_set}_gene_seqs/"

basename = os.path.splitext(os.path.basename(gff_file))[0]
end = basename.find("genomic") + len("genomic")
ref_name = basename[:end]

def overlap(a1, a2, b1, b2):
	"""Check if two numeric ranges overlap."""
	return a1 <= b2 and b1 <= a2

def parse_range(s):
	"""Parse a numeric range string into a tuple."""
	return tuple(map(int, s.split('-')))

blast_data = {}
with open(blast_file, 'r') as blast:
	for line in blast:
		fields = line.strip().split('\t')
		blast_data[(fields[0], int(fields[1]), int(fields[2]))] = fields[3]


if annotated == "annotated": 
	with open(gff_file, 'r') as gff:
		outstr = ""

		for line in gff:
			fields = line.strip().split('\t')
			if len(fields) == 9: #skip non-entry lines
				if fields[2] == 'gene': #extract only 'transcript' lines
					ranges = [int(fields[3]), int(fields[4])]
					# range1 = parse_range(fields[3])
					# print(range1)
					# range2 = parse_range(fields[4])

					for key, value in blast_data.items():
						# print(key)
						##match scaffold name and check for range overlaps between gene and blast hit 
						if (key[0] == fields[0]) and overlap(ranges[0], ranges[1], key[1], key[2]):
							##extract gene name from last column
							start = fields[8].find("ID=gene-") + len("ID=gene-")
							end = fields[8].find(";", start)
							gene = fields[8][start:end]
							outstr += f"{value}\t{gene}\n"

	with open(f"{out_dir}{ref_name}_Mphar_prot_matches.txt", "w") as out:
		out.write(outstr)


else: #unannotated
	with open(gff_file, 'r') as gff:
		outstr = ""

		for line in gff:
			fields = line.strip().split('\t')
			if len(fields) == 9: #skip non-entry lines
				if fields[2] == 'transcript': #extract only 'transcript' lines
					ranges = [int(fields[3]), int(fields[4])]
					# range1 = parse_range(fields[3])
					# print(range1)
					# range2 = parse_range(fields[4])

					for key, value in blast_data.items():
						# print(key)
						##match scaffold name and check for range overlaps between gene and blast hit 
						if (key[0] == fields[0]) and overlap(ranges[0], ranges[1], key[1], key[2]):
							##extract gene name from last column
							start = fields[8].find("ID=") + len("ID=")
							end = fields[8].find(";", start)
							gene = fields[8][start:end]
							outstr += f"{value}\t{gene}\n"

	with open(f"{out_dir}{ref_name}_Mphar_prot_matches.txt", "w") as out:
		out.write(outstr)
