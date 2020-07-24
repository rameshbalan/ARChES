## This is the add_OG_and_annotate.py

import argparse


parser = argparse.ArgumentParser(description='Using the reference transcriptome hits, this script adds chromosome number/name to each orthogroup')
parser.add_argument("GFF",type=str,help='Path to the reference gff file to parse.\
	Example: ~/Dosage_Comp/data/T_cast_ref.gff')
parser.add_argument("FILE",type=str,help='Path to BLASTed Orthogroup file to be annotated\
									Example: ~/Dosage_Comp/results/proteomes/Results_Jul05')

args = parser.parse_args()

## This chunk creates a dictionary with tcas cds regions numbers with their corresponding name.
prot_dict = {}
chromosome_dict = {}
with open(args.GFF, 'r') as gff:
	for line in gff:
		if not line.startswith("#"):
			split_line = line.rstrip().split("\t")
			if split_line[2] == "region":
				try:
					chrom = split_line[8].rstrip().split("Name=")[1].split(";")[0]
					chromosome_dict[split_line[0].rstrip()] = [chrom.rstrip()]
					# prot_dict[prot] = [chrom.rstrip()]
					if chrom == "LG2":
						chromosome_dict[split_line[0].rstrip()].append("Neo-X")
					elif chrom == "LGX":
						chromosome_dict[split_line[0].rstrip()].append("X")
					elif chrom == "Unknown" or chrom == "MT":
						chromosome_dict[split_line[0].rstrip()].append("Un")
					else:
						chromosome_dict[split_line[0].rstrip()].append("Autosome")
				except IndexError:
					chrom = split_line[8].rstrip().split("linkage-group=")[1].split(";")[0]
					chromosome_dict[split_line[0].rstrip()] = [chrom.rstrip()]
					if chrom == "LG2":
						chromosome_dict[split_line[0].rstrip()].append("Neo-X")
					elif chrom == "LGX":
						chromosome_dict[split_line[0].rstrip()].append("X")
					elif chrom == "Unknown" or chrom == "MT":
						chromosome_dict[split_line[0].rstrip()].append("Un")
					else:
						chromosome_dict[split_line[0].rstrip()].append("Autosome")
			elif split_line[2] == "CDS":
				prot = split_line[8].rstrip().split("Name=")[1].split(";")[0]
				prot_dict[prot] = chromosome_dict[split_line[0].rstrip()]

file_name = args.FILE
outF_name = file_name.split(".txt")[0]

## This chunk adds orthogroup, chromosome number and name to each blast output file..
with open(file_name,'r') as infile:
	with open(outF_name+"_out.txt",'w') as outfile:
		for line in infile:
			line_split = line.rstrip().split("\t")
			line_split.append(prot_dict[line_split[1]][0])
			line_split.append(prot_dict[line_split[1]][1])
			line_split.append(outF_name.split('blast_files/')[1].split(".txt")[0])
			outfile.write("\t".join(line_split))
			outfile.write("\n")
