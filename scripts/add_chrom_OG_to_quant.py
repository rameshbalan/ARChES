#!/usr/bin/python
import sys
import argparse


parser = argparse.ArgumentParser(description='Using the annotated_blast_file, this script adds chromosome number,name and orthogroup to expression values')
parser.add_argument("BLAST",type=str,help='Path to the annotated_blast_file file to parse.\
  Example: ~/Dosage_Comp/data/T_cast_ref.gff')
parser.add_argument("FILE",type=str,help='Path to TMM file to be annotated\
                  Example: results/G_corn_TMM.txt')
parser.add_argument("OUT",type=str,help='Path to annotated TMM file\
                  Example: results/G_corn_TMM.txt')
args = parser.parse_args()

annotate_dict = {}
with open(args.BLAST,'r') as annotateF:
  for line in annotateF:
    line_split = line.rstrip().split("\t")
    annotate_dict[line_split[0]] = [line_split[3],line_split[4],line_split[5]]

inputfile = args.FILE
outputfile = args.OUT
with open(inputfile, 'r') as infile:
  with open(outputfile, 'w') as outfile:
    next(infile)
    for line in infile:
      line_split = line.rstrip().split("\t")
      try:
        line_split.append(annotate_dict[line_split[0].rstrip('"').lstrip('"')][0])
        line_split.append(annotate_dict[line_split[0].rstrip('"').lstrip('"')][1])
        line_split.append(annotate_dict[line_split[0].rstrip('"').lstrip('"')][2])
        outfile.write("\t".join(line_split))
        outfile.write("\n")
      except KeyError:
        line_split.append("Unknown")
        line_split.append("Unknown")
        line_split.append("Unknown")
        outfile.write("\t".join(line_split))
        outfile.write("\n")