import argparse
import glob

parser = argparse.ArgumentParser(description='Pull out common orthogroups among all species')
parser.add_argument("N",type=str,help='Number of species in your orthofinder input')
parser.add_argument("P",type=str,help='Path to Orthogroup.tsv file in orthofinder results.\
									Example: ~/Dosage_Comp/results/proteomes/Results_Jul05')

args = parser.parse_args()

results_path = glob.glob(args.P+"/Results_*")[0]

with open(results_path+"/Orthogroups/Orthogroups.tsv","r") as infile:
	with open(results_path+"/Orthogroups/common_orthogroups.txt","w") as outfile:
		next(infile)
		for line in infile:
			counter = 0
			each_orthogroup = line.rstrip().split("\t")
			for orthologs in each_orthogroup:
				if orthologs.strip()=='' or len(each_orthogroup) != (int(args.N)+1):
					counter += 1
			if counter == 0:
				outfile.write(line)
