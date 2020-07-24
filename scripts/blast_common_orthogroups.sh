#!/bin/bash

##########################################
# $1 - Path to Reference protein
# $2 - Path to common_orthogroups.tsv file
# $3 - Number of threads/cores
##########################################

echo "========================================================="
echo "1st argument is the reference protein file"
echo "2nd argument is the base directory of orthofinder results"
echo "Ex: results/proteomes/OrthoFinder"
echo "3rd argument is the number of threads to use"
echo "========================================================="

mkdir results/blast_files

results_path=`ls $2`

common_orthogroups=`awk '{print $1}' "$2/$results_path/Orthogroups/common_orthogroups.txt" | tr "\n" " "`
for og in $common_orthogroups;
do
	echo "BLASTing $2/$results_path/Orthogroup_Sequences/${og}.fa against ${1}"
	blastp -query "$2/$results_path/Orthogroup_Sequences/${og}.fa" -subject "${1}" -outfmt "6 qseqid sseqid evalue" -max_target_seqs 1 -out "${og}.txt" -num_threads $3
	mv "${og}.txt" results/blast_files
done;
