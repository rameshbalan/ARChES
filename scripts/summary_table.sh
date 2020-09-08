#!/bin/bash

for value in "$@";
do
	echo "$value" > "$value".txt
	grep -c ">" results/"$value"/trinity_out_dir/"$value".fasta >> "$value".txt
	grep -c ">" "$value".fasta.transdecoder.cds >> "$value".txt
	grep -c ">" results/"$value"_nr95.fasta >> "$value".txt
	grep -c ">" "$value"_nr95.fasta.transdecoder.cds >> "$value".txt
	paste -sd, "$value".txt >> results/summary_table.csv
done;

