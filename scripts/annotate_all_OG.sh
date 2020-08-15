
# $1 is the gff file

for file in results/blast_files/*.txt;
do
	echo "${file}"
	python scripts/add_OG_and_annotate.py --GFF $1 --FILE "${file}" --neox $3 --x $4 --un $5
done;

cat results/blast_files/*_out.txt > $2
