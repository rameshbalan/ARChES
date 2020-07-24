
# $1 is the gff file

for file in results/blast_files/*.txt;
do
	echo "${file}"
	python scripts/add_OG_and_annotate.py $1 "${file}"
done;

cat results/blast_files/*_out.txt > $2
