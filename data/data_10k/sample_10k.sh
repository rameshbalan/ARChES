mkdir data_10k
for file in *.fq.gz;
do
  seqtk sample -s100 "${file}" 10000 > samp10k_"${file}"
  mv samp10k_"$file" data_10k/"${file}"
done;
cd data_10k
rename 's/.fq.gz/.fq/' *.fq.gz
for file in *.fq;
do
  gzip "${file}"
done;
