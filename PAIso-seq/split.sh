# bash split.sh sample_name
gunzip -c $1.fastq.gz >$1.ccs.fastq
fastq-splitter.pl --n-parts 100 $1.ccs.fastq

for i in {001..100}; do echo "python ./PAIso-seq/CCS_split_clean_end_extension_v1.py $1.ccs.part-${i}.fastq $1.ccs.pass.txt $1.barcode.fa 2 1>$1.${i}.out.txt 2>$1.${i}.err.txt"; done >demultiplex_run.sh

ParaFly -c demultiplex_run.sh -CPU 48 

cat $1.*.out.txt >$1.out.txt
cat $1.*.err.txt >$1.err.txt
cat $1.out.txt | awk '{print "@"$2"\n"$6"\n+\n"$7}' | gzip -nc >$1.clean.fastq.gz
wc -l $1.out.txt >> $1.counts.txt
wc -l $1.err.txt >> $1.counts.txt
count fq $1.clean.fastq.gz >> $1.counts.txt
sed -i '$ s/$/\t'"$1"'.clean.fastq.gz\n/' $1.counts.txt

gzip $1.err.txt $1.out.txt

rm $1.*.err.txt $1.*.out.txt *fastq.log.txt
