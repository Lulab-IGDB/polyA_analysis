# bash split.sh sample_name
gunzip -c $1.ccs.fastq.gz >$1.ccs.fastq
fastq-splitter.pl --n-parts 100 $1.ccs.fastq

for i in {001..100}; do echo "./PAIso-seq2/CCS_split_clean_UMI_V4.py $1.ccs.part-${i}.fastq $1.ccs.pass.txt barcode.fa 2 1>$1.${i}.out.txt 2>$1.${i}.err.txt"; done >demultiplex_run.sh

ParaFly -c demultiplex_run.sh -CPU 40

cat $1.*.out.txt >$1.out.txt
cat $1.*.err.txt >$1.err.txt
cat $1.out.txt | awk '{print "@"$2"\n"$6"\n+\n"$7}' | gzip -nc >$1.clean.fastq.gz
wc -l $1.out.txt >> $1.counts.txt
wc -l $1.err.txt >> $1.counts.txt
count fq $1.clean.fastq.gz >> $1.counts.txt
sed -i '$ s/$/\t'"$1"'.clean.fastq.gz\n/' $1.counts.txt

gzip $1.err.txt $1.out.txt

rm $1.*.err.txt $1.*.out.txt *fastq.log.txt

