### bash reads_classification_APAsites.sh sample_name PAL_cutoff APAsites.csv_file_name
bedtools bamtobed -i $1.clean.filter.bam > $1.clean.filter.bed

cat $1.chr_pc.polyA_note.pass10.txt | awk ' $4 >= 10 && $6+$10 >= '"$2"'' > $1.chr_pc.polyA_note.pass10.pal$2.txt

./bin/readsOnAPA_v4.py $1.clean.filter.bed $3 $1.chr_pc.polyA_note.pass10.pal$2.txt $1.V7.1.V4