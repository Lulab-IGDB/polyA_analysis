### bash APA_site_calling.sh sample_name PAL_cutoff
bedtools bamtobed -i $1.clean.filter.bam > $1.clean.filter.bed

cat $1.chr_pc.polyA_note.pass10.txt | awk ' $4 >= 10 && $6+$10 >= '"$2"'' > $1.chr_pc.polyA_note.pass10.pal$2.txt

./bin/findAPA_v7.1.py $1.chr_pc.polyA_note.pass10.pal$2.txt $1.clean.filter.bed $1.v7.1