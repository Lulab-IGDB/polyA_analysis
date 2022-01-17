# bash extract_and_filter.sh sample_name
fisher.pl --except -fc 2 $1.chr_rRNA_reads.txt $1.polyA_note.txt | fisher.pl -fc 3 ./annotations/human/chr_pc.txt - >$1.chr_pc.polyA_note.pass1.txt
fisher.pl --except -fc 2 $1.chr_rRNA_reads.txt $1.polyA_note.txt | fisher.pl -fc 3 ./annotations/human/chr_lnc.txt - >$1.chr_lnc.polyA_note.pass1.txt
fisher.pl --except -fc 2 $1.chr_rRNA_reads.txt $1.polyA_note.txt | fisher.pl -fc 3 ./annotations/human/mt_pc.txt - >$1.mt_pc.polyA_note.pass1.txt
fisher.pl --except -fc 2 $1.chr_rRNA_reads.txt $1.polyA_note.txt | fisher.pl -fc 3 ./annotations/human/mt_rRNA.txt - >$1.mt_rRNA.polyA_note.pass1.txt
fisher.pl --except -fc 2 $1.chr_rRNA_reads.txt $1.polyA_note.txt | fisher.pl -fc 3 ./annotations/human/hist.txt - >$1.hist.polyA_note.pass1.txt
fisher.pl --except -fc 2 $1.chr_rRNA_reads.txt $1.polyA_note.txt | fisher.pl -fc 3 ./annotations/human/hist_var.txt - >$1.hist_var.polyA_note.pass1.txt
fisher.pl -fc 2 $1.chr_rRNA_reads.txt $1.polyA_note.txt >$1.chr_rRNA.polyA_note.pass1.txt
fisher.pl --except -fc 2 $1.chr_rRNA_reads.txt $1.polyA_note.txt | awk '$3 == "Unknown"' >$1.unknown.polyA_note.pass1.txt

# pass >= 10
cat $1.chr_pc.polyA_note.pass1.txt | awk '$4 >= 10' >$1.chr_pc.polyA_note.pass10.txt
cat $1.chr_lnc.polyA_note.pass1.txt | awk '$4 >= 10' >$1.chr_lnc.polyA_note.pass10.txt
cat $1.mt_pc.polyA_note.pass1.txt | awk '$4 >= 10' >$1.mt_pc.polyA_note.pass10.txt
cat $1.mt_rRNA.polyA_note.pass1.txt | awk '$4 >= 10' >$1.mt_rRNA.polyA_note.pass10.txt
cat $1.hist.polyA_note.pass1.txt | awk '$4 >= 10' >$1.hist.polyA_note.pass10.txt
cat $1.hist_var.polyA_note.pass1.txt | awk '$4 >= 10' >$1.hist_var.polyA_note.pass10.txt
cat $1.chr_rRNA.polyA_note.pass1.txt | awk '$4 >= 10' >$1.chr_rRNA.polyA_note.pass10.txt
cat $1.unknown.polyA_note.pass1.txt | awk '$4 >= 10' >$1.unknown.polyA_note.pass10.txt
