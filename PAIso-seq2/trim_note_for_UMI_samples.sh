# bash trim_note_for_UMI_samples.sh sample_name
./bin/PolyA_trim_V5.4.1.py $1.clean.filter.bam >$1.polyA_trim.out.txt
./bin/PolyA_note_V2.1.py $1.polyA_trim.out.txt $1.clean.filter.bam.featureCounts 1>$1.polyA_note.txt 2> $1.polyA_note.err.txt

bedtools bamtobed -i $1.clean.filter.bam > $1.clean.filter.bed
./PAIso-seq2/pickOne_V2.py $1.clean.filter.bed > $1.UMI_uniq.out.txt
./bin/fisher.pl -fc 2 $1.UMI_uniq.out.txt $1.polyA_note.txt > $1.polyA_note.UMI_uniq.txt

rm $1.clean.filter.bed
