# bash trim_note_for_UMI_samples_human_last_exon.sh sample_name
./bin/PolyA_trim_V5.5.2.1.py ./annotation/human/gencode.v36.primary_assembly.annotation.transcript.last_exon.txt $1.clean.filter.bam.featureCounts $1.clean.filter.bam | awk '$1=="Last_TRUE"' > $1.PolyA_trim_V5.5.2.1.txt

./bin/PolyA_note_V2.py $1.PolyA_trim_V5.5.2.1.txt $1.clean.filter.bam.featureCounts 1>$1.polyA_note.txt 2> $1.polyA_note.err.txt

bedtools bamtobed -i $1.clean.filter.bam > $1.clean.filter.bed
./PAIso-seq2/pickOne_V2.py $1.clean.filter.bed > $1.UMI_uniq.out.txt
./bin/fisher.pl -fc 2 $1.UMI_uniq.out.txt $1.polyA_note.txt > $1.polyA_note.UMI_uniq.txt

rm $1.clean.filter.bed