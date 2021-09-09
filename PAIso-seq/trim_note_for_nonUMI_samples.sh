./bin/PolyA_trim_V5.4.1.py $1.clean.filter.bam >$1.polyA_trim.out.txt
./bin/PolyA_note_V2.1.py $1.polyA_trim.out.txt $1.clean.filter.bam.featureCounts 1>$1.polyA_note.txt 2> $1.polyA_note.err.txt