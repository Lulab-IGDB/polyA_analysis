# bash mapping.sh sample_name
minimap2 -ax splice -uf --secondary=no -t 15 -L --MD --cs --junc-bed ./annotations/mouse/gencode.vM25.primary_assembly.annotation.bed ./annotations/mouse/mouse.mmi $1.clean.fastq.gz 2>align.log | samtools view -bS >$1.clean.aligned.bam

samtools sort $1.clean.aligned.bam > $1.clean.aligned.sort.bam
samtools index $1.clean.aligned.sort.bam
samtools view -F 3844 -bS $1.clean.aligned.bam > $1.clean.filter.bam
samtools view -f 4 -bS $1.clean.aligned.bam > $1.clean.unmapped.bam

./script_soft/featureCounts -L -g gene_id -t exon -s 1 -R CORE -a ./annotations/mouse/gencode.vM25.primary_assembly.annotation.gtf -o $1.featureCounts $1.clean.filter.bam &>featureCounts.log

bedtools bamtobed -i $1.clean.filter.bam > $1.clean.filter.bed

bwa mem -x pacbio -t 15 ./annotations/mouse/Rn45s.fa $1.clean.fastq.gz 2>align.log | samtools view -bS -F 3844 >$1.chr_rRNA.align.bam
samtools view $1.chr_rRNA.align.bam | awk '{print $1}' >$1.chr_rRNA_reads.txt
rm $1.featureCounts $1.clean.filter.bed $1.chr_rRNA.align.bam