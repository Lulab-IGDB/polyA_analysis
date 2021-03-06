# bash mapping.sh sample_name

wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.primary_assembly.annotation.gtf.gz -P ./annotation/human/
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/GRCh38.primary_assembly.genome.fa.gz  -P ./annotation/human/

paftools.js gff2bed ./annotation/human/gencode.vM25.primary_assembly.annotation.gtf > ./annotation/human/gencode.vM25.primary_assembly.annotation.bed
minimap2 -x splice -t 20 -d ./annotation/human/human.mmi ./annotation/human/GRCh38.primary_assembly.genome.fa.gz

minimap2 -ax splice -uf --secondary=no -t 40 -L --MD --cs --junc-bed ./annotation/human/gencode.vM25.primary_assembly.annotation.bed ./annotations/human/human.mmi $1.clean.fastq.gz 2>align.log | samtools view -bS >$1.clean.aligned.bam

samtools sort $1.clean.aligned.bam > $1.clean.aligned.sort.bam
samtools index $1.clean.aligned.sort.bam
samtools view -F 3844 -bS $1.clean.aligned.bam > $1.clean.filter.bam
samtools view -f 4 -bS $1.clean.aligned.bam > $1.clean.unmapped.bam

featureCounts -L -g gene_id -t exon -s 1 -R CORE -a ./annotation/human/gencode.vM25.primary_assembly.annotation.gtf -o $1.featureCounts $1.clean.filter.bam &>featureCounts.log

bwa mem -x pacbio -t 15 ./annotation/human/human_RNA45SN1.fa $1.clean.fastq.gz 2>align.log | samtools view -bS -F 3844 >$1.chr_rRNA.align.bam
samtools view $1.chr_rRNA.align.bam | awk '{print $1}' >$1.chr_rRNA_reads.txt
rm $1.featureCounts $1.clean.filter.bed $1.chr_rRNA.align.bam

