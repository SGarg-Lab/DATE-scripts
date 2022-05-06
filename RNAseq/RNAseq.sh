#!/bin/bash

# example script for processing RNAseq

module load TrimGalore/0.6.5 fastxtoolkit/0.0.13
module load hisat2/2.2.0 samtools/1.10
module load python/2.7.13 r/3.6.0
module load subread/1.6.2

cd /path/to/fastq/files
mkdir fastq_trim fastq_trim_filtlen fastq_trim_filtlen_filtq align

# list of fastq files
F="
sample1
sample2
sample3
"
mm10_index=path/to/hisat2/mm10/genome/index
GTF=path/to/gencode.vM15.primary_assembly.annotation.gtf

for F in ${List}; do

  echo trimming adapter for $F:
  trim_galore --paired --illumina --dont_gzip --length 20 ${F}_1.fq.gz ${F}_2.fq.gz -o fastq_trim/

  gunzip fastq_trim/${F}_1.fq.gz
  gunzip fastq_trim/${F}_2.fq.gz

  echo filtering reads for length 20 for $F:
  awk 'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0&&length(b)>20{print a"\n"b"\n"c"\n"$0;}' fastq_trim/${F}_1.fq > fastq_trim_filtlen/${F}_1.fq
  awk 'NR%4==1{a=$0} NR%4==2{b=$0} NR%4==3{c=$0} NR%4==0&&length(b)>20{print a"\n"b"\n"c"\n"$0;}' fastq_trim/${F}_2.fq > fastq_trim_filtlen/${F}_2.fq

  echo filtering reads for read quality for $F:
  fastq_quality_filter -Q33 -q 30 -p 90 -i fastq_trim_filtlen/${F}_1.fastq -o fastq_trim_filtlen_filtq/${F}_1.fastq -v
  fastq_quality_filter -Q33 -q 30 -p 90 -i fastq_trim_filtlen/${F}_2.fastq -o fastq_trim_filtlen_filtq/${F}_2.fastq -v

  echo aligning to mm10 genome for $F:
  hisat2 -q -x $mm10_index -1 fastq_trim_filtlen_filtq/${F}_1.fastq \
    -2 fastq_trim_filtlen_filtq/${F}_2.fastq \
    -S align/$F.sam \
    --rna-strandness FR --summary-file align/summary_${F}.txt
  samtools sort align/$F.sam > align/$F.bam
  samtools index -b align/$F.bam
  mv *.bai align

done

# count reads
featureCounts --verbose -p -a $GTF -g gene_name -o RNAseq_count.txt align/*.bam

# differential expression analysis
Rscript /path/to/RNAseq_DE.R
