#!/bin/bash

# example script for processing ChIPseq for TF binding analysis

module load perl/5.24.1 sra/2.10.0 samtools/1.10 bowtie2/2.3.5.1 macs2/2.2.7.1

mm10_ref=path/to/bowtie2/mm10/genome/reference

sample="
IPsample_rep1
IPsample_rep2
"

control="
inputsample
"


mkdir fastq bam

for F in ${sample}; do

  prefetch -v $F
  fastq-dump --outdir fastq --split-files $F/${F}.sra
  bowtie2 -x $mm10_ref -U fastq/${F}_1.fastq | samtools sort > bam/$F.bam
  samtools index $F.bam

done


for F in ${control}; do

  prefetch -v $F
  fastq-dump --outdir fastq --split-files $F/${F}.sra
  bowtie2 -x $mm10_ref -U fastq/${F}_1.fastq | samtools sort > bam/$F.bam
  samtools index $F.bam

done

# call peaks
macs2 callpeak -t bam/sample1.bam bam/sample2.bam -c bam/input.bam -f BAM -g mm -B --outdir peaks
