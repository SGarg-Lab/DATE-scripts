#!/bin/bash

# example script for processing ATACseq

module load TrimGalore/0.6.5 fastxtoolkit/0.0.13
module load bowtie2/2.3.5.1 samtools/1.10
module load Genrich/0.6 r/3.6.0

cd /path/to/fastq/files
mkdir fastq_trim align peak

# list of fastq files
F="
sample1
sample2
sample3
"

mm10_index=path/to/bowtie2/mm10/genome/index
GTF=path/to/gencode.vM15.primary_assembly.annotation.gtf
blacklist=path/to/ENCFF547MET_ENCODE_mm10_blacklist.bed

for F in ${List}; do

  echo trimming adapter for $F:
  trim_galore --paired --fastqc --nextera --gzip ${F}_1.fq.gz ${F}_2.fq.gz -o fast_trim/

  echo aligning to mm10 genome for $F:
  bowtie2 --very-sensitive -k 10 -x $mm10_index -1 fast_trim/${F}_1.fq.gz -2 fast_trim/${F}_2.fq.gz -S align/$F.sam

  # only keep autosomal, properly-paired, mapped read pairs with mapping quality >= 30
  samtools view -bS align/${F}.sam -u -b -h -f 3 -F 4 -F 8 -F 256 -F 1024 -F 2048 -q 30 | samtools sort - -o align/${F}.bam
  samtools index align/${F}.bam
  mv *.bai align

  echo calling peaks for $F:
  Genrich -t align/${F}.bam -o peak/$F.narrowPeak -a 10 -j -y -r -e chrM,chrY -v -E $blacklist

done

# compile consensus peak list
cat peak/*.narrowPeak | bedtools sort -i - | bedtools merge -i - -d 500 > consensus_peak_list.bed

# differential accessibility analysis
Rscript /path/to/ATACseq_DA.R
