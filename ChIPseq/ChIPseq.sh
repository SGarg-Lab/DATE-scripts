#!/bin/bash

# example script for processing ChIPseq

module load bwa/0.7.17 samtools/1.10 bedtools/2.29.2 ucsc-tools/20170321 python/2.7.13 deeptools/3.0.1
mm10_ref=path/to/bwa/mm10/genome/reference
blacklist=path/to/ENCFF547MET_ENCODE_mm10_blacklist.bed

cd /path/to/fastq/files
mkdir bam bigwig bigwig_subtract peaks

# list of fastq files
List="
sample1
sample2
control1
control2
"

# align ChIP data
for F in ${List}; do

  echo aligning $F:
  bwa mem $mm10_ref ${F}1_sequence.fastq ${F}2_sequence.fastq  | samtools sort > bam/$F.bam
  samtools sort -n bam/$F.bam -o bam/${F}_namesorted.bam

  echo converting $F to bigwig:
  bamCoverage -b bam/${F}_namesorted.bam -o bigwig/${F}.bigwig --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --centerReads -p 6 --extendReads

done

# subtract control ChIP signal
bigwigCompare -b1 bigwig/sample.bigwig -b2 bigwig/control.bigwig -o bigwig_subtract/sample.bigwig

List2="
sample1
sample2
"

# call peaks
for F2 in ${List}2; do

  echo callping peaks for $F2:
  Genrich -t bam/${F2}_namesorted.bam -o peaks/${F2}.narrowPeak -y -a 1 -e chrM,chrY -v -E $blacklist

done
