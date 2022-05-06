#!/bin/bash

# example script for plotting profiles of normalized ChIP/mPRO/ATAC signal across regions

module load python/2.7.13 deeptools/3.0.1 r/4.0.4

regions=path/to/regions.bed
sample1=path/to/sample1.bigwig
sample2=path/to/sample2.bigwig

computeMatrix reference-point -R $regions -S $sample1 -o sample1_matrix.gz --referencePoint center -a 1500 -b 1500
computeMatrix reference-point -R $regions -S $sample2 -o sample2_matrix.gz --referencePoint center -a 1500 -b 1500

gunzip sample1_matrix.gz
gunzip sample2_matrix.gz

Rscript plot_profile.R
