library(universalmotif)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(devtools)
library(viridis)
library(usedist)
library(TFBSTools)

setwd("/path/to/PWM/files")
dicty <- c()
for (file1 in list.files(getwd())){
  for (file2 in list.files(getwd())){

    # read motif for first TF, convert to PWM
    tf1 <- basename(file1)
    if (grepl('h_', tf1, fixed = TRUE)){
      motif1 <- read_homer(file1)
    } else{
      motif1 <- read_meme(file1)[[1]]
    }
    pfm1 <- convert_motifs(motif1, class="TFBSTools-PFMatrix")
    pwm1 <- toPWM(pfm1, type="prob")

    # find reverse complement
    pfm1_RC <- convert_motifs(motif_rc(motif1), class="TFBSTools-PFMatrix")
    pwm1_RC <- toPWM(pfm1_RC, type="prob")

    # repeat for second TF
    tf2 <- basename(file2)
    if (grepl('h_', tf2, fixed = TRUE)){
      motif2 <- read_homer(file2)
    } else{
      motif2 <- read_meme(file2)[[1]]
    }
    pfm2 <- convert_motifs(motif2, class="TFBSTools-PFMatrix")
    pwm2 <- toPWM(pfm2, type="prob")
    pfm2_RC <- convert_motifs(motif_rc(motif2), class="TFBSTools-PFMatrix")
    pwm2_RC <- toPWM(pfm2_RC, type="prob")

    # calculate KL distance
    val_00 <- PWMSimilarity(pwm1, pwm2, method="KL")
    val_01 <- PWMSimilarity(pwm1, pwm2_RC, method="KL")
    val_10 <- PWMSimilarity(pwm1_RC, pwm2, method="KL")
    val_11 <- PWMSimilarity(pwm1_RC, pwm2_RC, method="KL")
    dicty[paste(tf1, tf2, "00", sep=" ")] <- val_00
    dicty[paste(tf1, tf2, "01", sep=" ")] <- val_01
    dicty[paste(tf1, tf2, "10", sep=" ")] <- val_10
    dicty[paste(tf1, tf2, "11", sep=" ")] <- val_11
  }
}

# output to text file
filename='pwmsimilarity.txt'
file.create(filename)
for (key in names(dicty)){
  cat(c(key, as.character(dicty[key])), "\n", file=filename, append=TRUE)
}
cat(readLines(filename), sep="\n")

# take minimum of every four rows to determine minimum KL distance for each TF pair
pwmsimilarity <- read.table("pwmsimilarity.txt", header = F)
kldistance <- data.frame()
for(i in 1:nrow(pwmsimilarity)/4){
  pwmsimilarity_subset <- pwmsimilarity[c(4*i-3:4*i), ]
  kldistance[i, ] <- c(
    pwmsimilarity_subset[1,1],
    pwmsimilarity_subset[1,2],
    minimum(pwmsimilarity_subset[ ,4])
  )
}
