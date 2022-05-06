# Differential accessibility analysis

library(GenomicRanges)
library(csaw)
library(limma)
library(edgeR)
library(ggplot2)

bampath <- "path/to/bam/files"
genotype <- c("Klf4KO","Klf4KO","Klf4KO","Zfp281KO","Zfp281KO","Zfp281KO","WT","WT","WT")
sampleid <- c("Klf4KO_1","Klf4KO_2","Klf4KO_3","Zfp281KO_1","Zfp281KO_2","Zfp281KO_3","WT_1","WT_2","WT_3")

### Specify input files
peaks <- read.table("consensus_peak_list.bed")
colnames(peaks) <- c("chr","start","end")
peaks <- GRanges(peaks)

# specify bam files
setwd(bampath)
bams <- list.files(bampath, pattern = "*_coordsort.bam$")

# specify blacklist region
blacklist <- read.table("ENCFF547MET_ENCODE_mm10_blacklist.bed")
colnames(blacklist) <- c("chr","start","end")
blacklist <- GRanges(blacklist)

# specify read paramters
standard.chr <- paste0("chr", c(1:19))
param <- readParam(pe = "both", discard = blacklist, restrict = standard.chr)

### Count reads in peaks
peak.counts <- regionCounts(bams, peaks, param=param)
colnames(peak.counts) <- sampleid
ranges <- data.frame(rowRanges(peak.counts))
ranges$id <- paste0(ranges$seqnames, "_", ranges$start, "_", ranges$end)
rownames(peak.counts) <- ranges$id

### Set up DGE
dge <- DGEList(counts = assay(peak.counts, "counts"), samples = genotype)

# Filter for lowly enriched peaks
keep <- rowSums(cpm(dge) > 1) >= 5
dge <- dge[keep, ]

# Normalize read counts
dge <- calcNormFactors(dge)
dge$samples$norm.factors
# 1/NormFactor can be used to scale bigWig file for visualization

### Differential accessibility analysis
# Set up design matrix
groups <- data.frame(file = bams, genotype = genotype)
condition <- factor(groups$genotype, levels = c("WT", "Klf4KO", "Zfp281KO"))
mydesign <- model.matrix(~0 + genotype)

# Model mean-variance relationship with voom
myvoom <- voom(counts = dge, design = mydesign)
plotMDS(myvoom, top = 100)

contrasts <- makeContrasts(
  "Klf4KO_vs_WT" = genotypeKlf4KO - genotypeWT,
  "Zfp281KO_vs_WT" = genotypeZfp281KO - genotypeWT,
  "Klf4KO_vs_Zfp281KO" = genotypeKlf4KO - genotypeZfp281KO,
  levels = myvoom$design
)

# Linear modelling
fit <- lmFit(myvoom, design = myvoom$design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
summary(decideTests(fit))

# Output results
results <- data.frame(
  FC_Klf4KO_vs_WT       = topTable(fit,sort="none",n=Inf, coef = "Klf4KO_vs_WT")$logFC,
  pv_Klf4KO_vs_WT       = topTable(fit,sort="none",n=Inf, coef = "Klf4KO_vs_WT")$adj.P.Val,
  FC_Zfp281KO_vs_WT     = topTable(fit,sort="none",n=Inf, coef = "Zfp281KO_vs_WT")$logFC,
  pv_Zfp281KO_vs_WT     = topTable(fit,sort="none",n=Inf, coef = "Zfp281KO_vs_WT")$adj.P.Val,
  FC_vs_Zfp281KO        = topTable(fit,sort="none",n=Inf, coef = "Klf4KO_vs_Zfp281KO")$logFC,
  pv_Klf4KO_vs_Zfp281KO = topTable(fit,sort="none",n=Inf, coef = "Klf4KO_vs_Zfp281KO")$adj.P.Val,
  pv_all                = topTable(fit,sort="none",n=Inf)$adj.P.Val,
  avg_expression        = topTable(fit,sort="none",n=Inf)$AveExpr
)
