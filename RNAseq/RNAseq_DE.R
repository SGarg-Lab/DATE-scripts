library(DESeq2)
library(limma)
library(Glimma)
library(edgeR)

### Load data
ct <- read.table("RNAseq_count.txt", header = T)
meta <- read.table("/path/to/rnaseq_metadata.txt", header = T)
meta$rep <- paste0("rep_", meta$rep)
meta$lane <- paste0("lane_", meta$lane)

geneinfo <- ct[ ,c(1:6)]
geneinfo$Chr <- sapply(1:nrow(geneinfo), function(x) substr(geneinfo$Chr[x], 1, 4))
ct <- ct[ ,c(7:33)]
colnames(ct) <- meta$sample
rownames(ct) <- geneinfo$Geneid

### Set up DGE
dge <- DGEList(counts = ct, samples = meta$sample)
dge$samples$clone <- as.factor(meta$clone)
dge$samples$genotype <- as.factor(meta$genotype)
dge$samples$lane <- as.factor(meta$lane)
dge$genes <- geneinfo

cpm <- cpm(dge)
lcpm <- cpm(dge, log=TRUE)
rpkm <- rpkm(dge, gene.length = dge$genes$Length)
L <- mean(dge$samples$lib.size) * 1e-6
M <- median(dge$samples$lib.size) * 1e-6

# Filter for lowly expressed genes
keep.exprs <- filterByExpr(dge, group = dge$samples$genotype)
dge <- dge[keep.exprs,, keep.lib.sizes=FALSE]

# Normalize read counts
dge <- calcNormFactors(dge, method = "TMM")
dge$samples$norm.factors
# 1/NormFactor can be used to scale bigWig file for visualization

### Differential expression analysis
# Set up design matrix
design <- model.matrix(~0+meta$genotype)
rownames(design) <- meta$sample
colnames(design) <- c("Klf4KO","WT","Zfp281KO")

# Model mean-variance relationship with voom
contrasts <- makeContrasts(
   WTvsKlf4KO = WT-Klf4KO,
   WTvsZfp281KO = WT-Zfp281KO,
   Klf4KOvsZfp281KO = Klf4KO-Zfp281KO,
   levels = colnames(design)
)

# Linear modelling
v <- voom(x, design, plot=TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contrasts)
efit <- eBayes(vfit)
summary(decideTests(efit))

# Output results
results <- data.frame(
  genes               = dge$genes,
  FC_WTvsKlf4KO       = topTable(efit,sort="none",n=Inf, coef = "WTvsKlf4KO")$logFC,
  pv_WTvsKlf4KO       = topTable(efit,sort="none",n=Inf, coef = "WTvsKlf4KO")$adj.P.Val,
  FC_WTvsZfp281KO     = topTable(efit,sort="none",n=Inf, coef = "WTvsZfp281KO")$logFC,
  pv_WTvsZfp281KO     = topTable(efit,sort="none",n=Inf, coef = "WTvsZfp281KO")$adj.P.Val,
  FC_Klf4KOvsZfp281KO = topTable(efit,sort="none",n=Inf, coef = "Klf4KOvsZfp281KO")$logFC,
  pv_Klf4KOvsZfp281KO = topTable(efit,sort="none",n=Inf, coef = "Klf4KOvsZfp281KO")$adj.P.Val,
  avg_expression      = topTable(efit,sort="none",n=Inf)$AveExpr
)


