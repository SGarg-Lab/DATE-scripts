normct <- read.table("normalized.counts.txt", header = T)
datebed <- read.table("DATEs.bed")
satebed <- read.table("SATEs.bed")

# fourth replicate did not contain control sgRNAs
normct[normct$Gene == "Control", c("R4_NH","R4_NL","R4_UN")] <- NA

# subset sgRNAs that are present in all replicates
normct_use <- normct[normct$Gene == "Control" | ((normct$R1_UN > 0 & normct$R2_UN > 0 & normct$R3_UN > 0 & normct$R4_UN > 0) & rowSums(normct[ ,c("R1_NH","R1_NL","R2_NH","R2_NL","R3_NH","R3_NL", "R4_NH","R4_NL")] > 100)), ]

# format data and calculate enrichment
df <- data.frame(sgRNA = normct_use$sgRNA, gene = normct_use$Gene,
                 NHmean = rowMeans(normct_use[ ,grep("NH", colnames(normct_use))], na.rm = T),
                 NLmean = rowMeans(normct_use[ ,grep("NL", colnames(normct_use))], na.rm = T),
                 UNmean = rowMeans(normct_use[ ,grep("UN", colnames(normct_use))], na.rm = T))
df$NHUN <- log2(df$NHmean / df$UNmean)
df$NLUN <- log2(df$NLmean / df$UNmean)
df$NHNL <- log2(df$NHmean / df$NLmean)

df$type <- "SATE"
df$type[df$gene %in% date] <- "DATE"
df$type[df$gene == "Control"] <- "Control"

df_complete <- df[is.finite(df$NHUN) & is.finite(df$NLUN) & is.finite(df$NHNL), ]
df_complete <- df_complete[complete.cases(df_complete), ]
rownames(df_complete) <- df_complete$sgRNA

df_complete <- df_complete[order(-df_complete$NHUN), ]
df_complete$rank_NHUN <- 1:nrow(df_complete)
df_complete <- df_complete[order(-df_complete$NLUN), ]
df_complete$rank_NLUN <- 1:nrow(df_complete)
df_complete <- df_complete[order(-df_complete$NHNL), ]
df_complete$rank_NHNL <- 1:nrow(df_complete)

r <- 150
top_hits <- df_complete[df_complete$rank_NHNL <= r & df_complete$NHUN > 0 & df_complete$NLUN < 0, ]
bot_hits <- df_complete[(df_complete$rank_NHNL > nrow(df_complete) - r) & df_complete$NHUN < 0 & df_complete$NLUN > 0, ]
