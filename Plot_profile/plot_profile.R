library(ggplot2)
library(data.table)
library(patchwork)
library(robustbase)
library(caTools)
theme_set(theme_classic())
theme_replace(panel.border = element_rect(fill=NA, size = 0.3),
              axis.title = element_blank(),
              axis.line = element_blank(),
              legend.position = "none",
              plot.title = element_text(hjust = 0.5, size = 5))

matformat <- function(mat, bed, label){
  bedid <- paste0(bed$V1, ":", bed$V2, "-", bed$V3)
  matbed <- as.matrix(mat[mat$V4 %in% bedid, -(1:6)])
  matbedmean <- colMeans(matbed, na.rm = T) - min(colMeans(matbed, na.rm = T))
  df <- data.frame(x = 1:length(matbedmean),y = matbedmean, id = label)
  return(df)
}


# read in sample matrices and regions
sample1 <- fread("sample1_matrix")
sample2 <- fread("sample2_matrix")
regions <- read.table("path/to/regions.bed")

# format for plotting 
plotdf <- rbind(
  matformat(sample1, regions, "sample1"),
  matformat(sample2, regions, "sample2"))

# plot
s = 0.2 # smoothing span
ggplot() +
  geom_smooth(data = plotdf[plotdf$id == "sample1",], aes(x = x, y = y), size = 0.3, color = "red", span = s, se = F) +
  geom_smooth(data = plotdf[plotdf$id == "sample2",], aes(x = x, y = y), size = 0.3, color = "blue", span = s, se = F) +
  scale_x_continuous(breaks = c(1,150,300), labels = c("-1.5kb", "center", "1.5kb"))
