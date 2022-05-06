library(Seurat)
library(monocle)
library(ggplot2)
library(dplyr)
library(tibble)
library(caTools)

# Set up Seurat object from WT + KO cells
ids <- c("WT_rep1","WT_rep2","Klf4KO_rep1","Klf4KO_rep2","Zfp281KO_rep1","Zfp281KO_rep2")
data10x <- lapply(ids, function(x) Read10X(paste0(,x,"_filtered_feature_bc_matrix")))
objlist <- lapply(data10x, function(x) CreateSeuratObject(counts = x))
WTKO <- merge(objlist[[1]], y = c(objlist[[2]],objlist[[3]],objlist[[4]],objlist[[5]],objlist[[6]]),
                add.cell.ids = c("WT-1","WT-2","Klf4KO-1","Klf4KO-2","Zfp281KO-1","Zfp281KO-2"),
                project = "KO_ESC")
WTKO@meta.data$sample <- sapply(1:nrow(WTKO@meta.data), function(x) sapply(strsplit(rownames(WTKO@meta.data[x,]),"_"),`[`, 1))
WTKO@meta.data$genotype <- sapply(1:nrow(WTKO@meta.data), function(x) sapply(strsplit(WTKO@meta.data$sample[x],"-"),`[`, 1))


# QC and normalization
WTKO <- PercentageFeatureSet(WTKO, pattern = "^mt-", col.name = "percent.mt")
WTKO <- PercentageFeatureSet(WTKO, pattern = "^Rps", col.name = "percent.rps")
WTKO <- PercentageFeatureSet(WTKO, pattern = "^Rpl", col.name = "percent.rpl")
WTKO <- subset(WTKO, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 10 & percent.rpl > 5 & percent.rps > 5)

s.genes <- stringr::str_to_title(cc.genes$s.genes)
g2m.genes <- stringr::str_to_title(cc.genes$g2m.genes)
WTKO <- CellCycleScoring(WTKO, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
WTKO <- SCTransform(WTKO, vars.to.regress = c("percent.mt", "nCount_RNA", "Phase"), return.only.var.genes = FALSE, verbose = TRUE)


# Dimensionaltiy reduction
WTKO <- RunPCA(WTKO, features = rownames(WTKO))
WTKO <- FindNeighbors(WTKO , dim = 1:20)
WTKO <- FindClusters(WTKO, resolution = 1)
WTKO <- RunUMAP(WTKO, dims = 1:20)
saveRDS(WTKO, "scRNAseq_redo.RDS")


# Assign cells state score
addscore <- function(signature, label, seuratobj){
  dat <- data.frame(seuratobj[["SCT"]]@scale.data)
  countsum <- colSums(dat[rownames(dat) %in% signature, ])
  zscore <- (countsum - min(countsum)) / (max(countsum) - min(countsum))
  seuratobj <- AddMetaData(seuratobj, metadata = zscore, col.name = label)
  return(seuratobj)
}

gmx <- read.table("State_signatures.gmx", skip = 2)
WTKO <- addscore(as.character(unlist(gmx$V1)), "s1", WTKO)
WTKO <- addscore(as.character(unlist(gmx$V2)), "s2", WTKO)
WTKO <- addscore(as.character(unlist(gmx$V3)), "s3", WTKO)
WTKO@meta.data$state <- colnames(WTKO@meta.data[ ,c("s1","s2","s3")])[apply(WTKO@meta.data[ ,c("s1","s2","s3")],1,which.max)]


# Determine Klf4 and Zfp281 program score ratio
targetK <- read.table("Klf4_target.txt")$V1
targetZ <- read.table("Zfp281_target.txt")$V1
WTKO <- addscore(targetK, "Klf4target", WTKO)
WTKO <- addscore(targetZ, "Zfp281target", WTKO)
WTKO$Klf4Zfp281targetratio <- log2(WTKO$Klf4target / WTKO$Zfp281target)



# Sample equal number of cells for trajectory analysis
n <- 10000
seurat_wt <- WTKO@meta.data[WTKO@meta.data$genotype == "WT", ]
seurat_klf4ko <- WTKO@meta.data[WTKO@meta.data$genotype == "Klf4KO", ]
seurat_zfp281ko <- WTKO@meta.data[WTKO@meta.data$genotype == "Zfp281KO", ]
sample_wt <- rownames(seurat_wt[sample(nrow(seurat_wt), n, replace = F), ])
sample_klf4ko <- rownames(seurat_klf4ko[sample(nrow(seurat_klf4ko), n, replace = F), ])
sample_zfp281ko <- rownames(seurat_zfp281ko[sample(nrow(seurat_zfp281ko), n, replace = F), ])
WTKO_subset_cells <- WTKO[ ,colnames(WTKO) %in% c(sample_wt, sample_klf4ko, sample_zfp281ko)]


# Set up monocle CDS for trajectory analysis
data <- as(as.matrix(WTKO_subset_cells[["SCT"]]@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = WTKO_subset_cells@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
WTKO_mon <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = uninormal())

# Run ordering algorithm
var_genes <- WTKO_subset_cells[["SCT"]]@var.features
WTKO_mon <- setOrderingFilter(WTKO_mon, var_genes)
WTKO_mon <- reduceDimension(WTKO_mon, norm_method="none", reduction_method="DDRTree")
WTKO_mon <- orderCells(WTKO_mon, num_paths = 2)

# Plot trajectory
sample_state <- pData(WTKO_mon)$State
lib_info_with_pseudo <- pData(WTKO_mon)
reduced_dim_coords <- reducedDimK(WTKO_mon)
ica_space_df <- Matrix::t(reduced_dim_coords) %>%
  as.data.frame() %>%
  select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
  mutate(sample_name = rownames(.), sample_state = rownames(.))
dp_mst <- minSpanningTree(WTKO_mon)
edge_df <- dp_mst %>%
  igraph::as_data_frame() %>%
  select_(source = "from", target = "to") %>%
  left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
  left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")
data_df <- t(monocle::reducedDimS(WTKO_mon)) %>%
  as.data.frame() %>%
  select_(data_dim_1 = 1, data_dim_2 = 2) %>%
  rownames_to_column("sample_name") %>%
  mutate(sample_state) %>%
  left_join(lib_info_with_pseudo %>% rownames_to_column("sample_name"), by = "sample_name")

return_rotation_mat <- function(theta) {
  theta <- theta / 180 * pi
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
}
rot_mat <- return_rotation_mat(270)
data_df[, cn1] <- as.matrix(data_df[, c("data_dim_1", "data_dim_2")]) %*% t(rot_mat)
edge_df[, cn2] <- as.matrix(edge_df[, c("source_prin_graph_dim_1", "source_prin_graph_dim_2")]) %*% t(rot_mat)
edge_df[, cn3] <- as.matrix(edge_df[,  c("target_prin_graph_dim_1", "target_prin_graph_dim_2")]) %*% t(rot_mat)

# Format data - use for plotting
plotdf <- data.frame(Embeddings(WTKO, "umap"), 
                     WTKO@meta.data,
                     trajectory_dim_1 = data_df$data_dim_1,
                     trajectory_dim_2 = data_df$data_dim_2,
                     progression = pData(WTKO_mon$Pseudotime),
                     state = pData(WTKO_mon$State))

# Example UMAP plot
ggplot(pltodf, aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(aes(color = genotype))

# Example trajectory plot 
ggplot(plotdf, aes(x = trajectory_dim_1, y = trajectory_dim_2)) + 
  geom_point(aes(color = genotype)) + 
  xlab("Trajectory Component 1") + ylab("Trajectory Component 2") +  
  geom_segment(data=edge_df, aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), linetype="solid", na.rm=TRUE)

# Example branch plot 
plotdf$branch <- "Root"
plotdf$branch[plotdf$State %in% c(4,5)] <- "Branch1"
plotdf$branch[plotdf$State %in% c(1,3)] <- "Branch2"
plotdf$branch[plotdf$State %in% c(7,8,2)] <- "Empty"
plotdf$Klf4KOpresence <- 0 
plotdf$Klf4KOpresence[plotdf$genotype == "Klf4KO"] <- 1
plotdf$ZfpKOpresence <- 0
plotdf$ZfpKOpresence[plotdf$genotype == "Zfp281KO"] <- 1
plotdf$WTpresence <- 0
plotdf$WTpresence[plotdf$genotype == "WT"] <- 1

calcbranch <- function(plotdf, branchselect){
  branchdf <- plotdf[plotdf$branch == branchselect, ]
  k <- 500
  
  branchdf_calc <- data.frame(
    S1 = unlist(runquantile(branchdf$s1, k, probs = 0.5)),
    S2 = unlist(runquantile(branchdf$s2, k, probs = 0.5)),
    S3 = unlist(runquantile(branchdf$s3, k, probs = 0.5)),
    S1q25 = unlist(runquantile(branchdf$top200.s1, k, probs = 0.25)),
    S2q25 = unlist(runquantile(branchdf$top200.s2, k, probs = 0.25)),
    S3q25 = unlist(runquantile(branchdf$top200.s3, k, probs = 0.25)),
    S1q75 = unlist(runquantile(branchdf$top200.s1, k, probs = 0.75)),
    S2q75 = unlist(runquantile(branchdf$top200.s2, k, probs = 0.75)),
    S3q75 = unlist(runquantile(branchdf$top200.s3, k, probs = 0.75)), 
    Pseudotime = unlist(runmean(branchdf$Pseudotime, k)), 
    Klf4KO = unlist(runmean(branchdf$Klf4KOpresence, k)),
    Zfp281KO = unlist(runmean(branchdf$ZfpKOpresence, k)),
    WT = unlist(runmean(branchdf$WTpresence,k)),
    KZ = unlist(runmean(branchdf$Klf4Zfp281targetratio, k)),
    KZq25 = unlist(runquantile(branchdf$Klf4Zfp281targetratio, k, probs = 0.25)),
    KZq75 = unlist(runquantile(branchdf$Klf4Zfp281targetratio, k, probs = 0.75)),
    sample = branchselect)
  
  return(branchdf_calc)
}

root <- calcbranch(plotdf, "Root")

ggplot(root) + 
  geom_line(aes(x = Pseudotime, y = KZ)) + 
  geom_ribbon(aes(x = Pseudotime, ymin = KZq25, ymax = KZq75), fill = "grey50", alpha = 0.3)


