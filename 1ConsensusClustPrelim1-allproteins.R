library(haven)
library(tidyverse)
library(dplyr)
library(mice)
library(rstatix)
library(tibble)
library(CorLevelPlot)
library(rstatix)
library(reshape2)
library(ggstatsplot)
library(gridExtra)
library(cluster)
library(factoextra)
library(ConsensusClusterPlus)



#-------------------------------------------------------------------------------
# Author: Dr. Ellen C Francis
# Topic: comparison of unsupervised data analyses methods for molecular biologists
# 
# The goal is to compare WGCNA, PCA, and consensus clustering (beta div is supervised) 
# .	We do this because we want to know whether there are other pathways working
#      indirectly with insulin growth factor signaling that we wouldn't have thought 
#      of using our biological knowledge alone. 
# .	To compare the utility of each multivariate framework, we used regression-based 
#      approaches with derived features from each method (principal components from PCA, 
#      module eigengenes from WGCNA, and ordination axes from beta diversity analyses) 
#      as predictors of continuous metabolic outcomes in childhood. We evaluated model 
#      performance using Akaike Information Criterion (AIC) to assess relative model fit 
#      while accounting for differences in model complexity.
# .	Then we pick the method with the best AIC and use it in association with childhood outcomes
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# 1. DATA PREP: selecting protein data, imputation, scaling
#-------------------------------------------------------------------------------

p <- read.csv("placenta_raw.csv")

#-------------------------------------------------------------------------------
#examines the pattern of missingness
md.pattern(p, rotate.names = TRUE)

# Total number of columns (variables)
n_cols <- ncol(p)

# Proportion of missingness for each row
row_missing_pct <- rowSums(is.na(p)) / n_cols

# Number of rows with >10% missing
n_rows_gt10pct <- sum(row_missing_pct > 0.1)

cat("Number of samples with >10% missing data:", n_rows_gt10pct, "\n")

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
## check for highly correlated variables. 
cor_matrix <- cor(p, use = "pairwise.complete.obs")
findHighlyCorrelated <- function(mat, threshold = 0.999) {
  which(abs(mat) > threshold & abs(mat) < 1, arr.ind = TRUE)
}
findHighlyCorrelated(cor_matrix)

#exclude ratios from prediction
init <- mice(p, maxit = 0)
pred <- init$predictorMatrix

# Optionally: exclude ratio variables from predicting others
#ratios <- grep("ratio", names(p), value = TRUE)
#pred[, ratios] <- 0  # prevent using ratio variables to impute others
#imp <- mice(p, 
#            method = "cart", 
#            predictorMatrix = pred, 
#            ridge = 1e-5, 
#            m = 5, 
#            seed = 123)

#compares PPM we only have continuous variables so we use predictive mean matching (PPM)
imp <- mice(data = p, 
            method = "cart", 
            m = 5, 
            seed = 123)

completed_data <- complete(imp, 1)
##### Plot of the distributions of the clinical outcomes after applying MICE
stripplot(imp, pch = 20, cex = 1.2)  # visualization of imputed vs. observed

#put ids back
completed_data$pid <- row.names(p)


#-------------------------------------------------------------------------------
#assign to new dataframe so i have a clean one to go back to if needed
bdf1 <- completed_data

#restrict to only total and ratios
my_proteins <- c(
  "AKT", "AKT.308_ratio", "AKT.473_ratio", "AMPK", "AMPK_ratio",
  "eiF2a", "eiF2a_ratio", "ERK.1.2.", "ERK_ratio",
  "GSK3_ratio", "GSK3beta", "IGF1r", "IL.1beta", "IRBeta",
  "JNK1", "JNK1_ratio", "JNK2", "JNK2_ratio", "OGT", "p38MAPK",
  "p70S6K_ratio", "PGC1.alpha", "PKCa", "S6", "S6_ratio_rerun",
  "STAT3", "STAT3_ratio", "t_p70S6k", "tPI3K_p85a",
  "X.Pro..Caspase.1", "X11b.HSD2", "X4E.BP1", "X4EBP1_ratio", "pid"
)


protein_labels <- c(
  "AKT" = "AKT",
  "AKT.308_ratio" = "AKT Thr308 (ratio)",
  "AKT.473_ratio" = "AKT Ser473 (ratio)",
  "AMPK" = "AMPK",
  "AMPK_ratio" = "AMPK (ratio)",
  "eiF2a" = "eIF2a",
  "eiF2a_ratio" = "eIF2a (ratio)",
  "ERK.1.2." = "ERK1/2",
  "ERK_ratio" = "ERK1/2 (ratio)",
  "GSK3_ratio" = "GSK3 (ratio)",
  "GSK3beta" = "GSK3??",
  "IGF1r" = "IGF1r",
  "IL.1beta" = "IL-1??",
  "IRBeta" = "IR??",
  "JNK1" = "JNK1",
  "JNK1_ratio" = "JNK1 (ratio)",
  "JNK2" = "JNK2",
  "JNK2_ratio" = "JNK2 (ratio)",
  "OGT" = "OGT",
  "p38MAPK" = "p38MAPK",
  "p70S6K_ratio" = "p70S6K (ratio)",
  "PGC1.alpha" = "PGC1a",
  "PKCa" = "PKCÎa",
  "S6" = "RPS6",
  "S6_ratio_rerun" = "RPS6 (ratio)",
  "STAT3" = "STAT3",
  "STAT3_ratio" = "STAT3 (ratio)",
  "t_p70S6k" = "p70S6K",
  "tPI3K_p85a" = "PI3K-p85a",
  "X.Pro..Caspase.1" = "Pro-Caspase 1",
  "X11b.HSD2" = "11??-HSD2",
  "X4E.BP1" = "4EBP1",
  "X4EBP1_ratio" = "4EBP1 (ratio)"
)


bdf1_flt <- bdf1[, my_proteins]


# 
bdf1_flt<- column_to_rownames(bdf1_flt, var="pid")

bdf1_scale <- scale(bdf1_flt) #scale only
bdf1_log <- (log(bdf1_flt + 1))#log only

# calculating correlation of placental protein data 
corr_mat <- cor(bdf1_scale, method = "spearman")
corr_vals <- corr_mat[lower.tri(corr_mat)]  # Extract unique correlations
summary(abs(corr_vals))
#individual protein correlations
summary(abs(cor(bdf1_scale)))

#variance
apply(bdf1_scale, 2, var)
summary(apply(bdf1_scale, 2, var))
#-------------------------------------------------------------------------------


#Rows = samples, columns = proteins.
#Use raw or imputed values, not scaled (unless using euclidean dist.)
#using consensus clustering on scaled protein data

results <- ConsensusClusterPlus(
  d = bdf1_scale,          # if t() then transpose so proteins are rows and ends up clustering ppl
  maxK = 6,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1.0,
  clusterAlg = "hc",
  distance = "spearman",        # or "pearson"
  seed = 123,
  plot = "pdf"
)

#NEW NOTE: k = 4 using cumulative distribution function to assess gains in stability, 
# k2-k4 show increases. K3 resulted in 1 large cluster and 2 small. so used k4 

#Delta area represents the relative increase in the cumulative distribution function (CDF) between successive values of K.
#After K = 4, the increase in stability drops off steeply, suggesting diminishing returns in partitioning the data further.

#visualization
#-------------------------------------------------------------------------------
# Extract consensus matrix and class labels for K =4
cons_matrix_p <- results[[4]]$consensusMatrix
class_labels_p <- as.character(results[[4]]$consensusClass)

#------------------------------------------------------------------------------
#extract clusters
protein_clusters <- results[[4]]$consensusClass
protein_clusters_df <- data.frame(
  Protein = colnames(bdf1_scale),  # assuming proteins are columns in bdf1_scale
  Cluster = factor(protein_clusters)
)
write.csv (protein_clusters_df, "consen-clust.csv")
protein_clusters_df <- read.csv("consen-clust.csv")
#examine clusters
split(protein_clusters_df$Protein, protein_clusters_df$Cluster)

# Create row and column annotations
library(ComplexHeatmap)
library(circlize)

## new names and order
cluster_map <- c(
  `2` = "mTOR/\nAMPK",
  `4` = "IGF/Mito.\nBiogenesis",
  `1` = "Placental Ins. Coord.",
  `3` = "Infl./Stress"
)
order <- c(
  "mTOR/\nAMPK",
  "IGF/Mito.\nBiogenesis",
  "Placental Ins. Coord.",
  "Infl./Stress"
)
protein_clusters_df <- protein_clusters_df %>%
  mutate(
    Cluster_num  = as.character(Cluster),
    Cluster_name = recode(Cluster_num, !!!cluster_map),
    Cluster_name = factor(Cluster_name, levels = order)
  )

# Sort proteins by cluster assignment
sorted_proteins <- protein_clusters_df %>%
  arrange(Cluster) %>%
  pull(Protein)

# assign names to consensus clustering matrix rows and columns
rownames(cons_matrix_p) <- colnames(cons_matrix_p) <- colnames(bdf1_scale)

# Subset and reorder the consensus matrix and labels
cons_matrix_sorted <- cons_matrix_p[sorted_proteins, sorted_proteins]
#index protein labels from nice labels
protein_labels_sorted <- protein_labels[sorted_proteins]

## helper: vector of cluster names aligned to sorted_proteins
sorted_cluster_names <- protein_clusters_df$Cluster_name[match(sorted_proteins, protein_clusters_df$Protein)]

##--- 4) Nice, named colors for the four biological clusters ------------
cluster_colors <- c(
  "mTOR/\nAMPK"              = "#B2182B",
  "IGF/Mito.\nBiogenesis"   = "#2166AC",
  "Placental Ins. Coord."  = "#4DAF4A",
  "Infl./Stress"           = "#984EA3"
)


ha_row <- rowAnnotation(
  Cluster = sorted_cluster_names,
  col = list(Cluster = cluster_colors),
  show_legend = FALSE
  )

ha_col <- HeatmapAnnotation(
  Cluster = sorted_cluster_names,
  col = list(Cluster = cluster_colors),
  show_legend = FALSE
)

##--- 5) Plot -----------------------------------------------------------
split_vec <- factor(sorted_cluster_names, levels = order)

Heatmap(
  cons_matrix_sorted,
  name = "Co-clustering\nproportion",
  col = colorRamp2(c(0, 0.5, 1), c("white", "#bdd7e7", "#08306b")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  ## >>> force block order and keep it <<<
  row_split = split_vec,
  column_split = split_vec,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_labels = protein_labels_sorted,
  column_labels = protein_labels_sorted,
  left_annotation = ha_row,
  top_annotation  = ha_col,
  row_title   = NULL,     # remove y-axis cluster titles
  #column_title_rot = 45,  # rotate top cluster titles 45°
  # optional: styling for top titles
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  rect_gp = gpar(col = "white"),
  heatmap_legend_param = list(
    legend_direction = "horizontal",
    title_position = "topcenter",
    color_bar = "continuous",
    legend_width = unit(4, "cm"),
    legend_height = unit(4, "mm")
  ),
  border = TRUE
)

## draw without merging legends (not required, but explicit)
raw(ht, heatmap_legend_side = "bottom")   # only the heatmap legend will show


#------------------------------------------------------------------------------

#### next step cluster centroid score of each cluster
#------------------------------------------------------------------------------
# Split proteins into lists by cluster
cluster_proteins <- split(protein_clusters_df$Protein, protein_clusters_df$Cluster)

# For each cluster, calculate mean expression per sample
cluster_centroids <- lapply(cluster_proteins, function(protein_list) {
  rowMeans(bdf1_scale[, protein_list, drop = FALSE], na.rm = TRUE)
})

# Combine into a dataframe
cluster_centroid_df <- as.data.frame(cluster_centroids)
colnames(cluster_centroid_df) <- paste0("CC", names(cluster_centroids))

#  add sample IDs
cluster_centroid_df$childID <- rownames(bdf1_scale)
write.csv(cluster_centroid_df, "cluster_centroid_df.csv")

# View the result
head(cluster_centroid_df)





#-------------------------------------------------------------------------------
##this results in clustering of samples if you run a transposed 
results_t <- ConsensusClusterPlus(
  d = t(bdf1_scale),          # if t() then transpose so proteins are rows and ends up clustering ppl
  maxK = 6,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1.0,
  clusterAlg = "hc",
  distance = "spearman",        # or "pearson"
  seed = 123,
  plot = "pdf"
)
cons_matrix <- results_t[[4]]$consensusMatrix
cluster_labels <- results_t[[4]]$consensusClass

# Assign row/column names as protein names if not already assigned
rownames(cons_matrix) <- colnames(cons_matrix) <- names(cluster_labels)

# Create annotation
row_anno <- HeatmapAnnotation(
  Cluster = factor(cluster_labels),
  col = list(Cluster = structure(RColorBrewer::brewer.pal(3, "Set1"),
                                 names = levels(factor(cluster_labels)))),
  which = "row"
)

col_anno <- HeatmapAnnotation(
  Cluster = factor(cluster_labels),
  col = list(Cluster = structure(RColorBrewer::brewer.pal(3, "Set1"),
                                 names = levels(factor(cluster_labels)))),
  which = "column"
)

# Plot
Heatmap(
  matrix = cons_matrix,
  name = "Consensus\nIndex",
  show_row_names = TRUE,
  show_column_names = FALSE,
  top_annotation = col_anno,
  left_annotation = row_anno,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 8)
)
print(hplot)

#------------------------------------------------------------------------------