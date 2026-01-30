library(tidyverse)
#BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
#devtools::install_github("kevinblighe/CorLevelPlot")
library(CorLevelPlot)
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)



#-------------------------------------------------------------------------------
# Author: Dr. Ellen C Francis
# Topic: unsupervised data analyses methods for identifying latent patterns
# We apply consensus clustering to derive latent variables that capture placental signaling pathways
#-------------------------------------------------------------------------------

# read in data prepared with script (0_dataPrep)
bdf1_scalelg<- read.csv(("bdf1_scalelg_rr.csv"))
bdf1_scalelg <- select(bdf1_scalelg, -"X")

#move id to dataframe
bdf1_scalelg <- column_to_rownames(bdf1_scalelg, var="ChildID")


# Convert  data frame to a numeric matrix so it runs in consensus clustering
bdf1_scalelg_matrix <- data.matrix(bdf1_scalelg)

#Rows = samples, columns = proteins.
#using consensus clustering on log transformed and scaled protein data

results <- ConsensusClusterPlus(
  d = bdf1_scalelg_matrix,          # if t() then transpose so proteins are rows and ends up clustering ppl
  maxK = 6,
  reps = 1000,
  pItem = 0.8,
  pFeature = 1.0,
  clusterAlg = "hc",
  distance = "spearman",        # or "pearson"
  title = "consensus_results_rr",
  seed = 123,
  plot = "pdf"
)

#The largest gain in stability is between k = 3 and 4.
#Beyond k = 4, the delta drops sharply and then plateaus, suggesting limited benefit from more clusters.

#Delta area represents the relative increase in the cumulative distribution function (CDF) between successive values of K.
#After K = 4, the increase in stability drops off steeply, suggesting diminishing returns in partitioning the data further.



#-------------------------------------------------------------------------------
# Extract consensus matrix and class labels for K=4
cons_matrix_p <- results[[4]]$consensusMatrix
class_labels_p <- as.character(results[[4]]$consensusClass)

#------------------------------------------------------------------------------
#extract clusters
protein_clusters <- results[[4]]$consensusClass
protein_clusters_df <- data.frame(
  Protein = colnames(bdf1_scalelg_matrix),  # assuming proteins are columns in bdf1_scale
  Cluster = factor(protein_clusters)
)

#save protein clusters
#write.csv (protein_clusters_df, "consen-clust-rr.csv")
#protein_clusters_df <- read.csv("consen-clust-rr.csv")


#examine clusters
split(protein_clusters_df$Protein, protein_clusters_df$Cluster)

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
cluster_centroid_df$childID <- rownames(bdf1_scalelg)
write.csv(cluster_centroid_df, "cluster_centroid_df_rr.csv")

# View the result
head(cluster_centroid_df)



#------------------------------------------------------------------------------
#visualizations
#------------------------------------------------------------------------------

#  add sample IDs
bdf1_scalelg <- rownames_to_column(bdf1_scalelg, "ChildID")


bdf_long <- reshape2::melt(bdf1_scalelg)
colnames(bdf_long) <- c("Sample", "Protein", "value")
# Add cluster assignment to each protein
bdf_long$Cluster <- protein_clusters[as.character(bdf_long$Protein)]

bdf_long$Sample <- factor(bdf_long$Sample)
bdf_long$Cluster <- factor(bdf_long$Cluster)


# Plot protein expression profiles grouped by cluster
ggplot(bdf_long, aes(x = Sample, y = value, group = Protein, color = Cluster)) +
  geom_line(alpha = 0.8) +
  facet_wrap(~ Cluster, scales = "free_y") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Protein Expression Profiles by Consensus Cluster",
    y = "Scaled Protein Expression",
    x = "Sample"
  )
#------------------------------------------------------------------------------
## HEAT MAP VISUAL
#------------------------------------------------------------------------------

# Create row and column annotations
##  names and order.
cluster_map <- c(
  `3` = "Cluster 3",
  `2` = "Cluster 2",
  `1` = "Cluster 1",
  `4` = "Cluster 4"
)
order <- c(
  "Cluster 3",
  "Cluster 2",
  "Cluster 1",
  "Cluster 4"
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
rownames(cons_matrix_p) <- colnames(cons_matrix_p) <- colnames(bdf1_scalelg_matrix)

# Subset and reorder the consensus matrix and labels
cons_matrix_sorted <- cons_matrix_p[sorted_proteins, sorted_proteins]
#index protein labels from nice labels
protein_labels_sorted <- protein_labels[sorted_proteins]

## helper: vector of cluster names aligned to sorted_proteins
sorted_cluster_names <- protein_clusters_df$Cluster_name[match(sorted_proteins, protein_clusters_df$Protein)]

##--- 4) Nice, named colors for the four biological clusters ------------
cluster_colors <- c(
  "mTOR/\nAMPK"              = "#B2182B",
  "Ins. Coord./\nMito. Biogenesis"   = "#2166AC",
  "Placental Sig Coord."  = "#4DAF4A",
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
  row_title   = NULL,     # remove y-axis cluster titles
  #column_title_rot = 45,  # rotate top cluster titles 450
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

