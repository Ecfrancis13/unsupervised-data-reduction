
library(tidyverse)
library(psych)
library(ggpubr)
library(pheatmap)
#-------------------------------------------------------------------------------
# Author: Dr. Ellen C Francis
# Topic: unsupervised data analyses methods for identifying latent patterns
# We apply principal component analysis to derive latent variables that capture placental signaling pathways
#-------------------------------------------------------------------------------

# read in data prepared with script (0_dataPrep)
bdf1_scalelg<- read.csv(("bdf1_scalelg_rr.csv"))
bdf1_scalelg <- select(bdf1_scalelg, -"X")

#move id to dataframe
bdf1_scalelg <- column_to_rownames(bdf1_scalelg, var="ChildID")
#-------------------------------------------------------------------------------

#Rows = samples, columns = proteins.

#log transformed and scaled data are best for PCA. WE DO NOT USE THE internal
#FUNCTISN BECAUSE THE DATA ARE ALREADY log and scaled and centered
#scale function

pca_result <- prcomp(bdf1_scalelg, center = FALSE, scale. = FALSE)

#-------------------------------------------------------------------------------
#Decide how many PCs to retain 
#-------------------------------------------------------------------------------

# Compute Eigenvalues
eigenvalues <- pca_result$sdev^2
eigenvalues

# Compute Cumulative Variance Explained
cumulative_variance <- cumsum(eigenvalues) / sum(eigenvalues)
cumulative_variance

# Kaiser's Rule (PCs with Eigenvalue > 1)
kaiser_rule <- sum(eigenvalues > 1)

# Parallel Analysis (PA) - using log and scaled data to match the prcomp function
pa_result <- fa.parallel(bdf1_scalelg, fa = "pc", n.iter = 100, show.legend = FALSE, main = "Parallel Analysis")
parallel_analysis_n <- sum(eigenvalues > pa_result$fa.values)

print(pa_result$fa.values)
###
#Blue X's = actual eigenvalues from data. Red dotted line = average eigenvalues from simulated random data (null model)
# Horizontal black line at 1 = Kaiser rule (eigenvalue > 1)

#Interpretation: retain components to the left of where the blue crosses fall below the red line.
# we retain the first 9 components because their eigenvalues exceed the null distribution.
#After PC9, the data is indistinguishable from noise.
#we selected the number of principal components based on parallel analysis, which compares observed eigenvalues to those expected under random data.
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
#Extract loadings and scores
#-------------------------------------------------------------------------------
# Extract scores (Here we extract PCs 1:9 PCs)
pc_scores <- as.data.frame(pca_result$x[, 1:9])

# add sample IDs 
pc_scores$childID <- rownames(pca_result$x)
#write.csv(pc_scores, "pc_scores_df.csv")


# Extract the loadings (aka eigenvectors)
protein_loadings <- pca_result$rotation  # rows = proteins, cols = PCs

# Convert to a data frame for easy export
protein_loadings_df <- as.data.frame(protein_loadings)

# Optionally round for readability
protein_loadings_df <- round(protein_loadings_df, 2)


# Save to CSV
#write.csv(protein_loadings_df, "pca_protein_loadings.csv", row.names = TRUE)


#-------------------------------------------------------------------------------
#Visualizations
#-------------------------------------------------------------------------------
# Get top 5 proteins for each of the first 3 PCs
top_loadings <- function(pc_num, n = 5) {
  loadings <- pca_result$rotation[, pc_num]
  sorted <- sort(abs(loadings), decreasing = TRUE)
  top_proteins <- names(sorted)[1:n]
  data.frame(Protein = top_proteins,
             Loading = loadings[top_proteins],
             PC = paste0("PC", pc_num))
}

# Combine top contributors from multiple PCs
top_df <- do.call(rbind, lapply(1:9, top_loadings, n = 5))
print(top_df)


#label proteins

protein_labels <- c(
  "AKT" = "AKT",
  "AKT.308_ratio" = "AKT Thr308 (ratio)",
  "AKT.473_ratio" = "AKT Ser473 (ratio)",
  "AMPK" = "AMPK",
  "AMPK_ratio" = "AMPK (ratio)",
  "eiF2a" = "eIF2α",
  "eiF2a_ratio" = "eIF2α (ratio)",
  "ERK.1.2." = "ERK1/2",
  "ERK_ratio" = "ERK1/2 (ratio)",
  "GSK3_ratio" = "GSK3β (ratio)",
  "GSK3beta" = "GSK3β",
  "IGF1r" = "IGF1r",
  "IL.1beta" = "IL-1β",
  "IRBeta" = "IRβ",
  "JNK1" = "JNK1",
  "JNK1_ratio" = "JNK1 (ratio)",
  "JNK2" = "JNK2",
  "JNK2_ratio" = "JNK2 (ratio)",
  "OGT" = "OGT",
  "p38MAPK" = "p38MAPK",
  "p70S6K_ratio" = "p70S6K (ratio)",
  "PGC1.alpha" = "PGC1α",
  "PKCa" = "PKCα",
  "S6" = "RPS6",
  "S6_ratio_rerun" = "RPS6 (ratio)",
  "STAT3" = "STAT3",
  "STAT3_ratio" = "STAT3 (ratio)",
  "t_p70S6k" = "p70S6K",
  "tPI3K_p85a" = "PI3K-p85α",
  "X.Pro..Caspase.1" = "Pro-Caspase 1",
  "X11b.HSD2" = "11β-HSD2",
  "X4E.BP1" = "4EBP1",
  "X4EBP1_ratio" = "4EBP1 (ratio)"
)

# Step 1: Extract absolute loadings to identify top proteins
loadings_mat <- abs(pca_result$rotation[, 1:9])
top_protein_idx <- order(rowSums(loadings_mat), decreasing = TRUE)[1:20]

# Step 2: Subset original (signed) loadings to preserve sign
top_loadings <- pca_result$rotation[top_protein_idx, 1:9]

# Step 3: Apply custom protein name labels to rows
rownames(top_loadings) <- protein_labels[rownames(top_loadings)]

# Order by PC with strongest loading
dominant_pc <- apply(abs(top_loadings), 1, which.max)
top_loadings_ordered <- top_loadings[order(dominant_pc), ]

pheatmap(top_loadings_ordered,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         main = "Top Protein Loadings by PC")


