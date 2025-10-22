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
pdf1 <- completed_data

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

pdf1_flt <- pdf1[, my_proteins]


# 
pdf1_flt<- column_to_rownames(pdf1_flt, var="pid")

pdf1_scale <- (scale(pdf1_flt)) #scale only
pdf1_log <- (log(pdf1_flt + 1))#log only

pdf1_scalelg <- (scale(pdf1_log))

#-------------------------------------------------------------------------------


#Rows = samples, columns = proteins.
#log transformed and scaled data are best for PCA. USE log data and internal
#scale function

pca_result <- prcomp(pdf1_log, center = TRUE, scale. = TRUE)

#Decide how many PCs to retain 
library(psych)
library(ggpubr)
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
pa_result <- fa.parallel(pdf1_scalelg, fa = "pc", n.iter = 100, show.legend = FALSE, main = "Parallel Analysis")
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
#visualization
#-------------------------------------------------------------------------------


# Variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

df <- data.frame(
  PC = 1:length(var_explained),
  VarianceExplained = var_explained
)

ggplot(df[1:10, ], aes(x = PC, y = VarianceExplained)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 9, linetype = "dashed", color = "red") +
  annotate("text", x = 9.5, y = var_explained[9],
           label = "PC 9", color = "red", hjust = 0) +
  labs(
    title = "Scree Plot with Parallel Analysis Cutoff",
    x = "Principal Component",
    y = "Proportion of Variance Explained"
  ) +
  theme_minimal()

#-------------------------------------------------------------------------------
#Extract loadings and scores
#-------------------------------------------------------------------------------
# Extract scores (samples × PCs)
pc_scores <- as.data.frame(pca_result$x[, 1:9])

# add sample IDs 
pc_scores$childID <- rownames(pca_result$x)
write.csv(pc_scores, "pc_scores_df.csv")


# Extract the loadings (aka eigenvectors)
protein_loadings <- pca_result$rotation  # rows = proteins, cols = PCs

# Convert to a data frame for easy export
protein_loadings_df <- as.data.frame(protein_loadings)

# Optionally round for readability
protein_loadings_df <- round(protein_loadings_df, 2)


# Save to CSV
write.csv(protein_loadings_df, "pca_protein_loadings.csv", row.names = TRUE)


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

library(pheatmap)
#label proteins
protein_labels <- c(
  "IRBeta" = "IR??",
  "PGC1.alpha" = "PGC1??",
  "X4E.BP1" = "4EBP1",
  "STAT3_ratio" = "STAT3 (ratio)",
  "eiF2a" = "eIF2??",
  "AKT" = "AKT",
  "t_p70S6k" = "p70S6K",
  "tPI3K_p85a" = "PI3K-p85??",
  "S6" = "S6",
  "JNK1" = "JNK1",
  "p38MAPK" = "p38MAPK",
  "AKT.473_ratio" = "AKT Ser473 (ratio)",
  "AKT.308_ratio" = "AKT Thr308 (ratio)",
  "OGT" = "OGT",
  "S6_ratio_rerun" = "S6 (ratio)",
  "AMPK_ratio" = "AMPK (ratio)",
  "X4EBP1_ratio" = "4EBP1 (ratio)",
  "JNK1_ratio" = "JNK1 (ratio)",
  "ERK_ratio" = "ERK1/2 (ratio)",
  "eiF2a_ratio" = "eIF2?? (ratio)"
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

library(ComplexHeatmap)
library(circlize)

Heatmap(top_loadings_ordered,
        name = "Loading",
        col = colorRamp2(c(-0.5, 0, 0.5), 
                         c("#4575b4", "#ffffff00", "#d73027")),  # white with 0 alpha at 0
        show_row_names = TRUE,
        show_column_names = TRUE)

