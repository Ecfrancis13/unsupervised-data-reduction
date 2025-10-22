library(haven)
library(tidyverse)
library(dplyr)
library(mice)
library(rstatix)
library(WGCNA)
library(tibble)
library(CorLevelPlot)
library(rstatix)
library(reshape2)
library(ggstatsplot)
library(gridExtra)
library(flashClust)




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
wdf1 <- completed_data

#restrict to only total and ratios
my_proteins <- c(
  "AKT", "AKT.308_ratio", "AKT.473_ratio", "AMPK", "AMPK_ratio",
  "eiF2a", "eiF2a_ratio", "ERK.1.2.", "ERK_ratio",
  "GSK3_ratio", "GSK3beta", "IGF1r", "IL.1beta", "IRBeta",
  "JNK1", "JNK1_ratio", "JNK2", "JNK2_ratio", "OGT", "p38MAPK",
  "p70S6K_ratio", "PGC1.alpha", "PKCa", "S6", "S6_ratio_rerun",
  "STAT3", "STAT3_ratio", "t_p70S6k", "tPI3K_p85a",
  "X.Pro..Caspase.1", "X11b.HSD2", "X4E.BP1", "X4EBP1_ratio","pid"
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
  "PGC1.alpha" = "PGC1??",
  "PKCa" = "PKC??",
  "S6" = "S6",
  "S6_ratio_rerun" = "S6 (ratio)",
  "STAT3" = "STAT3",
  "STAT3_ratio" = "STAT3 (ratio)",
  "t_p70S6k" = "p70S6K",
  "tPI3K_p85a" = "PI3K-p85??",
  "X.Pro..Caspase.1" = "Pro-Caspase 1",
  "X11b.HSD2" = "11??-HSD2",
  "X4E.BP1" = "4EBP1",
  "X4EBP1_ratio" = "4EBP1 (ratio)"
)

wdf1_flt <- wdf1[, my_proteins]

wdf1_flt<- column_to_rownames(wdf1_flt, var="pid")


# 
wdf1_scalelg <- wdf1_flt %>%
  mutate(across(1:17, ~ log(. + 1))) %>%   # log transform with log(x + 1) to handle zeros
  mutate(across(1:17, scale))              # scale the log-transformed columns

wdf1_scale <- (scale(wdf1_flt)) #scale only
wdf1_log <- (log(wdf1_flt + 1))#log only


#-------------------------------------------------------------------------------
datExpr<-wdf1_scale

# WGCNA expects samples in rows and network data, eg placenta proteins, in columns. 
# this conforms to a typical dataframe in a biological experiment.

gsg <- goodSamplesGenes(datExpr)
summary(gsg)
gsg$allOK
#look to see how many are outliers. All proteins and samples ok.
table(gsg$goodGenes)
table(gsg$goodSamples)

#detect outliers using hierarchical clustering 

htree <- hclust(dist(datExpr), method ="average")
plot(htree) # looks good.
#-------------------------------------------------------------------------------


#Mean expression across samples
#-------------------------------------------------------------------------------
meanExpressionByArray=apply(wdf1,1,mean, na.rm=T)

barplot(meanExpressionByArray,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples", cex.names = 0.7)
#-------------------------------------------------------------------------------


#Picking power
#-------------------------------------------------------------------------------
options(stringsAsFactors = FALSE);
# choose a set of soft thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by =2))
sft <- pickSoftThreshold(wdf1_scale,
                         powerVector = power, 
                         networkType = "signed",
                         verbose = 5)
#-------------------------------------------------------------------------------


#visualization to pick power. you want a power that gives the maximum Rsq (SFT.R.sq) and min mean connectivity (mean.k)
#-------------------------------------------------------------------------------
sft.data <- sft$fitIndices

a1<- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() + 
  geom_text(nudge_y = 0.1) +   #move label above
  geom_hline(yintercept = 0.7, color = "navy") + #add navy line at .8 power
  geom_vline(xintercept = 8, color = "red") + #add red line at chosen power
  labs(x= "Power", y="signed R^2") +
  theme_classic()

a2<- ggplot(sft.data, aes(Power, mean.k., label = Power)) + 
  geom_point() + 
  geom_text(nudge_y = 0.1) +   #move label above
  geom_vline(xintercept = 7, color = "red") + #add purple line at chosen power
  labs(x= "Power", y="Mean Connectivity") +
  theme_classic()

# you want to choose a plot that has a power that has an r^2 greater than .8 and has low mean connectivity 
# with ratios, raw normalized data OR scaled (zscore) data 5 is best power
#     WITH UNSIGNED NETWORK.
# with ratios AND total protein vars with zscore [wdf1_scaleflt], 7 is best power R2=.7 
#     WITH SIGNED NETWORK. 

grid.arrange(a1, a2, nrow = 2)
#-------------------------------------------------------------------------------

net <- blockwiseModules(
  wdf1_scale,
  power = 7,  # the power you choose
  networkType = "signed",  # or "unsigned" / "signed"
  TOMType = "signed",             # use "signed" or "unsigned" to match above
  corFnc = "cor",
  corOptions = list(method = "spearman"),
  minModuleSize = 3,              # small value since we only have 18 proteins
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  verbose = 3
)

#table of module size and assign colors back
#-------------------------------------------------------------------------------

# Convert numeric labels to color names (optional)
moduleColors <- labels2colors(net$colors)
table(moduleColors)

# Get protein names in dendrogram order
proteinNames <- colnames(wdf1_scale)[net$blockGenes[[1]]]

nice_proteinNames <- protein_labels[proteinNames]

#-------------------------------------------------------------------------------



# Plot the dendrogram with module colors
#-------------------------------------------------------------------------------
plotDendroAndColors(
  dendro = net$dendrograms[[1]],               # one block, so use [[1]]. If doing in larger datasets (genes) pick block run
  colors = moduleColors[net$blockGenes[[1]]],  # match module colors to dendrogram order
  groupLabels = "Module Colors",
  dendroLabels = nice_proteinNames,         # label proteins
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.2,
  cex.dendroLabels = 0.8               # Adjust font size
)

# Assessing deep splits. SENSITIVITY ANALYSIS
#-------------------------------------------------------------------------------
# Calculate adjacency
adjacency <- adjacency(wdf1_scale, power = 7, type = "signed")

# Calculate TOM
TOM <- TOMsimilarity(adjacency, TOMType = "signed")
dissTOM <- 1 - TOM

#extract and assign the dissimilarity matrix
geneTree <- hclust(as.dist(dissTOM), method = "average")

mColorh <- NULL
for (ds in 0:3) {
  tree <- cutreeHybrid(
    dendro = geneTree,
    distM = dissTOM,
    deepSplit = ds, #us diff deep splits
    pamStage = FALSE,
    cutHeight = 0.99,
    minClusterSize = 3
  )
  mColorh <- cbind(mColorh, labels2colors(tree$labels))
}
colnames(mColorh) <- paste0("ds_", 0:3)
plotDendroAndColors(geneTree, 
                    mColorh, 
                    paste("dpSplt =",0:3), 
                    main = "",
                    dendroLabels = proteinNames,         # label proteins
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    cex.dendroLabels = 0.6)
#If deepSplit = 2 and 3 produce the same module assignment, 
#this suggests we're already at the natural resolution limit of our data.
#-------------------------------------------------------------------------------


#------------------------------------------------------------------------
#We first will get the kME values, along with their associated p-values for 
#the network and will then output the resulting
#table to a file ("kMEtable1.csv").

#also will output ME scores for the samples for each module
#------------------------------------------------------------------------

MEList = moduleEigengenes(wdf1_scale, colors = moduleColors)
MEs = MEList$eigengenes

# Calculate module membership (kME) for each metabolite
geneModuleMembership = signedKME(wdf1_scale, MEs)
MMPvalue = corPvalueStudent(as.matrix(geneModuleMembership), nrow(wdf1_scale))
colnames(MMPvalue) = paste0("pval_", names(MEs))

# Extract gene (metabolite) names (use column names of wdf1_scale)
Gene = colnames(wdf1_scale)

# Create the full kME table
kMEtable = data.frame(Metabolite = Gene, 
                      Module = moduleColors,  # Correct module assignments
                      geneModuleMembership, 
                      MMPvalue)

# Save the full kME and p-value table to CSV
write.csv(kMEtable, "kMEtable.csv", row.names = FALSE)

# These are your sample-level module scores
module_scores <- MEs

#  Add sample IDs 
module_scores$childID <- rownames(wdf1_scale)
write.csv(module_scores, "module_scores_df.csv")
#-------------------------------------------------------------------------------


#assess module connectivity
#-------------------------------------------------------------------------------

intra_mod_conn <- intramodularConnectivity(adjacency, colors = moduleColors)
# Mean kWithin per module
mean_conn_by_module <- tapply(intra_mod_conn$kWithin, moduleColors, mean)
top_modules <- names(sort(mean_conn_by_module, decreasing = TRUE))[1:3]
print(top_modules) #"brown"     "turquoise" "blue"   

#HUBS
#-------------------------------------------------------------------------------
# Initialize list to store top hub genes per module
topGenesKME = list()

# Extract module names from the kME matrix
module_names = colnames(geneModuleMembership)

# Loop through each module to identify top hub genes
for (module in module_names) {
  
  # Rank genes within the module (higher kME is better for hubs)
  kMErank = rank(-geneModuleMembership[, module])
  
  # Extract the top 10 hub genes
  topGenes = Gene[kMErank <= 10]
  
  # Store in the list
  module_name_cleaned = sub("kME_", "", module)
  topGenesKME[[module_name_cleaned]] = topGenes
}

# Convert the list to a data frame for easier viewing
topGenesKME_df = as.data.frame(topGenesKME)

# Save the top hub genes to a CSV file
write.csv(topGenesKME_df, "TopGenesKME.csv", row.names = FALSE)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

#PROTEIN LEVEL/INTENSITY
#-------------------------------------------------------------------------------
# Loop through each module
for (module in unique(moduleColors)) {
  
  # Extract the ME for the current module
  ME_values = MEs[, paste0("ME", module)]
  
  # Order the samples by ME values
  ordered_indices = order(ME_values)
  ordered_ME_values = ME_values[ordered_indices]
  
  # Create the bar plot
  bp = barplot(ordered_ME_values, 
               col = module, 
               main = paste("Module:", module), 
               cex.main = 0.8, 
               ylab = "Eigengene Expression", 
               xlab = "Sample", 
               border = NA, 
               las = 2)
  
  # Calculate the 25th, 50th, and 75th percentiles
  n_samples = length(ordered_ME_values)
  quantiles = quantile(1:n_samples, probs = c(0.25, 0.5, 0.75))
  
  # Add vertical lines at the quantiles
  abline(v = bp[quantiles], col = "black", lty = 2, lwd = 2)
  
  # Add percentile labels
  text(x = bp[quantiles], y = max(ordered_ME_values) * 0.95, 
       labels = c("25%", "50%", "75%"), 
       pos = 3, cex = 0.7, col = "black")
}