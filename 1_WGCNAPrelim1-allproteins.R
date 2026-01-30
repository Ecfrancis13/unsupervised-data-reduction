
library(tidyverse)
#BiocManager::install("WGCNA")
#BiocManager::install("GO.db")
#BiocManager::install("AnnotationDbi")
library(WGCNA)


#-------------------------------------------------------------------------------
# Author: Dr. Ellen C Francis
# Topic: unsupervised data analyses methods for identifying latent patterns
# We apply WGCNA to derive latent variables that capture placental signaling pathways
#-------------------------------------------------------------------------------


# read in data prepared with script (0_dataPrep)
bdf1_scalelg<- read.csv(("bdf1_scalelg_rr.csv"))
bdf1_scalelg <- select(bdf1_scalelg, -"X")

#move id to dataframe
bdf1_scalelg <- column_to_rownames(bdf1_scalelg, var="ChildID")

###############################################################################
# there are a few approaches to deriving a network. This script uses two, and at the end selects 
# what we end up using is on line 145
###############################################################################
#-------------------------------------------------------------------------------
datExpr<-bdf1_scalelg

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
meanExpressionByArray=apply(bdf1_scalelg,1,mean, na.rm=T)

barplot(meanExpressionByArray,
        xlab = "Sample", ylab = "Mean expression",
        main ="Mean expression across samples", cex.names = 0.7)
#-------------------------------------------------------------------------------


#Picking power
#-------------------------------------------------------------------------------
options(stringsAsFactors = FALSE);
# choose a set of soft thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by =2))
sft <- pickSoftThreshold(bdf1_scalelg,
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
  geom_vline(xintercept = 7, color = "red") + #add red line at chosen power
  labs(x= "Power", y="signed R^2") +
  theme_classic()

a2<- ggplot(sft.data, aes(Power, mean.k., label = Power)) + 
  geom_point() + 
  geom_text(nudge_y = 0.1) +   #move label above
  geom_vline(xintercept = 7, color = "red") + #add purple line at chosen power
  labs(x= "Power", y="Mean Connectivity") +
  theme_classic()

# you want to choose a plot that has a power that has an r^2 greater than .8 and has low mean connectivity 
#with log data the R^2 doesn't achieve scale free topology eg>8, but with just scaled data (no log) there
#is beter topology but interestingly no matter if the data are log scaled or just scaled the modules and network is the same.

grid.arrange(a1, a2, nrow = 2)
#-------------------------------------------------------------------------------

net <- blockwiseModules(
  bdf1_scalelg,
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
proteinNames <- colnames(bdf1_scalelg)[net$blockGenes[[1]]]

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
###############################################################################
# Assessing deep splits. SENSITIVITY ANALYSIS. THIS IS WHAT the study USED 
###############################################################################
# Calculate adjacency
adjacency <- adjacency(bdf1_scalelg, power = 7, type = "signed")

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
#If deepSplit = 2 and 1 produce the same module assignment, 
#this suggests wsame resolution limit of our data.

# we use deep split 3
ds3 <- cutreeHybrid(
  dendro = geneTree,
  distM = dissTOM,
  deepSplit = 3,
  pamStage = FALSE,
  cutHeight = 0.99,
  minClusterSize = 3
)

moduleColors_ds3 <- labels2colors(ds3$labels)
#-------------------------------------------------------------------------------
plotDendroAndColors(
  geneTree, 
  moduleColors_ds3,
  "Module Colors",
  dendroLabels = nice_proteinNames, # label proteins
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.02,
  cex.dendroLabels = 0.8
)
#-------------------------------------------------------------------------------


#------------------------------------------------------------------------
#We first will get the kME values, along with their associated p-values for 
#the network and will then output the resulting
#table to a file ("kMEtable1.csv").

#also will output ME scores for the samples for each module
#------------------------------------------------------------------------

MEList = moduleEigengenes(bdf1_scalelg, colors = moduleColors_ds3)
MEs = MEList$eigengenes

# Calculate module membership (kME) for each metabolite
geneModuleMembership = signedKME(bdf1_scalelg, MEs)
MMPvalue = corPvalueStudent(as.matrix(geneModuleMembership), nrow(bdf1_scalelg))
colnames(MMPvalue) = paste0("pval_", names(MEs))

# Extract gene (metabolite) names (use column names of bdf1_scalelg)
Gene = colnames(bdf1_scalelg)

# Create the full kME table
kMEtable = data.frame(Metabolite = Gene, 
                      Module = moduleColors_ds3,  # Correct module assignments
                      geneModuleMembership, 
                      MMPvalue)

# Save the full kME and p-value table to CSV
write.csv(kMEtable, "kMEtable_rr.csv", row.names = FALSE)

# These are your sample-level module scores
module_scores <- MEs

#  Add sample IDs 
module_scores$childID <- rownames(bdf1_scalelg)
write.csv(module_scores, "module_scores_df_rr.csv")
#-------------------------------------------------------------------------------


#assess module connectivity
#-------------------------------------------------------------------------------

intra_mod_conn <- intramodularConnectivity(adjacency, colors = moduleColors_ds3)
# Mean kWithin per module
mean_conn_by_module <- tapply(intra_mod_conn$kWithin, moduleColors_ds3, mean)
top_modules <- names(sort(mean_conn_by_module, decreasing = TRUE))[1:3]
print(top_modules) #   

#HUBS
#-------------------------------------------------------------------------------
# Initialize list to store top hub genes per module
chooseTopHubInEachModule(bdf1_scalelg, moduleColors_ds3,omitColors = "grey", 
  power = 7, 
  type = "signed")


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
write.csv(topGenesKME_df, "TopGenesKME_rr.csv", row.names = FALSE)
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