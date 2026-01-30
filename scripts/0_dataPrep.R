
#-------------------------------------------------------------------------------
# Author: Dr. Ellen C Francis
# Topic: unsupervised data analyses methods for identifying latent patterns
#-------------------------------------------------------------------------------

#libraries
library(naniar)
library(VIM)
library(finalfit)
library(tidyverse)

setwd('C:/Users/ef447/OneDrive - Rutgers University/pmomic/PJ1/Rfiles/Data')


#-------------------------------------------------------------------------------
# 1. DATA PREP: selecting protein data, investigate missing patterns, imputation, scaling
#-------------------------------------------------------------------------------
#create list of proteins that will be used, and a clean label

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


#read in data
p <- read.csv("placenta_raw.csv")
#move ids to row names, drop infant sex and other id vars
p_id <- column_to_rownames(p, var="PID")
p <- p %>% select(-PID) #have to use a dataframe without row names for missing investigation
#------------------------------------------------------------------------------
#missing
#------------------------------------------------------------------------------
### Investigate placentaRaw dataset. We do not select ratios because if total or phospyrlated protein is missing then raito is also going to be missing. 
placentaRawSub <- subset(p, select = -c(STAT3_ratio,JNK1_ratio,JNK2_ratio,GSK3_ratio,X4EBP1_ratio,p_p70S6K_389,AKT.473_ratio,AKT.308_ratio,eiF2a_ratio,p70S6K_ratio,AMPK_ratio,S6_ratio_rerun,pS6_rerun,ERK_ratio))
#LittleMCAR(placentaRawSub) #0.08
aggr(placentaRawSub[1:16], numbers = TRUE, prop = FALSE, sortVar = TRUE,
     cex.axis = 0.3) 
aggr(placentaRawSub[17:32], numbers = TRUE, prop = FALSE, sortVar = TRUE,
     cex.axis = 0.3) 
aggr(placentaRawSub, numbers = TRUE, prop = FALSE, sortVar = TRUE,
     cex.axis = 0.3)
mcar_test(placentaRawSub) #0.07
#- does not violate mcar assumption   
#- does not prove mcar, but more comfortable using complete case and not expecting it to be biased


#examines the pattern of missingness
md.pattern(p, rotate.names = TRUE)
#UPDATED THIS FIG WITH NICE NAMES USING "PROTEIN LABEL"
p_df<-p
names(p_df) <- protein_labels[names(p_df)]
md.pattern(p_df, rotate.names = TRUE)


# Total number of columns (variables)
n_cols <- ncol(p)

# Proportion of missingness for each row
row_missing_pct <- rowSums(is.na(p)) / n_cols

# Number of rows with >10% missing
n_rows_gt10pct <- sum(row_missing_pct > 0.1)

cat("Number of samples with >10% missing data:", n_rows_gt10pct, "\n")
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#Correlation among variables
#-------------------------------------------------------------------------------
## check for highly correlated variables. All are just correlted with themselves
cor_matrix <- cor(p, use = "pairwise.complete.obs")
findHighlyCorrelated <- function(mat, threshold = 0.999) {
  which(abs(mat) > threshold & abs(mat) < 1, arr.ind = TRUE)
}
findHighlyCorrelated(cor_matrix)

#-------------------------------------------------------------------------------
#make placenta protein data set based on complete case
#-------------------------------------------------------------------------------
#put ids back
completed_data <- as.data.frame(na.omit(p_id))
completed_data<-rownames_to_column(completed_data,var = "pid")


#write.csv(completed_data, "completed_data_rr.csv")
#completed_data<- read.csv(("completed_data_rr.csv"))
#completed_data <- select(completed_data, -"X")
#------------------------------------------------
#assign to new dataframe so there is a back up if needed
bdf1 <- completed_data

#select only the proteins of interest based on "my protein list"
bdf1_flt <- bdf1[, my_proteins]

#put id back back to row name for scaling
bdf1_flt<- column_to_rownames(bdf1_flt, var="pid")

bdf1_scale <- scale(bdf1_flt) #scale only
bdf1_log <- (log(bdf1_flt + 1))#log(x + 1) to handle zeros
#scales and centers the log transformed data
bdf1_scalelg <- bdf1_log %>% scale()

#-------------------------------------------------------------------------------
#Correlation among variables to be used after transformation and scaling
#-------------------------------------------------------------------------------
# calculating correlation of placental protein data 
corr_mat <- cor(bdf1_scalelg, method = "spearman")
corr_vals <- corr_mat[lower.tri(corr_mat)]  # Extract unique correlations
summary(abs(corr_vals))
#individual protein correlations
summary(abs(cor(bdf1_scalelg)))

#variance
apply(bdf1_scalelg, 2, var)
summary(apply(bdf1_scalelg, 2, var))
bdf1_scalelg <- as.data.frame(bdf1_scalelg)

#  add sample IDs
bdf1_scalelg <- rownames_to_column(bdf1_scalelg, "ChildID")

write.csv(bdf1_scalelg, "bdf1_scalelg_rr.csv")
bdf1_scalelg<- read.csv(("bdf1_scalelg_rr.csv"))
bdf1_scalelg <- select(bdf1_scalelg, -"X")
#--------------------------------------
