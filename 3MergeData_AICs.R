#-------------------------------------------------------------------------------
# Author: Dr. Ellen C Francis
# 
# Here we merge the data sets from consensus clustering, wgcna, pca, and meta
# then we build our models and use AIC to pick the data reduction method?
#-------------------------------------------------------------------------------

library(tidyverse)

metaSel<- read.csv("metaSelect.csv") %>% select(-"X") #meta data

module_scores<- read.csv("module_scores_df.csv") %>% select(-"X")
pca_scores<- read.csv("pc_scores_df.csv") %>% select(-"X")
clust_scores<- read.csv("cluster_centroid_df.csv") %>% select(-"X")

df1 <- merge(module_scores, pca_scores, by="childID")
df2 <- merge(df1, clust_scores, by="childID")
df_merged<- right_join(metaSel, df2, by="childID")


#drop if missing birthweight because this person dropped out of study
#-------------------------------------------------------------------------------

df_merged<- filter(df_merged, !is.na(wtkgIPV3))#108 observations
#-------------------------------------------------------------------------------

#compare AIC for birthweight, fmp at birth, and childweight using  clusters, PCs, and Modules
#-------------------------------------------------------------------------------
################################################################################
#BIRTH  DATA
#
#
#1. PCA-based model (with PC1 to PC9)
bmodel_pca <- lm(fmpIPV3 ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + infant_sex + GA_at_birth_in_days, data = df_merged)

#2. WGCNA-based model (with module eigengenes)
bmodel_wgcna <- lm(fmpIPV3 ~ MEblue + MEturquoise + MEbrown + MEyellow + MEgrey + infant_sex + GA_at_birth_in_days, data = df_merged)

#3. Consensus clustering model (with cluster centroid scores)
bmodel_cc <- lm(fmpIPV3 ~ CC1 + CC2 + CC3 + CC4 + infant_sex + GA_at_birth_in_days, data = df_merged)

#4. Compare AIC across models
bwt_aic <-AIC(bmodel_pca, bmodel_wgcna, bmodel_cc)

#5. Compare R2 across models
# List your models
bmodel_list <- list(
  PCA = bmodel_pca,
  WGCNA = bmodel_wgcna,
  ConsensusClustering = bmodel_cc
)

# Extract R² and adjusted R²
bmodel_r2_df <- lapply(bmodel_list, function(mod) {
  summ <- summary(mod)
  data.frame(
    R2 = summ$r.squared,
    Adj_R2 = summ$adj.r.squared
  )
}) %>% bind_rows(.id = "Method")

# Calculate RMSE and MAE for birth models
birth_rmse_mae <- lapply(bmodel_list, function(mod) {
  obs <- mod$model$fmpIPV3
  pred <- predict(mod)
  rmse <- sqrt(mean((obs - pred)^2))
  mae <- mean(abs(obs - pred))
  data.frame(RMSE = rmse, MAE = mae)
}) %>% bind_rows(.id = "Method")

# Combine R2 and RMSE/MAE into one table
birth_model_metrics <- left_join(bmodel_r2_df, birth_rmse_mae, by = "Method")
print(birth_model_metrics)


#testing differences in models. AIC cant test directly unless a nested model
#this gives a probabilistic interpretation of which is better
aic_df <- AIC(bmodel_pca, bmodel_wgcna, bmodel_cc)
aic_df$delta_AIC <- aic_df$AIC - min(aic_df$AIC)
aic_df$akaike_weight <- exp(-0.5 * aic_df$delta_AIC) / sum(exp(-0.5 * aic_df$delta_AIC))
print(aic_df)


## bootstrap to see if there is a difference in RSME 
set.seed(123)
n_boot <- 1000

boot_rmse <- function(mod) {
  obs <- mod$model$fmpIPV3
  pred <- predict(mod)
  replicate(n_boot, {
    idx <- sample(seq_along(obs), replace = TRUE)
    sqrt(mean((obs[idx] - pred[idx])^2))
  })
}

booted_rmses <- lapply(bmodel_list, boot_rmse)

# Convert to data.frame for t-tests or plotting
rmse_df <- as.data.frame(booted_rmses)
names(rmse_df) <- names(bmodel_list)

# paired t-test: CC vs WGCNA
t.test(rmse_df$ConsensusClustering, rmse_df$WGCNA, paired = TRUE)
t.test(rmse_df$ConsensusClustering, rmse_df$PCA, paired = TRUE)



################################################################################
#CHILD DATA
#
#

#1. PCA-based model (with PC1 to PC9)
model_pca <- lm(fmpHSII ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + infant_sex + ageyrsHSII, data = df_merged)

#2. WGCNA-based model (with module eigengenes)
model_wgcna <- lm(fmpHSII ~ MEblue + MEturquoise + MEbrown + MEyellow  + MEgrey + infant_sex + ageyrsHSII, data = df_merged)

#3. Consensus clustering model (with cluster centroid scores)
model_cc <- lm(fmpHSII ~ CC1 + CC2 + CC3 + CC4 + infant_sex + ageyrsHSII, data = df_merged)

#4. Compare AIC across models
cwt_aic <-AIC(model_pca, model_wgcna, model_cc)

#5. Compare R2 across models
# List your models
model_list <- list(
  PCA = model_pca,
  WGCNA = model_wgcna,
  ConsensusClustering = model_cc
)

# Extract R² and adjusted R²
model_r2_df <- lapply(model_list, function(mod) {
  summ <- summary(mod)
  data.frame(
    R2 = summ$r.squared,
    Adj_R2 = summ$adj.r.squared
  )
}) %>% bind_rows(.id = "Method")

# Calculate RMSE and MAE for childhood models
child_rmse_mae <- lapply(model_list, function(mod) {
  obs <- mod$model$fmpHSII
  pred <- predict(mod)
  rmse <- sqrt(mean((obs - pred)^2))
  mae <- mean(abs(obs - pred))
  data.frame(RMSE = rmse, MAE = mae)
}) %>% bind_rows(.id = "Method")

# Combine R2 and RMSE/MAE into one table
child_model_metrics <- left_join(model_r2_df, child_rmse_mae, by = "Method")
print(child_model_metrics)


#testing differences in models. AIC cant test directly unless a nested model
#this gives a probabilistic interpretaiton of which is better
aic_df <- AIC(model_pca, model_wgcna, model_cc)
aic_df$delta_AIC <- aic_df$AIC - min(aic_df$AIC)
aic_df$akaike_weight <- exp(-0.5 * aic_df$delta_AIC) / sum(exp(-0.5 * aic_df$delta_AIC))
print(aic_df)


## bootstrap to see if there is a difference in RSME 
set.seed(123)
n_boot <- 1000

boot_rmse <- function(mod) {
  obs <- mod$model$fmpHSII
  pred <- predict(mod)
  replicate(n_boot, {
    idx <- sample(seq_along(obs), replace = TRUE)
    sqrt(mean((obs[idx] - pred[idx])^2))
  })
}

booted_rmses <- lapply(model_list, boot_rmse)

# Convert to data.frame for t-tests or plotting
rmse_df <- as.data.frame(booted_rmses)
names(rmse_df) <- names(model_list)

# paired t-test: CC vs WGCNA
t.test(rmse_df$ConsensusClustering, rmse_df$WGCNA, paired = TRUE)
t.test(rmse_df$ConsensusClustering, rmse_df$PCA, paired = TRUE)




#compare AIC for birthweight, fmp at birth, and childweight using top 3 clusters, PCs, and Modules
#top_modules) WGCNA modules top based on mean connectivity are #"brown"     "turquoise" "blue" 
#-------------------------------------------------------------------------------
################################################################################
#BIRTH OUTCOMES
#
#

#1. PCA-based model (with PC1 to PC3)
bmodel_pca <- lm(fmpIPV3 ~ PC1 + PC2 + PC3 + infant_sex + GA_at_birth_in_days, data = df_merged)

#2. WGCNA-based model (with module eigengenes)
bmodel_wgcna <- lm(fmpIPV3 ~ MEblue + MEturquoise + MEbrown + infant_sex + GA_at_birth_in_days, data = df_merged)

#3. Consensus clustering model (with cluster centroid scores)
bmodel_cc <- lm(fmpIPV3 ~ CC1 + CC2 + CC3 + infant_sex + GA_at_birth_in_days, data = df_merged)

#4. Compare AIC across models
bwt_aic <-AIC(bmodel_pca, bmodel_wgcna, bmodel_cc)

#5. Compare R2 across models
# List your models
bmodel_list <- list(
  PCA = bmodel_pca,
  WGCNA = bmodel_wgcna,
  ConsensusClustering = bmodel_cc
)

# Extract R² and adjusted R²
bmodel_r2_df <- lapply(bmodel_list, function(mod) {
  summ <- summary(mod)
  data.frame(
    R2 = summ$r.squared,
    Adj_R2 = summ$adj.r.squared
  )
}) %>% bind_rows(.id = "Method")

# Calculate RMSE and MAE for birth models
birth_rmse_mae <- lapply(bmodel_list, function(mod) {
  obs <- mod$model$fmpIPV3
  pred <- predict(mod)
  rmse <- sqrt(mean((obs - pred)^2))
  mae <- mean(abs(obs - pred))
  data.frame(RMSE = rmse, MAE = mae)
}) %>% bind_rows(.id = "Method")

# Combine R2 and RMSE/MAE into one table
birth_model_metrics <- left_join(bmodel_r2_df, birth_rmse_mae, by = "Method")
print(birth_model_metrics)


#testing differences in models. AIC cant test directly unless a nested model
#this gives a probabilistic interpretaiton of which is better
aic_df <- AIC(bmodel_pca, bmodel_wgcna, bmodel_cc)
aic_df$delta_AIC <- aic_df$AIC - min(aic_df$AIC)
aic_df$akaike_weight <- exp(-0.5 * aic_df$delta_AIC) / sum(exp(-0.5 * aic_df$delta_AIC))
print(aic_df)


## bootstrap to see if there is a difference in RSME 
set.seed(123)
n_boot <- 1000

boot_rmse <- function(mod) {
  obs <- mod$model$fmpIPV3
  pred <- predict(mod)
  replicate(n_boot, {
    idx <- sample(seq_along(obs), replace = TRUE)
    sqrt(mean((obs[idx] - pred[idx])^2))
  })
}

booted_rmses <- lapply(bmodel_list, boot_rmse)

# Convert to data.frame for t-tests or plotting
rmse_df <- as.data.frame(booted_rmses)
names(rmse_df) <- names(bmodel_list)

# paired t-test: CC vs WGCNA
t.test(rmse_df$ConsensusClustering, rmse_df$WGCNA, paired = TRUE)
t.test(rmse_df$ConsensusClustering, rmse_df$PCA, paired = TRUE)



################################################################################
#CHILD OUTCOMES
#
#

#1. PCA-based model (with PC1 to PC3)
model_pca <- lm(fmpHSII ~ PC1 + PC2 + PC3 + infant_sex + ageyrsHSII, data = df_merged)

#2. WGCNA-based model (with module eigengenes)
model_wgcna <- lm(fmpHSII ~ MEblue + MEturquoise + MEbrown + infant_sex + ageyrsHSII, data = df_merged)

#3. Consensus clustering model (with cluster centroid scores)
model_cc <- lm(fmpHSII ~ CC1 + CC2 + CC3 + infant_sex + ageyrsHSII, data = df_merged)

#4. Compare AIC across models
cwt_aic <-AIC(model_pca, model_wgcna, model_cc)

#5. Compare R2 across models
# List your models
model_list <- list(
  PCA = model_pca,
  WGCNA = model_wgcna,
  ConsensusClustering = model_cc
)

# Extract R² and adjusted R²
model_r2_df <- lapply(model_list, function(mod) {
  summ <- summary(mod)
  data.frame(
    R2 = summ$r.squared,
    Adj_R2 = summ$adj.r.squared
  )
}) %>% bind_rows(.id = "Method")

# Calculate RMSE and MAE for childhood models
child_rmse_mae <- lapply(model_list, function(mod) {
  obs <- mod$model$fmpHSII
  pred <- predict(mod)
  rmse <- sqrt(mean((obs - pred)^2))
  mae <- mean(abs(obs - pred))
  data.frame(RMSE = rmse, MAE = mae)
}) %>% bind_rows(.id = "Method")

# Combine R2 and RMSE/MAE into one table
child_model_metrics <- left_join(model_r2_df, child_rmse_mae, by = "Method")
print(child_model_metrics)


aic_df <- AIC(model_pca, model_wgcna, model_cc)
aic_df$delta_AIC <- aic_df$AIC - min(aic_df$AIC)
aic_df$akaike_weight <- exp(-0.5 * aic_df$delta_AIC) / sum(exp(-0.5 * aic_df$delta_AIC))
print(aic_df)


## bootstrap to see if there is a difference in RSME 
set.seed(123)
n_boot <- 1000

boot_rmse <- function(mod) {
  obs <- mod$model$fmpHSII
  pred <- predict(mod)
  replicate(n_boot, {
    idx <- sample(seq_along(obs), replace = TRUE)
    sqrt(mean((obs[idx] - pred[idx])^2))
  })
}

booted_rmses <- lapply(model_list, boot_rmse)

# Convert to data.frame for t-tests or plotting
rmse_df <- as.data.frame(booted_rmses)
names(rmse_df) <- names(model_list)

# paired t-test: CC vs WGCNA
t.test(rmse_df$ConsensusClustering, rmse_df$WGCNA, paired = TRUE)
t.test(rmse_df$ConsensusClustering, rmse_df$PCA, paired = TRUE)









