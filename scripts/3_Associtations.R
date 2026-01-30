library(tidyverse)
library(glmnet)
library(purrr)
library(broom)

#metabolic markers
bldresults<-metaSel<- read.csv("metaSelect.csv") %>% select(-"X") %>% select(childID, ResultGlucose,ResultInsulin,homair, ResultAdiponectin, ResultLeptin, ResultCholesterol,
                                                                             ResultTriglyceride,ResultHDL,ResultLDL)
metaSel<- read.csv("metaSelect_use_rr.csv") %>% select(-"X")
#93 variables

df_merged <- metaSel
#62 observations after merge

module_scores<- read.csv("module_scores_df_rr.csv") %>% select(-"X")
pca_scores<- read.csv("pc_scores_df_RR.csv") %>% select(-"X")
clust_scores<- read.csv("cluster_centroid_df_rr.csv") %>% select(-"X")

df1 <- merge(module_scores, pca_scores, by="childID")
df2 <- merge(df1, clust_scores, by="childID")
df_merged<- right_join(metaSel, df2, by="childID")

df_m1 <- df_merged %>%
  select(CC1, CC2, CC3, CC4, 
         PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9,
         MEblue, MEturquoise, MEbrown, MEyellow,
         infant_sex, ageyrsHSII, fmpHSII) %>%
  drop_na()


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#CLUSTERS
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
X <- model.matrix(~ CC1 + CC2 + CC3 + CC4 + infant_sex + ageyrsHSII, data = df_m1)[, -1]
y <- df_m1$fmpHSII

# Penalize clusters, not covariates
penalty <- c(CC1=1, CC2=1, CC3=1, CC4=1, infant_sex=0, ageyrsHSII=0)

set.seed(1234)
cc_fit_cv <- cv.glmnet(
  X, y,
  alpha = 1,              # LASSO
  family = "gaussian",
  penalty.factor = penalty
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# WGCNA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
X <- model.matrix(~ MEblue + MEturquoise + MEbrown + MEyellow  + infant_sex + ageyrsHSII, data = df_m1)[, -1]
y <- df_m1$fmpHSII

# Penalize clusters, not covariates
penalty <- c(MEblue=1, MEturquoise=1, MEbrown=1, MEyellow=1, infant_sex=0, ageyrsHSII=0)
set.seed(1234)

wgcna_fit_cv <- cv.glmnet(
  X, y,
  alpha = 1,              # LASSO
  family = "gaussian",
  penalty.factor = penalty
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# PCA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
X <- model.matrix(~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + infant_sex + ageyrsHSII, data = df_m1)[, -1]
y <- df_m1$fmpHSII

# Penalize clusters, not covariates
penalty <- c(PC1=1, PC2=1, PC3=1, PC4=1, PC5=1,  PC6=1,  PC7=1,  PC8=1, 
             PC9=1,  infant_sex=0, ageyrsHSII=0)
set.seed(1234)

pca_fit_cv <- cv.glmnet(
  X, y,
  alpha = 1,              # LASSO
  family = "gaussian",
  penalty.factor = penalty
)

################################################################################
#compare patterns that weren't penalized to zero
coef(cc_fit_cv, s = "lambda.min")
coef(wgcna_fit_cv, s = "lambda.min")
coef(pca_fit_cv, s = "lambda.min")

################################################################################

#now run the association tests
# ----------------------------
# Outcome label-variable mapping
# ----------------------------
outcome_map <- tibble::tibble(
  label = c(
    "Glucose, mg/dL",
    "Insulin, μIU/mL",
    "HOMA-IR, %",
    "Adiponectin, ug/ml",
    "Leptin, μg/l",
    "Cholesterol, mg/dL",
    "Triglycerides, mg/dL",
    "HDL, mg/dL",
    "LDL, mg/dL",
    "%Fat mass",
    "BMI percentile"
  ),
  var = c(
    "ResultGlucose",
    "ResultInsulin",
    "homair",
    "ResultAdiponectin",
    "ResultLeptin",
    "ResultCholesterol",
    "ResultTriglyceride",
    "ResultHDL",
    "ResultLDL",
    "fmpHSII",
    "childbmipercentile"
  )
)
##############################################################################
# ----------------------------
# Setup Predictors and Model, these are the patterns that were retained based on LASSo
# ----------------------------
pred_set1c3 <- c("CC3")

pred_set2p6 <- c("PC6")
pred_set2p7 <- c("PC7")
pred_set2p8 <- c("PC8")

model1_covars <- c("infant_sex", "ageyrsHSII")  # Adjust for your setup! 
model2_covars <- c("infant_sex", "ageyrsHSII", "pbmi")  # Adjust for your setup! 

pred_sets <- list(
  CC3 = pred_set1c3,
  PC7 = pred_set2p7,
  PC6 = pred_set2p6,
  PC8 = pred_set2p8
)

adj_sets <- list(
  Model1 = model1_covars,
  Model2 = model2_covars
)

df_merged2 <- left_join(df_merged, bldresults, by="childID")

# ----------------------------
# Function to loop over outcomes and clusters
# ----------------------------
get_results <- function(outcome_var, outcome_label) {
  map_dfr(names(pred_sets), function(pset_name) {
    predictors <- pred_sets[[pset_name]]
    
    map_dfr(names(adj_sets), function(adjust_name) {
      covars <- adj_sets[[adjust_name]]
      
      map_dfr(predictors, function(pred) {
        
        # build formula
        fml <- as.formula(
          paste(outcome_var, "~", pred, "+", paste(covars, collapse = "+"))
        )
        
        fit <- lm(fml, data = df_merged2)
        tidy_res <- tidy(fit, conf.int = TRUE, conf.level = 0.95) %>%
          filter(term == pred)
        
        tibble(
          Outcome     = outcome_label,
          OutcomeVar  = outcome_var,
          PredSet     = pset_name,
          Predictor   = pred,
          Model       = adjust_name,
          N           = nobs(fit),
          Beta        = tidy_res$estimate,
          Lower95CI   = tidy_res$conf.low,
          Upper95CI   = tidy_res$conf.high,
          Pvalue      = tidy_res$p.value,
          BIC         = BIC(fit)
        )
      })
    })
  })
}
# ----------------------------
# Apply across all outcomes
# ----------------------------
all_results <- map2_dfr(outcome_map$var, outcome_map$label, get_results)

#Because childhood fat mass percent (fmpHSII) was the primary outcome of interest, 
#false discovery rate correction was applied to secondary outcomes only; results for fmpHSII are presented with nominal p-values.
# split
results_no_fmp <- all_results %>%
  filter(OutcomeVar != "fmpHSII")

results_fmp <- all_results %>%
  filter(OutcomeVar == "fmpHSII") %>%
  mutate(Qvalue = NA_real_)

# apply FDR
results_no_fmp <- results_no_fmp %>%
  group_by(OutcomeVar, Model) %>%
  mutate(Qvalue = p.adjust(Pvalue, method = "BH")) %>%
  ungroup()

# reassemble
all_results <- bind_rows(results_no_fmp, results_fmp)
#################################################################################

#write.csv(all_results, "all_results_rr95v2.csv")

