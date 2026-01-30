#-------------------------------------------------------------------------------
# Author: Dr. Ellen C Francis
# Topic: comparison of omic analyses for molecular biologists
# 
# Here we merge the data sets from consensus clustering, wgcna, pca, and meta
# then we build our models and use AIC to pick the data reduction method?
#-------------------------------------------------------------------------------

library(tidyverse)
library(knitr)
library(VIM)
library(finalfit)
library(naniar)
#meta data
metaSel<- read.csv("metaSelect.csv") %>% select(-"X")

#unsupervised latent pattern variables
module_scores<- read.csv("module_scores_df_rr.csv") %>% select(-"X")
pca_scores<- read.csv("pc_scores_df_RR.csv") %>% select(-"X")
clust_scores<- read.csv("cluster_centroid_df_rr.csv") %>% select(-"X")

df1 <- merge(module_scores, pca_scores, by="childID")
df2 <- merge(df1, clust_scores, by="childID")
df_merged<- right_join(metaSel, df2, by="childID")

#drop if missing birth weight because this person dropped out of study
#-------------------------------------------------------------------------------
df_merged<- filter(df_merged, !is.na(wtkgIPV3))#108 observations
#-------------------------------------------------------------------------------

################################################################################

#
#how many completed an in person visit?
df_merged<- filter(df_merged, !is.na(BMIHSII))#108 observations


### Number missing child fat mass
sum(is.na(df_merged$fmpHSII))

### Complete-case dataset (dropping missing BMI)
df_merged_complete <- df_merged |> 
  filter(!is.na(fmpHSII))
nrow(df_merged_complete)



# 1) Overall BMI missingness
df_merged |>
  summarise(pct_missing_fmphsii = mean(is.na(fmpHSII)) * 100) |>
  mutate(pct_missing_fmphsii = round(pct_missing_fmphsii, 1)) |>
  kable(caption = "Overall fmp missingness (%)")

# we include as predictors other data, and a few collected at followup visit at age 5


#Good predictors:  Enrollment maternal demographics: Baseline clinical:
#maternal BMI at enrollment;Perinatal exposures:GDM, smoking in pregnancy;
#breastfeeding duration, physical activity, diet, bmi
#early life anthropometrics (0–12 months)

predvars <- c("bfcat",
              "GA_at_birth_in_days",
              "gdm",
              "infant_sex",
              "maternal_age",
              "momedu",
              "gwg_iom",
              "Parity",
              "pbmi",
              "pred_kcal",
              "race_cat",
              "wtkgIPV3",
              "gestsmoking",
              "child_hei_score",
              "childbmipercentile",
              "pahsii"

)

metaSelectSub <- select(df_merged, fmpHSII)
metaSelectSub2 <- df_merged[,c(predvars)]
metaSelectSub3 <- cbind(metaSelectSub,metaSelectSub2)
metaSelectSub3$momedu <- as.factor(metaSelectSub3$momedu)


covar_names <- c(
  bfcat               = "Infant feeding",
  GA_at_birth_in_days = "Gestational age",
  gdm                 = "GDM exposure",
  infant_sex          = "Offspring sex",
  maternal_age        = "Maternal age",
  momedu              = "Education",
  gwg_iom             = "Gest. kg gain",
  Parity              = "Parity",
  pbmi                = "pBMI",
  pred_kcal           = "Maternal kcal",
  race_cat            = "Race/ethnicity",
  wtkgIPV3            = "Birthweight",
  gestsmoking         = "Maternal smoking",
  fmpHSII             = "Fat mass (age 4)",
  child_hei_score     = "Child diet",
  childbmipercentile  = "Child BMI",
  pahsii              = "Child activity"
)

metaSelectSub3_renamed <- metaSelectSub3
names(metaSelectSub3_renamed) <- covar_names[names(metaSelectSub3)]


aggr(metaSelectSub3_renamed, numbers = TRUE, prop = FALSE, sortVar = TRUE,
     cex.axis = 0.6) 
mcar_test(metaSelectSub3)
#evidence supports mcar

# Does fat mass missingness vary by observed variables? (MAR hint)

dependent = "fmpHSII"
metaSelectSub3 %>% 
  missing_compare(dependent, predvars) %>% 
  knitr::kable(row.names=FALSE, align = c("l", "l", "r", "r", "r"), 
               caption = "Mean comparisons between values of responders (Not missing) and 
        non-responders (Missing) on the fat mass variable.")

mod <- glm(is.na(fmpHSII) ~pred_kcal_child+ maternal_age + gestsmoking + gdm + edlevel+bfcat+Parity, family=binomial,data=metaSelect)
summary(mod)
#VIF(mod)
#- complete case analysis, mean and single regression imputation are still very often applied (Eekhout et al. 2012) but generally not recommended because they decreases statistical power or lead to an incorrect estimation of standard errors when the data is MCAR, MAR and MNAR     

##- MCAR means missingness is completely unrelated to anything else.   
#If missingness is significantly associated with observed variables → MCAR is ruled out.   
#- MNAR not testable, but could consider if likely to be dependent on itself or other unobserved variables   
# missingness in Visit 2 outcomes was largely consistent with MCAR,
######################################################################################
pid<- select(df_merged, childID, ageyrsHSII)
metaSelectSub4 <- cbind(metaSelectSub3,pid)
#write.csv(metaSelectSub4, "metaSelect_use_rr.csv")







