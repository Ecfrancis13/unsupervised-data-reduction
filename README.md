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
