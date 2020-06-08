## boot.glmm.R: R code for estimating P-values by applying the bootstrap to a GLMM likelihood ratio statistic.

## David I. Warton and Francis K. C. Hui
## School of Mathematics and Statistics
## The University of New South Wales
## Last modified: 27/7/10

## This code could be modified to calculate any other statistic by modifying the definition of “p.obs” (the observed P-value) and “ps[i]” (the P-value for the ith resample) – the key point is that the same code needs to be used to calculate each, except “p.obs” uses the original dataset (“dataset”), and “ps[i]” uses the resampled dataset (“data.i”).

## The code has been written so that it can be applied to data with the same structure as that of the supplementary file glmmeg.R:
## “dataset” is a dataframe containing:
## - “success” (the binomial counts)
## - “sample.size” (the number of trials for each count, which need not be constant)
## - “location” (the explanatory variable being used in testing, whether a continuous variable, a factor describing treatment levels, ...)
## - “species” (indicating the different species used in analysis, for which different random effects will be fitted).

## Bootstrapping GLMM's is computationally intensive - if we calculate 1000 GLMM’s, for a dataset of the size of that generated in the glmm.R supplementary file, it might take a couple of minutes to run.  For a smaller dataset (for which resampling is more likely to actually be needed), it will take less time.

####################################################################

#First calculate a test statistic for the original dataset:
n.bootstrap=10
fit  = glmer(cbind(success, sample.size - success) ~ location + (1|species), family = binomial, data=dataset)
fit2 = glmer(cbind(success, sample.size - success) ~ 1 + (1|species), family = binomial, data=dataset)
p.obs = anova(fit2, fit)$"Pr(>Chisq)"[2] # LR test
n.obs = dim(dataset)[1]

#Now generate references for bootstrap samples:
boot.ref = sample(n.obs, n.bootstrap*n.obs, replace=T)
boot.ref = matrix(boot.ref, n.obs, n.bootstrap)
p.count  = 0
eps      = 1.e-8
data.i=dataset

ps = rep(NA, n.bootstrap)

#Recalculate the test stat for each bootstrap sample:
for( i in 1:n.bootstrap )
{
	print(i)
  i.ref              = boot.ref[,i]
	data.i$success     = dataset$success[i.ref]
	data.i$sample.size = dataset$sample.size[i.ref]

	fit  = glmer(cbind(success, sample.size - success) ~ location + (1|species), binomial, data = data.i)
	fit2 = glmer(cbind(success, sample.size - success) ~ 1 + (1|species), binomial, data = data.i)
	ps[i] = anova(fit2, fit)$"Pr(>Chisq)"[2] # LR test
}
ps[is.na(ps)]=1
p.boot = mean(ps<(p.obs+eps))
print(paste("Bootstrap P-value:",p.boot))


