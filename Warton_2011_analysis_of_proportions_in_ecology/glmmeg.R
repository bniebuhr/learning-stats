## glmmeg.R: R code demonstrating how to fit a logistic regression model, with a random intercept term, to randomly generated overdispersed binomial data. 

## David I. Warton and Francis K. C. Hui
## School of Mathematics and Statistics
## The University of New South Wales
## Last modified: 27/7/10

## Some brief details: GLMM’s are fitted in the following code using the lme4 package on R, which you will need to have installed from CRAN. This package fits GLMM’s using Laplace quadrature, which usually provides a good approximation, particularly when fitting a model with one or two random effects terms. If you get a warning that convergence has not been achieved, try using the nAGQ argument (e.g. add ", nAGQ=4" to the line where you call the glmer function) to fit the model using a more accurate but more computationally intensive approach known as adaptive quadrature.

## REFERENCES ##

## For a general introduction to GLMM:
## Benjamin M. Bolker, Mollie E. Brooks, Connie J. Clark, Shane W. Geange, John R. Poulsen, M. Henry H. Stevens and Jada-Simone S. White (2009) Generalized linear mixed models: a practical guide for ecology and evolution. Trends in Ecology & Evolution, 24 (3), 127-135.

## For further details on implementing GLMM in R or S-Plus:
## José C. Pinheiro, Douglas M. Bates (2009) Mixed-Effects Models in S and S-PLUS, Second edition.  Springer-Verlag, New York, USA.

## For details on implementing GLMM in SAS:
## Ramon C. Littell, George A. Milliken, Walter W. Stroup, Russell D. Wolfinger, Oliver Schabenberber (2006) SAS for Mixed Models, Second Edition. SAS Institute, Cary, USA.

## NOTE - the below code currently does not run when using R 2.9.1 or a later version, instead returning the error "Number of levels of a grouping factor for the random effects must be less than the number of observations". This error message should not appear, and if it does appear the problem can be avoided by expanding the dataset out into a Bernoulli response (see end of code), or by downloading and installing an older version of the package:
## Matrix Package:- R package version 0.999375-24
## lme4 Package: - R package version 0.999375-31

## Enjoy massaging your data!

#########################################################################################

## GENERATING AN EXAMPLE DATASET FOR ANALYSIS ##

## In order to illustrate the used of GLMM, over-dispersed binomial data are generated here according to a balanced one-way ANOVA design, with 15 “species” at each of four levels of the factor “location”.

species     = 1:60 # We assume that the 60 rows of the dataset correspond to 60 different species.
location    = c(rep("australia",15), rep("canada",15), rep("argentina",15), rep("china",15)) 
sample.size = 12
p           = c(rep(0.3,15), rep(0.4,15), rep(0.5,15), rep(0.6,15))
eta         = log( p/(1-p) ) + rnorm(60)
p           = exp(eta) / ( exp(eta) + 1 )

success     = rbinom(60, size=sample.size, prob=p) 	
failure     = sample.size - success
location    = factor(location)
dataset     = data.frame(location, species, sample.size, success, failure)
rm(location, species, success, failure)


#########################################################################################

## ANALYSIS OF EXAMPLE DATASET ##

attach(dataset)

## Plot the sample proportions against location
plot(success/sample.size ~ location)

## Logistic regression (fixed effects)
fit.glm     = glm(cbind(success, failure) ~ location, family = binomial)#, dataset=dataset)
anova(fit.glm, test = "Chisq")
summary(fit.glm)
(fit.p.hat <- exp(coef(fit.glm)) / (exp(coef(fit.glm)) + 1)) # expected probabilities

## Check to see if residual deviance is large relative to residual df.
## Note that for the data generated above, the residual deviance is over twice as large as the residual df, so there is clear evidence that the data are overdispersed.
## This means that a GLMM should be fitted, with a random term for species (row):

## GLMM
library(lme4)
fit.glmm           = glmer(cbind(success, failure) ~ location + (1|species), family = binomial, data=dataset)
fit.glmm.intercept = glmer(cbind(success, failure) ~ 1 + (1|species), family = binomial, data=dataset)
anova(fit.glmm.intercept, fit.glmm)
## Note the significant evidence of a location effect.
## If you got the error message "Number of levels of a grouping factor for the random effects must be less than the number of observations", see code at the end.


#########################################################################################

## PRODUCING DIAGNOSTIC PLOTS ##

## Logistic regression
plot(fit.glm$fit, residuals(fit.glm), pch = 19, las = 1, cex = 1.4)
abline(0,0,lwd = 1.5)
## check for no pattern

## GLMM
par(mfrow=c(1,2))
## first we plot random effects against the predicted values from the fixed effect component of the model and check for no trend:
m        = model.matrix(fit.glmm)
ft.fix   = m %*% fixef(fit.glmm)
plot(ft.fix, ranef(fit.glmm, drop = T)$species, pch = 19, las = 1, cex = 1.4)
abline(0,0,lwd = 1.5)

## now check for approximate normality of random effects:
qqnorm(ranef(fit.glmm, drop = T)$species, pch = 19, las = 1, cex = 1.4)


#########################################################################################

## WORK_AROUND IF YOU GOT AN ERROR RUNNING GLMM ##

## If you got the error message "Number of levels of a grouping factor for the random effects must be less than the number of observations": this is because lme4 has a bug in R version 2.9.1 (which we hope will be fixed soon)! You can either use an older version of lme4 for analysis or as a work-around you can expand your dataset out and fit the glmm as in the below:

detach(dataset)
## Create an expanded dataset based on dataset (this can readily be used on your own data, it just requires a data.frame called "dataset" containing successes and failures labelled as "success" and "failure"):

dataset.expanded = dataset[0,]
for (i in 1:length(dataset$success))
{
	if(dataset$success[i]>0)
{
dataset.add.succ = dataset[rep(i,dataset$success[i]),]
		dataset.add.succ$success=1
		dataset.add.succ$failure=0
		dataset.expanded=rbind(dataset.expanded, dataset.add.succ)
	}
	if(dataset$failure[i]>0)
{
dataset.add.fail = dataset[rep(i,dataset$failure[i]),]
		dataset.add.fail$success=0
		dataset.add.fail$failure=1
		dataset.expanded=rbind(dataset.expanded, dataset.add.fail)
	}
}

## Fit the GLMM’s
fit.glmm = glmer(success ~ location + (1|species), family = binomial, data=dataset.expanded)
fit.glmm.intercept = glmer(success ~ 1 + (1|species), family = binomial, data=dataset.expanded)
anova(fit.glmm.intercept, fit.glmm)

## Construct residual plots:
par(mfrow=c(1,2))
## first plot random effects against the predicted values and check for no trend:
m        = model.matrix(fit.glmm)
ft.fix   = m%*%fixef(fit.glmm)
rans           = t(as.matrix(fit.glmm@Zt)) %*% ranef(fit.glmm)$species[[1]]
plot(ft.fix, rans, pch = 19, las = 1, cex = 1.4)
abline(0,0,lwd = 1.5)

## now check for approximate normality of random effects:
qqnorm(ranef(fit.glmm, drop = T)$species, pch = 19, las = 1, cex = 1.4)


