pkgname <- "semTools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('semTools')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("EFA-class")
### * EFA-class

flush(stderr()); flush(stdout())

### Name: EFA-class
### Title: Class For Rotated Results from EFA
### Aliases: EFA-class show,EFA-method summary,EFA-method

### ** Examples

unrotated <- efaUnrotate(HolzingerSwineford1939, nf=3, varList=paste0("x", 1:9), estimator="mlr")
summary(unrotated, std=TRUE)
inspect(unrotated, "standardized")

# Rotated by Quartimin
rotated <- oblqRotate(unrotated, method="quartimin")
summary(rotated)



cleanEx()
nameEx("FitDiff-class")
### * FitDiff-class

flush(stderr()); flush(stdout())

### Name: FitDiff-class
### Title: Class For Representing A Template of Model Fit Comparisons
### Aliases: FitDiff-class show,FitDiff-method summary,FitDiff-method

### ** Examples

HW.model <- ' visual =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed =~ x7 + x8 + x9 '

out <- measurementInvariance(HW.model, data=HolzingerSwineford1939, group="school", quiet=TRUE)
modelDiff <- compareFit(out)
summary(modelDiff)
summary(modelDiff, fit.measures="all")
summary(modelDiff, fit.measures=c("aic", "bic"))

## Not run: 
##D # Save results to a file 
##D saveFile(modelDiff, file="modelDiff.txt")
##D 
##D # Copy to a clipboard
##D clipboard(modelDiff)
## End(Not run)



cleanEx()
nameEx("auxiliary")
### * auxiliary

flush(stderr()); flush(stdout())

### Name: auxiliary
### Title: Analyzing data with full-information maximum likelihood with
###   auxiliary variables
### Aliases: auxiliary cfa.auxiliary sem.auxiliary growth.auxiliary
###   lavaan.auxiliary

### ** Examples

# Example of confirmatory factor analysis

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
			  
dat <- data.frame(HolzingerSwineford1939, z=rnorm(nrow(HolzingerSwineford1939), 0, 1))
			  
fit <- cfa(HS.model, data=dat, meanstructure=TRUE) 
fitaux <- auxiliary(HS.model, aux="z", data=dat, fun="cfa") # Use lavaan script
fitaux <- cfa.auxiliary(fit, aux="z", data=dat) # Use lavaan output

# Example of multiple groups confirmatory factor analysis

fitgroup <- cfa(HS.model, data=dat, group="school", meanstructure=TRUE)
fitgroupaux <- cfa.auxiliary(fitgroup, aux="z", data=dat, group="school")

# Example of path analysis

mod <- ' x5 ~ x4
x4 ~ x3
x3 ~ x1 + x2'

fitpath <- sem(mod, data=dat, fixed.x=FALSE, meanstructure=TRUE) # fixed.x must be FALSE
fitpathaux <- sem.auxiliary(fitpath, aux="z", data=dat)

# Example of full structural equation modeling

dat2 <- data.frame(PoliticalDemocracy, z=rnorm(nrow(PoliticalDemocracy), 0, 1))
model <- ' 
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

    dem60 ~ ind60
    dem65 ~ ind60 + dem60

    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
fitsem <- sem(model, data=dat2, meanstructure=TRUE)
fitsemaux <- sem.auxiliary(fitsem, aux="z", data=dat2, meanstructure=TRUE)

# Example of covariate at the factor level

HS.model.cov <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
             speed   =~ x7 + x8 + x9 
			  visual ~ sex
			  textual ~ sex
			  speed ~ sex'
	  
fitcov <- cfa(HS.model.cov, data=dat, fixed.x=FALSE, meanstructure=TRUE) 
fitcovaux <- cfa.auxiliary(fitcov, aux="z", data=dat)

# Example of  Endogenous variable with single indicator 
HS.model.cov2 <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              x7 ~ visual + textual'
 	  
fitcov2 <- sem(HS.model.cov2, data=dat, fixed.x=FALSE, meanstructure=TRUE) 
fitcov2aux <- sem.auxiliary(fitcov2, aux="z", data=dat)

# Multiple auxiliary variables
HS.model2 <- ' visual  =~ x1 + x2 + x3
              speed   =~ x7 + x8 + x9'
fit <- cfa(HS.model2, data=HolzingerSwineford1939)
fitaux <- cfa.auxiliary(HS.model2, data=HolzingerSwineford1939, aux=c("x4", "x5")) 



cleanEx()
nameEx("clipboard")
### * clipboard

flush(stderr()); flush(stdout())

### Name: clipboard_saveFile
### Title: Copy or save the result of 'lavaan' or 'FitDiff' objects into a
###   clipboard or a file
### Aliases: clipboard saveFile

### ** Examples

## Not run: 
##D library(lavaan)
##D HW.model <- ' visual  =~ x1 + c1*x2 + x3
##D               textual =~ x4 + c1*x5 + x6
##D                speed   =~ x7 + x8 + x9 '
##D 
##D fit <- cfa(HW.model, data=HolzingerSwineford1939, group="school", meanstructure=TRUE)
##D 
##D # Copy the summary of the lavaan object
##D clipboard(fit)
##D 
##D # Copy the modification indices and the model fit from the miPowerFit function
##D clipboard(fit, "mifit")
##D 
##D # Copy the parameter estimates
##D clipboard(fit, "coef")
##D 
##D # Copy the standard errors
##D clipboard(fit, "se")
##D 
##D # Copy the sample statistics
##D clipboard(fit, "samp")
##D 
##D # Copy the fit measures
##D clipboard(fit, "fit")
##D 
##D # Save the summary of the lavaan object
##D saveFile(fit, "out.txt")
##D 
##D # Save the modification indices and the model fit from the miPowerFit function
##D saveFile(fit, "out.txt", "mifit")
##D 
##D # Save the parameter estimates
##D saveFile(fit, "out.txt", "coef")
##D 
##D # Save the standard errors
##D saveFile(fit, "out.txt", "se")
##D 
##D # Save the sample statistics
##D saveFile(fit, "out.txt", "samp")
##D 
##D # Save the fit measures
##D saveFile(fit, "out.txt", "fit")
## End(Not run)



cleanEx()
nameEx("compareFit")
### * compareFit

flush(stderr()); flush(stdout())

### Name: compareFit
### Title: Build an object summarizing fit indices across multiple models
### Aliases: compareFit

### ** Examples

m1 <- ' visual  =~ x1 + x2 + x3
        textual =~ x4 + x5 + x6
        speed   =~ x7 + x8 + x9 '

fit1 <- cfa(m1, data=HolzingerSwineford1939)

m2 <- ' f1  =~ x1 + x2 + x3 + x4 
        f2 =~ x5 + x6 + x7 + x8 + x9 '
fit2 <- cfa(m2, data=HolzingerSwineford1939)
compareFit(fit1, fit2, nested=FALSE)

HW.model <- ' visual =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed =~ x7 + x8 + x9 '

out <- measurementInvariance(HW.model, data=HolzingerSwineford1939, group="school", quiet=TRUE)
compareFit(out)



cleanEx()
nameEx("dat2way")
### * dat2way

flush(stderr()); flush(stdout())

### Name: dat2way
### Title: Simulated Dataset to Demonstrate Two-way Latent Interaction
### Aliases: dat2way

### ** Examples

head(dat2way)



cleanEx()
nameEx("dat3way")
### * dat3way

flush(stderr()); flush(stdout())

### Name: dat3way
### Title: Simulated Dataset to Demonstrate Three-way Latent Interaction
### Aliases: dat3way

### ** Examples

head(dat3way)



cleanEx()
nameEx("efaUnrotate")
### * efaUnrotate

flush(stderr()); flush(stdout())

### Name: efaUnrotate
### Title: Analyze Unrotated Exploratory Factor Analysis Model
### Aliases: efaUnrotate

### ** Examples

unrotated <- efaUnrotate(HolzingerSwineford1939, nf=3, varList=paste0("x", 1:9), estimator="mlr")
summary(unrotated, std=TRUE)
inspect(unrotated, "standardized")

dat <- data.frame(HolzingerSwineford1939, z=rnorm(nrow(HolzingerSwineford1939), 0, 1))
unrotated2 <- efaUnrotate(dat, nf=2, varList=paste0("x", 1:9), aux="z")



cleanEx()
nameEx("exLong")
### * exLong

flush(stderr()); flush(stdout())

### Name: exLong
### Title: Simulated Data set to Demonstrate Longitudinal Measurement
###   Invariance
### Aliases: exLong

### ** Examples

head(exLong)



cleanEx()
nameEx("findRMSEApower")
### * findRMSEApower

flush(stderr()); flush(stdout())

### Name: findRMSEApower
### Title: Find the statistical power based on population RMSEA
### Aliases: findRMSEApower

### ** Examples

findRMSEApower(rmsea0=.05, rmseaA=.08, df=20, n=200)



cleanEx()
nameEx("findRMSEApowernested")
### * findRMSEApowernested

flush(stderr()); flush(stdout())

### Name: findRMSEApowernested
### Title: Find power given a sample size in nested model comparison
### Aliases: findRMSEApowernested

### ** Examples

findRMSEApowernested(rmsea0A = 0.06, rmsea0B = 0.05, rmsea1A = 0.08, 
rmsea1B = 0.05, dfA = 22, dfB = 20, n = 200, alpha = 0.05, group = 1)



cleanEx()
nameEx("findRMSEAsamplesize")
### * findRMSEAsamplesize

flush(stderr()); flush(stdout())

### Name: findRMSEAsamplesize
### Title: Find the minimum sample size for a given statistical power based
###   on population RMSEA
### Aliases: findRMSEAsamplesize

### ** Examples

findRMSEAsamplesize(rmsea0=.05, rmseaA=.08, df=20, power=0.80)



cleanEx()
nameEx("findRMSEAsamplesizenested")
### * findRMSEAsamplesizenested

flush(stderr()); flush(stdout())

### Name: findRMSEAsamplesizenested
### Title: Find sample size given a power in nested model comparison
### Aliases: findRMSEAsamplesizenested

### ** Examples

findRMSEAsamplesizenested(rmsea0A = 0, rmsea0B = 0, rmsea1A = 0.06, 
rmsea1B = 0.05, dfA = 22, dfB = 20, power=0.80, alpha=.05, group=1) 



cleanEx()
nameEx("fitMeasuresMx")
### * fitMeasuresMx

flush(stderr()); flush(stdout())

### Name: fitMeasuresMx
### Title: Find fit measures from an MxModel result
### Aliases: fitMeasuresMx

### ** Examples

## Not run: 
##D library(OpenMx)
##D data(demoOneFactor)
##D manifests <- names(demoOneFactor)
##D latents <- c("G")
##D factorModel <- mxModel("One Factor", 
##D     type="RAM",
##D     manifestVars=manifests, 
##D     latentVars=latents,
##D     mxPath(from=latents, to=manifests),
##D     mxPath(from=manifests, arrows=2),
##D     mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
##D     mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
##D )
##D factorFit <- mxRun(factorModel)
##D round(fitMeasuresMx(factorFit), 3)
##D 
##D # Compare with lavaan
##D library(lavaan)
##D script <- "f1 =~ x1 + x2 + x3 + x4 + x5"
##D fitMeasures(cfa(script, sample.cov = cov(demoOneFactor), sample.nobs = 500, std.lv = TRUE))
## End(Not run)



cleanEx()
nameEx("fmi")
### * fmi

flush(stderr()); flush(stdout())

### Name: fmi
### Title: Fraction of Missing Information.
### Aliases: fmi

### ** Examples

library(Amelia)
library(lavaan)

modsim <- '
f1 =~ 0.7*y1+0.7*y2+0.7*y3
f2 =~ 0.7*y4+0.7*y5+0.7*y6
f3 =~ 0.7*y7+0.7*y8+0.7*y9'

datsim <- simulateData(modsim,model.type="cfa", meanstructure=TRUE, 
                       std.lv=TRUE, sample.nobs=c(200,200))
randomMiss2 <- rbinom(prod(dim(datsim)), 1, 0.1)
randomMiss2 <- matrix(as.logical(randomMiss2), nrow=nrow(datsim))
randomMiss2[,10] <- FALSE
datsim[randomMiss2] <- NA
datsimMI <- amelia(datsim,m=3,idvars="group")

out1 <- fmi(datsimMI$imputations, exclude="group")
out1
                       
out2 <- fmi(datsimMI$imputations, exclude="group", method="null")
out2
                       
out3 <- fmi(datsimMI$imputations, varnames=c("y1","y2","y3","y4"))
out3

out4 <- fmi(datsimMI$imputations, group="group")
out4




cleanEx()
nameEx("impliedFactorStat")
### * impliedFactorStat

flush(stderr()); flush(stdout())

### Name: impliedFactorStat
### Title: Calculate the model-implied factor means and covariance matrix.
### Aliases: impliedFactorMean impliedFactorCov impliedFactorStat

### ** Examples

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939, group="school")
impliedFactorStat(fit)



cleanEx()
nameEx("imposeStart")
### * imposeStart

flush(stderr()); flush(stdout())

### Name: imposeStart
### Title: Specify starting values from a lavaan output
### Aliases: imposeStart

### ** Examples

# The following example show that the longitudinal weak invariance model
# using effect coding was not convergent with three time points but convergent
# with two time points. Thus, the parameter estimates from the model with
# two time points are used as starting values of the three time points.
# The model with new starting values is convergent properly.

weak2time <- ' 
	# Loadings
	f1t1 =~ LOAD1*y1t1 + LOAD2*y2t1 + LOAD3*y3t1
    f1t2 =~ LOAD1*y1t2 + LOAD2*y2t2 + LOAD3*y3t2
	
	# Factor Variances
	f1t1 ~~ f1t1
	f1t2 ~~ f1t2
	
	# Factor Covariances
	f1t1 ~~ f1t2 
	
	# Error Variances
	y1t1 ~~ y1t1
	y2t1 ~~ y2t1
	y3t1 ~~ y3t1
	y1t2 ~~ y1t2
	y2t2 ~~ y2t2
	y3t2 ~~ y3t2
	
	# Error Covariances
	y1t1 ~~ y1t2 
	y2t1 ~~ y2t2 
	y3t1 ~~ y3t2
	
	# Factor Means
	f1t1 ~ NA*1
	f1t2 ~ NA*1
	
	# Measurement Intercepts
	y1t1 ~ INT1*1
	y2t1 ~ INT2*1
	y3t1 ~ INT3*1
	y1t2 ~ INT4*1
	y2t2 ~ INT5*1
	y3t2 ~ INT6*1
	
	# Constraints for Effect-coding Identification
	LOAD1 == 3 - LOAD2 - LOAD3
	INT1 == 0 - INT2 - INT3
	INT4 == 0 - INT5 - INT6
'
model2time <- lavaan(weak2time, data = exLong)

weak3time <- ' 
	# Loadings
	f1t1 =~ LOAD1*y1t1 + LOAD2*y2t1 + LOAD3*y3t1
    f1t2 =~ LOAD1*y1t2 + LOAD2*y2t2 + LOAD3*y3t2
    f1t3 =~ LOAD1*y1t3 + LOAD2*y2t3 + LOAD3*y3t3
	
	# Factor Variances
	f1t1 ~~ f1t1
	f1t2 ~~ f1t2
	f1t3 ~~ f1t3
	
	# Factor Covariances
	f1t1 ~~ f1t2 + f1t3
	f1t2 ~~ f1t3 
	
	# Error Variances
	y1t1 ~~ y1t1
	y2t1 ~~ y2t1
	y3t1 ~~ y3t1
	y1t2 ~~ y1t2
	y2t2 ~~ y2t2
	y3t2 ~~ y3t2
	y1t3 ~~ y1t3
	y2t3 ~~ y2t3
	y3t3 ~~ y3t3
	
	# Error Covariances
	y1t1 ~~ y1t2 
	y2t1 ~~ y2t2 
	y3t1 ~~ y3t2
	y1t1 ~~ y1t3 
	y2t1 ~~ y2t3 
	y3t1 ~~ y3t3
	y1t2 ~~ y1t3
	y2t2 ~~ y2t3 
	y3t2 ~~ y3t3
	
	# Factor Means
	f1t1 ~ NA*1
	f1t2 ~ NA*1
	f1t3 ~ NA*1
	
	# Measurement Intercepts
	y1t1 ~ INT1*1
	y2t1 ~ INT2*1
	y3t1 ~ INT3*1
	y1t2 ~ INT4*1
	y2t2 ~ INT5*1
	y3t2 ~ INT6*1
	y1t3 ~ INT7*1
	y2t3 ~ INT8*1
	y3t3 ~ INT9*1
	
	# Constraints for Effect-coding Identification
	LOAD1 == 3 - LOAD2 - LOAD3
	INT1 == 0 - INT2 - INT3
	INT4 == 0 - INT5 - INT6
	INT7 == 0 - INT8 - INT9
'
### The following command does not provide convergent result
# model3time <- lavaan(weak3time, data = exLong)

### Use starting values from the model with two time points
model3time <- imposeStart(model2time, lavaan(weak3time, data = exLong))
summary(model3time)



cleanEx()
nameEx("indProd")
### * indProd

flush(stderr()); flush(stdout())

### Name: indProd
### Title: Make products of indicators using no centering, mean centering,
###   double-mean centering, or residual centering
### Aliases: indProd orthogonalize

### ** Examples

# Mean centering / two-way interaction / match-paired
dat <- indProd(attitude[,-1], var1=1:3, var2=4:6)

# Residual centering / two-way interaction / match-paired
dat2 <- indProd(attitude[,-1], var1=1:3, var2=4:6, match=FALSE, meanC=FALSE, residualC=TRUE, doubleMC=FALSE)

# Double-mean centering / two-way interaction / match-paired
dat3 <- indProd(attitude[,-1], var1=1:3, var2=4:6, match=FALSE, meanC=TRUE, residualC=FALSE, doubleMC=TRUE)

# Mean centering / three-way interaction / match-paired
dat4 <- indProd(attitude[,-1], var1=1:2, var2=3:4, var3=5:6)

# Residual centering / three-way interaction / match-paired
dat5 <- indProd(attitude[,-1], var1=1:2, var2=3:4, var3=5:6, match=FALSE, meanC=FALSE, residualC=TRUE, doubleMC=FALSE)

# Double-mean centering / three-way interaction / match-paired
dat6 <- indProd(attitude[,-1], var1=1:2, var2=3:4, var3=5:6, match=FALSE, meanC=TRUE, residualC=TRUE, doubleMC=TRUE)



cleanEx()
nameEx("kd")
### * kd

flush(stderr()); flush(stdout())

### Name: kd
### Title: Generate data via the Kaiser-Dickman (1962) algorithm.
### Aliases: kd

### ** Examples

#### First Example

## Get data
dat <- HolzingerSwineford1939[,7:15]
hs.n <- nrow(dat)

## Covariance matrix divided by n
hscov <- ((hs.n-1)/hs.n) * cov(dat)

## Generate new, raw data corresponding to hscov
newdat <- kd(hscov, hs.n)

## Difference between new covariance matrix and hscov is minimal
newcov <- (hs.n-1)/hs.n * cov(newdat)
summary(as.numeric(hscov - newcov))

## Generate sample data, treating hscov as population matrix
newdat2 <- kd(hscov, hs.n, type="sample")

#### Another example

## Define a covariance matrix
covmat <- matrix(0, 3, 3); diag(covmat) <- 1.5; covmat[2:3,1] <- c(1.3, 1.7); covmat[3,2] <- 2.1
covmat <- covmat + t(covmat)

## Generate data of size 300 that have this covariance matrix
rawdat <- kd(covmat, 300)

## Covariances are exact if we compute sample covariance matrix by
## dividing by n (vs by n-1)
summary(as.numeric((299/300)*cov(rawdat) - covmat))

## Generate data of size 300 where covmat is the population covariance matrix
rawdat2 <- kd(covmat, 300)



cleanEx()
nameEx("kurtosis")
### * kurtosis

flush(stderr()); flush(stdout())

### Name: kurtosis
### Title: Finding excessive kurtosis
### Aliases: kurtosis

### ** Examples

kurtosis(1:5)



cleanEx()
nameEx("lavaanStar-class")
### * lavaanStar-class

flush(stderr()); flush(stdout())

### Name: lavaanStar-class
### Title: Class For Representing A (Fitted) Latent Variable Model with
###   Additional Elements
### Aliases: lavaanStar-class inspect,lavaanStar-method
###   summary,lavaanStar-method

### ** Examples

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
			  
dat <- data.frame(HolzingerSwineford1939, z=rnorm(nrow(HolzingerSwineford1939), 0, 1))
			  
fit <- cfa(HS.model, data=dat) 
fitaux <- auxiliary(HS.model, aux="z", data=dat, fun="cfa")



cleanEx()
nameEx("lisrel2lavaan")
### * lisrel2lavaan

flush(stderr()); flush(stdout())

### Name: lisrel2lavaan
### Title: Latent variable modeling in 'lavaan' using LISREL syntax
### Aliases: lisrel2lavaan

### ** Examples

## Not run: 
##D 	## calling lisrel2lavaan without specifying the filename argument will  
##D 	## open a file browser window with which LISREL syntax can be selected. 
##D 	
##D 	## any additional arguments to be passed to lavaan for data analysis can
##D 	## be specified normally. 
##D 	
##D 	lisrel2lavaan(se="standard")
##D 	## lavaan output summary printed to screen
##D 	## lavaan fit object returned silently
##D 	
##D 	## manual file specification 
##D 	
##D 	lisrel2lavaan(filename="myFile.LS8", se="standard")
##D 	## lavaan output summary printed to screen
##D 	## lavaan fit object returned silently
## End(Not run)



cleanEx()
nameEx("loadingFromAlpha")
### * loadingFromAlpha

flush(stderr()); flush(stdout())

### Name: loadingFromAlpha
### Title: Find standardized factor loading from coefficient alpha
### Aliases: loadingFromAlpha

### ** Examples

    loadingFromAlpha(0.8, 4)



cleanEx()
nameEx("longInvariance")
### * longInvariance

flush(stderr()); flush(stdout())

### Name: longInvariance
### Title: Measurement Invariance Tests Within Person
### Aliases: longInvariance longInvariance

### ** Examples

model <- ' f1t1 =~ y1t1 + y2t1 + y3t1
              f1t2 =~ y1t2 + y2t2 + y3t2
			  f1t3 =~ y1t3 + y2t3 + y3t3'

# Create list of variables
var1 <- c("y1t1", "y2t1", "y3t1")
var2 <- c("y1t2", "y2t2", "y3t2")
var3 <- c("y1t3", "y2t3", "y3t3")
constrainedVar <- list(var1, var2, var3)

# Invariance of the same factor across timepoints
longInvariance(model, auto=1, constrainAuto=TRUE, varList=constrainedVar, data=exLong)

# Invariance of the same factor across timepoints and groups
longInvariance(model, auto=1, constrainAuto=TRUE, varList=constrainedVar, data=exLong, group="sex", group.equal=c("loadings", "intercepts"))



cleanEx()
nameEx("mardiaKurtosis")
### * mardiaKurtosis

flush(stderr()); flush(stdout())

### Name: mardiaKurtosis
### Title: Finding Mardia's multivariate kurtosis
### Aliases: mardiaKurtosis

### ** Examples

library(lavaan)
mardiaKurtosis(HolzingerSwineford1939[,paste("x", 1:9, sep="")])



cleanEx()
nameEx("mardiaSkew")
### * mardiaSkew

flush(stderr()); flush(stdout())

### Name: mardiaSkew
### Title: Finding Mardia's multivariate skewness
### Aliases: mardiaSkew

### ** Examples

library(lavaan)
mardiaSkew(HolzingerSwineford1939[,paste("x", 1:9, sep="")])



cleanEx()
nameEx("measurementInvariance")
### * measurementInvariance

flush(stderr()); flush(stdout())

### Name: measurementInvariance
### Title: Measurement Invariance Tests
### Aliases: measurementInvariance measurementinvariance

### ** Examples

HW.model <- ' visual =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed =~ x7 + x8 + x9 '

measurementInvariance(HW.model, data=HolzingerSwineford1939, group="school")



cleanEx()
nameEx("miPowerFit")
### * miPowerFit

flush(stderr()); flush(stdout())

### Name: miPowerFit
### Title: Modification indices and their power approach for model fit
###   evaluation
### Aliases: miPowerFit miPowerFit

### ** Examples

library(lavaan)

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939, group="sex", meanstructure=TRUE)
miPowerFit(fit)

model <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
fit2 <- sem(model, data=PoliticalDemocracy, meanstructure=TRUE)
miPowerFit(fit2, stdLoad=0.3, cor=0.2, stdBeta=0.2, intcept=0.5)



cleanEx()
nameEx("monteCarloMed")
### * monteCarloMed

flush(stderr()); flush(stdout())

### Name: monteCarloMed
### Title: Monte Carlo Confidence Intervals to Test Complex Indirect
###   Effects
### Aliases: monteCarloMed

### ** Examples

#Simple two path mediation
#Write expression of indirect effect
med <- 'a*b'
#Paramter values from analyses
aparam <- 1
bparam<-2
#Asymptotic covariance matrix from analyses
AC <- matrix(c(.01,.00002,
               .00002,.02), nrow=2, byrow=TRUE)
#Compute CI, include a plot
monteCarloMed(med, coef1=aparam, coef2=bparam, outputValues=FALSE, plot=TRUE, ACM=AC)

#Use a matrix of parameter estimates as input
aparam<-c(1,2)
monteCarloMed(med, coef1=aparam, outputValues=FALSE, plot=TRUE, ACM=AC)



#complex mediation with two paths for the indirect effect
#Write expression of indirect effect
med <- 'a1*b1 + a1*b2'
#Paramter values and standard errors from analyses
aparam <- 1
b1param<-2
b2param<-1
#Asymptotic covariance matrix from analyses
AC <- matrix(c(1,.00002, .00003,
                    .00002,1, .00002,
					.00003, .00002, 1), nrow=3, byrow=TRUE)
#Compute CI do not include a plot
monteCarloMed(med, coef1=aparam, coef2=b1param, coef3=b2param, ACM=AC)



cleanEx()
nameEx("moreFitIndices")
### * moreFitIndices

flush(stderr()); flush(stdout())

### Name: moreFitIndices
### Title: Calculate more fit indices
### Aliases: moreFitIndices

### ** Examples

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
moreFitIndices(fit)

fit2 <- cfa(HS.model, data=HolzingerSwineford1939, estimator="mlr")
moreFitIndices(fit2)



cleanEx()
nameEx("nullMx")
### * nullMx

flush(stderr()); flush(stdout())

### Name: nullMx
### Title: Analyzing data using a null model
### Aliases: nullMx

### ** Examples

## Not run: 
##D library(OpenMx)
##D data(demoOneFactor)
##D nullModel <- nullMx(demoOneFactor)
## End(Not run)



cleanEx()
nameEx("nullRMSEA")
### * nullRMSEA

flush(stderr()); flush(stdout())

### Name: nullRMSEA
### Title: Calculate the RMSEA of the null model
### Aliases: nullRMSEA

### ** Examples

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
nullRMSEA(fit)



cleanEx()
nameEx("parcelAllocation")
### * parcelAllocation

flush(stderr()); flush(stdout())

### Name: parcelAllocation
### Title: Random Allocation of Items to Parcels in a Structural Equation
###   Model
### Aliases: parcelAllocation

### ** Examples

#Fit 3 factor CFA to simulated data.
#Each factor has 9 indicators that are randomly parceled into 3 parcels
#Lavaan syntax for the model to be fit to parceled data
syntax <- 'La =~ V1 + V2 + V3 
           Lb =~ V4 + V5 + V6
'
#Parcel and fit data 20 times. The actual parcel number should be higher than 20 times.
name1 <- colnames(simParcel)[1:9]
name2 <- colnames(simParcel)[10:18]
parcelAllocation(list(c(3,3,3),c(3,3,3)), list(name1, name2), nAlloc=20, syntax=syntax, dataset=simParcel)



cleanEx()
nameEx("plotProbe")
### * plotProbe

flush(stderr()); flush(stdout())

### Name: plotProbe
### Title: Plot the graphs for probing latent interaction
### Aliases: plotProbe

### ** Examples

library(lavaan) 

dat2wayMC <- indProd(dat2way, 1:3, 4:6)

model1 <- "
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
f12 =~ x1.x4 + x2.x5 + x3.x6
f3 =~ x7 + x8 + x9
f3 ~ f1 + f2 + f12
f12 ~~0*f1
f12 ~~ 0*f2
x1 ~ 0*1
x4 ~ 0*1
x1.x4 ~ 0*1
x7 ~ 0*1
f1 ~ NA*1
f2 ~ NA*1
f12 ~ NA*1
f3 ~ NA*1
"

fitMC2way <- sem(model1, data=dat2wayMC, meanstructure=TRUE, std.lv=FALSE)
result2wayMC <- probe2WayMC(fitMC2way, c("f1", "f2", "f12"), "f3", "f2", c(-1, 0, 1))
plotProbe(result2wayMC, xlim=c(-2, 2))


dat3wayMC <- indProd(dat3way, 1:3, 4:6, 7:9)

model3 <- "
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
f3 =~ x7 + x8 + x9
f12 =~ x1.x4 + x2.x5 + x3.x6
f13 =~ x1.x7 + x2.x8 + x3.x9
f23 =~ x4.x7 + x5.x8 + x6.x9
f123 =~ x1.x4.x7 + x2.x5.x8 + x3.x6.x9
f4 =~ x10 + x11 + x12
f4 ~ f1 + f2 + f3 + f12 + f13 + f23 + f123
f1 ~~ 0*f12
f1 ~~ 0*f13
f1 ~~ 0*f123
f2 ~~ 0*f12
f2 ~~ 0*f23
f2 ~~ 0*f123
f3 ~~ 0*f13
f3 ~~ 0*f23
f3 ~~ 0*f123
f12 ~~ 0*f123
f13 ~~ 0*f123
f23 ~~ 0*f123
x1 ~ 0*1
x4 ~ 0*1
x7 ~ 0*1
x10 ~ 0*1
x1.x4 ~ 0*1
x1.x7 ~ 0*1
x4.x7 ~ 0*1
x1.x4.x7 ~ 0*1
f1 ~ NA*1
f2 ~ NA*1
f3 ~ NA*1
f12 ~ NA*1
f13 ~ NA*1
f23 ~ NA*1
f123 ~ NA*1
f4 ~ NA*1
" 

fitMC3way <- sem(model3, data=dat3wayMC, meanstructure=TRUE, std.lv=FALSE)
result3wayMC <- probe3WayMC(fitMC3way, c("f1", "f2", "f3", "f12", "f13", "f23", "f123"), "f4", c("f1", "f2"), c(-1, 0, 1), c(-1, 0, 1))
plotProbe(result3wayMC, xlim=c(-2, 2))



cleanEx()
nameEx("plotRMSEAdist")
### * plotRMSEAdist

flush(stderr()); flush(stdout())

### Name: plotRMSEAdist
### Title: Plot the sampling distributions of RMSEA
### Aliases: plotRMSEAdist

### ** Examples

plotRMSEAdist(rmsea=c(.05, .08), n=200, df=20, ptile=0.95, rmseaScale = TRUE)
plotRMSEAdist(rmsea=c(.05, .01), n=200, df=20, ptile=0.05, rmseaScale = FALSE)



cleanEx()
nameEx("plotRMSEApower")
### * plotRMSEApower

flush(stderr()); flush(stdout())

### Name: plotRMSEApower
### Title: Plot power curves for RMSEA
### Aliases: plotRMSEApower

### ** Examples

plotRMSEApower(.025, .075, 23, 100, 500, 10)



cleanEx()
nameEx("plotRMSEApowernested")
### * plotRMSEApowernested

flush(stderr()); flush(stdout())

### Name: plotRMSEApowernested
### Title: Plot power of nested model RMSEA
### Aliases: plotRMSEApowernested

### ** Examples

plotRMSEApowernested(rmsea0A = 0, rmsea0B = 0, rmsea1A = 0.06, rmsea1B = 0.05, 
dfA=22, dfB=20, nlow=50, nhigh=500, steps=1, alpha=.05, group=1)  



cleanEx()
nameEx("probe2WayMC")
### * probe2WayMC

flush(stderr()); flush(stdout())

### Name: probe2WayMC
### Title: Probing two-way interaction on the residual-centered latent
###   interaction
### Aliases: probe2WayMC

### ** Examples

library(lavaan) 

dat2wayMC <- indProd(dat2way, 1:3, 4:6)

model1 <- "
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
f12 =~ x1.x4 + x2.x5 + x3.x6
f3 =~ x7 + x8 + x9
f3 ~ f1 + f2 + f12
f12 ~~0*f1
f12 ~~ 0*f2
x1 ~ 0*1
x4 ~ 0*1
x1.x4 ~ 0*1
x7 ~ 0*1
f1 ~ NA*1
f2 ~ NA*1
f12 ~ NA*1
f3 ~ NA*1
"

fitMC2way <- sem(model1, data=dat2wayMC, meanstructure=TRUE, std.lv=FALSE)
summary(fitMC2way)

result2wayMC <- probe2WayMC(fitMC2way, c("f1", "f2", "f12"), "f3", "f2", c(-1, 0, 1))
result2wayMC



cleanEx()
nameEx("probe2WayRC")
### * probe2WayRC

flush(stderr()); flush(stdout())

### Name: probe2WayRC
### Title: Probing two-way interaction on the residual-centered latent
###   interaction
### Aliases: probe2WayRC

### ** Examples

library(lavaan) 

dat2wayRC <- orthogonalize(dat2way, 1:3, 4:6)

model1 <- "
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
f12 =~ x1.x4 + x2.x5 + x3.x6
f3 =~ x7 + x8 + x9
f3 ~ f1 + f2 + f12
f12 ~~0*f1
f12 ~~ 0*f2
x1 ~ 0*1
x4 ~ 0*1
x1.x4 ~ 0*1
x7 ~ 0*1
f1 ~ NA*1
f2 ~ NA*1
f12 ~ NA*1
f3 ~ NA*1
"

fitRC2way <- sem(model1, data=dat2wayRC, meanstructure=TRUE, std.lv=FALSE)
summary(fitRC2way)

result2wayRC <- probe2WayRC(fitRC2way, c("f1", "f2", "f12"), "f3", "f2", c(-1, 0, 1))
result2wayRC



cleanEx()
nameEx("probe3WayMC")
### * probe3WayMC

flush(stderr()); flush(stdout())

### Name: probe3WayMC
### Title: Probing two-way interaction on the residual-centered latent
###   interaction
### Aliases: probe3WayMC

### ** Examples

library(lavaan)

dat3wayMC <- indProd(dat3way, 1:3, 4:6, 7:9)

model3 <- "
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
f3 =~ x7 + x8 + x9
f12 =~ x1.x4 + x2.x5 + x3.x6
f13 =~ x1.x7 + x2.x8 + x3.x9
f23 =~ x4.x7 + x5.x8 + x6.x9
f123 =~ x1.x4.x7 + x2.x5.x8 + x3.x6.x9
f4 =~ x10 + x11 + x12
f4 ~ f1 + f2 + f3 + f12 + f13 + f23 + f123
f1 ~~ 0*f12
f1 ~~ 0*f13
f1 ~~ 0*f123
f2 ~~ 0*f12
f2 ~~ 0*f23
f2 ~~ 0*f123
f3 ~~ 0*f13
f3 ~~ 0*f23
f3 ~~ 0*f123
f12 ~~ 0*f123
f13 ~~ 0*f123
f23 ~~ 0*f123
x1 ~ 0*1
x4 ~ 0*1
x7 ~ 0*1
x10 ~ 0*1
x1.x4 ~ 0*1
x1.x7 ~ 0*1
x4.x7 ~ 0*1
x1.x4.x7 ~ 0*1
f1 ~ NA*1
f2 ~ NA*1
f3 ~ NA*1
f12 ~ NA*1
f13 ~ NA*1
f23 ~ NA*1
f123 ~ NA*1
f4 ~ NA*1
" 

fitMC3way <- sem(model3, data=dat3wayMC, meanstructure=TRUE, std.lv=FALSE)
summary(fitMC3way)

result3wayMC <- probe3WayMC(fitMC3way, c("f1", "f2", "f3", "f12", "f13", "f23", "f123"), "f4", c("f1", "f2"), c(-1, 0, 1), c(-1, 0, 1))
result3wayMC



cleanEx()
nameEx("probe3WayRC")
### * probe3WayRC

flush(stderr()); flush(stdout())

### Name: probe3WayRC
### Title: Probing three-way interaction on the residual-centered latent
###   interaction
### Aliases: probe3WayRC

### ** Examples

library(lavaan)

dat3wayRC <- orthogonalize(dat3way, 1:3, 4:6, 7:9)

model3 <- "
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
f3 =~ x7 + x8 + x9
f12 =~ x1.x4 + x2.x5 + x3.x6
f13 =~ x1.x7 + x2.x8 + x3.x9
f23 =~ x4.x7 + x5.x8 + x6.x9
f123 =~ x1.x4.x7 + x2.x5.x8 + x3.x6.x9
f4 =~ x10 + x11 + x12
f4 ~ f1 + f2 + f3 + f12 + f13 + f23 + f123
f1 ~~ 0*f12
f1 ~~ 0*f13
f1 ~~ 0*f123
f2 ~~ 0*f12
f2 ~~ 0*f23
f2 ~~ 0*f123
f3 ~~ 0*f13
f3 ~~ 0*f23
f3 ~~ 0*f123
f12 ~~ 0*f123
f13 ~~ 0*f123
f23 ~~ 0*f123
x1 ~ 0*1
x4 ~ 0*1
x7 ~ 0*1
x10 ~ 0*1
x1.x4 ~ 0*1
x1.x7 ~ 0*1
x4.x7 ~ 0*1
x1.x4.x7 ~ 0*1
f1 ~ NA*1
f2 ~ NA*1
f3 ~ NA*1
f12 ~ NA*1
f13 ~ NA*1
f23 ~ NA*1
f123 ~ NA*1
f4 ~ NA*1
" 

fitRC3way <- sem(model3, data=dat3wayRC, meanstructure=TRUE, std.lv=FALSE)
summary(fitRC3way)

result3wayRC <- probe3WayRC(fitRC3way, c("f1", "f2", "f3", "f12", "f13", "f23", "f123"), "f4", c("f1", "f2"), c(-1, 0, 1), c(-1, 0, 1))
result3wayRC



cleanEx()
nameEx("reliability")
### * reliability

flush(stderr()); flush(stdout())

### Name: reliability
### Title: Calculate reliability values of factors
### Aliases: reliability

### ** Examples

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
reliability(fit)



cleanEx()
nameEx("reliabilityL2")
### * reliabilityL2

flush(stderr()); flush(stdout())

### Name: reliabilityL2
### Title: Calculate the reliability values of a second-order factor
### Aliases: reliabilityL2

### ** Examples

HS.model3 <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 
			  higher =~ visual + textual + speed'

fit6 <- cfa(HS.model3, data=HolzingerSwineford1939)
reliability(fit6) # Should provide a warning for the endogenous variable
reliabilityL2(fit6, "higher")



cleanEx()
nameEx("residualCovariate")
### * residualCovariate

flush(stderr()); flush(stdout())

### Name: residualCovariate
### Title: Residual centered all target indicators by covariates
### Aliases: residualCovariate

### ** Examples

dat <- residualCovariate(attitude, 2:7, 1)



cleanEx()
nameEx("rotate")
### * rotate

flush(stderr()); flush(stdout())

### Name: rotate
### Title: Implement orthogonal or oblique rotation
### Aliases: orthRotate oblqRotate funRotate

### ** Examples

unrotated <- efaUnrotate(HolzingerSwineford1939, nf=3, varList=paste0("x", 1:9), estimator="mlr")

# Orthogonal varimax
out.varimax <- orthRotate(unrotated, method="varimax")
summary(out.varimax, sort=FALSE, suppress=0.3)

# Orthogonal Quartimin
orthRotate(unrotated, method="quartimin")

# Oblique Quartimin
oblqRotate(unrotated, method="quartimin")

# Geomin
oblqRotate(unrotated, method="geomin")

## Not run: 
##D # Target rotation
##D library(GPArotation)
##D target <- matrix(0, 9, 3)
##D target[1:3, 1] <- NA
##D target[4:6, 2] <- NA
##D target[7:9, 3] <- NA
##D colnames(target) <- c("factor1", "factor2", "factor3")
##D # This function works with GPArotation version 2012.3-1
##D funRotate(unrotated, fun="targetQ", Target=target) 
## End(Not run)



cleanEx()
nameEx("runMI")
### * runMI

flush(stderr()); flush(stdout())

### Name: runMI
### Title: Multiply impute and analyze data using lavaan
### Aliases: runMI cfa.mi sem.mi growth.mi lavaan.mi

### ** Examples

library(lavaan)

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

HSMiss <- HolzingerSwineford1939[,paste("x", 1:9, sep="")]
randomMiss <- rbinom(prod(dim(HSMiss)), 1, 0.1)
randomMiss <- matrix(as.logical(randomMiss), nrow=nrow(HSMiss))
HSMiss[randomMiss] <- NA

out <- cfa.mi(HS.model, data=HSMiss, m = 3, chi="all")
summary(out)
inspect(out, "fit")
inspect(out, "impute")

## Not run: 
##D ##Multiple group example
##D HSMiss2 <- cbind(HSMiss, school = HolzingerSwineford1939[,"school"])
##D out2 <- cfa.mi(HS.model, data=HSMiss2, m = 3, miArgs=list(noms="school"), chi="MR", group="school")
##D summary(out2)
##D inspect(out2, "fit")
##D inspect(out2, "impute")
##D 
##D ##Example using previously imputed data with runMI
##D library(Amelia)
##D 
##D modsim <- '
##D f1 =~ 0.7*y1+0.7*y2+0.7*y3
##D f2 =~ 0.7*y4+0.7*y5+0.7*y6
##D f3 =~ 0.7*y7+0.7*y8+0.7*y9'
##D 
##D mod <- '
##D f1 =~ y1+y2+y3
##D f2 =~ y4+y5+y6
##D f3 =~ y7+y8+y9'
##D 
##D datsim <- simulateData(modsim,model.type="cfa", meanstructure=TRUE, 
##D 	std.lv=TRUE, sample.nobs=c(200,200))
##D randomMiss2 <- rbinom(prod(dim(datsim)), 1, 0.1)
##D randomMiss2 <- matrix(as.logical(randomMiss2), nrow=nrow(datsim))
##D datsim[randomMiss2] <- NA
##D datsimMI <- amelia(datsim,m=3, noms="group")
##D 
##D out3 <- runMI(mod, data=datsimMI$imputations, chi="LMRR", group="group", fun="cfa")
##D summary(out3)
##D inspect(out3, "fit")
##D inspect(out3, "impute")
##D 
##D # Categorical variables
##D dat <- simulateData(popModel, sample.nobs  = 200L)
##D miss.pat <- matrix(as.logical(rbinom(prod(dim(dat)), 1, 0.2)), nrow(dat), ncol(dat))
##D dat[miss.pat] <- NA
##D out5 <- cfa.mi(analyzeModel, data=dat, ordered=paste0("y", 1:4), m = 3, miArgs=list(ords = c("y1", "y2", "y3", "y4")))
##D summary(out5)
##D inspect(out5, "fit")
##D inspect(out5, "impute")
##D 
## End(Not run)



cleanEx()
nameEx("saturateMx")
### * saturateMx

flush(stderr()); flush(stdout())

### Name: saturateMx
### Title: Analyzing data using a saturate model
### Aliases: saturateMx

### ** Examples

## Not run: 
##D library(OpenMx)
##D data(demoOneFactor)
##D satModel <- saturateMx(demoOneFactor)
## End(Not run)



cleanEx()
nameEx("simParcel")
### * simParcel

flush(stderr()); flush(stdout())

### Name: simParcel
### Title: Simulated Data set to Demonstrate Random Allocations of Parcels
### Aliases: simParcel

### ** Examples

head(simParcel)



cleanEx()
nameEx("skew")
### * skew

flush(stderr()); flush(stdout())

### Name: skew
### Title: Finding skewness
### Aliases: skew

### ** Examples

skew(1:5)



cleanEx()
nameEx("splitSample")
### * splitSample

flush(stderr()); flush(stdout())

### Name: splitSample
### Title: Randomly Split a Data Set into Halves
### Aliases: splitSample

### ** Examples

#### Input is .dat file
#splitSample("C:/Users/Default/Desktop/MYDATA.dat")
#### Output saved to "C:/Users/Default/Desktop/" in .dat format
#### Names are "MYDATA_s1.dat" and "MYDATA_s2.dat"

#### Input is R object
##Split C02 dataset from the datasets package
library(datasets)
splitMyData <- splitSample(CO2, path="object")
summary(splitMyData[[1]])
summary(splitMyData[[2]])
#### Output object splitMyData becomes list of output data sets

#### Input is .dat file in "C:/" folder
#splitSample("C:/testdata.dat", path = "C:/Users/Default/Desktop/", type = "csv")
#### Output saved to "C:/Users/Default/Desktop/" in .csv format
#### Names are "testdata_s1.csv" and "testdata_s2.csv"

#### Input is R object
#splitSample(myData, path = "C:/Users/Default/Desktop/", name = "splitdata")
#### Output saved to "C:/Users/Default/Desktop/" in .dat format
#### Names are "splitdata_s1.dat" and "splitdata_s2.dat"



cleanEx()
nameEx("standardizeMx")
### * standardizeMx

flush(stderr()); flush(stdout())

### Name: standardizeMx
### Title: Find standardized estimates for OpenMx output
### Aliases: standardizeMx

### ** Examples

## Not run: 
##D library(OpenMx)
##D data(myFADataRaw)
##D myFADataRaw <- myFADataRaw[,c("x1","x2","x3","x4","x5","x6")]
##D oneFactorModel <- mxModel("Common Factor Model Path Specification", 
##D 	type="RAM",
##D 	mxData(
##D 		observed=myFADataRaw, 
##D 		type="raw"
##D 	),
##D 	manifestVars=c("x1","x2","x3","x4","x5","x6"),
##D 	latentVars="F1",
##D 	mxPath(from=c("x1","x2","x3","x4","x5","x6"),
##D 		arrows=2,
##D 		free=TRUE,
##D 		values=c(1,1,1,1,1,1),
##D 		labels=c("e1","e2","e3","e4","e5","e6")
##D 	), 
##D 	# residual variances
##D 	# -------------------------------------
##D 	mxPath(from="F1",
##D 		arrows=2,
##D 		free=TRUE,
##D 		values=1,
##D 		labels ="varF1"
##D 	), 
##D 	# latent variance
##D 	# -------------------------------------
##D 	mxPath(from="F1",
##D 		to=c("x1","x2","x3","x4","x5","x6"),
##D 		arrows=1,
##D 		free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
##D 		values=c(1,1,1,1,1,1),
##D 		labels =c("l1","l2","l3","l4","l5","l6")
##D 	), 
##D 	# factor loadings
##D 	# -------------------------------------
##D 	mxPath(from="one",
##D 		to=c("x1","x2","x3","x4","x5","x6","F1"),
##D 		arrows=1,
##D 		free=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE),
##D 		values=c(1,1,1,1,1,1,0),
##D 		labels =c("meanx1","meanx2","meanx3","meanx4","meanx5","meanx6",NA)
##D 	) 
##D 	# means
##D 	# -------------------------------------
##D ) # close model
##D # Create an MxModel object
##D # -----------------------------------------------------------------------------
##D oneFactorFit <- mxRun(oneFactorModel)      
##D standardizeMx(oneFactorFit)
##D 
##D # Compare with lavaan
##D library(lavaan)
##D script <- "f1 =~ x1 + x2 + x3 + x4 + x5 + x6"
##D fit <- cfa(script, data=myFADataRaw, meanstructure=TRUE)
##D standardize(fit)
## End(Not run)



cleanEx()
nameEx("tukeySEM")
### * tukeySEM

flush(stderr()); flush(stdout())

### Name: tukeySEM
### Title: Tukey's WSD post-hoc test of means for unequal variance and
###   sample size
### Aliases: tukeySEM

### ** Examples

##For a case where three groups have been compared:
##Group 1: mean = 3.91, var = 0.46, n = 246
##Group 2: mean = 3.96, var = 0.62, n = 465
##Group 3: mean = 2.94, var = 1.07, n = 64

#compare group 1 and group 2
tukeySEM(3.91, 3.96, 0.46, 0.62, 246, 425, 3)

#compare group 1 and group 3
tukeySEM(3.91, 2.94, 0.46, 1.07, 246, 64, 3)

#compare group 2 and group 3
tukeySEM(3.96, 2.94, 0.62, 1.07, 465, 64, 3)



cleanEx()
nameEx("wald")
### * wald

flush(stderr()); flush(stdout())

### Name: wald
### Title: Calculate multivariate Wald statistics
### Aliases: wald

### ** Examples

# Test the difference in factor loadings
HS.model <- ' visual  =~ x1 + con1*x2 + con1*x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + con2*x8 + con2*x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
wald(fit, "con2 - con1")

# Simultaneously test the difference in the influences 
# of x1 and x2 on intercept and slope
model.syntax <- '
    i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
    s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
    i ~ x1 + x2
    s ~ x1 + x2
    t1 ~ c1
    t2 ~ c2
    t3 ~ c3
    t4 ~ c4
'

fit2 <- growth(model.syntax, data=Demo.growth)
wald.syntax <- '
	i~x1 - i~x2
	1/2*s~x1 - 1/2*s~x2
'
wald(fit2, wald.syntax)

# Mplus example of MODEL TEST
model3 <- ' f1  =~ x1 + p2*x2 + p3*x3 + p4*x4 + p5*x5 + p6*x6
			p4 == 2*p2'

fit3 <- cfa(model3, data=HolzingerSwineford1939)
wald(fit3, "p3; p6 - 0.5*p5")



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
