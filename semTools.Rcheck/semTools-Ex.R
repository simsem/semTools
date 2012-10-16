pkgname <- "semTools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('semTools')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("auxiliary")
### * auxiliary

flush(stderr()); flush(stdout())

### Name: auxiliary
### Title: Analyzing data with full-information maximum likelihood with
###   auxiliary variables
### Aliases: auxiliary

### ** Examples

# Example of confirmatory factor analysis

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
			  
dat <- data.frame(HolzingerSwineford1939, z=rnorm(nrow(HolzingerSwineford1939), 0, 1))
			  
fit <- cfa(HS.model, data=dat) #, group="sex", meanstructure=TRUE)
fitaux <- auxiliary(fit, aux="z", data=dat)

# Example of multiple groups confirmatory factor analysis

fitgroup <- cfa(HS.model, data=dat, group="school")
fitgroupaux <- auxiliary(fitgroup, aux="z", data=dat, group="school")

# Example of path analysis

mod <- ' x5 ~ x4
x4 ~ x3
x3 ~ x1 + x2'

fitpath <- sem(mod, data=dat)
fitpathaux <- auxiliary(fitpath, aux="z", data=dat)

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
fitsemaux <- auxiliary(fitsem, aux="z", data=dat2, meanstructure=TRUE)

#########################################################
## These following codes show the models that do not work with the current function

##### 1. covariate at the factor level
## HS.model.cov <- ' visual  =~ x1 + x2 + x3
##              textual =~ x4 + x5 + x6
##             speed   =~ x7 + x8 + x9 
##			  visual ~ sex
##			  textual ~ sex
##			  speed ~ sex'
	  
## fitcov <- cfa(HS.model.cov, data=dat) 
## fitcovaux <- auxiliary(fitcov, aux="z", data=dat)

### The auxiliary code does not work when specifying manually.
## HS.model.covxx <- ' visual  =~ x1 + x2 + x3
##              textual =~ x4 + x5 + x6
##              speed   =~ x7 + x8 + x9 
##			  visual ~ sex
##			  textual ~ sex
##			  speed ~ sex
##			  z ~~ z
##			  z ~~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
##			  z ~ sex'
	  
## fitcovxx <- cfa(HS.model.covxx, data=dat) 
## fitcovaux <- auxiliary(fitcov, aux="z", data=dat)

##### 2. Endogenous variable with single indicator 
## HS.model.cov2 <- ' visual  =~ x1 + x2 + x3
##               textual =~ x4 + x5 + x6
##               x7 ~ visual + textual'
## 	  
## fitcov2 <- sem (HS.model.cov2, data=dat, fixed.x=FALSE) #, group="sex", meanstructure=TRUE)
## fitcov2aux <- auxiliary(fitcov2, aux="z", data=dat)


### The auxiliary code does not work when specifying manually.
## HS.model.covyy <- ' visual  =~ x1 + x2 + x3
##               textual =~ x4 + x5 + x6
##               x7 ~ visual + textual
## 			  z ~~ x1 + x2 + x3 + x4 + x5 + x6 + x7'
## fitcovyy <- sem(HS.model.covyy, data=dat) #, group="sex", meanstructure=TRUE)
		  



cleanEx()
nameEx("clipboard")
### * clipboard

flush(stderr()); flush(stdout())

### Name: clipboard_saveFile
### Title: Copy or save the result of 'lavaan' object into a clipboard or a
###   file
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
			  
fit <- cfa(HS.model, data=dat) #, group="sex", meanstructure=TRUE)
fitaux <- auxiliary(fit, aux="z", data=dat)



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
### Aliases: moreFitIndices moreFitIndices

### ** Examples

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
moreFitIndices(fit)

fit2 <- cfa(HS.model, data=HolzingerSwineford1939, estimator="mlr")
moreFitIndices(fit2)



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
nameEx("residualCovariate")
### * residualCovariate

flush(stderr()); flush(stdout())

### Name: residualCovariate
### Title: Residual centered all target indicators by covariates
### Aliases: residualCovariate

### ** Examples

dat <- residualCovariate(attitude, 2:7, 1)



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
