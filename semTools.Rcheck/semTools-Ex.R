pkgname <- "semTools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('semTools')

assign(".oldSearch", search(), pos = 'CheckExEnv')
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



cleanEx()
nameEx("orthogonalize")
### * orthogonalize

flush(stderr()); flush(stdout())

### Name: orthogonalize
### Title: Orthogonalize data for 2-way interaction in SEM
### Aliases: orthogonalize

### ** Examples


library(MASS)

n <- 500
means <- c(0,0)
covmat <- matrix(c(1, 0.3, 0.3, 1),nrow=2)

data <- mvrnorm(n,means,covmat)

x<-as.vector(data[,1])
z<-as.vector(data[,2])

y<-rnorm(n,0,1)+.4*x+.4*z+.2*x*z

x1<-rnorm(n,0.2,.2)+.7*x
x2<-rnorm(n,0.2,.2)+.7*x
x3<-rnorm(n,0.2,.2)+.7*x
z1<-rnorm(n,0.2,.2)+.7*z
z2<-rnorm(n,0.2,.2)+.7*z
z3<-rnorm(n,0.2,.2)+.7*z
y1<-rnorm(n,0.2,.2)+.7*y
y2<-rnorm(n,0.2,.2)+.7*y
y3<-rnorm(n,0.2,.2)+.7*y

dat<-data.frame(cbind(x1,x2,x3,z1,z2,z3,y1,y2,y3))

datOrth <-orthogonalize(dat,(1:3), (4:6))

#Fit model in Lavaan
library(lavaan)

syntax <- ' 
x =~ x1 + x2 +x3
z =~ z1 + z2 + z3
xz =~ x1z1 + x1z2 + x1z3 + x2z1 + x2z2 + x2z3 + x3z1 + x3z2 + x3z3
y =~ y1 + y2 + y3
x ~~ z
x ~~ 0*xz
z ~~ 0*xz
y ~ x + z +xz
'

fit <- sem(model = syntax, data=datOrth, std.lv=TRUE)
summary(fit, fit.measures=TRUE)



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
nameEx("runMI")
### * runMI

flush(stderr()); flush(stdout())

### Name: runMI
### Title: Multiply impute and analyze data using lavaan
### Aliases: runMI

### ** Examples

## No Example




cleanEx()
nameEx("skew")
### * skew

flush(stderr()); flush(stdout())

### Name: skew
### Title: Finding skewness
### Aliases: skew

### ** Examples

skew(1:5)



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
