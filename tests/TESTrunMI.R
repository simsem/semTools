### Terrence D. Jorgensen
### Last updated: 28 August 2019
### test runMI

library(lavaan)
library(Amelia)

# source("../FuturePlans/NEWrunMI.R")
library(semTools)


## --------
## Examples
## --------

set.seed(12345)
HSMiss <- HolzingerSwineford1939[ , paste("x", 1:9, sep = "")]
HSMiss$x5 <- ifelse(HSMiss$x1 <= quantile(HSMiss$x1, .3), NA, HSMiss$x5)
HSMiss$x9 <- ifelse(is.na(HSMiss$x5), NA, HSMiss$x9)
HSMiss$school <- HolzingerSwineford1939$school
HS.amelia <- amelia(HSMiss, m = 20, noms = "school")
imps <- HS.amelia$imputations

HS.model <- '
visual  =~ x1 + a*x2 + b*x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9
ab := a*b
'
## fit single-group model
fit1 <- cfa.mi(HS.model, data = imps, std.lv = TRUE, meanstructure = TRUE)
fit0 <- cfa.mi(HS.model, data = fit1, std.lv = TRUE, meanstructure = TRUE,
               #se = "none",
               #estimator = "mlm",
               constraints = '.p2. == .p3. ; .p5. == .p6. ; .p8. == .p9.')
## use methods
summary(fit1, stand = TRUE, rsq = TRUE) # mimics parameterEstimates()
coef(fit1)           # pooled coefs
vcov(fit1)[1:4, 1:4] # pooled sampling covariance matrix
fitted(fit1)         # model-implied moments evaluated at pooled coefficients
resid(fit1, type = "cor.bentler") # (average sampstats) - fitted
nobs(fit1)

lavTestLRT.mi(fit1)                     # pooled LRT by default (D3 statistic)
lavTestLRT.mi(fit1, test = "D2")        # or pooled Wald test (D2 statistic)
lavTestLRT.mi(fit1, asymptotic = TRUE)  # as chisq == F * df1
lavTestLRT.mi(fit0, h1 = fit1, asymptotic = TRUE)  # compare fit
lavTestLRT.mi(fit0, h1 = fit1, test = "D2")  # compare fit using D2
fitMeasures(fit1)
fitMeasures(fit1, output = "text")

lavTestScore.mi(fit0, epc = TRUE)
modindices.mi(fit1)



## fit multigroup model
mgfit1 <- cfa.mi(HS.model, data = imps, group = "school", estimator = "mlm")
mgfit0 <- cfa.mi(HS.model, data = imps, group = "school", estimator = "mlm",
                 group.equal = c("loadings","intercepts"))
## use methods
summary(mgfit0, standardized = "std.all") # can also request "std.lv"
summary(mgfit0, ci = FALSE, fmi = TRUE, asymptotic = TRUE)
coef(mgfit0)           # pooled coefs
vcov(mgfit0)[1:4, 1:4] # pooled sampling covariance matrix
fitted(mgfit0)         # model-implied moments evaluated at pooled coefficients
resid(mgfit0, type = "cor.bentler") # (average sampstats) - fitted
nobs(mgfit0)
nobs(mgfit0, total = FALSE) # N per group

lavTestLRT.mi(mgfit0)              # pooled LRT by default (D3 statistic)
lavTestLRT.mi(mgfit0, test = 'D2') # use D2 method (necessary for categorical data)
lavTestLRT.mi(mgfit0, test = 'D2', pool.robust = TRUE)
lavTestLRT.mi(mgfit0, test = 'D2', pool.robust = TRUE, asymptotic = TRUE)
lavTestLRT.mi(mgfit0, h1 = mgfit1)           # compare fit
lavTestLRT.mi(mgfit0, h1 = mgfit1, test = 'D2') # robustifies pooled naive statistic
lavTestLRT.mi(mgfit0, h1 = mgfit1, test = 'D2', pool.robust = TRUE) # pools robust statistic
fitMeasures(mgfit1, output = "text")
print(fitMeasures(mgfit0, output = "text", test = "D2"), add.h0 = TRUE)

## use D1 to test a parametrically nested model (whether latent means are ==)
lavTestWald.mi(mgfit0, test = "D1", asymptotic = TRUE, constraints = '
      .p70. == 0
      .p71. == 0
      .p72. == 0')
## compare to complete data and FIML
mgfit <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
             estimator = "mlr", group.equal = c("loadings","intercepts"))
mgfitm <- cfa(HS.model, data = HSMiss, group = "school", missing = "fiml",
              estimator = "mlr", group.equal = c("loadings","intercepts"))
lavTestWald(mgfit, constraints = '
            .p70. == 0
            .p71. == 0
            .p72. == 0')
lavTestWald(mgfitm, constraints = '
            .p70. == 0
            .p71. == 0
            .p72. == 0')


## Score test for metric invariance
lavTestScore.mi(mgfit0, release = 1:6, test = "D1", scale.W = TRUE)
lavTestScore.mi(mgfit0, release = 1:6, test = "D1", scale.W = FALSE)
## compare univariate score test to modification indices
lavTestScore.mi(mgfit1, add = 'x1 ~~ x9 ; x4 ~~ x7', epc = TRUE,
                test = "D1", asymptotic = TRUE)
modindices.mi(mgfit1, op = '~~', test = "D1")
## match perfectly! (not if scale.W = TRUE)
## also comparable with D2 pooled results:
lavTestScore.mi(mgfit1, add = 'x1 ~~ x9 ; x4 ~~ x7', asymptotic = TRUE)
modindices.mi(mgfit1, op = '~~')



## -------------------------
## fixed.x and conditional.x
## -------------------------

## Using help-page example
data(datCat)
datCat$u2 <- as.integer(datCat$u2) # mixture of ordered/continuous indicators
datCat$u3 <- as.integer(datCat$u3)
datCat$u5 <- as.integer(datCat$u5) # exogenous predictors must be numeric
datCat$u6 <- as.integer(datCat$u6)

## impose missing values
set.seed(456)
for (i in 1:8) datCat[sample(1:nrow(datCat), size = .1*nrow(datCat)), i] <- NA
library(Amelia)
catimps <- amelia(datCat, m = 20, ords = paste0("u", 1:8), noms = "g")$imputations

catmod <- '
f =~ 1*u1 + 1*u2 + 1*u3 + 1*u4
u1 + u2 ~ u5 + u6
'
fitx <- cfa.mi(catmod, data = catimps, fixed.x = TRUE, conditional.x = FALSE)
fitxg <- cfa.mi(catmod, data = catimps, fixed.x = TRUE, conditional.x = FALSE,
                group = "g")
fit.x <- cfa.mi(catmod, data = catimps, fixed.x = TRUE, conditional.x = TRUE)
fit.xg <- cfa.mi(catmod, data = catimps, fixed.x = TRUE, conditional.x = TRUE,
                 group = "g")

fitted(fitx)
resid(fitx)
resid(fitx, type = "crmr")
resid(fitx, type = "srmr")
fitMeasures(fitx, fit.measures = "rmr")

fitted(fitxg)
resid(fitxg)
resid(fitxg, type = "crmr")
resid(fitxg, type = "srmr")
fitMeasures(fitxg, fit.measures = "rmr")

fitted(fit.x)
resid(fit.x)
resid(fit.x, type = "crmr")
resid(fit.x, type = "srmr")
fitMeasures(fit.x, fit.measures = "rmr")

fitted(fit.xg)
resid(fit.xg)
resid(fit.xg, type = "crmr")
resid(fit.xg, type = "srmr")
fitMeasures(fit.xg, fit.measures = "rmr")

summary(fitx, stand=TRUE) # only problem left: standardizing requires cov.x
summary(fitxg, stand=TRUE) # only problem left: standardizing requires cov.x
summary(fit.x, stand=TRUE) # only problem left: standardizing requires cov.x
summary(fit.xg, stand=TRUE) # only problem left: standardizing requires cov.x


## ------------------------------------
## check lavaan.mi with multilevel data
## ------------------------------------

library(lavaan)

data(Demo.twolevel)
Demo.twolevel$id <-  paste0("id", Demo.twolevel$cluster) # character IDs
## assign clusters to arbitrary groups
Demo.twolevel$g <- ifelse(Demo.twolevel$cluster %% 2L, "foo", "bar")
## randomize rows so cluster IDs are out of order
set.seed(123)
Demo.twolevel <- Demo.twolevel[sample(nrow(Demo.twolevel)), ]
## create missing data
set.seed(456)
Demo.twolevel$y1[sample(nrow(Demo.twolevel), size = .1*nrow(Demo.twolevel)) ] <- NA

model <- '
level: within
  fw =~ y1 + L2w*y2 + L3w*y3
  fw ~ x1 + x2 + x3
level: between
  fb =~ y1 + L2b*y2 + L3b*y3
  fb ~ b1*w1 + b2*w2

## constraints
L2w == L2b
L3w == L3b
'
model2 <- ' group: foo
level: within
  fw =~ y1 + L2*y2 + L3*y3
  fw ~ x1 + x2 + x3
level: between
  fb =~ y1 + L2*y2 + L3*y3
  fb ~ w1 + w2

group: bar

level: within
  fw =~ y1 + L2*y2 + L3*y3
  fw ~ x1 + x2 + x3
level: between
  fb =~ y1 + L2*y2 + L3*y3
  fb ~ w1 + w2
'

fit2 <- sem(model, data = Demo.twolevel, cluster = "cluster")
mgfit2 <- sem(model2, data = Demo.twolevel, cluster = "id", group = "g")

summary(fit2, fit.measures = TRUE, standardized = TRUE)
lavInspect(fit2, "ngroups")
lavInspect(fit2, "nlevels")
lavInspect(fit2, "cluster")
lavInspect(fit2, "nclusters")
lavInspect(fit2, "ncluster.size")
lavInspect(fit2, "cluster.size")
lavInspect(fit2, "cluster.sizes")
lavInspect(fit2, "average.cluster.size")
lavInspect(fit2, "cluster.id")
lavInspect(fit2, "cluster.idx")
lavInspect(fit2, "cluster.label")

lavInspect(fit2, "cov.lv")
lavInspect(fit2, "theta")

## impute data
library(mice)
m <- 5
mice.out <- mice(Demo.twolevel, m = m, diagnostics = TRUE)
imputedData <- list()
for (i in 1:m) {
  imputedData[[i]] <- complete(data = mice.out, action = i, include = FALSE)
}

library(semTools)
fit2 <- sem(model, data = imputedData[[1]], cluster = "id")
fitList2 <- semList(model, dataList = imputedData, cluster = "id",
                    store.slots = "baseline")
fit2.mi <- sem.mi(model, data = imputedData,
                  group = "g",
                  cluster = "id")
## check methods
fit2.mi
summary(fit2.mi, ci = TRUE, standardized = TRUE, rsquare = TRUE, fmi = TRUE)
coef(fit2.mi)
vcov(fit2.mi)
nobs(fit2.mi)
fitted(fit2.mi)
resid(fit2.mi, type = "cor")
modindices.mi(fit2.mi)
lavTestScore.mi(fit2.mi, test = "D2")
lavTestScore.mi(fit2.mi, test = "D1") # very different results
lavTestWald.mi(fit2.mi, constraints = 'b1 == b2')
lavTestWald.mi(fit2.mi, constraints = 'b1 == b2', test = "D2")
lavTestLRT.mi(fit2.mi, test = "D2")
anova(fit2.mi)
anova(fit2.mi, asymptotic = TRUE) # comparable to anova(fit)
fitMeasures(fit2.mi)
## custom baseline model
basemod <- '
level: within
  y1 ~~ y1
  y2 ~~ y2
  y3 ~~ y3
  y1 + y2 + y3 ~ x1 + x2 + x3
level: between
  y1 ~~ y1
  y2 ~~ y2
  y3 ~~ y3
  y1 + y2 + y3 ~ 1 + w1 + w2
'
base.mi <- lavaan.mi(basemod, data = imputedData, cluster = "id")
fitMeasures(fit2.mi, baseline.model = base.mi, test = "D2") # works now



## -------------------------------------
## check lavaan.mi with sampling weights
## -------------------------------------

library(lavaan)
example(cfa)
## use agemo as weights (internally rescaled)
fit.w <- update(fit, sampling.weights = "agemo")
## coef method
rbind(unweighted = coef(fit), weighted = coef(fit.w))
## parTable
rbind(unweighted = parTable(fit)$est, weighted = parTable(fit.w)$est)
## GLIST format
lavInspect(fit, "est")$lambda
lavInspect(fit.w, "est")$lambda


## bootstrap with lavaanList
bootFun <- function() {
  N <- nrow(lavaan::HolzingerSwineford1939)
  idx <- sample(1:N, size = N, replace = TRUE)
  lavaan::HolzingerSwineford1939[idx, ]
}
set.seed(123)
lst <- cfaList(HS.model, dataFunction = bootFun, ndat = 5)
set.seed(123)
lst.w <- cfaList(HS.model, dataFunction = bootFun, ndat = 5,
                 sampling.weights = "agemo")
## coef method
data.frame(unweighted = coef(lst), weighted = coef(lst.w))
## parTable
data.frame(lst@ParTableList[[1]][c("lhs","op","rhs","est")],
           weighted = lst.w@ParTableList[[1]]$est)
## same.  Problem:
fit.w@Options$sampling.weights


## multiple imputations
HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
                                      "ageyr","agemo","school")]
set.seed(12345)
HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
age <- HSMiss$ageyr + HSMiss$agemo/12
HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
## impute missing data first
library(Amelia)
set.seed(12345)
HS.amelia <- amelia(HSMiss, m = 20, noms = "school", p2s = FALSE)
imps <- HS.amelia$imputations
## fit model both ways
library(semTools)
mi <- cfa.mi(HS.model, data = imps)
mi.w <- cfa.mi(HS.model, data = imps, sampling.weights = "agemo")
data.frame(unweighted = coef(mi), weighted = coef(mi.w))



## ----------------------------------------
## Compare to old runMI output and to Mplus
## ----------------------------------------

library(semTools)
mgfit0.old <- semTools::runMI(HS.model, data = imps, group = "school", fun = "cfa",
                               group.equal = c("loadings","intercepts"))
inspect(mgfit0.old, "impute")



library(MplusAutomation)

setwd("FuturePlans/MplusTest/")

## save imputed data in Mplus format
for (i in names(imps)) prepareMplusData(imps[[i]], filename = paste0(i, ".dat"))
## save list of imputed data files
cat(paste0(names(imps), ".dat\n"), file = "imps.dat")

## write single-group model to file
cat('TITLE: single-group 3-factor model
DATA:
  FILE = "imps.dat";
  TYPE IS IMPUTATION;
VARIABLE:
  NAMES = x1 x2 x3 x4 x5 x6 x7 x8 x9 school;
  USEVAR = x1-x9;
MODEL:
  f1 BY x1* x2 x3;
  f2 BY x4* x5 x6;
  f3 BY x7* x8 x9;
  f1@1 f2@1 f3@1;
', file = "single.inp")

## write multigroup model to file
cat('TITLE: multigroup scalar invariance model
DATA:
  FILE = "imps.dat";
  TYPE IS IMPUTATION;
VARIABLE:
  NAMES = x1 x2 x3 x4 x5 x6 x7 x8 x9 school;
  USEVAR = x1-x9;
  GROUPING = school (1 = Grant_White   2 = Pasteur);
MODEL:
  f1 BY x1-x3;
  f2 BY x4-x6;
  f3 BY x7-x9;
MODEL Grant_White:
  [f1 f2 f3];
MODEL Pasteur:
  [f1@0 f2@0 f3@0];
', file = "scalar.inp")

## run the model, extract fit measures
runModels(logFile = NULL)



## -----------------------------------------------------
## Investigate problem with CFI when test is MV-adjusted
## https://github.com/simsem/semTools/issues/29
## -----------------------------------------------------

library(mice)

# Data obtained from http://openpsychometrics.org/
# Load Rosenberg self-esteem scale data from http://openpsychometrics.org/_rawdata/RSE.zip
# Only first 10000 participants for faster computation
rosenberg.data <- read.table("tests/data/data.txt", header = TRUE)[1:10000,1:10]

## impose 50% missing values
# set.seed(123)
# for (COL in 1:10) {
#   ROW <- sample(1:nrow(rosenberg.data), size = .5*nrow(rosenberg.data))
#   rosenberg.data[ROW, COL] <- NA
# }
#
# # Impute with mice
nImputations <- 5
# imp <- mice(rosenberg.data, m = nImputations)
# # Create list of imputated data sets
# impList <- list()
# for (i in 1:nImputations) impList[[i]] <- complete(imp, action = i)
# saveRDS(impList, "tests/data/imps.rds")
impList <- readRDS("tests/data/imps.rds")

## estimate the fraction of missing information for summary stats
fmi(impList, ords = paste0("Q", 1:10))


# CFA model with all items loading on one factor
rosenberg.model <- '
A =~ Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + Q8 + Q9 + Q10
'
# Fit CFA across all imputed datasets and pool results
# All items are treated as ordinal
rosenberg.fit <- cfa.mi(rosenberg.model, data = impList, estimator = "WLSMV",
                        ordered = paste0("Q", 1:10), FUN = fitMeasures)
# Obtain pooled fit indices
pooled.indices <- anova(rosenberg.fit, indices = TRUE)
# Loop over all imputed datasets and extract fit indices
indices <- do.call(rbind, lapply(rosenberg.fit@funList, function(x) x$userFUN1))

### Compare pooled fit indices to average fit indices

# Average fit indices across imputations (cfi.scaled = .948, rmsea.scaled = .175)
round(colMeans(indices), 3)[grep("cfi|rmsea|chisq", colnames(indices))]
# Fit indices from pooled statistics (cfi.scaled = .052, rmsea.scaled = .083)
round(pooled.indices[grep("cfi|rmsea|chisq", names(pooled.indices))], 3)


## check between-imputation variability of chi-squareds
round(apply(indices[ , grep("chisq", colnames(indices))], 2, sd), 3)

## average relative increase in variance (ARIV)
apply(indices[ , grep("chisq", colnames(indices))], 2,
      function(w) (1 + 1/nImputations) * var(sqrt(w)))

## calculate D2 for baseline.chisq manually from the quantities above
X2 <- 404005.7773054838 # print(colMeans(indices)["baseline.chisq"], 16)
DF <- 45 # degrees of freedom for baseline model
M <- (nImputations + 1) / (nImputations - 1)
ARIV <- 77.65339227 # average relative increase in variance
((X2/DF - M*ARIV) / (1 + ARIV)) * DF # this D2 matches pooled.indices output
pooled.indices["baseline.chisq"]
