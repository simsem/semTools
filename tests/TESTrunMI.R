### Terrence D. Jorgensen
### Last updated: 16 March 2019
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

lavTestLRT.mi(mgfit0)          # pooled LRT by default (D3 statistic)
lavTestLRT.mi(mgfit0, asymptotic = TRUE)  # reported by Mplus
lavTestLRT.mi(mgfit0, test = 'D2')        # use D2 method (necessary for categorical data)
lavTestLRT.mi(mgfit0, h1 = mgfit1)           # compare fit
lavTestLRT.mi(mgfit0, h1 = mgfit1, test = 'D2') # robustifies pooled naive statistic
lavTestLRT.mi(mgfit0, h1 = mgfit1, test = 'D2', pool.robust = TRUE) # pools robust statistic


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



## ------------------------------------
## check lavaan.mi with multilevel data
## ------------------------------------

Galo <- read.table("https://users.ugent.be/~yrosseel/lavaan/kiel2018/Galo.dat")
names(Galo) <- c("school", "sex", "galo", "advice", "feduc", "meduc",
                 "focc", "denom")

# missing data
Galo[Galo == 999] <- NA
Galo$denom1 <- ifelse(Galo$denom == 1, 1, 0)
Galo$denom2 <- ifelse(Galo$denom == 2, 1, 0)
Galo$g <- ifelse(Galo$school < 30, "foo", "bar") # randomly assign schools to groups
Galo$sch <- paste0("sch", Galo$school) # character IDs
Galo <- Galo[sample(nrow(Galo)), ] # randomize rows so cluster IDs are out of order

model <- '
level: within
  wses =~ a*focc + b*meduc + c*feduc
  advice ~ wc*wses + wb*galo
  galo   ~ wa*wses
  # residual correlation
  focc ~~ feduc
level: between
  bses =~ a*focc + b*meduc + c*feduc
  advice ~ bc*bses + bb*galo
  galo   ~ ba*bses + denom1 + denom2
  feduc ~~ 0*feduc
# defined parameters
  wi := wa * wb
  bi := ba * bb
'
G.in.L <- '
level: within

group: foo
  galo ~ focc
group: bar
  galo ~ focc

level: between

group: foo
  galo ~ focc
group: bar
  galo ~ focc
'
L.in.G <- ' group: foo
level: within
  galo ~ focc
level: between
  galo ~ focc

group: bar

level: within
  galo ~ focc
level: between
  galo ~ focc
'

library(lavaan)
fit <- sem(model, data = Galo, cluster = "school", fixed.x = FALSE,
           std.lv = TRUE, h1 = TRUE)
fitmg <- sem(L.in.G, data = Galo, cluster = "sch", fixed.x = FALSE,
             std.lv = TRUE, h1 = TRUE, group = "g")

summary(fit, fit.measures = TRUE, standardized = TRUE)
lavInspect(fit, "ngroups")
lavInspect(fit, "nlevels")
fit@Data@nlevels
lavInspect(fit, "cov.lv")
lavInspect(fit, "theta")

## impute data
library(mice)
m <- 5
mice.out <- mice(Galo, m = m, diagnostics = TRUE)
imputedData <- list()
for (i in 1:m) {
  imputedData[[i]] <- complete(data = mice.out, action = i, include = FALSE)
}

library(semTools)
fit.mi <- sem.mi(model, data = imputedData, cluster = "school", fixed.x = FALSE,
                 std.lv = TRUE)
## runs fine, check methods
fit.mi
summary(fit.mi, ci = TRUE, standardized = TRUE, rsquare = TRUE, fmi = TRUE)
coef(fit.mi)
vcov(fit.mi)
nobs(fit.mi)
fitted(fit.mi)
resid(fit.mi, type = "cor") # works, but resid(fit) generates an error
modindices.mi(fit.mi) # modindices(fit) also generates an error
lavTestScore.mi(fit.mi, test = "D2") # lavTestScore(fit, epc = T) generates an error, expects a "group" column
lavTestScore.mi(fit.mi, test = "rubin") # very different results
lavTestWald.mi(fit.mi, constraints = 'wa == ba ; wb == bb')
lavTestLRT.mi(fit.mi, test = "D2")
anova(fit.mi) # D3 takes an obscenely long time
anova(fit.mi, asymptotic = TRUE) # comparable to anova(fit)
fitMeasures(fit.mi) # fails to fit the baseline.model
fitMeasures(fit.mi, fit.measures = c("rmsea","rmr","aic"))



## -------------------------------------
## check lavaan.mi with sampling weights
## -------------------------------------

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
