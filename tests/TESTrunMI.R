### Terrence D. Jorgensen
### Last updated: 18 April 2018
### try new runMI

library(lavaan)
library(Amelia)

# source("../FuturePlans/NEWrunMI.R")
library(semTools)


##############
## Examples ##
##############

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

anova(fit1)                   # pooled LRT by default (D3 statistic)
anova(fit1, test = "D2")      # or pooled Wald test (D2 statistic)
anova(fit1, asymptotic = TRUE, indices = TRUE)  # as chisq == F * df1, add indices
anova(fit0, h1 = fit1, asymptotic = TRUE)  # compare fit
anova(fit0, h1 = fit1, test = "D2")  # compare fit



## fit multigroup model
mgfit1 <- cfa.mi(HS.model, data = imps, group = "school", estimator = "mlmv")
mgfit0 <- cfa.mi(HS.model, data = imps, group = "school", estimator = "mlmv",
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

anova(mgfit0)          # pooled LRT by default (D3 statistic)
anova(mgfit0, asymptotic = TRUE, indices = TRUE)  # reported by Mplus
anova(mgfit0, test = 'D2', indices = "all")        # add fit indices
anova(mgfit0, h1 = mgfit1)           # compare fit
anova(mgfit0, h1 = mgfit1, test = 'D2') # robustifies pooled naive statistic
anova(mgfit0, h1 = mgfit1, test = 'D2', pool.robust = TRUE) # pools robust statistic

## use D1 to test a parametrically nested model (whether latent means are ==)
anova(mgfit0, test = "D1", constraints = '
      .p70. == 0
      .p71. == 0
      .p72. == 0')

## Score test for metric invariance
lavTestScore.mi(mgfit1, add = 'x1 ~~ x9 ; x4 ~~ x7',
                type = "Rubin", scale.W = TRUE)
lavTestScore.mi(mgfit1, add = 'x1 ~~ x9 ; x4 ~~ x7',
                type = "Rubin", scale.W = FALSE)
lavTestScore.mi(mgfit0, release = 1:6, type = "Rubin", scale.W = TRUE)
lavTestScore.mi(mgfit0, release = 1:6, type = "Rubin", scale.W = FALSE)


##############################################
## Compare to old runMI output and to Mplus ##
##############################################

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
## Investigate probelm with CFI when test is MV-adjusted
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
