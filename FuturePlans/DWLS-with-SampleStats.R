### Terrence D. Jorgensen
### Last updated: 1 November 2022
### sandbox for manual (D)WLS estimation with sample stats
### Eventual uses:
### - Fit SRM in lavaan
### - pool imputed polychorics before fitting model in "single" step:
###     https://doi.org/10.1080/00273171.2018.1523000
### - MASEM: https://groups.google.com/d/msg/lavaan/aLZJ1Fahrnk/Bv2BqE1uAQAJ


library(lavaan)

## make some ordinal data
HS9 <- HolzingerSwineford1939[,c("x1","x2","x3","x4","x5",
                                 "x6","x7","x8","x9")]
Data <- as.data.frame( lapply(HS9, cut, 3, labels=FALSE) )
N <- nrow(Data)
Data$school <- HolzingerSwineford1939$school



## ---------------------------
## Estimate Summary Statistics
## ---------------------------

## estimate polychorics (DWLS and PML for comparison)
cor.wls <- lavCor(Data[paste0("x", 1:6)], ordered = paste0("x", 1:6),
                  output = "fit",
                  se = "robust.sem", h1 = TRUE, estimator = "DWLS")
cor.pml <- lavCor(Data[paste0("x", 1:6)], ordered = paste0("x", 1:6),
                  output = "fit",
                  se = "robust.huber.white", h1 = TRUE, estimator = "PML")

## estimate polyserials (only treat 2:4 as ordered)
cov.wls <- lavCor(Data[paste0("x", 1:6)], ordered = paste0("x", 2:4),
                  output = "fit",
                  se = "robust.sem", h1 = TRUE, estimator = "DWLS")
cov.pml <- lavCor(Data[paste0("x", 1:6)], ordered = paste0("x", 2:4),
                  output = "fit",
                  se = "robust.huber.white", h1 = TRUE, estimator = "PML")



## ----------------------------------
## Fit a one-factor model to raw data
## ----------------------------------

model <- ' trait =~ x1 + x2 + x3 + x4 + x5 + x6 '
fit3.wls <- sem(model, data = Data, ordered = paste0("x", 2:4), std.lv = TRUE)
fit6.wls <- sem(model, data = Data, ordered = paste0("x", 1:6), std.lv = TRUE)



## -------------------------------------
## Fit same model to summary statistics:
##   ALL 6 variables are ordered
## -------------------------------------

## save data summary from fitted model
## NOTE: Equivalently, use cor.wls (saturated model)
sample.cov  <- lavInspect(fit6.wls, "sampstat")$cov
sample.mean <- lavInspect(fit6.wls, "sampstat")$mean
sample.th   <- lavInspect(fit6.wls, "sampstat")$th
attr(sample.th, "th.idx") <- lavInspect(fit6.wls, "th.idx")
sample.nobs <- lavInspect(fit6.wls, "nobs")
WLS.V <- lavInspect(fit6.wls, "WLS.V")
NACOV <- lavInspect(fit6.wls, "gamma")

# using sample statistics FROM A FITTED lavaan MODEL
fit6.sum1 <- cfa(model, std.lv = TRUE,
                 sample.cov = sample.cov, sample.mean = sample.mean,
                 sample.th = sample.th, sample.nobs = sample.nobs,
                 WLS.V = WLS.V, NACOV = NACOV)
fit6.wls # identical



## Compare information stored in fitted model to saturated model from lavCor()

## vcov(cor.wls) in different order
idx <- c(16:27, 1:15) # if all are polychorics

NACOV.wls <- (N-1)*vcov(cor.wls)[idx, idx]
NACOV.pml <- (N-1)*vcov(cor.pml)[idx, idx]
all.equal(unclass(NACOV), NACOV.wls, check.attributes = FALSE)
all.equal(unclass(NACOV), NACOV.pml, check.attributes = FALSE) # not quite

WLS.V.wls <- solve(diag(diag(NACOV.wls)))
WLS.V.pml <- solve(diag(diag(NACOV.pml)))
all.equal(unclass(WLS.V), WLS.V.wls)
all.equal(unclass(WLS.V), WLS.V.pml) # not quite


# using sample statistics FROM SATURATED MODEL
Thr <- lavInspect(cor.wls, "th")
attr(Thr, "th.idx") <- lavInspect(cor.wls, "th.idx")
fit6.sum2 <- cfa(model, std.lv = TRUE,
                 sample.cov = lavInspect(cor.wls, "cov.ov"),
                 sample.mean = lavInspect(cor.wls, "mean.ov"),
                 sample.th = Thr, sample.nobs = N,
                 WLS.V = WLS.V.wls, NACOV = NACOV.wls)
fit6.wls # identical



## ----------------------------------------
## DWLS on correlations, missing thresholds
## (mimic a meta-analysis)
## ----------------------------------------

NACOV.hack <- matrix(0, nrow = 21, ncol = 21)
NACOV.hack[7:21, 7:21] <- NACOV.wls[13:27, 13:27]
rownames(NACOV.hack)[7:21] <- rownames(NACOV.wls)[13:27]
colnames(NACOV.hack)[7:21] <- colnames(NACOV.wls)[13:27]
diag(NACOV.hack)[1:6] <- 1
WLS.V.hack <- MASS::ginv(diag(diag(NACOV.hack)))
Thr.hack <- rep(0, 6)
attr(Thr.hack, "th.idx") <- 1:6
names(Thr.hack) <- names(attr(Thr.hack, "th.idx")) <- paste0("x", 1:6, "|t1")
fit6.sum3 <- cfa(model, std.lv = TRUE,
                 sample.cov = lavInspect(cor.wls, "cov.ov"),
                 sample.mean = lavInspect(cor.wls, "mean.ov"),
                 sample.th = Thr.hack, sample.nobs = N,
                 WLS.V = WLS.V.hack, NACOV = NACOV.hack)
summary(fit6.sum3)
summary(fit6.wls)



## -------------------------------------------
## Full WLS on correlations, missing variances
## (mimic a meta-analysis)
## -------------------------------------------

## WLS expects order to be vech(cov, diag = TRUE), not c(var, vech(cor, diag = FALSE))
idx <- c(1, 7:11, 2, 12:15, 3, 16:18, 4, 19:20, 5, 21, 6)

NACOV.hack2 <- NACOV.hack[idx, idx]
WLS.V.hack2 <- MASS::ginv(NACOV.hack2)
fit6.sum4 <- cfa(model, std.lv = TRUE,
                 sample.cov = lavInspect(cor.wls, "cov.ov"), # * (N-1) / N,
                 sample.nobs = N, estimator = "wls",
                 WLS.V = WLS.V.hack2, NACOV = NACOV.hack2)
summary(fit6.sum4) # similar, but not identical (worse fit)
summary(fit6.wls)



## -------------------------------------
## Fit same model to summary statistics:
##   ONLY 3 variables are ordered
## -------------------------------------

## save data summary from fitted model
sample.cov  <- lavInspect(fit3.wls, "sampstat")$cov
sample.mean <- lavInspect(fit3.wls, "sampstat")$mean
sample.th   <- lavInspect(fit3.wls, "sampstat")$th
attr(sample.th, "th.idx") <- lavInspect(fit3.wls, "th.idx")
sample.nobs <- lavInspect(fit3.wls, "nobs")
WLS.V <- lavInspect(fit3.wls, "WLS.V")
NACOV <- lavInspect(fit3.wls, "gamma")

# using sample statistics FROM A FITTED lavaan MODEL
fit3.sum1 <- cfa(model, std.lv = TRUE,
                 sample.cov = sample.cov, sample.mean = sample.mean,
                 sample.th = sample.th, sample.nobs = sample.nobs,
                 WLS.V = WLS.V, NACOV = NACOV)
fit3.wls # identical


## Equivalently, use cov.wls (saturated-model estimates)
sample.cov  <- lavInspect(cov.wls, "cov.ov")
sample.mean <- lavInspect(cov.wls, "mean.ov")
sample.th   <- lavInspect(cov.wls, "thresholds")
attr(sample.th, "th.idx") <- lavInspect(cov.wls, "th.idx")
sample.nobs <- lavInspect(cov.wls, "nobs")
WLS.V <- lavInspect(cov.wls, "WLS.V")
NACOV <- lavInspect(cov.wls, "gamma")

# using sample statistics FROM A FITTED lavaan MODEL
fit3.sum2 <- cfa(model, std.lv = TRUE,
                 sample.cov = sample.cov, sample.mean = sample.mean,
                 sample.th = sample.th, sample.nobs = sample.nobs,
                 WLS.V = WLS.V, NACOV = NACOV)
fit3.wls # identical




## Compare information stored in fitted model to saturated model from lavCor()

## vcov(cov.wls) in different order: interleaved means + thresholds,
##                                   variances, polychorics/serials
idx <- c(25, 19:24, 26:27, 1:18)
## Another issue: sampling error for means is in the opposite direction
(N-1)*vcov(cov.wls)[c(25, 19:24, 26:27), c(25, 19:24, 26:27)]
NACOV[1:9, 1:9]
## apply a sign-switch for rows/columns of means with other rows/columns
neg1 <- matrix(1, nrow = nrow(NACOV), ncol = ncol(NACOV))
neg1[25:27, -(25:27)] <- -1
neg1[-(25:27), 25:27] <- -1


NACOV.wls <- (N-1)*(neg1*vcov(cov.wls))[idx, idx]
NACOV.pml <- (N-1)*(neg1*vcov(cov.pml))[idx, idx]
all.equal(unclass(NACOV), NACOV.wls, check.attributes = FALSE)
all.equal(unclass(NACOV), NACOV.pml, check.attributes = FALSE) # not quite

WLS.V.wls <- solve(diag(diag(NACOV.wls)))
WLS.V.pml <- solve(diag(diag(NACOV.pml)))
all.equal(unclass(WLS.V), WLS.V.wls)
all.equal(unclass(WLS.V), WLS.V.pml) # not quite


# using sample statistics FROM SATURATED MODEL
Thr <- lavInspect(cov.wls, "th")
attr(Thr, "th.idx") <- lavInspect(cov.wls, "th.idx")
fit3.sum2 <- cfa(model, std.lv = TRUE,
                 sample.cov = lavInspect(cov.wls, "cov.ov"),
                 sample.mean = lavInspect(cov.wls, "mean.ov"),
                 sample.th = Thr, sample.nobs = N,
                 WLS.V = WLS.V.wls, NACOV = NACOV.wls)
fit3.wls # identical



## ------------------------------------------------
## try capturing order from lavInspect(fit, "free")
## ------------------------------------------------

free.idx <- lavInspect(cov.wls, "free")
nn0    <- lavNames(cov.wls)
nn.num <- lavNames(cov.wls, "ov.num")
nn.ord <- lavNames(cov.wls, "ov.ord")

mth.idx <- integer()
for (i in nn0) {
  if (i %in% nn.num) {
    mth.idx <- c(mth.idx, free.idx$nu[i, 1])
  }
  if (i %in% nn.ord) {
    nn.thr <- rownames(free.idx$tau)
    row.thr <- grep(pattern = paste0(i, "|"), nn.thr, fixed = TRUE)
    mth.idx <- c(mth.idx, free.idx$tau[row.thr, 1])
  }
  #TODO: any other possibilities?
}

## or, simpler?
mth.idx <- lavInspect(cov.wls, "th.idx")
num.idx <- setdiff(seq_along(nn0), mth.idx)
mth.idx[mth.idx != 0] <- free.idx$tau[ , 1]
mth.idx[mth.idx == 0] <- free.idx$nu[num.idx, 1]

wls.idx <- c(mth.idx, # interleaved thresholds + (negative) means
             diag(free.idx$theta)[num.idx], # variances
             lav_matrix_vech(free.idx$theta, diagonal = FALSE)) # covariances
all(idx == wls.idx) # same: cbind(idx, wls.idx)


