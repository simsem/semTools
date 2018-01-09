### Terrence D. Jorgensen
### Last updated: 6 February 2017
### compare semTools::twostage() to lavaan(missing = "two.stage")

# devtools::install_github("simsem/semTools/semTools")
library(semTools)

# impose missing data for example
HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
                                      "ageyr","agemo","school")]
set.seed(12345)
HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
age <- HSMiss$ageyr + HSMiss$agemo/12
HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)

## specify CFA model from lavaan's ?cfa help page
HS.model <- '
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
'


# relevant default options are:\
# information = "observed" (versus "expected")
# h1.information = "structured" (versus "unstructured")
# observed.information = "h1" (versus "hessian") ############## lavOptions() always sets to "h1"



## semTools
out <- cfa.2stage(model = model, data = dat1, observed.information = "hessian")
## robust
out.r <- cfa.2stage(model = model, data = dat1, se = "robust.huber.white")

## lavaan
fit <- cfa(model = model, data = dat1, missing = "two.stage")
## robust
fit.r <- cfa(model = model, data = dat1, missing = "robust.two.stage")

## check defaults
lavInspect(out@target, "options")[c("information","h1.information","observed.information")]
lavInspect(out.r@target, "options")[c("information","h1.information","observed.information")]
lavInspect(fit, "options")[c("information","h1.information","observed.information")]
lavInspect(fit.r, "options")[c("information","h1.information","observed.information")]


## compare standard errors
cbind(mine = coef(out), my.r = coef(out.r), Yves = coef(fit), Yv.r = coef(fit.r))
cbind(mine = summary(out)$se, my.r = summary(out.r)$se,
      Yves = parameterEstimates(fit)$se, Yv.r = parameterEstimates(fit.r)$se)

## compare fit statistics
anova(out)
anova(out.r)
fitMeasures(fit, c("chisq","chisq.scaled"))
fitMeasures(fit.r, c("chisq","chisq.scaled"))




## FIXME: add these options to above
un1 <- cfa(model = HS.model, data = HSMiss, missing = "two.stage",
           h1.information = "unstructured", observed.information = "h1")
unHes <- cfa(model = HS.model, data = HSMiss, missing = "two.stage",
             h1.information = "unstructured", observed.information = "hessian")
str1 <- cfa(model = HS.model, data = HSMiss, missing = "two.stage",
            h1.information = "structured", observed.information = "h1")
strHes <- cfa(model = HS.model, data = HSMiss, missing = "two.stage",
              h1.information = "structured", observed.information = "hessian")
## check coefs
data.frame(Terry = coef(out.t),
           un1 = coef(un1),
           unHes = coef(unHes),
           str1 = coef(un1),
           strHes = coef(strHes))
## check SEs
data.frame(Terry = sqrt(diag(vcov(out.t))),
           un1 = sqrt(diag(vcov(un1))),
           unHes = sqrt(diag(vcov(unHes))),
           str1 = sqrt(diag(vcov(un1))),
           strHes = sqrt(diag(vcov(strHes))))





# 1) names.gw  (scalar)
# 2) names.th  (thresholds and numeric intercepts are interweaved)
# 3) names.mu  (vector)
# 4) names.pi  (as vec(PI))
# 5) names.cov (as vech(cov, diag = TRUE)
# 6) names.var (numeric variables only)
# 7) names.cor (as vech(cor, diag = FALSE)
#
# If data is continuous, names.th, names.var and names.cor will be empty;
# If data is categorical, names.mu and names.cov will be empty;
# If conditional.x = TRUE, and there are some exogenous covariates,
# names.pi will be non-empty, and cov/cor matrices only contain non-exogenous variables

