### Terrence D. Jorgensen
### Last updated: 14 May 2018
### test categorical invariance in various ways:
###  - std.lv = TRUE or FALSE
###  - various parameterizations: "delta", "theta", or Millsap & Tein (2004)
###  - method = "satorra.2000" or "satorra.bentler.2001" -- why?
###  - A.method = "exact" or "delta"
###  -


library(semTools) # for example data


## ----------
## fixed to 1
## ----------

model <- ' f1 =~ u1 + u2 + u3 + u4'
fit.config <- cfa(model, data = datCat, group = "g",
                  parameterization = "theta", ordered = c("u1","u2","u3","u4"))
fit.metric <- cfa(model, data = datCat, group = "g", group.equal = c("loadings"),
                  parameterization = "theta", ordered = c("u1","u2","u3","u4"))
lavTestLRT(fit.metric, fit.config) # no problem
fit.scalar <- cfa(model, data = datCat, group = "g",
                  group.equal = c("loadings","thresholds"),
                  parameterization = "theta", ordered = c("u1","u2","u3","u4"))
lavTestLRT(fit.metric, fit.scalar) # no problem



## ----------------
## constrained to 1
## ----------------

model <- ' f1 =~ u1 + u2 + u3 + u4
u1 ~~ c(theta1.1, theta1.2)*u1
u2 ~~ c(theta2.1, theta2.2)*u2
u3 ~~ c(theta3.1, theta3.2)*u3
u4 ~~ c(theta4.1, theta4.2)*u4
f1 ~  c(alpha1.1, alpha1.2)*1

theta1.1 == 1
theta2.1 == 1
theta3.1 == 1
theta4.1 == 1
alpha1.1 == 0
'
con.config <- '
theta1.2 == 1
theta2.2 == 1
theta3.2 == 1
theta4.2 == 1
alpha1.2 == 0
'
fit.config <- cfa(model, data = datCat, constraints = con.config,
                  group = "g",
                  parameterization = "theta", ordered = c("u1","u2","u3","u4"))
fit.metric <- cfa(model, data = datCat, constraints = con.config,
                  group = "g", group.equal = c("loadings"),
                  parameterization = "theta", ordered = c("u1","u2","u3","u4"))
lavTestLRT(fit.metric, fit.config) # no problem
fit.scalar <- cfa(model, data = datCat,
                  group = "g", group.equal = c("loadings","thresholds"),
                  parameterization = "theta", ordered = c("u1","u2","u3","u4"))
lavTestLRT(fit.config, fit.metric, fit.scalar) # no problem



## -------------------------------------
## compare to measurementInvarianceCat()
## -------------------------------------

syntax <- ' f1 =~ u1 + u2 + u3 + u4'

measurementInvarianceCat(model = syntax, data = datCat, group = "g",
                         parameterization = "theta", estimator = "wlsmv",
                         ordered = c("u1", "u2", "u3", "u4"))





