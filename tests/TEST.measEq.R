### Terrence D. Jorgensen
### Last updated: 7 January 2019
### check measEq syntax-writing engine

# devtools::install_github("simsem/semTools/semTools")
library(semTools)


## TODO:  Arguments to vary for quality control (check statistical equivalence):
##        - ID.fac
##        - ID.cat (what happens when 1, 2, or 3 thresholds are constrained)
##        - parameterization
##        - group.equal and long.equal together and separate
##        - auto
##        - longIndNames empty or not
##        - second-order factors
##        - orthogonal (try longitudinal bifactor model)
##        - single-indicator constructs
##        - check sample.cov when is.null(sample.mean) omits mean-structure
##        - lavaan.mi



## --------------------------
## Example 1: Continuous data
## --------------------------

## example population: 2 groups, 2 occations, 2 constructs (partial overlap)
ex.pop <- '
## Lambda
math1 =~ c(.5, .5)*add1 + c(.6, .6)*subtract1 + c(.7, .7)*multiply1 + c(.8, .8)*longDiv1 + c(.9, .9)*fractions1
math2 =~ c(.7, .7)*multiply2 + c(.7, .7)*longDiv2 + c(.8, .8)*fractions2 + c(.8, .8)*algebra2 + c(.8, .8)*geometry2
verbal1 =~ c(.5, .5)*spell1 + c(.6, .6)*match1 + c(.7, .7)*read1 + c(.8, .8)*write1
verbal2 =~ c(.7, .7)*read2 + c(.7, .7)*write2 + c(.7, .7)*synonyms2 + c(.7, .7)*analogy2
success =~ c(.7, .7)*self + c(.7, .7)*peer + c(.7, .7)*teacher + c(.7, .7)*parent

## Psi
math1   ~~ c(  1, 1  )*math1   + c(.5, .5)*math2   + c(.2, .3)*verbal1 + c(.3, .2)*verbal2 + c(.3, .3)*success
math2   ~~ c(1.5, 2  )*math2   + c(.2, .3)*verbal1 + c(.3, .2)*verbal2 + c(.5, .5)*success
verbal1 ~~ c(  1, 1  )*verbal1 + c(.4, .5)*verbal2 + c(.3, .3)*success
verbal2 ~~ c(  2, 1.5)*verbal2 + c(.5, .5)*success
success ~~ c(  1, 1  )*success

## (off-diagonal) Theta
# multiply1  ~~ c(.1, .1)*multiply2
# longDiv1   ~~ c(.1, .1)*longDiv2
# fractions1 ~~ c(.1, .1)*fractions2
# read1      ~~ c(.1, .1)*read2
# write1     ~~ c(.1, .1)*write2

## nu (all zero) and alpha
math1   ~ c(.0, -.1)*1
math2   ~ c(.3,  .2)*1
verbal1 ~ c(.0,  .1)*1
verbal2 ~ c(.2,  .5)*1
success ~ c( 3,   3)*1
'
ex.dat <- simulateData(ex.pop, model.type = "cfa", ov.var = 1, seed = 12345,
                       sample.nobs = c(male = 200L, female = 200L),
                       group.label = c("male","female"))

## generic configural model, all that users should have to specify
ex.mod <- '
math1 =~ add1 + subtract1 + multiply1 + longDiv1 + fractions1
math2 =~ multiply2 + longDiv2 + fractions2 + algebra2 + geometry2
verbal1 =~ spell1 + match1 + read1 + write1
verbal2 =~ read2 + write2 + synonyms2 + analogy2
success =~ self + peer + teacher #+ parent
Parent     =~ parent             # try a single-indicator construct
'
ex.fit <- cfa(ex.mod, data = ex.dat, std.lv = TRUE, meanstructure = TRUE, group = "group")
# summary(ex.fit, fit = TRUE, std = TRUE)

## parse group.partial
group <- "group"
group.equal <- c("loadings","intercepts","residuals")
group.partial <- "math1 =~ longDiv1 ; longDiv1 ~ 1 ; longDiv1 ~~ longDiv1
math2 =~ longDiv2 ; longDiv2 ~ 1 ; longDiv2 ~~ longDiv2 "
## parse long.partial
long.equal <- c("loadings","intercepts","residuals")
long.partial <- "math =~ div ; div ~ 1 ; div ~~ div"
longFacNames <- list(math = c("math1","math2"),
                     verbal = c("verbal1","verbal2"))
longIndNames <- list(mult = c("multiply1","multiply2"),
                     div = c("longDiv1","longDiv2"),
                     frac = c("fractions1","fractions2"),
                     read = c("read1","read2"),
                     write = c("write1","write2"))

# lavTemplate <- cfa(ex.mod, data = ex.dat, std.lv = TRUE, meanstructure = TRUE,
#                    group = "group", do.fit = FALSE)

## names of ordinal indicators, number of categories
# allOrdNames <- lavaan::lavNames(lavTemplate, type = "ov.ord")

# nG <- lavInspect(lavTemplate, "ngroups")
# parameterization <- lavInspect(lavTemplate, "options")$parameterization
# orthogonal <- lavInspect(lavTemplate, "options")$orthogonal

syntax.ex <- measEq.syntax(configural.model = ex.mod, data = ex.dat,
                           ID.fac = "fx", group = group,
                           group.equal = group.equal, group.partial = group.partial,
                           longFacNames = longFacNames, longIndNames = longIndNames,
                           long.equal = long.equal, long.partial = long.partial)
cat(as.character(syntax.ex))
summary(syntax.ex)



## ---------------------------
## Example 2: Categorical data
## ---------------------------

## another example, using categorical data (pretend f1 and f2 are longitudinal)
cat.mod <- ' FU1 =~ u1 + u2 + u3 + u4
             FU2 =~ u5 + u6 + u7 + u8 '
cat.fit <- cfa(cat.mod, data = semTools::datCat, std.lv = TRUE, group = "g",
               parameterization = "theta")
# summary(cat.fit, fit = TRUE, std = TRUE)

## parse group.partial
group <- "g"
group.equal <- c("thresholds","loadings","intercepts","residuals")
group.partial <- "FU1 =~ u3 ; u3 ~ 1 ; u3 ~~ u3 ; u3 | t1
                  FU2 =~ u7 ; u7 ~ 1 ; u7 ~~ u7 ; u7 | t1"
## parse long.partial
long.equal <- c("thresholds","loadings","intercepts","residuals")
long.partial <- "FU =~ ._FU_.ind.3           ;  ._FU_.ind.3 ~ 1 ;
                 ._FU_.ind.3 ~~ ._FU_.ind.3  ;  ._FU_.ind.1 | t2" # automatic names
longFacNames <- list(FU = c("FU1","FU2"))
# long.partial <- "FU =~ v3 ; v3 ~ 1 ; v3 ~~ v3 ; v1 | t2" # manual names
# longIndNames <- lapply(1:4, function(i) paste0("u", c(i, i + 4)))
# names(longIndNames) <- paste0("v", 1:4)

# lavTemplate <- lavaan::cfa(syntax, data = datCat, std.lv = TRUE, do.fit = FALSE,
#                         group = "g", parameterization = "theta")

syntax.cat <- measEq.syntax(configural.model = cat.mod, parameterization = "theta",
                            data = semTools::datCat, # without this, assumes continuous, even if ordered=paste0("u", 1:8)
                            ID.fac = "std.lv", ID.cat = "Wu.Estabrook.2016",
                            ID.thr = c(1L, 2L), # only for ID.cat == "Millsap.Tein.2004"
                            group = group, group.equal = group.equal, group.partial = group.partial,
                            longFacNames = longFacNames, #longIndNames = longIndNames,
                            long.equal = long.equal, long.partial = long.partial)
cat(as.character(syntax.cat))
summary(syntax.cat)


long.equal <- group.equal <- "thresholds"
## constrain 0, 1, or 2 thresholds for first indicator, check ID.cat constraints
group.partial <- "         u1 | t1 + t2 + t3 #+ t4
                           u5 | t1 + t2 + t3 #+ t4 "
long.partial <- " ._FU_.ind.1 | t1 + t2 + t3 #+ t4 "

syntax.cat <- measEq.syntax(configural.model = cat.mod, parameterization = "theta",
                            data = semTools::datCat, # without this, assumes continuous, even if ordered=paste0("u", 1:8)
                            ID.fac = "std.lv", ID.cat = "wu",
                            group = group, group.equal = group.equal, group.partial = group.partial,
                            longFacNames = longFacNames, #longIndNames = longIndNames,
                            long.equal = long.equal, long.partial = long.partial)
cat(as.character(syntax.cat))
summary(syntax.cat)



## -----------------------
## Example 2b: Binary data
## -----------------------

## Unique issues with choosing defaults when following Wu & Estabrook (2016)
##  - more models are equivalent to configural
##  - 1. constrain threshold, why free intercept instead of variance?
##      - either way, equivalent to configural
##  - 2. constrain c("thresholds","loadings","intercepts"), less restricted
##       model than (1.) because latent means and residual variances BOTH free

myData <- read.table("http://www.statmodel.com/usersguide/chap5/ex5.16.dat")
names(myData) <- c("u1","u2","u3","u4","u5","u6","x1","x2","x3","g")
## can pretend longitudinal, use longFacNames from Ex. 2

bin.mod <- '
FU1 =~ u1 + u2 + u3
FU2 =~ u4 + u5 + u6
'
bin.fit <- cfa(bin.mod, data = myData, ordered = paste0("u", 1:6),
               std.lv = TRUE, group = "g", parameterization = "theta")

mod.config <- measEq.syntax(configural.model = bin.fit, group = "g")
mod.thresh <- measEq.syntax(configural.model = bin.fit,group = "g",
                            group.equal = c("thresholds"))
mod.metric <- measEq.syntax(configural.model = bin.fit, group = "g",
                            group.equal = c("thresholds","loadings"))
mod.scalar <- measEq.syntax(configural.model = bin.fit, group = "g",
                            group.equal = c("thresholds","loadings","intercepts"))
mod.strict <- measEq.syntax(configural.model = bin.fit, group = "g",
                            group.equal = c("thresholds","loadings","intercepts","residuals"))

cat(as.character(mod.config))
cat(as.character(mod.thresh))
cat(as.character(mod.metric))
cat(as.character(mod.scalar))
cat(as.character(mod.strict))


## better order
test.seq <- list(strong = c("thresholds","loadings","intercepts"),
                 means = "means", strict = "residuals", homo = "lv.variances")
meq.list <- list()
for (i in 0:length(test.seq)) {
  if (i == 0L) {
    meq.label <- "configural"
    group.equal <- ""
    long.equal <- ""
  } else {
    meq.label <- names(test.seq)[i]
    group.equal <- unlist(test.seq[1:i])
    # long.equal <- unlist(test.seq[1:i])
  }
  meq.list[[meq.label]] <- measEq.syntax(configural.model = bin.mod,
                                         data = myData,
                                         ordered = paste0("u", 1:6),
                                         parameterization = "theta",
                                         ID.fac = "std.lv",
                                         ID.cat = "Wu.Estabrook.2016",
                                         group = "g",
                                         group.equal = group.equal,
                                         #longFacNames = longFacNames,
                                         #long.equal = long.equal,
                                         return.fit = TRUE)
}

compareFit(meq.list)
#TODO: add this example to the help page




## ------------------------------
## Example 3: Higher-Order Factor
## ------------------------------

HS.higher <- '
visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9
intel   =~ visual + textual + speed    # + ageyr   # mix of latent & manifest
'
## try with sample moments instead of data
sample.cov <- lapply(unique(HolzingerSwineford1939$school), function(g) {
  cov(HolzingerSwineford1939[HolzingerSwineford1939$school == g,
                             paste0("x", 1:9)])
})
sample.mean <- lapply(unique(HolzingerSwineford1939$school), function(g) {
  colMeans(HolzingerSwineford1939[HolzingerSwineford1939$school == g,
                                  paste0("x", 1:9)])
})
sample.nobs <- table(HolzingerSwineford1939$school)

## try missing data with auxiliaries
# dat1 <- lavaan::HolzingerSwineford1939
# set.seed(12345)
# dat1$z <- rnorm(nrow(dat1))
# dat1$x5 <- ifelse(dat1$z < quantile(dat1$z, .3), NA, dat1$x5)
# dat1$x9 <- ifelse(dat1$z > quantile(dat1$z, .8), NA, dat1$x9)
# fit.aux <- cfa.auxiliary(HS.higher, data = dat1, group = "school",
#                          aux = c("z","ageyr","grade"), std.lv = TRUE)

## try multiple imputations
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
## fit multigroup model without invariance constraints
mgfit1 <- cfa.mi(HS.higher, data = imps, std.lv = TRUE, group = "school")



## parse group.partial
group <- "school"
group.equal <- c("loadings",
                 #"regressions",
                 "intercepts", "means")
group.partial <- " intel ~ 1 "
##TODO: Notes for help page:
##  - "regressions" only constraints beta, so a higher-order factor with both
##    manifest and latent indicators will be constrained in separate steps

syntax.hi <- measEq.syntax(configural.model = mgfit1, #HS.higher,
                           ID.fac = "std.lv", # resets to "ul"
                           #data = HolzingerSwineford1939,
                           #data = dat1, # intercepts & (co)variances get constrained! use aux() after
                           #sample.cov = sample.cov, sample.nobs = sample.nobs,
                           #sample.mean = sample.mean,
                           group = group, group.equal = group.equal,
                           group.partial = group.partial)
cat(as.character(syntax.hi))
summary(syntax.hi)



## -------------------------
## Example 4: Bifactor Model
## -------------------------

## add noise to indicators twice, to simulate longitudinal measurement
HS.orig <- HolzingerSwineford1939[paste0("x", 1:9)]
set.seed(123)
HS.1 <- lapply(HS.orig, function(x) x + rnorm(length(x), sd = .5))
HS.2 <- lapply(HS.orig, function(x) x + rnorm(length(x), mean = 1, sd = .5))
names(HS.2) <- paste0("y", 1:9)
HS.bi <- data.frame(HS.1, HS.2)

HS.bifac <- '
## Time 1
visual1  =~ x1 + x2 + x3
textual1 =~ x4 + x5 + x6
speed1   =~ x7 + x8 + x9
intel1   =~ x2 + x1 + x3 + x4 + x5 + x6 + x7 + x8 + x9

## Time 2
visual2  =~ y1 + y2 + y3
textual2 =~ y4 + y5 + y6
speed2   =~ y7 + y8 + y9
intel2   =~ y2 + y1 + y3 + y4 + y5 + y6 + y7 + y8 + y9
'
long.equal <- c("loadings","intercepts")
##TODO: Notes for help page:
##  - without a unique indicator per construct, constraining intercepts does
##    not free latent means, even if they can be identified.
##  - Must free manually if (ID.fac == "uv") OR set ID.fac == "ul" (not "fx")
longFacNames <- list(visual = c("visual1","visual2"),
                     textual = c("textual1","textual2"),
                     speed = c("speed1","speed2"),
                     intel = c("intel1","intel2"))
longIndNames <- lapply(1:9, function(i) paste0(c("x","y"), i))
names(longIndNames) <- paste0("v", 1:9)


## set orthogonal (bifactor model), but autocovariances should still be free
syntax.bi <- measEq.syntax(configural.model = HS.bifac, data = HS.bi,
                           ID.fac = "ul", long.equal = long.equal,
                           longIndNames = longIndNames,
                           longFacNames = longFacNames, orthogonal = TRUE)
cat(as.character(syntax.bi))
summary(syntax.bi)



## ------------------------------
## Example 5: Multiple Imputation
## ------------------------------


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

mgfit1 <- cfa.mi(HS.model, data = imps, std.lv = TRUE, group = "school")

syntax.mi <- measEq.syntax(configural.model = mgfit1, group = "school",
                           group.equal = c("loadings","intercepts"))
cat(as.character(syntax.mi))
summary(syntax.mi)

## return fit
fit.mi <- measEq.syntax(configural.model = mgfit1, group = "school",
                        group.equal = c("loadings","intercepts"),
                        return.fit = TRUE)
summary(fit.mi)





