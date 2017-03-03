\name{runMI}
\alias{runMI}
\alias{lavaan.mi}
\alias{cfa.mi}
\alias{sem.mi}
\alias{growth.mi}
\title{
Fit a lavaan Model to Multiple Imputed Data Sets
}
\description{
This function fits a lavaan model to a list of imputed data sets, and can also implement multiple imputation for a single \code{data.frame} with missing observations, using either the Amelia package or the mice package.
}
\usage{
runMI(model, data, fun = "lavaan", ...,
      m, miArgs = list(), miPackage = "Amelia", seed = 12345)
lavaan.mi(model, data, ...,
          m, miArgs = list(), miPackage = "Amelia", seed = 12345)
cfa.mi(model, data, ...,
       m, miArgs = list(), miPackage = "Amelia", seed = 12345)
sem.mi(model, data, ...,
       m, miArgs = list(), miPackage = "Amelia", seed = 12345)
growth.mi(model, data, ...,
          m, miArgs = list(), miPackage = "Amelia", seed = 12345)
}
\arguments{
  \item{model}{
	The analysis model can be specified using lavaan \code{\link[lavaan]{model.syntax}} or a parameter table (as returned by \code{\link[lavaan]{parTable}}).
}
  \item{data}{
A \code{data.frame} with missing observations, or a \code{list} of imputed data sets (if data are imputed already). If \code{runMI} has already been called, then imputed data sets are stored in the \code{@DataList} slot, so \code{data} can also be a \code{lavaan.mi} object from which the same imputed data will be used for additional analyses.
}
 \item{fun}{
	\code{character}. Name of a specific lavaan function used to fit \code{model} to \code{data} (i.e., \code{"lavaan"}, \code{"cfa"}, \code{"sem"}, or \code{"growth"}). Only required for \code{runMI}.
}
  \item{\dots}{
	additional arguments to pass to \code{\link[lavaan]{lavaan}} or \code{\link[lavaan]{lavaanList}}. See also \code{\link[lavaan]{lavOptions}}. Note that \code{lavaanList} provides parallel computing options, as well as a \code{FUN} argument so the user can extract custom output after the model is fitted to each imputed data set (see \strong{Examples}).  TIP: If a custom \code{FUN} is used \emph{and} \code{parallel = "snow"} is requested, the user-supplied function should explicitly call \code{library} or use \code{\link[base]{"::"}} for any functions not part of the base distribution.
}
  \item{m}{
\code{integer}. Request the number of imputations. Ignored if \code{data} is already a \code{list} of imputed data sets or a \code{lavaan.mi} object.
}
  \item{miArgs}{
Addition arguments for the multiple-imputation function (\code{miPackage}). The arguments should be put in a list (see example below). Ignored if \code{data} is already a \code{list} of imputed data sets or a \code{lavaan.mi} object.
}
  \item{miPackage}{
Package to be used for imputation. Currently these functions only support \code{"Amelia"} or \code{"mice"} for imputation. Ignored if \code{data} is already a \code{list} of imputed data sets or a \code{lavaan.mi} object.
}
  \item{seed}{
	Random number seed to be set before imputing the data. Ignored if \code{data} is already a \code{list} of imputed data sets or a \code{lavaan.mi} object.
}
}
\value{
	A \code{\linkS4class{lavaan.mi}} object
}
\references{
Enders, C. K. (2010). \emph{Applied missing data analysis}. New York, NY: Guilford.

Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}. New York, NY: Wiley.
}
\author{
  Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\examples{ \dontrun{
## impose missing data for example
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

## impute data within runMI...
out1 <- cfa.mi(HS.model, data = HSMiss, m = 20, seed = 12345,
               miArgs = list(noms = "school"))

## ... or impute missing data first
library(Amelia)
set.seed(12345)
HS.amelia <- amelia(HSMiss, m = 20, noms = "school", p2s = FALSE)
imps <- HS.amelia$imputations
out2 <- cfa.mi(HS.model, data = imps)

## same results (using the same seed results in the same imputations)
cbind(impute.within = coef(out1), impute.first = coef(out2))

summary(out1)
summary(out1, ci = FALSE, fmi = TRUE, add.attributes = FALSE)
summary(out1, ci = FALSE, stand = TRUE, rsq = TRUE)

## model fit. D3 includes information criteria
anova(out1)
anova(out1, test = "D2", indices = TRUE) # request D2 and fit indices



## fit multigroup model without invariance constraints
mgfit1 <- cfa.mi(HS.model, data = imps, estimator = "mlm", group = "school")
## add invariance constraints, and use previous fit as "data"
mgfit0 <- cfa.mi(HS.model, data = mgfit1, estimator = "mlm", group = "school",
                 group.equal = c("loadings","intercepts"))

## compare fit (scaled likelihood ratio test)
anova(mgfit0, h1 = mgfit1)

## correlation residuals
resid(mgfit0, type = "cor.bentler")


## use D1 to test a parametrically nested model (whether latent means are ==)
anova(mgfit0, test = "D1", constraints = '
      .p70. == 0
      .p71. == 0
      .p72. == 0')



## ordered-categorical data
data(datCat)
lapply(datCat, class)
## impose missing values
set.seed(123)
for (i in 1:8) datCat[sample(1:nrow(datCat), size = .1*nrow(datCat)), i] <- NA

catout <- cfa.mi(' f =~ u1 + u2 + u3 + u4 ', data = datCat,
                 m = 3, seed = 456,
                 miArgs = list(ords = paste0("u", 1:8), noms = "g"),
                 FUN = function(fit) {
                   list(wrmr = lavaan::fitMeasures(fit, "wrmr"),
                        zeroCells = lavaan::lavInspect(fit, "zero.cell.tables"))
                 })
summary(catout)
anova(catout, indices = "all") # note the scaled versions of indices, too

## extract custom output
sapply(catout@funList, function(x) x$wrmr) # WRMR for each imputation
catout@funList[[1]]$zeroCells # zero-cell tables for first imputation
catout@funList[[2]]$zeroCells # zero-cell tables for second imputation ...

} }
