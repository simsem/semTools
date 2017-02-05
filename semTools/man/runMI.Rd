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
	The analysis model can be specified using lavaan \code{\link[lavaan]{model.syntax}} or parameter table (as returned by \code{\link[lavaan]{parTable}}).
}
  \item{data}{
A \code{data.frame} with missing observations, or a \code{list} of imputed data sets (if data are imputed already). If \code{runMI} has already been called, then imputed data sets are stored in the \code{@DataList} slot, so \code{data} can also be a \code{lavaan.mi} object from which the same imputed data will be used for additional analyses.
}
 \item{fun}{
	\code{character}. Name of a specific lavaan function used to fit \code{model} to \code{data} (i.e., \code{"lavaan"}, \code{"cfa"}, \code{"sem"}, or \code{"growth"}). Only required for \code{runMI}.
}
  \item{\dots}{
	additional arguments to pass to \code{\link[lavaan]{lavaan}}
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
\examples{
\dontrun{

## impose missing data for example
set.seed(12345)
HSMiss <- HolzingerSwineford1939[ , paste("x", 1:9, sep = "")]
HSMiss$x5 <- ifelse(HSMiss$x1 <= quantile(HSMiss$x1, .3), NA, HSMiss$x5)
HSMiss$x9 <- ifelse(is.na(HSMiss$x5), NA, HSMiss$x9)
HSMiss$school <- HolzingerSwineford1939$school

## specify model
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
HS.amelia <- amelia(HSMiss, m = 20, noms = "school")
imps <- HS.amelia$imputations
out2 <- cfa.mi(HS.model, data = imps, m = 20, seed = 12345,
               miArgs = list(noms = "school"))

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



## Categorical variables
popModel <- "
  f1 =~ 0.6*y1 + 0.6*y2 + 0.6*y3 + 0.6*y4
  y1 ~*~ 1*y1
  y2 ~*~ 1*y2
  y3 ~*~ 1*y3
  y4 ~*~ 1*y4
  f1 ~~ 1*f1
  y1 | 0.5*t1
  y2 | 0.25*t1
  y3 | 0*t1
  y4 | -0.5*t1
"
analyzeModel <- "
  f1 =~ y1 + y2 + y3 + y4
  y1 ~*~ 1*y1
  y2 ~*~ 1*y2
  y3 ~*~ 1*y3
  y4 ~*~ 1*y4
"
set.seed(12345)
dat <- simulateData(popModel, sample.nobs  = 200L)
miss.pat <- matrix(as.logical(rbinom(prod(dim(dat)), 1, 0.2)),
                   nrow(dat), ncol(dat))
dat[miss.pat] <- NA

catout <- cfa.mi(analyzeModel, data = dat, ordered = paste0("y", 1:4),
                 m = 3,	miArgs = list(ords = c("y1", "y2", "y3", "y4")))

summary(catout) # standardized/R-squared not yet available for categorical data
anova(catout, indices = "all") # note the scaled versions of indices, too

}
}
