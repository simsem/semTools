\name{fmi}
\alias{fmi}
\title{
Fraction of Missing Information.
}
\description{
	This function estimates the Fraction of Missing Information (FMI) for summary statistics of each variable, using either an incomplete data set or a list of imputed data sets.
}
\usage{
fmi(data, method = "saturated", group = NULL, ords = NULL,
    varnames = NULL, exclude = NULL, fewImps = TRUE)
}
\arguments{
  \item{data}{
	Either a single \code{data.frame} with incomplete observations, or a \code{list} of imputed data sets.
}
  \item{method}{
	character. If \code{"saturated"} or \code{"sat"} (default), the model used to estimate FMI is a freely estimated covariance matrix and mean vector for numeric variables, and/or polychoric correlations and thresholds for ordered categorical variables, for each group (if applicable). If \code{"null"}, only means and variances are estimated for numeric variables, and/or thresholds for ordered categorical variables (i.e., covariances and/or polychoric correlations are constrained to zero). See Details for more information.
}
 \item{group}{
  character. The optional name of a grouping variable, to request FMI in each group.
}
 \item{ords}{
  character. Optional vector of names of ordered-categorical variables, which are not already stored as class \code{ordered} in \code{data}.
}
 \item{varnames}{
	character. Optional vector of variable names, to calculate FMI for a subset of variables in \code{data}. By default, all numeric and ordered variables will be included, unless \code{data} is a single incomplete \code{data.frame}, in which case only numeric variables can be used with FIML estimation. Other variable types will be removed.
}
 \item{exclude}{
  character. Optional vector of variable names to exclude from the analysis.
}
 \item{fewImps}{
	logical. If \code{TRUE} (default), use the estimate of FMI that applies a correction to the estimated between-imputation variance. Recommended when there are few imputations; makes little difference when there are many imputations.
}
}
\details{
The function estimates a saturated model with \code{\link[lavaan]{lavaan}} for a single incomplete data set using FIML, or with \code{\link{lavaan.mi}} for a list of imputed data sets. If method = \code{"saturated"}, FMI will be estiamted for all summary statistics, which could take a lot of time with big data sets. If method = \code{"null"}, FMI will only be estimated for univariate statistics (e.g., means, variances, thresholds). The saturated model gives more reliable estimates, so it could also help to request a subset of variables from a large data set.
}
\value{
\code{fmi} returns a list with at least 2 of the following:
\item{Covariances}{A list of symmetric matrices: (1) the estimated/pooled covariance matrix, or a list of group-specific matrices (if applicable) and (2) a matrix of FMI, or a list of group-specific matrices (if applicable). Only available if \code{method = "saturated"}.}
\item{Variances}{The estimated/pooled variance for each numeric variable. Only available if \code{method = "null"} (otherwise, it is on the diagonal of Covariances).}
\item{Means}{The estimated/pooled mean for each numeric variable.}
\item{Thresholds}{The estimated/pooled threshold(s) for each ordered-categorical variable.}
\item{message}{A message indicating caution when the null model is used.}
}
\references{
Rubin, D.B. (1987) \emph{Multiple imputation for nonresponse in surveys.} New York, NY: Wiley.

Savalei, V. & Rhemtulla, M. (2012) On obtaining estimates of the fraction
of missing information from full information maximum likelihood, \emph{Structural Equation Modeling, 19}(3), 477-494. doi:10.1080/10705511.2012.687669

Wagner, J. (2010) The fraction of missing information as a tool for monitoring the quality of survey data, \emph{Public Opinion Quarterly, 74}(2), 223-243. doi:10.1093/poq/nfq007
}
\author{
  Mauricio Garnier Villarreal (University of Kansas; \email{mauricio.garniervillarreal@marquette.edu})
  Terrence Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\examples{
library(Amelia)
library(lavaan)

## specify population and simulate a data set
modsim <- '
f1 =~ 0.7*y1+0.7*y2+0.7*y3
f2 =~ 0.7*y4+0.7*y5+0.7*y6
f3 =~ 0.7*y7+0.7*y8+0.7*y9
'
datsim <- simulateData(modsim,model.type="cfa", meanstructure=TRUE,
                       std.lv=TRUE, sample.nobs=c(200,200))

## impose missing values
randomMiss2 <- rbinom(prod(dim(datsim)), 1, 0.1)
randomMiss2 <- matrix(as.logical(randomMiss2), nrow=nrow(datsim))
randomMiss2[,10] <- FALSE
datsim[randomMiss2] <- NA

## calculate FMI (using FIML, provide partially observed data set)
(out1 <- fmi(datsim, exclude = "group"))
(out2 <- fmi(datsim, exclude = "group", method = "null"))
(out3 <- fmi(datsim, varnames = c("y1","y2","y3","y4")))
(out4 <- fmi(datsim, group = "group"))

## ordered-categorical data
data(datCat)
lapply(datCat, class)
## impose missing values
for (i in 1:8) datCat[sample(1:nrow(datCat), size = .2*nrow(datCat)), i] <- NA
## impute data m = 3 times
impout <- amelia(datCat, m = 3, noms = "g", ords = paste0("u", 1:8))
imps <- impout$imputations
## calculate FMI, using list of imputed data sets
fmi(imps, group = "g")
}
