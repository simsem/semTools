### Terrence D. Jorgensen
### Last updated: 5 May 2022
### package documentation, along with convenience documentation (e.g., imports)


##' semTools: Useful Tools for Structural Equation Modeling
##'
##' The \pkg{semTools} package provides many miscellaneous functions that are
##' useful for statistical analysis involving SEM in R.  Many functions extend
##' the funtionality of the \pkg{lavaan} package.  Some sets of functions in
##' \pkg{semTools} correspond to the same theme. We call such a collection of
##' functions a \emph{suite}. Our suites include:
##' \itemize{
##' \item{Model Fit Evaluation:
##'   \code{\link{moreFitIndices}},
##'   \code{\link{nullRMSEA}},
##'   \code{\link{singleParamTest}},
##'   \code{\link{miPowerFit}}, and
##'   \code{\link{chisqSmallN}}}
##' \item{Measurement Invariance:
##'   \code{\link{measEq.syntax}},
##'   \code{\link{partialInvariance}},
##'   \code{\link{partialInvarianceCat}}, and
##'   \code{\link{permuteMeasEq}}}
##' \item{Power Analysis:
##'   \code{\link{SSpower}},
##'   \code{\link{findRMSEApower}},
##'   \code{\link{plotRMSEApower}},
##'   \code{\link{plotRMSEAdist}},
##'   \code{\link{findRMSEAsamplesize}},
##'   \code{\link{findRMSEApowernested}},
##'   \code{\link{plotRMSEApowernested}}, and
##'   \code{\link{findRMSEAsamplesizenested}}}
##' \item{Missing Data Analysis:
##'   \code{\link{runMI}},
##'   \code{\link{auxiliary}},
##'   \code{\link{twostage}},
##'   \code{\link{fmi}},
##'   \code{\link{bsBootMiss}},
##'   \code{\link{quark}}, and
##'   \code{\link{combinequark}}}
##' \item{Latent Interactions:
##'   \code{\link{indProd}},
##'   \code{\link{orthogonalize}},
##'   \code{\link{probe2WayMC}},
##'   \code{\link{probe3WayMC}},
##'   \code{\link{probe2WayRC}},
##'   \code{\link{probe3WayRC}}, and
##'   \code{\link{plotProbe}}}
##' \item{Exploratory Factor Analysis (EFA):
##'   \code{\link{efa.ekc}},
##'   \code{\link{efaUnrotate}},
##'   \code{\link{orthRotate}},
##'   \code{\link{oblqRotate}}, and
##'   \code{\link{funRotate}}}
##' \item{Reliability Estimation:
##'   \code{\link{compRelSEM}} and
##'   \code{\link{maximalRelia}}
##'   (see also \code{\link{AVE}})}
##' \item{Parceling:
##'   \code{\link{parcelAllocation}},
##'   \code{\link{PAVranking}}, and
##'   \code{\link{poolMAlloc}}}
##' \item{Non-Normality:
##'   \code{\link{skew}},
##'   \code{\link{kurtosis}},
##'   \code{\link{mardiaSkew}},
##'   \code{\link{mardiaKurtosis}}, and
##'   \code{\link{mvrnonnorm}}}
##' }
##' All users of R (or SEM) are invited to submit functions or ideas for
##' functions by contacting the maintainer, Terrence Jorgensen
##' (\email{TJorgensen314@gmail.com}). Contributors are encouraged to use
##' \code{Roxygen} comments to document their contributed code, which is
##' consistent with the rest of \pkg{semTools}.  Read the vignette from the
##' \pkg{roxygen2} package for details:
##' \code{vignette("rd", package = "roxygen2")}
##'
##' @docType package
##' @name semTools
NULL



##' @importFrom methods setClass setMethod getMethod show is new slot as hasArg
NULL


##' @importFrom graphics hist plot par abline lines legend
NULL


##' @importFrom stats nobs residuals resid fitted fitted.values coef anova vcov
NULL


