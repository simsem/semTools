### Terrence D. Jorgensen
### Last updated: 10 June 2024
### package documentation, along with convenience documentation (e.g., imports)


##' semTools: Useful Tools for Structural Equation Modeling
##'
##' The \pkg{semTools} package provides many miscellaneous functions that are
##' useful for statistical analysis involving SEM in R.  Many functions extend
##' the funtionality of the \pkg{lavaan} package.  Some sets of functions in
##' \pkg{semTools} correspond to the same theme. We call such a collection of
##' functions a *suite*. Our suites include:
##' \itemize{
##' \item{Model Fit Evaluation:
##'   [moreFitIndices()],
##'   [nullRMSEA()],
##'   [singleParamTest()],
##'   [miPowerFit()], and
##'   [chisqSmallN()]}
##' \item{Measurement Invariance:
##'   [measEq.syntax()],
##'   [partialInvariance()],
##'   [partialInvarianceCat()], and
##'   [permuteMeasEq()]}
##' \item{Power Analysis:
##'   [SSpower()],
##'   [findRMSEApower()],
##'   [plotRMSEApower()],
##'   [plotRMSEAdist()],
##'   [findRMSEAsamplesize()],
##'   [findRMSEApowernested()],
##'   [plotRMSEApowernested()], and
##'   [findRMSEAsamplesizenested()]}
##' \item{Missing Data Analysis:
##                     [runMI()],   # DEPRECATE
##'   [auxiliary()],
##'   [twostage()],
##'   [fmi()],
##'   [bsBootMiss()],
##'   [quark()], and
##'   [combinequark()]}
##' \item{Latent Interactions:
##'   [indProd()],
##'   [orthogonalize()],
##'   [probe2WayMC()],
##'   [probe3WayMC()],
##'   [probe2WayRC()],
##'   [probe3WayRC()], and
##'   [plotProbe()]}
##' \item{Exploratory Factor Analysis (EFA):
##'   [efa.ekc()]}
##' \item{Reliability Estimation:
##'   [compRelSEM()] and
##'   [maximalRelia()]
##'   (see also [AVE()])}
##' \item{Parceling:
##'   [parcelAllocation()],
##'   [PAVranking()], and
##'   [poolMAlloc()]}
##' \item{Non-Normality:
##'   [skew()],
##'   [kurtosis()],
##'   [mardiaSkew()],
##'   [mardiaKurtosis()], and
##'   [mvrnonnorm()]}
##' }
##' All users of R (or SEM) are invited to submit functions or ideas for
##' functions by contacting the maintainer, Terrence Jorgensen
##' (\email{TJorgensen314@gmail.com}). Contributors are encouraged to use
##' `Roxygen` comments to document their contributed code, which is
##' consistent with the rest of \pkg{semTools}.  Read the vignette from the
##' \pkg{roxygen2} package for details:
##' `vignette("rd", package = "roxygen2")`
##'
##' @name semTools
##' @aliases semTools-package
NULL



##' @importFrom methods setClass setMethod getMethod show is new slot as hasArg
NULL


##' @importFrom graphics hist plot par abline lines legend
NULL


##' @importFrom stats nobs residuals resid fitted fitted.values coef anova vcov
NULL


