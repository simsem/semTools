### Terrence D. Jorgensen
### Last updated: 8 April 2017
### package documentation, along with convenience documentation (e.g., imports)


#' semTools: Useful Tools for Structural Equation Modeling
#'
#' The \code{semTools} package provides many miscellaneous functions that are
#' useful for statistical analysis involving SEM in R.  Many functions extend
#' the funtionality of the \code{lavaan} package.  There are some sets of
#' functions in \code{semTools} that fit into certain themes, called "suites",
#' such as:
#' \item{Model Fit Evaluation}{\code{moreFitIndices}, \code{nullRMSEA},
#'    \code{singleParamTest}, \code{miPowerFit}, and \code{chisqSmallN}}
#' \item{Measurement Invariance}{\code{measurementInvariance},
#'    \code{measurementInvarianceCat}, \code{longInvariance},
#'    \code{partialInvariance}, \code{partialInvarianceCat}, and
#'    \code{permuteMeasEq}}
#' \item{Power Analysis}{\code{SSpower}, \code{findRMSEApower},
#'    \code{plotRMSEApower}, \code{plotRMSEAdist}, \code{findRMSEAsamplesize},
#'    \code{findRMSEApowernested}, \code{plotRMSEApowernested}, and
#'    \code{findRMSEAsamplesizenested}}
#' \item{Missing Data Analysis}{\code{auxiliary}, \code{runMI}, \code{twostage},
#'    \code{fmi}, \code{bsBootMiss}, \code{quark}, and \code{combinequark}}
#' \item{Latent Interactions}{\code{indProd}, \code{orthogonalize},
#'    \code{probe2WayMC}, \code{probe3WayMC}, \code{probe2WayRC},
#'    \code{probe3WayRC}, and \code{plotProbe}}
#' \item{Exploratory Factor Analysis (EFA)}{\code{efaUnrotate}, \code{efa.ekc},
#'    \code{orthRotate}, \code{oblqRotate}, and \code{funRotate}}
#' \item{Reliability of a Composite Score}{\code{reliability},
#'    \code{reliabiltyL2}, and \code{maximalRelia}}
#' \item{Parcelling}{\code{parcelAllocation}, \code{PAVranking}, and
#'    \code{poolMAlloc}}
#' \item{Non-Normality}{\code{skew}, \code{kurtosis}, \code{mardiaSkew},
#'    \code{mardiaKurtosis}, and \code{mvrnonnorm}}
#' 
#' @docType package
#' @name semTools
NULL



#' @importFrom methods setClass setMethod getMethod show is new slot as hasArg
NULL


#' @importFrom graphics hist plot par abline lines legend
NULL


#' @importFrom stats nobs residuals resid fitted fitted.values coef anova vcov
NULL


