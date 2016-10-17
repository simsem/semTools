\name{twostage-class}
\docType{class}
\alias{twostage-class}
\alias{show,twostage-method}
\alias{summary,twostage-method}
\alias{anova,twostage-method}
\alias{vcov,twostage-method}
\alias{coef,twostage-method}
\alias{fitted.values,twostage-method}
\alias{fitted,twostage-method}
\alias{residuals,twostage-method}
\alias{resid,twostage-method}
\alias{nobs,twostage-method}
\title{
	Class for the Results of 2-Stage Maximum Likelihood (TSML) Estimation for Missing Data
}
\description{
	This class contains the results of 2-Stage Maximum Likelihood (TSML) estimation for missing data.  The \code{summary}, \code{anova}, \code{vcov} methods return corrected \emph{SE}s and test statistics.  Other methods are simply wrappers around the corresponding \code{\linkS4class{lavaan}} methods.
}
\section{Objects from the Class}{
	Objects can be created via the \code{\link{twostage}} function.
}
\section{Slots}{
  \describe{
    \item{\code{saturated}:}{A fitted \code{\linkS4class{lavaan}} object containing the saturated model results.}
    \item{\code{target}:}{A fitted \code{\linkS4class{lavaan}} object containing the target/hypothesized model results.}
    \item{\code{baseline}:}{A fitted \code{\linkS4class{lavaan}} object containing the baseline/null model results.}
    \item{\code{auxNames}:}{A character string (potentially of \code{length == 0}) of any auxiliary variable names, if used.}
  }
}
\section{methods}{
  \describe{
    \item{anova}{\code{signature(object = "twostage", h1 = NULL, baseline = FALSE:} The \code{anova} function returns the residual-based chi-squared test statistic result, as well as the scaled chi-squared test statistic result, for the model in the \code{target} slot, or for the model in the \code{baseline} slot if \code{baseline = TRUE}.  The user can also provide a single additional \code{twostage} object to the \code{h1} argument, in which case \code{anova} returns residual-based and scaled chi-squared difference test results, under the assumption that the models are nested.  The models will be automatically sorted according their degrees of freedom.}
    \item{show}{\code{signature(object = "twostage"):} The \code{show} function is used to display the results of the \code{anova} method, as well as the header of the (uncorrected) target model results.}
    \item{summary}{\code{signature(object = "twostage", ...):} The summary function prints the same information from the \code{show} method, but also provides (and returns) the output of \code{\link[lavaan]{parameterEstimates}(object@target, ...)} with corrected \emph{SE}s, test statistics, and confidence intervals.  Additional arguments can be passed to \code{\link[lavaan]{parameterEstimates}}, including \code{fmi = TRUE} to provide an estimate of the fraction of missing information.}
    \item{vcov}{\code{signature(object = "twostage", baseline = FALSE:} Returns the asymptotic covariance matrix of the estimated parameters (corrected for additional uncertainty due to missing data) for the model in the \code{target} slot, or for the model in the \code{baseline} slot if \code{baseline = TRUE}.}
    \item{nobs}{\code{signature(object = "twostage", type = c("ntotal", "ngroups", "n.per.group", "norig", "patterns", "coverage")):} The \code{nobs} function will return the total sample sized used in the analysis by default.  Also available are the number of groups or the sample size per group, the original sample size (if any rows were deleted because all variables were missing), the missing data patterns, and the matrix of coverage (diagonal is the proportion of sample observed on each variable, and off-diagonal is the proportion observed for both of each pair of variables).}
    \item{coef}{\code{signature(object = "twostage", type = c("free", "user"):} This is simply a wrapper around the corresponding \code{\linkS4class{lavaan}} method, providing point estimates from the \code{target} slot.}
    \item{fitted.values}{\code{signature(object = "twostage", model = c("target", "saturated", "baseline"):} This is simply a wrapper around the corresponding \code{\linkS4class{lavaan}} method, providing model-implied sample moments from the slot specified in the \code{model} argument.}
    \item{fitted}{\code{signature(object = "twostage", model = c("target", "saturated", "baseline"):} an alias for \code{fitted.values}.}
    \item{residuals}{\code{signature(object = "twostage", type = c("raw", "cor", "normalized", "standardized"):} This is simply a wrapper around the corresponding \code{\linkS4class{lavaan}} method, providing residuals of the specified \code{type} from the \code{target} slot.}
    \item{resid}{\code{signature(object = "twostage", model = c("raw", "cor", "normalized", "standardized"):} an alias for \code{residuals}.}
  }
}
\author{
    Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\seealso{
\code{\link{twostage}}
}
\examples{
# See the example from the twostage function
}
