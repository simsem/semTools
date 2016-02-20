\name{longInvariance}
\alias{longInvariance}
\alias{longInvariance}
\title{
Measurement Invariance Tests Within Person
}
\description{
Testing measurement invariance across timepoints (longitudinal) or any context involving the use of the same scale in one case (e.g., a dyad case with husband and wife answering the same scale). The measurement invariance uses a typical sequence of model comparison tests. This function currently works with only one scale.
}
\usage{
longInvariance(model, varList, auto = "all", constrainAuto = FALSE, 
fixed.x = TRUE, std.lv = FALSE, group=NULL, group.equal="", 
group.partial="", warn=TRUE, debug=FALSE, strict = FALSE, quiet = FALSE, 
fit.measures = "default", method = "satorra.bentler.2001", ...)
}
\arguments{
	\item{model}{lavaan syntax or parameter table}
	\item{varList}{A list containing indicator names of factors used in the invariance testing, such as the list that the first element is the vector of indicator names in the first timepoint and the second element is the vector of indicator names in the second timepoint. The order of indicator names should be the same (but measured in different times or different units).}
	\item{auto}{The order of autocorrelation on the measurement errors on the similar items across factor (e.g., Item 1 in Time 1 and Time 2). If 0 is specified, the autocorrelation will be not imposed. If 1 is specified, the autocorrelation will imposed for the adjacent factor listed in \code{varList}. The maximum number can be specified is the number of factors specified minus 1. If \code{"all"} is specified, the maximum number of order will be used.}
	\item{constrainAuto}{If \code{TRUE}, the function will equate the auto-\emph{covariance} to be equal within the same item across factors. For example, the covariance of item 1 in time 1 and time 2 is equal to the covariance of item 1 in time 2 and time 3.}
	\item{fixed.x}{See \code{\link[lavaan]{lavaan}.}}
	\item{std.lv}{See \code{\link[lavaan]{lavaan}.}}
	\item{group}{See \code{\link[lavaan]{lavaan}.}}
	\item{group.equal}{See \code{\link[lavaan]{lavaan}.}}
	\item{group.partial}{See \code{\link[lavaan]{lavaan}.}}
	\item{warn}{See \code{\link[lavaan]{lavaan}.}}
	\item{debug}{See \code{\link[lavaan]{lavaan}.}}
	\item{strict}{If \code{TRUE}, the sequence requires `strict' invariance. See details for more information.}
	\item{quiet}{If \code{TRUE}, a summary is printed out containing an overview of the different models that are fitted, together with some model comparison tests.}
	\item{fit.measures}{Fit measures used to calculate the differences between nested models.}
	\item{method}{The method used to calculate likelihood ratio test. See \code{\link[lavaan]{lavTestLRT}} for available options}
	\item{...}{Additional arguments in the \code{\link[lavaan]{lavaan}} function.}
}
\details{
If \code{strict = FALSE}, the following four models are tested in order:
\enumerate{
  \item{Model 1: configural invariance. The same factor structure is imposed
        on all units.}
  \item{Model 2: weak invariance. The factor loadings are constrained to be
        equal across units.}
  \item{Model 3: strong invariance. The factor loadings and intercepts are
        constrained to be equal across units.}
  \item{Model 4: The factor loadings, intercepts and means are constrained to
        be equal across units.}
}
Each time a more restricted model is fitted, a chi-square difference test
is reported, comparing the current model with the previous one, and comparing
the current model to the baseline model (Model 1). In addition, the difference
in cfi is also reported (delta.cfi).

If \code{strict = TRUE}, the following five models are tested in order:
\enumerate{
  \item{Model 1: configural invariance. The same factor structure is imposed
        on all units.}
  \item{Model 2: weak invariance. The factor loadings are constrained to be
        equal across units.}
  \item{Model 3: strong invariance. The factor loadings and intercepts are
        constrained to be equal across units.}
  \item{Model 4: strict invariance. The factor loadings, intercepts and
        residual variances are constrained to be equal across units.}
  \item{Model 5: The factor loadings, intercepts, residual variances and means
        are constrained to be equal across units.}
}

Note that if the chi-square test statistic is scaled (eg. a
Satorra-Bentler or Yuan-Bentler test statistic), a special version of the
chi-square difference test is used as described in
\url{http://www.statmodel.com/chidiff.shtml}
}

\value{
Invisibly, all model fits in the sequence are returned as a list.
}
\references{
Vandenberg, R. J., and Lance, C. E. (2000). A review and synthesis of the measurement invariance literature: Suggestions, practices, and recommendations for organizational research. \emph{Organizational Research Methods, 3,} 4-70.
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com}); Yves Rosseel (Ghent University; \email{Yves.Rosseel@UGent.be})
}
\seealso{
	\code{\link{measurementinvariance}} For the measurement invariance test between groups
}
\examples{
model <- ' f1t1 =~ y1t1 + y2t1 + y3t1
              f1t2 =~ y1t2 + y2t2 + y3t2
			  f1t3 =~ y1t3 + y2t3 + y3t3'

# Create list of variables
var1 <- c("y1t1", "y2t1", "y3t1")
var2 <- c("y1t2", "y2t2", "y3t2")
var3 <- c("y1t3", "y2t3", "y3t3")
constrainedVar <- list(var1, var2, var3)

# Invariance of the same factor across timepoints
longInvariance(model, auto=1, constrainAuto=TRUE, varList=constrainedVar, data=exLong)

# Invariance of the same factor across timepoints and groups
longInvariance(model, auto=1, constrainAuto=TRUE, varList=constrainedVar, data=exLong, group="sex", 
	group.equal=c("loadings", "intercepts"))
}
