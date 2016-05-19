\name{measurementInvariance}
\alias{measurementInvariance}
\alias{measurementinvariance}
\title{
Measurement Invariance Tests
}
\description{
Testing measurement invariance across groups using a typical sequence of
model comparison tests.
}
\usage{
measurementInvariance(..., std.lv = FALSE, strict = FALSE, quiet = FALSE,
fit.measures = "default", method = "satorra.bentler.2001")
}
\arguments{
  \item{...}{The same arguments as for any lavaan model.
    See \code{\link{cfa}} for more information.}
  \item{std.lv}{If \code{TRUE}, the fixed-factor method of scale identification is used. If \code{FALSE}, the first variable for each factor is used as marker variable.}
  \item{strict}{If \code{TRUE}, the sequence requires `strict' invariance.
    See details for more information.}
  \item{quiet}{If \code{FALSE} (default), a summary is printed out containing an
    overview of the different models that are fitted, together with some
    model comparison tests. If \code{TRUE}, no summary is printed.}
  \item{fit.measures}{Fit measures used to calculate the differences between nested models.}
  \item{method}{The method used to calculate likelihood ratio test. See \code{\link[lavaan]{lavTestLRT}} for available options}
}
\details{
If \code{strict = FALSE}, the following four models are tested in order:
\enumerate{
  \item{Model 1: configural invariance. The same factor structure is imposed
        on all groups.}
  \item{Model 2: weak invariance. The factor loadings are constrained to be
        equal across groups.}
  \item{Model 3: strong invariance. The factor loadings and intercepts are
        constrained to be equal across groups.}
  \item{Model 4: The factor loadings, intercepts and means are constrained to
        be equal across groups.}
}
Each time a more restricted model is fitted, a chi-square difference test
is reported, comparing the current model with the previous one, and comparing
the current model to the baseline model (Model 1). In addition, the difference
in cfi is also reported (delta.cfi).

If \code{strict = TRUE}, the following five models are tested in order:
\enumerate{
  \item{Model 1: configural invariance. The same factor structure is imposed
        on all groups.}
  \item{Model 2: weak invariance. The factor loadings are constrained to be
        equal across groups.}
  \item{Model 3: strong invariance. The factor loadings and intercepts are
        constrained to be equal across groups.}
  \item{Model 4: strict invariance. The factor loadings, intercepts and
        residual variances are constrained to be equal across groups.}
  \item{Model 5: The factor loadings, intercepts, residual variances and means
        are constrained to be equal across groups.}
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
    Yves Rosseel (Ghent University; \email{Yves.Rosseel@UGent.be});
	Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\seealso{
	\code{\link{longInvariance}} for the measurement invariance test within person; \code{partialInvariance} for the automated function for finding partial invariance models
}
\examples{
HW.model <- ' visual =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed =~ x7 + x8 + x9 '

measurementInvariance(HW.model, data=HolzingerSwineford1939, group="school")
}
