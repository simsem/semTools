% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerAnalysisNested.R
\name{findRMSEApowernested}
\alias{findRMSEApowernested}
\title{Find power given a sample size in nested model comparison}
\usage{
findRMSEApowernested(rmsea0A = NULL, rmsea0B = NULL, rmsea1A,
  rmsea1B = NULL, dfA, dfB, n, alpha = 0.05, group = 1)
}
\arguments{
\item{rmsea0A}{The \eqn{H_0} baseline RMSEA}

\item{rmsea0B}{The \eqn{H_0} alternative RMSEA (trivial misfit)}

\item{rmsea1A}{The \eqn{H_1} baseline RMSEA}

\item{rmsea1B}{The \eqn{H_1} alternative RMSEA (target misfit to be rejected)}

\item{dfA}{degree of freedom of the more-restricted model}

\item{dfB}{degree of freedom of the less-restricted model}

\item{n}{Sample size}

\item{alpha}{The alpha level}

\item{group}{The number of group in calculating RMSEA}
}
\description{
Find the sample size that the power in rejection the samples from the
alternative pair of RMSEA is just over the specified power.
}
\examples{

findRMSEApowernested(rmsea0A = 0.06, rmsea0B = 0.05, rmsea1A = 0.08,
                     rmsea1B = 0.05, dfA = 22, dfB = 20, n = 200,
                     alpha = 0.05, group = 1)

}
\references{
MacCallum, R. C., Browne, M. W., & Cai, L. (2006). Testing
differences between nested covariance structure models: Power analysis and
null hypotheses. \emph{Psychological Methods, 11}(1), 19--35.
\doi{10.1037/1082-989X.11.1.19}
}
\seealso{
\itemize{
\item \code{\link[=plotRMSEApowernested]{plotRMSEApowernested()}} to plot the statistical power for
nested model comparison based on population RMSEA given the sample size
\item \code{\link[=findRMSEAsamplesizenested]{findRMSEAsamplesizenested()}} to find the minium sample
size for a given statistical power in nested model comparison based on
population RMSEA
}
}
\author{
Bell Clinton

Pavel Panko (Texas Tech University; \email{pavel.panko@ttu.edu})

Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
