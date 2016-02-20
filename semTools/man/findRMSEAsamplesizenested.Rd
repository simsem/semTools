\name{findRMSEAsamplesizenested}
\alias{findRMSEAsamplesizenested}
\title{Find sample size given a power in nested model comparison}
\description{
Find the sample size that the power in rejection the samples from the alternative pair of RMSEA is just over the specified power.
}
\usage{
findRMSEAsamplesizenested(rmsea0A = NULL, rmsea0B = NULL, rmsea1A, 
rmsea1B = NULL, dfA, dfB, power=0.80, alpha=.05, group=1)
}
\arguments{
  \item{rmsea0A}{The H0 baseline RMSEA.}
  \item{rmsea0B}{The H0 alternative RMSEA (trivial misfit).}
  \item{rmsea1A}{The H1 baseline RMSEA.}
  \item{rmsea1B}{The H1 alternative RMSEA (target misfit to be rejected).}
  \item{dfA}{degree of freedom of the more-restricted model.}
  \item{dfB}{degree of freedom of the less-restricted model.}
  \item{power}{The desired statistical power.}
  \item{alpha}{The alpha level.}
  \item{group}{The number of group in calculating RMSEA.}
}
\references{
MacCallum, R. C., Browne, M. W., & Cai, L. (2006). Testing differences between nested covariance structure models: Power analysis and null hypotheses. \emph{Psychological Methods, 11}, 19-35.
}
\author{
	Bell Clinton; Pavel Panko (Texas Tech University; \email{pavel.panko@ttu.edu}); Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\seealso{
	\itemize{
		\item \code{\link{plotRMSEApowernested}} to plot the statistical power for nested model comparison based on population RMSEA given the sample size
		\item \code{\link{findRMSEApowernested}} to find the power for a given sample size in nested model comparison based on population RMSEA
	}
}
\examples{
findRMSEAsamplesizenested(rmsea0A = 0, rmsea0B = 0, rmsea1A = 0.06, 
rmsea1B = 0.05, dfA = 22, dfB = 20, power=0.80, alpha=.05, group=1) 
}
