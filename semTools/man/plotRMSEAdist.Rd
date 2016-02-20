\name{plotRMSEAdist}
\alias{plotRMSEAdist}
\title{
Plot the sampling distributions of RMSEA
}
\description{
Plots the sampling distributions of RMSEA based on the noncentral chi-square distributions
}
\usage{
plotRMSEAdist(rmsea, n, df, ptile=NULL, caption=NULL, rmseaScale = TRUE, group=1)
}
\arguments{
  \item{rmsea}{The vector of RMSEA values to be plotted}
  \item{n}{Sample size of a dataset}
  \item{df}{Model degrees of freedom}
  \item{ptile}{The percentile rank of the distribution of the first RMSEA that users wish to plot a vertical line in the resulting graph}
  \item{caption}{The name vector of each element of \code{rmsea}}
  \item{rmseaScale}{If \code{TRUE}, the RMSEA scale is used in the x-axis. If \code{FALSE}, the chi-square scale is used in the x-axis.}
  \item{group}{The number of group that is used to calculate RMSEA.}
  }
\details{
This function creates overlappling plots of the sampling distribution of RMSEA based on noncentral chi-square distribution (MacCallum, Browne, & Suguwara, 1996). First, the noncentrality parameter (\eqn{\lambda}) is calculated from RMSEA (Steiger, 1998; Dudgeon, 2004) by
	\deqn{\lambda = (N - 1)d\varepsilon^2 / K,}
where \eqn{N} is sample size, \eqn{d} is the model degree of freedom, \eqn{K} is the number of groupand \eqn{\varepsilon} is the population RMSEA. Next, the noncentral chi-square distribution with a specified degree of freedom and noncentrality parameter is plotted. Thus, the x-axis represent the sample chi-square value. The sample chi-square value can be transformed to the sample RMSEA scale (\eqn{\hat{\varepsilon}}) by
	\deqn{\hat{\varepsilon} = \sqrt{K}\sqrt{\frac{\chi^2 - d}{(N - 1)d}},}
where \eqn{\chi^2} is the chi-square value obtained from the noncentral chi-square distribution.
}
\references{
Dudgeon, P. (2004). A note on extending Steiger's (1998) multiple sample RMSEA adjustment to other noncentrality parameter-based statistic. \emph{Structural Equation Modeling, 11}, 305-319.

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis and determination of sample size for covariance structure modeling. \emph{Psychological Methods, 1,} 130-149.

Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA fit index. \emph{Structural Equation Modeling, 5}, 411-419.
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\seealso{
	\itemize{
		\item \code{\link{plotRMSEApower}} to plot the statistical power based on population RMSEA given the sample size
		\item \code{\link{findRMSEApower}} to find the statistical power based on population RMSEA given a sample size
		\item \code{\link{findRMSEAsamplesize}} to find the minium sample size for a given statistical power based on population RMSEA
	}
}
\examples{
plotRMSEAdist(rmsea=c(.05, .08), n=200, df=20, ptile=0.95, rmseaScale = TRUE)
plotRMSEAdist(rmsea=c(.05, .01), n=200, df=20, ptile=0.05, rmseaScale = FALSE)
}