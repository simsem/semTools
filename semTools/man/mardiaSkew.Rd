\name{mardiaSkew}
\alias{mardiaSkew}
\title{
	Finding Mardia's multivariate skewness
}
\description{
	Finding Mardia's multivariate skewness of multiple variables
}
\usage{
mardiaSkew(dat, use = "everything")
}
\arguments{
  \item{dat}{
	The target matrix or data frame with multiple variables
}
  \item{use}{
	Missing data handling method from the \code{\link[stats]{cov}} function.
}
}
\value{
	A value of a Mardia's multivariate skewness with a test statistic 
}
\details{
	The Mardia's multivariate skewness formula (Mardia, 1970) is
	\deqn{
		b_{1, d} = \frac{1}{n^2}\sum^n_{i=1}\sum^n_{j=1}\left[ \left(\bold{X}_i - \bold{\bar{X}} \right)^{'} \bold{S}^{-1} \left(\bold{X}_j - \bold{\bar{X}} \right) \right]^3,
	}
	where \eqn{d} is the number of variables, \eqn{X} is the target dataset with multiple variables, \eqn{n} is the sample size, \eqn{\bold{S}} is the sample covariance matrix of the target dataset, and \eqn{\bold{\bar{X}}} is the mean vectors of the target dataset binded in \eqn{n} rows. When the population multivariate skewness is normal, the \eqn{\frac{n}{6}b_{1,d}} is asymptotically distributed as chi-square distribution with \eqn{d(d + 1)(d + 2)/6} degrees of freedom.
}
\references{
Mardia, K. V. (1970). Measures of multivariate skewness and kurtosis with applications. \emph{Biometrika, 57}, 519-530.
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\seealso{
	\itemize{
		\item \code{\link{skew}} Find the univariate skewness of a variable
		\item \code{\link{kurtosis}} Find the univariate excessive kurtosis of a variable
		\item \code{\link{mardiaKurtosis}} Find the Mardia's multivariate kurtosis of a set of variables
	}
}
\examples{
library(lavaan)
mardiaSkew(HolzingerSwineford1939[,paste("x", 1:9, sep="")])
}
