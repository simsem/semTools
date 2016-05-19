\name{mardiaKurtosis}
\alias{mardiaKurtosis}
\title{
	Finding Mardia's multivariate kurtosis
}
\description{
	Finding Mardia's multivariate kurtosis of multiple variables
}
\usage{
mardiaKurtosis(dat, use = "everything")
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
	A value of a Mardia's multivariate kurtosis with a test statistic 
}
\details{
	The Mardia's multivariate kurtosis formula (Mardia, 1970) is
	\deqn{
		b_{2, d} = \frac{1}{n}\sum^n_{i=1}\left[ \left(\bold{X}_i - \bold{\bar{X}} \right)^{'} \bold{S}^{-1} \left(\bold{X}_i - \bold{\bar{X}} \right) \right]^2,
	}
	where \eqn{d} is the number of variables, \eqn{X} is the target dataset with multiple variables, \eqn{n} is the sample size, \eqn{\bold{S}} is the sample covariance matrix of the target dataset, and \eqn{\bold{\bar{X}}} is the mean vectors of the target dataset binded in \eqn{n} rows. When the population multivariate kurtosis is normal, the \eqn{b_{2,d}} is asymptotically distributed as normal distribution with the mean of \eqn{d(d + 2)} and variance of \eqn{8d(d + 2)/n}.
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
		\item \code{\link{mardiaSkew}} Find the Mardia's multivariate skewness of a set of variables
	}
}
\examples{
library(lavaan)
mardiaKurtosis(HolzingerSwineford1939[,paste("x", 1:9, sep="")])
}
