\name{skew}
\alias{skew}
\title{
	Finding skewness
}
\description{
	Finding skewness (g1) of an object
}
\usage{
skew(object, population=FALSE)
}
\arguments{
  \item{object}{
	A vector used to find a skewness
}
  \item{population}{
	\code{TRUE} to compute the parameter formula. \code{FALSE} to compute the sample statistic formula.
  }
}
\value{
	A value of a skewness with a test statistic if the population is specified as \code{FALSE}
}
\details{
	The skewness computed is g1. The parameter skewness \eqn{\gamma_{2}} formula is
	
	\deqn{\gamma_{2} = \frac{\mu_{3}}{\mu^{3/2}_{2}},}
	
	where \eqn{\mu_{i}} denotes the \eqn{i} order central moment.
	
	The excessive kurtosis formula for sample statistic \eqn{g_{2}} is
	
	\deqn{g_{2} = \frac{k_{3}}{k^{2}_{2}},}
	
	where \eqn{k_{i}} are the \eqn{i} order \emph{k}-statistic.
	
	The standard error of the skewness is 
	
	\deqn{Var(\hat{g}_2) = \frac{6}{N}}
	
	where \eqn{N} is the sample size.
}
\references{
Weisstein, Eric W. (n.d.). \emph{Skewness.} Retrived from MathWorld--A Wolfram Web Resource \url{http://mathworld.wolfram.com/Skewness.html} 
}
\seealso{
	\itemize{
		\item \code{\link{kurtosis}} Find the univariate excessive kurtosis of a variable
		\item \code{\link{mardiaSkew}} Find the Mardia's multivariate skewness of a set of variables
		\item \code{\link{mardiaKurtosis}} Find the Mardia's multivariate kurtosis of a set of variables
	}
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
skew(1:5)
}
