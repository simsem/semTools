\name{efaUnrotate}
\alias{efaUnrotate}
\title{
	Analyze Unrotated Exploratory Factor Analysis Model
}
\description{
This function will analyze unrotated exploratory factor analysis model. The unrotated solution can be rotated by the \code{\link{orthRotate}} and \code{\link{oblqRotate}} functions.
}
\usage{
efaUnrotate(data, nf, varList=NULL, start=TRUE, aux=NULL, ...)
}
\arguments{
  \item{data}{
	A target data frame. 
}
  \item{nf}{
	The desired number of factors
}
  \item{varList}{
	Target observed variables. If not specified, all variables in the target data frame will be used.
}
  \item{start}{
	Use starting values in the analysis from the \code{\link{factanal}} function. If \code{FALSE}, the starting values from the \code{lavaan} package will be used.
}
  \item{aux}{
	The list of auxiliary variables. These variables will be included in the model by the saturated-correlates approach to account for missing information.
}
  \item{\dots}{
	Other arguments in the \code{\link[lavaan]{cfa}} function in the \code{lavaan} package, such as \code{ordered}, \code{se}, or \code{estimator}
}
}
\details{
	This function will generate a lavaan script for unrotated exploratory factor analysis model such that 1) all factor loadings are estimated, 2) factor variances are fixed to 1, 3) factor covariances are fixed to 0, and 4) the dot products of any pairs of columns in the factor loading matrix are fixed to zero (Johnson and Wichern, 2002). The reason for creating this function in addition to the \code{\link{factanal}} function is that users can enjoy some advanced features from the \code{lavaan} package such as scaled chi-square, diagonal weighted least square for ordinal indicators, or full-information maximum likelihood.
}
\value{
	A \code{lavaan} output of unrotated exploratory factor analysis solution.
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
unrotated <- efaUnrotate(HolzingerSwineford1939, nf=3, varList=paste0("x", 1:9), estimator="mlr")
summary(unrotated, std=TRUE)
inspect(unrotated, "std")

dat <- data.frame(HolzingerSwineford1939, z=rnorm(nrow(HolzingerSwineford1939), 0, 1))
unrotated2 <- efaUnrotate(dat, nf=2, varList=paste0("x", 1:9), aux="z")
}
