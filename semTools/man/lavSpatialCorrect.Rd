\name{spatialCorrect}
\alias{spatialCorrect}
\title{
	Correcting Standard Errors for Spatial Autocorrelation
}
\description{
This function addresses the problem of correcting sample sizes and standard errors in the presence of spatial autocorrelation in Structural Equation Models with with spatial structure in the autocorrelation of endogenous variables.
}
\usage{
spatialCorrect(obj, xvar, yvar, alpha=0.05)
}
\arguments{
  \item{obj}{
	A target \code{lavaan} object
}
  \item{xvar}{
	A variable indicating X-coordinate
}
  \item{yvar}{
	A variable indicating Y-coordinate
}
  \item{alpha}{
	The alpha level using for 
}
}
\details{
	This function will generate a lavaan script for unrotated exploratory factor analysis model such that 1) all factor loadings are estimated, 2) factor variances are fixed to 1, 3) factor covariances are fixed to 0, and 4) the dot products of any pairs of columns in the factor loading matrix are fixed to zero (Johnson and Wichern, 2002). The reason for creating this function in addition to the \code{\link{factanal}} function is that users can enjoy some advanced features from the \code{lavaan} package such as scaled chi-square, diagonal weighted least square for ordinal indicators, or full-information maximum likelihood.
}
\value{
	A \code{lavaan} output of unrotated exploratory factor analysis solution.
}
\author{
    Sunthud Pornprasertmanit (Texas Tech University; \email{sunthud.pornprasertmanit@ttu.edu})
}
\examples{
unrotated <- efaUnrotate(HolzingerSwineford1939, nf=3, varList=paste0("x", 1:9), estimator="mlr")
summary(unrotated, std=TRUE)
inspect(unrotated, "std")

dat <- data.frame(HolzingerSwineford1939, z=rnorm(nrow(HolzingerSwineford1939), 0, 1))
unrotated2 <- efaUnrotate(dat, nf=2, varList=paste0("x", 1:9), aux="z")
}
