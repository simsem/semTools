\name{spatialCorrect}
\alias{spatialCorrect}
\title{
Calculate reliability values of factors
}
\description{
Correcting sample sizes and standard errors in the presence of spatial autocorrelation in Structural Equation Models with with spatial structure in the autocorrelation of endogenous variables.
}
\usage{
spatialCorrect(obj, xvar, yvar, alpha=0.05)
}
\arguments{
  \item{obj}{The lavaan model object provided after running the \code{cfa}, \code{sem}, \code{growth}, or \code{lavaan} functions.}
  \item{xvar}{X-coordinate.}
  \item{yvar}{Y-coordinate.}
  \item{alpha}{The alpha level used to decide whether the effective sample size is used. If Moran's I is not significant, the original sample size is used for test statistic. Otherwise, the calculated effective sample size is used.}
}
\details{
This function implements the standard error correction for all endogenous variables, using Moran's I and an approximation of an effective sample size for large sample sizes. This function does not support the data with missing observations.
}
\value{
	\itemize{
		\item \emph{Morans_I} Moran's I statistics and effective sample sizes for all endogeneous variables.
		\item \emph{parameters} Adjusted standard errors and test statistics for each parameter estimate for each endogeneous variables.
	}
}
\author{
    Jarrett Byrnes (University of Massachusetts, Boston; \email{Jarrett.Byrnes@umb.edu})
}
\examples{
library(lavaan)
borModel <- '
  NDVI ~ nTot + T61 + Wet 
  nTot ~ T61
'
#note meanstructure = TRUE to obtain intercepts
borFit <- sem(borModel, data = boreal, meanstructure = TRUE)
spatialCorrect(borFit, boreal$x, boreal$y)
}
