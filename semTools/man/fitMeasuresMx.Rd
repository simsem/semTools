\name{fitMeasuresMx}
\alias{fitMeasuresMx}
\title{
	Find fit measures from an MxModel result
}
\description{
	Find fit measures from MxModel result. The saturate and null models are analyzed in the function and fit measures are calculated based on the comparison with the null and saturate models. The function is adjusted from the \code{fitMeasures} function in the lavaan package.
}
\usage{
fitMeasuresMx(object, fit.measures="all")
}
\arguments{
  \item{object}{
	The target \code{MxModel} object
}
  \item{fit.measures}{
	Target fit measures
  }
}
\value{
	A vector of fit measures
}
\seealso{
	\code{\link{nullMx}}, \code{\link{saturateMx}}, \code{\link{standardizeMx}}
}
\author{
    The original function is the \code{fitMeasures} function written by Yves Rosseel in the \code{lavaan} package. The function is adjusted for an \code{MxModel} object by Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
\dontrun{
library(OpenMx)
data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("One Factor", 
    type="RAM",
    manifestVars=manifests, 
    latentVars=latents,
    mxPath(from=latents, to=manifests),
    mxPath(from=manifests, arrows=2),
    mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
    mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
)
factorFit <- mxRun(factorModel)
round(fitMeasuresMx(factorFit), 3)

# Compare with lavaan
library(lavaan)
script <- "f1 =~ x1 + x2 + x3 + x4 + x5"
fitMeasures(cfa(script, sample.cov = cov(demoOneFactor), sample.nobs = 500, std.lv = TRUE))
}
}
