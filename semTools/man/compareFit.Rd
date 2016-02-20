\name{compareFit}
\alias{compareFit}
\title{
	Build an object summarizing fit indices across multiple models
}
\description{
This function will create the template that compare fit indices across multiple lavaan outputs. The results can be exported to a clipboard or a file later. 
}
\usage{
compareFit(..., nested = TRUE)
}
\arguments{
  \item{...}{
	\code{lavaan} outputs or lists of \code{lavaan} outputs
}
  \item{nested}{
	Logical whether the specified models are nested
}
}
\value{
	A \code{\linkS4class{FitDiff}} object that saves model fit comparisons across multiple models. If the output is not assigned as an object, the output is printed in two parts: 1) nested model comparison (if models are nested) and 2) fit indices summaries. In the fit indices summaries, daggers are tagged to the model with the best fit for each fit index.
}
\seealso{
	\code{\linkS4class{FitDiff}}, \code{\link{clipboard}}
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
m1 <- ' visual  =~ x1 + x2 + x3
        textual =~ x4 + x5 + x6
        speed   =~ x7 + x8 + x9 '

fit1 <- cfa(m1, data=HolzingerSwineford1939)

m2 <- ' f1  =~ x1 + x2 + x3 + x4 
        f2 =~ x5 + x6 + x7 + x8 + x9 '
fit2 <- cfa(m2, data=HolzingerSwineford1939)
compareFit(fit1, fit2, nested=FALSE)

HW.model <- ' visual =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed =~ x7 + x8 + x9 '

out <- measurementInvariance(HW.model, data=HolzingerSwineford1939, group="school", quiet=TRUE)
compareFit(out)
}
