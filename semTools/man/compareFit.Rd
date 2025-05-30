% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compareFit.R
\name{compareFit}
\alias{compareFit}
\title{Build an object summarizing fit indices across multiple models}
\usage{
compareFit(..., nested = TRUE, argsLRT = list(), indices = TRUE,
  moreIndices = FALSE, baseline.model = NULL, nPrior = 1)
}
\arguments{
\item{...}{fitted \code{lavaan} models or list(s) of \code{lavaan} objects.
\link[lavaan.mi:lavaan.mi-class]{lavaan.mi::lavaan.mi} objects are also accepted, but all models
must belong to the same class.}

\item{nested}{\code{logical} indicating whether the models in \code{...} are
nested. See \code{\link[=net]{net()}} for an empirical test of nesting.}

\item{argsLRT}{\code{list} of arguments to pass to
\code{\link[lavaan:lavTestLRT]{lavaan::lavTestLRT()}}, as well as to
\code{\link[lavaan.mi:lavTestLRT.mi]{lavaan.mi::lavTestLRT.mi()}} and \code{\link[lavaan:fitMeasures]{lavaan::fitMeasures()}} when
comparing \link[lavaan.mi:lavaan.mi-class]{lavaan.mi::lavaan.mi} models.}

\item{indices}{\code{logical} indicating whether to return fit indices from
the \code{\link[lavaan:fitMeasures]{lavaan::fitMeasures()}} function. Selecting particular
indices is controlled in the \code{summary} method; see
\linkS4class{FitDiff}.}

\item{moreIndices}{\code{logical} indicating whether to return fit indices
from the \code{\link[=moreFitIndices]{moreFitIndices()}} function. Selecting particular
indices is controlled in the \code{summary} method; see
\linkS4class{FitDiff}.}

\item{baseline.model}{optional fitted \link[lavaan:lavaan-class]{lavaan::lavaan} model passed to
\code{\link[lavaan:fitMeasures]{lavaan::fitMeasures()}} to calculate incremental fit indices.}

\item{nPrior}{passed to \code{\link[=moreFitIndices]{moreFitIndices()}}, if relevant}
}
\value{
A \linkS4class{FitDiff} object that saves model fit
comparisons across multiple models. If the models are not nested, only
fit indices for each model are returned. If the models are nested, the
differences in fit indices are additionally returned, as well as test
statistics comparing each sequential pair of models (ordered by their
degrees of freedom).
}
\description{
This function will create the template to compare fit indices across
multiple fitted lavaan objects. The results can be exported to a clipboard
or a file later.
}
\examples{

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

## non-nested models
fit1 <- cfa(HS.model, data = HolzingerSwineford1939)

m2 <- ' f1 =~ x1 + x2 + x3 + x4
        f2 =~ x5 + x6 + x7 + x8 + x9 '
fit2 <- cfa(m2, data = HolzingerSwineford1939)

(out1 <- compareFit(fit1, fit2, nested = FALSE))
summary(out1)


## nested model comparisons: measurement equivalence/invariance
fit.config <- cfa(HS.model, data = HolzingerSwineford1939, group = "school")
fit.metric <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
                  group.equal = "loadings")
fit.scalar <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
                  group.equal = c("loadings","intercepts"))
fit.strict <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
                  group.equal = c("loadings","intercepts","residuals"))

measEqOut <- compareFit(fit.config, fit.metric, fit.scalar, fit.strict,
                        moreIndices = TRUE) # include moreFitIndices()
summary(measEqOut)
summary(measEqOut, fit.measures = "all")
summary(measEqOut, fit.measures = c("aic", "bic", "sic", "ibic"))



}
\seealso{
\linkS4class{FitDiff}, \code{\link[=clipboard]{clipboard()}}
}
\author{
Terrence D. Jorgensen (University of Amsterdam;
\email{TJorgensen314@gmail.com})

Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
