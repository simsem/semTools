\name{chisqSmallN}
\alias{chisqSmallN}
\title{
  \emph{k}-factor correction for chi-squared test statistic
}
\description{
  Calculate \emph{k}-factor correction for chi-squared model-fit test statistic to adjust for small sample size.
}
\usage{
  chisqSmallN(fit0, fit1 = NULL, ...)
}
\arguments{
  \item{fit0}{The lavaan model object provided after running the \code{cfa}, \code{sem}, \code{growth}, or \code{lavaan} functions.}
  \item{fit1}{Optional additional \linkS4class{lavaan} model, in which \code{fit0} is nested.  If \code{fit0} has fewer \emph{df} than \code{fit1}, the models will be swapped, still on the assumption that they are nested.}
  \item{\dots}{Additional arguments to the \code{\link[lavaan]{lavTestLRT}} function.}
}
\details{
The \emph{k}-factor correction (Nevitt & Hancock, 2004) is a global fit index which can be computed by:

\deqn{ kc = 1 - \frac{2 \times P + 4 \times K + 5}{6 \times N}}

where \eqn{N} is the sample size when using normal likelihood, or \eqn{N - 1} when using \code{likelihood = 'wishart'}.
}

\value{
  A numeric vector including the unadjusted (naive) chi-squared test statistic, the \emph{k}-factor correction, the corrected test statistic, the \emph{df} for the test, and the \emph{p} value for the test under the null hypothesis that the model fits perfectly (or that the 2 models have equivalent fit).
}
\references{
  Nevitt, J., & Hancock, G. R. (2004). Evaluating small sample approaches for model test statistics in structural equation modeling. \emph{Multivariate Behavioral Research, 39}(3), 439-478. doi:10.1207/S15327906MBR3903_3
}
\author{
    Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\examples{
HS.model <- '
    visual  =~ x1 + b1*x2 + x3
    textual =~ x4 + b2*x5 + x6
    speed   =~ x7 + b3*x8 + x9
'
fit1 <- cfa(HS.model, data = HolzingerSwineford1939)
## test a single model (implicitly compared to a saturated model)
chisqSmallN(fit1)

## fit a more constrained model
fit0 <- cfa(HS.model, data = HolzingerSwineford1939,
            orthogonal = TRUE)
## compare 2 models
chisqSmallN(fit1, fit0)
}
