% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerAnalysisSS.R
\name{SSpower}
\alias{SSpower}
\title{Power for model parameters}
\usage{
SSpower(powerModel, n, nparam, popModel, mu, Sigma, fun = "sem",
  alpha = 0.05, ...)
}
\arguments{
\item{powerModel}{lavaan \code{\link[lavaan:model.syntax]{lavaan::model.syntax()}} for the model to
be analyzed. This syntax should constrain at least one nonzero parameter
to 0 (or another number).}

\item{n}{\code{integer}. Sample size used in power calculation, or a vector
of sample sizes if analyzing a multigroup model. If
\code{length(n) < length(Sigma)} when \code{Sigma} is a list, \code{n} will
be recycled. If \code{popModel} is used instead of \code{Sigma}, \code{n}
must specify a sample size for each group, because that is used to infer
the number of groups.}

\item{nparam}{\code{integer}. Number of invalid constraints in \code{powerModel}.}

\item{popModel}{lavaan \code{\link[lavaan:model.syntax]{lavaan::model.syntax()}} specifying the
data-generating model. This syntax should specify values for all nonzero
parameters in the model. If \code{length(n) > 1}, the same population
values will be used for each group, unless different population values are
specified per group, either in the lavaan \code{\link[lavaan:model.syntax]{lavaan::model.syntax()}}
or by utilizing a list of \code{Sigma} (and optionally \code{mu}).}

\item{mu}{\code{numeric} or \code{list}. For a single-group model, a vector
of population means. For a multigroup model, a list of vectors (one per
group). If \code{mu} and \code{popModel} are missing, mean structure will
be excluded from the analysis.}

\item{Sigma}{\code{matrix} or \code{list}. For a single-group model,
a population covariance matrix. For a multigroup model, a list of matrices
(one per group). If missing, \code{popModel} will be used to generate a
model-implied Sigma.}

\item{fun}{character. Name of \code{lavaan} function used to fit
\code{powerModel} (i.e., \code{"cfa"}, \code{"sem"}, \code{"growth"}, or
\code{"lavaan"}).}

\item{alpha}{Type I error rate used to set a criterion for rejecting H0.}

\item{...}{additional arguments to pass to \code{\link[lavaan:lavaan]{lavaan::lavaan()}}.
See also \code{\link[lavaan:lavOptions]{lavaan::lavOptions()}}.}
}
\description{
Apply Satorra & Saris (1985) method for chi-squared power analysis.
}
\details{
Specify all non-zero parameters in a population model, either by using
lavaan syntax (\code{popModel}) or by submitting a population covariance
matrix (\code{Sigma}) and optional mean vector (\code{mu}) implied by the
population model. Then specify an analysis model that places at least
one invalid constraint (note the number in the \code{nparam} argument).

There is also a Shiny app called "power4SEM" that provides a graphical user
interface for this functionality (Jak et al., in press).  It can be accessed
at \url{https://sjak.shinyapps.io/power4SEM/}.
}
\examples{
## Specify population values. Note every parameter has a fixed value.
modelP <- '
  f1 =~ .7*V1 + .7*V2 + .7*V3 + .7*V4
  f2 =~ .7*V5 + .7*V6 + .7*V7 + .7*V8
  f1 ~~ .3*f2
  f1 ~~ 1*f1
  f2 ~~ 1*f2
  V1 ~~ .51*V1
  V2 ~~ .51*V2
  V3 ~~ .51*V3
  V4 ~~ .51*V4
  V5 ~~ .51*V5
  V6 ~~ .51*V6
  V7 ~~ .51*V7
  V8 ~~ .51*V8
'
## Specify analysis model. Note parameter of interest f1~~f2 is fixed to 0.
modelA <- '
  f1 =~ V1 + V2 + V3 + V4
  f2 =~ V5 + V6 + V7 + V8
  f1 ~~ 0*f2
'
## Calculate power
SSpower(powerModel = modelA, popModel = modelP, n = 150, nparam = 1,
        std.lv = TRUE)

## Get power for a range of sample sizes
Ns <- seq(100, 500, 40)
Power <- rep(NA, length(Ns))
for(i in 1:length(Ns)) {
  Power[i] <- SSpower(powerModel = modelA, popModel = modelP,
                      n = Ns[i], nparam = 1, std.lv = TRUE)
}
plot(x = Ns, y = Power, type = "l", xlab = "Sample Size")


## Optionally specify different values for multiple populations

modelP2 <- '
  f1 =~ .7*V1 + .7*V2 + .7*V3 + .7*V4
  f2 =~ .7*V5 + .7*V6 + .7*V7 + .7*V8
  f1 ~~ c(-.3, .3)*f2                  # DIFFERENT ACROSS GROUPS
  f1 ~~ 1*f1
  f2 ~~ 1*f2
  V1 ~~ .51*V1
  V2 ~~ .51*V2
  V3 ~~ .51*V3
  V4 ~~ .51*V4
  V5 ~~ .51*V5
  V6 ~~ .51*V6
  V7 ~~ .51*V7
  V8 ~~ .51*V8
'
modelA2 <- '
  f1 =~ V1 + V2 + V3 + V4
  f2 =~ V5 + V6 + V7 + V8
  f1 ~~ c(psi21, psi21)*f2        # EQUALITY CONSTRAINT ACROSS GROUPS
'
## Calculate power
SSpower(powerModel = modelA2, popModel = modelP2, n = c(100, 100), nparam = 1,
        std.lv = TRUE)
## Get power for a range of sample sizes
Ns2 <- cbind(Group1 = seq(10, 100, 10), Group2 = seq(10, 100, 10))
Power2 <- apply(Ns2, MARGIN = 1, FUN = function(nn) {
  SSpower(powerModel = modelA2, popModel = modelP2, n = nn,
          nparam = 1, std.lv = TRUE)
})
plot(x = rowSums(Ns2), y = Power2, type = "l", xlab = "Total Sample Size",
     ylim = 0:1)
abline(h = c(.8, .9), lty = c("dotted","dashed"))
legend("bottomright", c("80\% Power","90\% Power"), lty = c("dotted","dashed"))

}
\references{
Satorra, A., & Saris, W. E. (1985). Power of the likelihood ratio
test in covariance structure analysis. \emph{Psychometrika, 50}(1), 83--90.
\doi{10.1007/BF02294150}

Jak, S., Jorgensen, T. D., Verdam, M. G., Oort, F. J., & Elffers, L.
(2021). Analytical power calculations for structural equation modeling:
A tutorial and Shiny app. \emph{Behavior Research Methods, 53}, 1385--1406.
\doi{10.3758/s13428-020-01479-0}
}
\author{
Alexander M. Schoemann (East Carolina University; \email{schoemanna@ecu.edu})

Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
