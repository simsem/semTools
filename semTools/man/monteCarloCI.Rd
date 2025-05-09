% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/monteCarloCI.R
\name{monteCarloCI}
\alias{monteCarloCI}
\alias{monteCarloMed}
\title{Monte Carlo Confidence Intervals to Test Functions of Parameter Estimates}
\usage{
monteCarloCI(object = NULL, expr, coefs, ACM, nRep = 20000,
  standardized = FALSE, fast = TRUE, level = 0.95, na.rm = TRUE,
  append.samples = FALSE, plot = FALSE,
  ask = getOption("device.ask.default"), ...)
}
\arguments{
\item{object}{A object of class \link[lavaan:lavaan-class]{lavaan::lavaan} in which
functions of parameters have already been defined using the \verb{:=}
operator in \code{lavaan}'s \code{\link[lavaan:model.syntax]{lavaan::model.syntax()}}. When
\code{NULL}, users must specify \code{expr}, \code{coefs}, and \code{ACM}.}

\item{expr}{Optional \code{character} vector specifying functions of model
parameters (e.g., an indirect effect). Ideally, the vector should have
names, which is necessary if any user-defined parameters refer to other
user-defined parameters defined earlier in the vector (order matters!).
All parameters appearing in the vector must be provided in \code{coefs},
or defined (as functions of \code{coefs}) earlier in \code{expr}. If
\code{length(expr) > 1L}, \code{nRep} samples will be drawn
simultaneously from a single multivariate distribution; thus,
\code{ACM} must include all parameters in \code{coefs}.}

\item{coefs}{\code{numeric} vector of parameter estimates used in
\code{expr}. Ignored when \code{object} is used.}

\item{ACM}{Symmetric \code{matrix} representing the asymptotic sampling
covariance matrix (ACOV) of the parameter estimates in \code{coefs}.
Ignored when \code{object} is used. Information on how to obtain the ACOV
in popular SEM software is described in \strong{Details}.}

\item{nRep}{\code{integer}. The number of samples to draw, to obtain an
empirical sampling distribution of model parameters. Many thousand are
recommended to minimize Monte Carlo error of the estimated CIs.}

\item{standardized}{\code{logical} indicating whether to obtain CIs for the
fully standardized (\code{"std.all"}) estimates, using their asymptotic
sampling covariance matrix.}

\item{fast}{\code{logical} indicating whether to use a fast algorithm that
assumes all functions of parameters (in \code{object} or \code{expr}) use
standard operations. Set to \code{FALSE} if using (e.g.) \code{\link[=c]{c()}}
to concatenate parameters in the definition, which would have unintended
consequences when vectorizing functions in \code{expr} across sampled
parameters.}

\item{level}{\code{numeric} confidence level, between 0--1}

\item{na.rm}{\code{logical} passed to \code{\link[stats:quantile]{stats::quantile()}}}

\item{append.samples}{\code{logical} indicating whether to return the
simulated empirical sampling distribution of parameters (in \code{coefs})
and functions (in \code{expr}) in a \code{list} with the results. This
could be useful to calculate more precise highest-density intervals (see
examples).}

\item{plot}{\code{logical} indicating whether to plot the empirical sampling
distribution of each function in \code{expr}}

\item{ask}{whether to prompt user before printing each plot}

\item{\dots}{arguments passed to \code{\link[graphics:hist]{graphics::hist()}} when
\code{plot = TRUE}.}
}
\value{
A \code{lavaan.data.frame} (to use lavaan's \code{print} method)
with point estimates and confidence limits of each requested function of
parameters in \code{expr} is returned. If \code{append.samples = TRUE},
output will be a \code{list} with the same \verb{$Results} along with a
second \code{data.frame} with the \verb{$Samples} (in rows) of each
parameter (in columns), and an additional column for each requested
function of those parameters.
}
\description{
Robust confidence intervals for functions of parameter estimates,
based on empirical sampling distributions of estimated model parameters.
}
\details{
This function implements the Monte Carlo method of obtaining an empirical
sampling distribution of estimated model parameters, as described by
MacKinnon et al. (2004) for testing indirect effects in mediation models.
This is essentially a parametric bootstrap method, which (re)samples
parameters (rather than raw data) from a multivariate-normal distribution
with mean vector equal to estimates in \code{coef()} and covariance matrix
equal to the asymptotic covariance matrix \code{vcov()} of estimated parameters.

The easiest way to use the function is to fit a SEM to data with
\code{\link[lavaan:lavaan]{lavaan::lavaan()}}, using the \verb{:=} operator in the
\code{\link[lavaan:model.syntax]{lavaan::model.syntax()}} to specify user-defined parameters.
All information is then available in the resulting
\link[lavaan:lavaan-class]{lavaan::lavaan} object.  Alternatively (especially when using
external SEM software to fit the model), the expression(s) can be explicitly
passed to the function, along with the vector of estimated model parameters
and their associated asymptotic sampling covariance matrix (ACOV).
For further information on the Monte Carlo method, see MacKinnon et al.
(2004) and Preacher & Selig (2012).

The asymptotic covariance matrix can be obtained easily from many popular
SEM software packages.
\itemize{
\item{LISREL: Including the EC option on the OU line will print the ACM
to a seperate file. The file contains the lower triangular elements of
the ACM in free format and scientific notation.}
\item{M\emph{plus}: Include the command TECH3; in the OUTPUT section.
The ACM will be printed in the output.}
\item{\code{lavaan}: Use the \code{\link[=vcov]{vcov()}} method on the fitted \link[lavaan:lavaan-class]{lavaan::lavaan}
object to return the ACM.}
}
}
\examples{

## From the mediation tutorial:
## http://lavaan.ugent.be/tutorial/mediation.html

set.seed(1234)
X <- rnorm(100)
M <- 0.5*X + rnorm(100)
Y <- 0.7*M + rnorm(100)
dat <- data.frame(X = X, Y = Y, M = M)

mod <- ' # direct effect
  Y ~ c*X
  # mediator
  M ~ a*X
  Y ~ b*M
  # indirect effect (a*b)
  ind := a*b
  # total effect
  total := ind + c
'
fit <- sem(mod, data = dat)
summary(fit, ci = TRUE) # print delta-method CIs

## Automatically extract information from lavaan object
set.seed(1234)
monteCarloCI(fit) # CIs more robust than delta method in smaller samples

## delta method for standardized solution
standardizedSolution(fit)
## compare to Monte Carlo CIs:
set.seed(1234)
monteCarloCI(fit, standardized = TRUE)

\donttest{
## save samples to calculate more precise intervals:
set.seed(1234)
foo <- monteCarloCI(fit, append.samples = TRUE)
# library(HDInterval) # not a dependency; must be installed
# hdi(foo$Samples)
}
## Parameters can also be obtained from an external analysis
myParams <- c("a","b","c")
(coefs <- coef(fit)[myParams]) # names must match those in the "expression"
## Asymptotic covariance matrix from an external analysis
(AsyCovMat <- vcov(fit)[myParams, myParams])
## Compute CI, include a plot
set.seed(1234)
monteCarloCI(expr = c(ind = 'a*b', total = 'ind + c',
                      ## other arbitrary functions are also possible
                      meaningless = 'sqrt(a)^b / log(abs(c))'),
             coefs = coefs, ACM = AsyCovMat,
             plot = TRUE, ask = TRUE) # print a plot for each

}
\references{
MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence limits
for the indirect effect: Distribution of the product and resampling methods.
\emph{Multivariate Behavioral Research, 39}(1) 99--128.
\doi{10.1207/s15327906mbr3901_4}

Preacher, K. J., & Selig, J. P. (2010, July). Monte Carlo method
for assessing multilevel mediation: An interactive tool for creating
confidence intervals for indirect effects in 1-1-1 multilevel models.
Computer software available from \url{http://quantpsy.org/}.

Preacher, K. J., & Selig, J. P. (2012). Advantages of Monte Carlo confidence
intervals for indirect effects. \emph{Communication Methods and Measures, 6}(2),
77--98. \doi{10.1080/19312458.2012.679848}

Selig, J. P., & Preacher, K. J. (2008, June). Monte Carlo method for
assessing mediation: An interactive tool for creating confidence intervals
for indirect effects. Computer software available from
\url{http://quantpsy.org/}.
}
\author{
Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
