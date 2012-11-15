\name{reliability}
\alias{reliability}
\title{
Calculate reliability values of factors
}
\description{
Calculate reliability values of factors by coefficient omega
}
\usage{
reliability(object)
}
\arguments{
  \item{object}{The lavaan model object provided after running the \code{cfa}, \code{sem}, \code{growth}, or \code{lavaan} functions.}
}
\details{
The coefficient omega (Raykov, 2001; Equation 5) can be calculated by 

\deqn{ \omega =\frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2} Var\left( \psi \right)}{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2} Var\left( \psi \right) + \sum^{k}_{i = 1} \theta_{ii} + 2\sum_{i < j} \theta_{ij} }, }

where \eqn{k} is the number of items in a factor, \eqn{\lambda_i} is the factor loading of item \emph{i}, \eqn{\psi} is the factor variance, \eqn{\theta_{ii}} is the variance of measurement errors of item \emph{i}, and \eqn{\theta_{ij}} is the covariance of measurement errors from item \emph{i} and \emph{j}. If measurement errors are not correlated, the last term in the denominator is dropped. For the categorical items, the reliability is calculated from the latent variable underlying the categorical items.
}
\value{
	Reliability values (coefficient omega) of each factor
}
\references{
Raykov, T. (2001). Estimation of congeneric scale reliability using covariance structure analysis with nonlinear constraints \emph{British Journal of Mathematical and Statistical Psychology, 54}, 315-323.
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
reliability(fit)
}
