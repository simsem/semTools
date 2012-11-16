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
The coefficient alpha (Cronbach, 1951) can be calculated by 

\deqn{ \alpha = \frac{k}{k - 1}\left[ 1 - \frac{\sum^{k}_{i = 1} \sigma_{ii}}{\sum^{k}_{i = 1} \sigma_{ii} + 2\sum_{i < j} \sigma_{ij}} \right],}

where \eqn{k} is the number of items in a factor, \eqn{\sigma_{ii}} is the item \emph{i} observed variances, \eqn{\sigma_{ij}} is the observed covariance of items \emph{i} and \emph{j}.

The coefficient omega (Raykov, 2001) can be calculated by 

\deqn{ \omega =\frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2} Var\left( \psi \right)}{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2} Var\left( \psi \right) + \sum^{k}_{i = 1} \theta_{ii} + 2\sum_{i < j} \theta_{ij} }, }

where \eqn{\lambda_i} is the factor loading of item \emph{i}, \eqn{\psi} is the factor variance, \eqn{\theta_{ii}} is the variance of measurement errors of item \emph{i}, and \eqn{\theta_{ij}} is the covariance of measurement errors from item \emph{i} and \emph{j}. 

The second coefficient omega (Bentler, 1972, 2009) can be calculated by 

\deqn{ \omega_2 = 1 - \frac{\bold{1}\prime \Theta \bold{1}}{\bold{1}^\prime \hat{\Sigma} \bold{1}}, }

where \eqn{\Theta} is the measurement error covariance matrix, \eqn{\hat{\Sigma}} is the model-implied covariance matrix, and \eqn{\bold{1}} is the \eqn{k}-dimensional vector of 1. The first and the second coefficients omega will have different values if there are dual loadings (or the existence of method factors). The first coefficient omega can be viewed as the reliability controlling for the other factors. The second coefficient omega can be viewed as the unconditional reliability. 

The third coefficient omega (McDonald, 1999) can be calculated by

\deqn{ \omega_3 = 1 - \frac{\bold{1}\prime \Theta \bold{1}}{\bold{1}^\prime \Sigma \bold{1}}, }

where \eqn{\Sigma} is the observed covariance matrix. If the model fits the data well, the third coefficient omega will be similar to the other two. Note that if there is a directional effect in the model, all coefficients omega will use the total factor variances, which is calculated by the \code{\link{impliedFactorCov}} function.

If measurement errors are not correlated, the last term in the denominator of the middle expression is dropped. For the categorical items, the reliability indices are calculated from the latent variable underlying the categorical items (i.e., using polychoric/polyserial correlations). Therefore, the coefficient alpha from this function may be not the same as the standard alpha calculation for categorical items. Researchers may check the \code{alpha} function in the \code{psych} package for the standard coefficient alpha calculation. 
}
\value{
	Reliability values (coefficient alpha, coefficients omega) of each factor in each group
}
\references{
Bentler, P. M. (1972). A lower-bound method for the dimension-free measurement of internal consistency. \emph{Social Science Research, 1}, 343-357.

Bentler, P. M. (2009). Alpha, dimension-free, and model-based internal consistency reliability. \emph{Psychometrika, 74}, 137-143.

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. \emph{Psychometrika, 16}, 297-334. 

McDonald, R. P. (1999). Test theory: A unified treatment. Mahwah, NJ: Erlbaum.

Raykov, T. (2001). Estimation of congeneric scale reliability using covariance structure analysis with nonlinear constraints \emph{British Journal of Mathematical and Statistical Psychology, 54}, 315-323.
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu}); Yves Rosseel (Ghent University; \email{Yves.Rosseel@UGent.be})
}
\seealso{
	\code{\link{reliabilityL2}} for reliability value of a desired second-order factor
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
reliability(fit)
}
