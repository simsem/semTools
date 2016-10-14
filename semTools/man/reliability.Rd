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

\deqn{ \omega_1 =\frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2} Var\left( \psi \right)}{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2} Var\left( \psi \right) + \sum^{k}_{i = 1} \theta_{ii} + 2\sum_{i < j} \theta_{ij} }, }

where \eqn{\lambda_i} is the factor loading of item \emph{i}, \eqn{\psi} is the factor variance, \eqn{\theta_{ii}} is the variance of measurement errors of item \emph{i}, and \eqn{\theta_{ij}} is the covariance of measurement errors from item \emph{i} and \emph{j}. 

The second coefficient omega (Bentler, 1972, 2009) can be calculated by 

\deqn{ \omega_2 = \frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2} Var\left( \psi \right)}{\bold{1}^\prime \hat{\Sigma} \bold{1}}, }

where \eqn{\hat{\Sigma}} is the model-implied covariance matrix, and \eqn{\bold{1}} is the \eqn{k}-dimensional vector of 1. The first and the second coefficients omega will have different values if there are dual loadings (or the existence of method factors). The first coefficient omega can be viewed as the reliability controlling for the other factors (like partial eta-squared in ANOVA). The second coefficient omega can be viewed as the unconditional reliability (like eta-squared in ANOVA). 

The third coefficient omega (McDonald, 1999), which is sometimes referred to hierarchical omega, can be calculated by

\deqn{ \omega_3 =\frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2} Var\left( \psi \right)}{\bold{1}^\prime \Sigma \bold{1}}, }

where \eqn{\Sigma} is the observed covariance matrix. If the model fits the data well, the third coefficient omega will be similar to the \eqn{\omega_2}. Note that if there is a directional effect in the model, all coefficients omega will use the total factor variances, which is calculated by \code{\link[lavaan]{lavInspect}(object, "cov.lv")}.

In conclusion, \eqn{\omega_1}, \eqn{\omega_2}, and \eqn{\omega_3} are different in the denominator. The denominator of the first formula assumes that a model is congeneric factor model where measurement errors are not correlated. The second formula is accounted for correlated measurement errors. However, these two formulas assume that the model-implied covariance matrix explains item relationships perfectly. The residuals are subject to sampling error. The third formula use observed covariance matrix instead of model-implied covariance matrix to calculate the observed total variance. This formula is the most conservative method in calculating coefficient omega.

The average variance extracted (AVE) can be calculated by 

\deqn{ AVE = \frac{\bold{1}^\prime \textrm{diag}\left(\Lambda\Psi\Lambda^\prime\right)\bold{1}}{\bold{1}^\prime \textrm{diag}\left(\hat{\Sigma}\right) \bold{1}}, }

Note that this formula is modified from Fornell & Larcker (1981) in the case that factor variances are not 1. The proposed formula from Fornell & Larcker (1981) assumes that the factor variances are 1. Note that AVE will not be provided for factors consisting of items with dual loadings. AVE is the property of items but not the property of factors.

Regarding to categorical items, coefficient alpha and AVE are calculated based on polychoric correlations. The coefficient alpha from this function may be not the same as the standard alpha calculation for categorical items. Researchers may check the \code{alpha} function in the \code{psych} package for the standard coefficient alpha calculation. 

Item thresholds are not accounted for. Coefficient omega for categorical items, however, is calculated by accounting for both item covariances and item thresholds using Green and Yang's (2009, formula 21) approach. Three types of coefficient omega indicate different methods to calculate item total variances. The original formula from Green and Yang is equivalent to \eqn{\omega_3} in this function.
}
\value{
	Reliability values (coefficient alpha, coefficients omega, average variance extracted) of each factor in each group
}
\references{
Bentler, P. M. (1972). A lower-bound method for the dimension-free measurement of internal consistency. \emph{Social Science Research, 1}, 343-357.

Bentler, P. M. (2009). Alpha, dimension-free, and model-based internal consistency reliability. \emph{Psychometrika, 74}, 137-143.

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of tests. \emph{Psychometrika, 16}, 297-334. 

Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models with unobservable variables and measurement errors. \emph{Journal of Marketing Research, 18}, 39-50.

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using structural equation modeling: An alternative to coefficient alpha. \emph{Psychometrika, 74}, 155-167.

McDonald, R. P. (1999). Test theory: A unified treatment. Mahwah, NJ: Erlbaum.

Raykov, T. (2001). Estimation of congeneric scale reliability using covariance structure analysis with nonlinear constraints \emph{British Journal of Mathematical and Statistical Psychology, 54}, 315-323.
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com}); Yves Rosseel (Ghent University; \email{Yves.Rosseel@UGent.be})
}
\seealso{
	\code{\link{reliabilityL2}} for reliability value of a desired second-order factor, \code{\link{maximalRelia}} for the maximal reliability of weighted composite
}
\examples{
library(lavaan)

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
reliability(fit)
}
