% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reliability.R
\name{reliability-deprecated}
\alias{reliability-deprecated}
\title{Composite Reliability using SEM}
\usage{
reliability(object, what = c("alpha", "omega", "omega2", "omega3", "ave"),
            return.total = FALSE, dropSingle = TRUE, omit.factors = character(0),
            omit.indicators = character(0), omit.imps = c("no.conv", "no.se"))
}
\arguments{
\item{object}{A \link[lavaan:lavaan-class]{lavaan::lavaan} or \link[lavaan.mi:lavaan.mi-class]{lavaan.mi::lavaan.mi} object,
expected to contain only exogenous common factors (i.e., a CFA model).}

\item{what}{\code{character} vector naming any reliability indices to
calculate. All are returned by default. When indicators are ordinal,
both traditional \code{"alpha"} and Zumbo et al.'s (2007) so-called
"ordinal alpha" (\code{"alpha.ord"}) are returned, though the latter is
arguably of dubious value (Chalmers, 2018).}

\item{return.total}{\code{logical} indicating whether to return a final
column containing the reliability of a composite of all indicators (not
listed in \code{omit.indicators}) of factors not listed in
\code{omit.factors}.  Ignored in 1-factor models, and should only be set
\code{TRUE} if all factors represent scale dimensions that could be
meaningfully collapsed to a single composite (scale sum or scale mean).}

\item{dropSingle}{\code{logical} indicating whether to exclude factors
defined by a single indicator from the returned results. If \code{TRUE}
(default), single indicators will still be included in the \code{total}
column when \code{return.total = TRUE}.}

\item{omit.factors}{\code{character} vector naming any common factors
modeled in \code{object} whose composite reliability is not of
interest. For example, higher-order or method factors. Note that
\code{\link[=reliabilityL2]{reliabilityL2()}} should be used to calculate composite
reliability of a higher-order factor.}

\item{omit.indicators}{\code{character} vector naming any observed variables
that should be ignored when calculating composite reliability. This can
be useful, for example, to estimate reliability when an indicator is
removed.}

\item{omit.imps}{\code{character} vector specifying criteria for omitting
imputations from pooled results.  Can include any of
\code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
default setting, which excludes any imputations that did not
converge or for which standard errors could not be computed.  The
last option (\code{"no.npd"}) would exclude any imputations which
yielded a nonpositive definite covariance matrix for observed or
latent variables, which would include any "improper solutions" such
as Heywood cases.  NPD solutions are not excluded by default because
they are likely to occur due to sampling error, especially in small
samples.  However, gross model misspecification could also cause
NPD solutions, users can compare pooled results with and without
this setting as a sensitivity analysis to see whether some
imputations warrant further investigation.}
}
\value{
Reliability values (coefficient alpha, coefficients omega, average
variance extracted) of each factor in each group. If there are multiple
factors, a \code{total} column can optionally be included.
}
\description{
Calculate composite reliability from estimated factor-model parameters
}
\details{
The coefficient alpha (Cronbach, 1951) can be calculated by

\deqn{ \alpha = \frac{k}{k - 1}\left[ 1 - \frac{\sum^{k}_{i = 1}
\sigma_{ii}}{\sum^{k}_{i = 1} \sigma_{ii} + 2\sum_{i < j} \sigma_{ij}}
\right],}

where \eqn{k} is the number of items in a factor, \eqn{\sigma_{ii}} is the
item \emph{i} observed variances, \eqn{\sigma_{ij}} is the observed
covariance of items \emph{i} and \emph{j}.

Several coefficients for factor-analysis reliability have been termed
"omega", which Cho (2021) argues is a misleading misnomer and argues for
using \eqn{\rho} to represent them all, differentiated by descriptive
subscripts.  In our package, we number \eqn{\omega} based on commonly
applied calculations.  Bentler (1968) first introduced factor-analysis
reliability for a unidimensional factor model with congeneric indicators.
However, assuming there are no cross-loadings in a multidimensional CFA,
this reliability coefficient can be calculated for each factor in the model.

\deqn{ \omega_1 =\frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
Var\left( \psi \right)}{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
Var\left( \psi \right) + \sum^{k}_{i = 1} \theta_{ii} + 2\sum_{i < j}
\theta_{ij} }, }

where \eqn{\lambda_i} is the factor loading of item \emph{i}, \eqn{\psi} is
the factor variance, \eqn{\theta_{ii}} is the variance of measurement errors
of item \emph{i}, and \eqn{\theta_{ij}} is the covariance of measurement
errors from item \emph{i} and \emph{j}. McDonald (1999) later referred to
this \emph{and other reliability coefficients} as "omega", which is a source
of confusion when reporting coefficients (Cho, 2021).

The additional coefficients generalize the first formula by accounting for
multidimenisionality (possibly with cross-loadings) and correlated errors.
By setting \code{return.total=TRUE}, one can estimate reliability for a
single composite calculated using all indicators in the multidimensional
CFA (Bentler, 1972, 2009).  \code{"omega2"} is calculated by

\deqn{ \omega_2 = \frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
Var\left( \psi \right)}{\bold{1}^\prime \hat{\Sigma} \bold{1}}, }

where \eqn{\hat{\Sigma}} is the model-implied covariance matrix, and
\eqn{\bold{1}} is the \eqn{k}-dimensional vector of 1. The first and the
second coefficients omega will have the same value per factor in models with
simple structure, but they differ when there are (e.g.) cross-loadings
or method factors. The first coefficient omega can be viewed as the
reliability controlling for the other factors (like \eqn{\eta^2_{partial}} in
ANOVA). The second coefficient omega can be viewed as the unconditional
reliability (like \eqn{\eta^2} in ANOVA).

The \code{"omega3"} coefficient (McDonald, 1999), sometimes referred to as
hierarchical omega, can be calculated by

\deqn{ \omega_3 =\frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
Var\left( \psi \right)}{\bold{1}^\prime \Sigma \bold{1}}, }

where \eqn{\Sigma} is the observed covariance matrix. If the model fits the
data well, \eqn{\omega_3} will be similar to \eqn{\omega_2}. Note that if
there is a directional effect in the model, all coefficients are calculated
from total factor variances: \code{lavInspect(object, "cov.lv")}.

In conclusion, \eqn{\omega_1}, \eqn{\omega_2}, and \eqn{\omega_3} are
different in the denominator. The denominator of the first formula assumes
that a model is congeneric factor model where measurement errors are not
correlated. The second formula accounts for correlated measurement errors.
However, these two formulas assume that the model-implied covariance matrix
explains item relationships perfectly. The residuals are subject to sampling
error. The third formula use observed covariance matrix instead of
model-implied covariance matrix to calculate the observed total variance.
This formula is the most conservative method in calculating coefficient
omega.

The average variance extracted (AVE) can be calculated by

\deqn{ AVE = \frac{\bold{1}^\prime
\textrm{diag}\left(\Lambda\Psi\Lambda^\prime\right)\bold{1}}{\bold{1}^\prime
\textrm{diag}\left(\hat{\Sigma}\right) \bold{1}}, }

Note that this formula is modified from Fornell & Larcker (1981) in the case
that factor variances are not 1. The proposed formula from Fornell & Larcker
(1981) assumes that the factor variances are 1. Note that AVE will not be
provided for factors consisting of items with dual loadings. AVE is the
property of items but not the property of factors. AVE is calculated with
polychoric correlations when ordinal indicators are used.

Coefficient alpha is by definition applied by treating indicators as numeric
(see Chalmers, 2018), which is consistent with the \code{alpha} function in
the \code{psych} package. When indicators are ordinal, \code{reliability}
additionally applies the standard alpha calculation to the polychoric
correlation matrix to return Zumbo et al.'s (2007) "ordinal alpha".

Coefficient omega for categorical items is calculated using Green and Yang's
(2009, formula 21) approach. Three types of coefficient omega indicate
different methods to calculate item total variances. The original formula
from Green and Yang is equivalent to \eqn{\omega_3} in this function.
Green and Yang did not propose a method for
calculating reliability with a mixture of categorical and continuous
indicators, and we are currently unaware of an appropriate method.
Therefore, when \code{reliability} detects both categorical and continuous
indicators of a factor, an error is returned. If the categorical indicators
load on a different factor(s) than continuous indicators, then reliability
will still be calculated separately for those factors, but
\code{return.total} must be \code{FALSE} (unless \code{omit.factors} is used
to isolate factors with indicators of the same type).
}
\examples{

data(HolzingerSwineford1939)
HS9 <- HolzingerSwineford1939[ , c("x7","x8","x9")]
HSbinary <- as.data.frame( lapply(HS9, cut, 2, labels=FALSE) )
names(HSbinary) <- c("y7","y8","y9")
HS <- cbind(HolzingerSwineford1939, HSbinary)

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ y7 + y8 + y9 '

fit <- cfa(HS.model, data = HS, ordered = c("y7","y8","y9"), std.lv = TRUE)

## works for factors with exclusively continuous OR categorical indicators
reliability(fit)

## reliability for ALL indicators only available when they are
## all continuous or all categorical
reliability(fit, omit.factors = "speed", return.total = TRUE)


## loop over visual indicators to calculate alpha if one indicator is removed
for (i in paste0("x", 1:3)) {
  cat("Drop x", i, ":\n")
  print(reliability(fit, omit.factors = c("textual","speed"),
                    omit.indicators = i, what = "alpha"))
}


## works for multigroup models and for multilevel models (and both)
data(Demo.twolevel)
## assign clusters to arbitrary groups
Demo.twolevel$g <- ifelse(Demo.twolevel$cluster \%\% 2L, "type1", "type2")
model2 <- ' group: type1
  level: within
    fac =~ y1 + L2*y2 + L3*y3
  level: between
    fac =~ y1 + L2*y2 + L3*y3

group: type2
  level: within
    fac =~ y1 + L2*y2 + L3*y3
  level: between
    fac =~ y1 + L2*y2 + L3*y3
'
fit2 <- sem(model2, data = Demo.twolevel, cluster = "cluster", group = "g")
reliability(fit2, what = c("alpha","omega3"))

}
\references{
Bentler, P. M. (1972). A lower-bound method for the dimension-free
measurement of internal consistency. \emph{Social Science Research, 1}(4),
343--357. \doi{10.1016/0049-089X(72)90082-8}

Bentler, P. M. (2009). Alpha, dimension-free, and model-based internal
consistency reliability. \emph{Psychometrika, 74}(1), 137--143.
\doi{10.1007/s11336-008-9100-1}

Chalmers, R. P. (2018). On misconceptions and the limited usefulness of
ordinal alpha. \emph{Educational and Psychological Measurement, 78}(6),
1056--1071. \doi{10.1177/0013164417727036}

Cho, E. (2021) Neither Cronbach’s alpha nor McDonald’s omega: A commentary
on Sijtsma and Pfadt. \emph{Psychometrika, 86}(4), 877--886.
\doi{10.1007/s11336-021-09801-1}

Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
tests. \emph{Psychometrika, 16}(3), 297--334. \doi{10.1007/BF02310555}

Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models
with unobservable variables and measurement errors. \emph{Journal of
Marketing Research, 18}(1), 39--50. \doi{10.2307/3151312}

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
structural equation modeling: An alternative to coefficient alpha.
\emph{Psychometrika, 74}(1), 155--167. \doi{10.1007/s11336-008-9099-3}

McDonald, R. P. (1999). \emph{Test theory: A unified treatment}. Mahwah, NJ:
Erlbaum.

Raykov, T. (2001). Estimation of congeneric scale reliability using
covariance structure analysis with nonlinear constraints \emph{British
Journal of Mathematical and Statistical Psychology, 54}(2), 315--323.
\doi{10.1348/000711001159582}

Zumbo, B. D., Gadermann, A. M., & Zeisser, C. (2007). Ordinal versions of
coefficients alpha and theta for Likert rating scales.
\emph{Journal of Modern Applied Statistical Methods, 6}(1), 21--29.
\doi{10.22237/jmasm/1177992180}

Zumbo, B. D., & Kroc, E. (2019). A measurement is a choice and Stevens’
scales of measurement do not help make it: A response to Chalmers.
\emph{Educational and Psychological Measurement, 79}(6), 1184--1197.
\doi{10.1177/0013164419844305}
}
\seealso{
\code{\link[=semTools-deprecated]{semTools-deprecated()}}
}
\author{
Sunthud Pornprasertmanit (\email{psunthud@gmail.com})

Yves Rosseel (Ghent University; \email{Yves.Rosseel@UGent.be})

Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\keyword{internal}
