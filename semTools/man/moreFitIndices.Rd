\name{moreFitIndices}
\alias{moreFitIndices}
\alias{moreFitIndices}
\title{
Calculate more fit indices
}
\description{
Calculate more fit indices that are not already provided in lavaan. 
}
\usage{
moreFitIndices(object, nPrior = 1)
}
\arguments{
  \item{object}{The lavaan model object provided after running the \code{cfa} or the \code{sem} functions.}
  \item{nPrior}{The sample size on which prior is based. This argument is used to compute BIC*.}
}
\details{
Normed Fit Index (nfi; West, Taylor, & Wu, 2012) is one of the relative fit indices which can be computed by

\deqn{ nfi =\frac{\chi^{2}_{0} - \chi^{2}_{k}}{\chi^{2}_{0}},}

where \eqn{\chi^{2}_{k}} is the chi-square test statistic value of the target model, \eqn{\chi^{2}_{0}} is the chi-square test statistic value of the null model.  

Incremental Fit Index (ifi; West, Taylor, & Wu, 2012) is one of the relative fit indices which can be computed by

\deqn{ ifi =\frac{\chi^{2}_{0} - \chi^{2}_{k}}{\chi^{2}_{0} - df_{k}},}

where \eqn{df_{k}} is the degree of freedom when fitting the target model  

Gamma Hat (gfi*; West, Taylor, & Wu, 2012) is one of the relative fit indices which can be computed by

\deqn{ gfi* =\frac{p}{p + 2 \times \frac{\chi^{2}_{k} - df_{k}}{N - 1}},}

where \eqn{N} is the sample size, \eqn{p} is the number of variables in the model. This formula assumes equal number of indicators across groups.

Adjusted Gamma Hat (agfi*; West, Taylor, & Wu, 2012) is one of the relative fit indices which can be computed by

\deqn{ agfi* = \left(1 - \frac{K \times p \times (p + 1)}{2 \times df_{k}} \right) \times \left( 1 - gfi* \right) ,}

where \eqn{K} is the number of group (please refer to Dudgeon, 2004 for the multiple-group adjustment for agfi*).

Corrected Akaike Information Criterion (AICc; Burnham & Anderson, 2003) is the corrected version of aic for small sample size:

\deqn{ aicc = f + \frac{2k(k + 1)}{N - k - 1},}

where \eqn{f} is the minimized discrepancy function, which is the product of the log likelihood and -2, and \eqn{k} is the number of parameters in the target model.

Expected Value of Cross-Validation Index (ECVI; West, Taylor, & Wu, 2012) is the average discrepancy in the fitted covariance matrices between two samples of equal sample size across all possible combinations of two samples from the same population:

\deqn{ ecvi = f + \frac{2 \times k}{N},}

Stochastic information criterion (sic; Preacher, 2006) is similar to aic or bic. This index will account for model complexity in the model's function form, in addition to the number of free parameters. sic can be computed by

\deqn{ sic = \frac{1}{2}\left(f - \log{\det{I(\hat{\theta})}}\right),}

where \eqn{I(\hat{\theta})} is the information matrix of the parameters.

Corrected Bayesian Information Criterion (BIC*; Kuha, 2004) is similar to bic but explicitly specifying the sample size on which the prior is based (\eqn{N_{prior}}).

\deqn{ bicc = f + k\log{(1 + N/N_{prior})},}

Hannan-Quinn Information Criterion (hqc; Hannan & Quinn, 1979) is used for model selection similar to aic or bic.

\deqn{ hqc = f + 2k\log{(\log{N})},}
}

\value{
\enumerate{
  \item{nfi} Normed Fit Index
  \item{ifi} Incremental Fit Index
  \item{gfi*} Gamma Hat
  \item{agfi*} Adjusted Gamma Hat
  \item{aicc} Corrected Akaike Information Criterion
  \item{ecvi} Expected Value of Cross-Validation Index 
  \item{sic} Stochastic Information Criterion
  \item{bic*} Bayesian Information Criterion with specifying the prior sample size
  \item{hqc} Hannan-Quinn Information Criterion
}
}
\references{
Burnham, K., & Anderson, D. (2003). \emph{Model selection and multimodel inference: A practical-theoretic approach.} New York, NY: Springer-Verlag.

Dudgeon, P. (2004). A note on extending Steiger's (1998) multiple sample RMSEA adjustment to other noncentrality parameter-based statistic. \emph{Structural Equation Modeling, 11}, 305-319.

Kuha, J. (2004). AIC and BIC: Comparisons of assumptions and performance. \emph{Sociological Methods Research, 33,} 188-229.

Preacher, K. J. (2006). Quantifying parsimony in structural equation modeling. \emph{Multivariate Behavioral Research, 43,} 227-259.

West, S. G., Taylor, A. B., & Wu, W. (2012). Model fit and model selection in structural equation modeling. In R. H. Hoyle (Ed.), \emph{Handbook of Structural Equation Modeling.} New York: Guilford.
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
	Aaron Boulton (University of Kansas; \email{aboulton@ku.edu})
}
\seealso{
	\code{\link{miPowerFit}} For the modification indices and their power approach for model fit evaluation
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
moreFitIndices(fit)
}
