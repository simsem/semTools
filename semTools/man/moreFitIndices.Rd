\name{moreFitIndices}
\alias{moreFitIndices}
\title{
Calculate more fit indices
}
\description{
Calculate more fit indices that are not already provided in lavaan.
}
\usage{
moreFitIndices(object, fit.measures = "all", nPrior = 1)
}
\arguments{
  \item{object}{The lavaan model object provided after running the \code{cfa}, \code{sem}, \code{growth}, or \code{lavaan} functions.}
  \item{fit.measures}{Additional fit measures to be calculated. All additional fit measures are calculated by default}
  \item{nPrior}{The sample size on which prior is based. This argument is used to compute BIC*.}
}
\details{
Gamma Hat (gammaHat; West, Taylor, & Wu, 2012) is a global fit index which can be computed by

\deqn{ gammaHat =\frac{p}{p + 2 \times \frac{\chi^{2}_{k} - df_{k}}{N - 1}},}

where \eqn{p} is the number of variables in the model, \eqn{\chi^{2}_{k}} is the chi-square test statistic value of the target model, \eqn{df_{k}} is the degree of freedom when fitting the target model, and \eqn{N} is the sample size. This formula assumes equal number of indicators across groups.

Adjusted Gamma Hat (adjGammaHat; West, Taylor, & Wu, 2012) is a global fit index which can be computed by

\deqn{ adjGammaHat = \left(1 - \frac{K \times p \times (p + 1)}{2 \times df_{k}} \right) \times \left( 1 - gammaHat \right) ,}

where \eqn{K} is the number of groups (please refer to Dudgeon, 2004 for the multiple-group adjustment for agfi*).

Corrected Akaike Information Criterion (aic.smallN; Burnham & Anderson, 2003) is the corrected version of aic for small sample size:

\deqn{ aic.smallN = f + \frac{2k(k + 1)}{N - k - 1},}

where \eqn{f} is the minimized discrepancy function, which is the product of the log likelihood and -2, and \eqn{k} is the number of parameters in the target model.

Corrected Bayesian Information Criterion (bic.priorN; Kuha, 2004) is similar to bic but explicitly specifying the sample size on which the prior is based (\eqn{N_{prior}}).

\deqn{ bic.priorN = f + k\log{(1 + N/N_{prior})},}

Stochastic information criterion (sic; Preacher, 2006) is similar to aic or bic. This index will account for model complexity in the model's function form, in addition to the number of free parameters. This index will be provided only when the chi-squared value is not scaled. The sic can be computed by

\deqn{ sic = \frac{1}{2}\left(f - \log{\det{I(\hat{\theta})}}\right),}

where \eqn{I(\hat{\theta})} is the information matrix of the parameters.

Hannan-Quinn Information Criterion (hqc; Hannan & Quinn, 1979) is used for model selection similar to aic or bic.

\deqn{ hqc = f + 2k\log{(\log{N})},}

Note that if Satorra-Bentler or Yuan-Bentler's method is used, the fit indices using the scaled chi-square values are also provided. 

See \code{\link{nullRMSEA}} for the further details of the computation of RMSEA of the null model.
}

\value{
\enumerate{
  \item{gammaHat} Gamma Hat
  \item{adjGammaHat} Adjusted Gamma Hat
  \item{baseline.rmsea} RMSEA of the Baseline (Null) Model
  \item{aic.smallN} Corrected (for small sample size) Akaike Information Criterion
  \item{bic.priorN} Bayesian Information Criterion with specifying the prior sample size
  \item{sic} Stochastic Information Criterion
  \item{hqc} Hannan-Quinn Information Criterion
  \item{gammaHat.scaled} Gamma Hat using Scaled Chi-square
  \item{adjGammaHat.scaled} Adjusted Gamma Hat using Scaled Chi-square
  \item{baseline.rmsea.scaled} RMSEA of the Baseline (Null) Model using Scaled Chi-square
}
}
\references{
Burnham, K., & Anderson, D. (2003). \emph{Model selection and multimodel inference: A practical-theoretic approach.} New York, NY: Springer-Verlag.

Dudgeon, P. (2004). A note on extending Steiger's (1998) multiple sample RMSEA adjustment to other noncentrality parameter-based statistic. \emph{Structural Equation Modeling, 11}, 305-319.

Kuha, J. (2004). AIC and BIC: Comparisons of assumptions and performance. \emph{Sociological Methods Research, 33}, 188-229.

Preacher, K. J. (2006). Quantifying parsimony in structural equation modeling. \emph{Multivariate Behavioral Research, 43}, 227-259.

West, S. G., Taylor, A. B., & Wu, W. (2012). Model fit and model selection in structural equation modeling. In R. H. Hoyle (Ed.), \emph{Handbook of Structural Equation Modeling.} New York: Guilford.
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
    Terrence Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
    Aaron Boulton (University of North Carolina, Chapel Hill; \email{aboulton@email.unc.edu})
    Ruben Arslan (Humboldt-University of Berlin, \email{rubenarslan@gmail.com})
    Yves Rosseel (Ghent University; \email{Yves.Rosseel@UGent.be})
}
\seealso{
	\itemize{
		\item \code{\link{miPowerFit}} For the modification indices and their power approach for model fit evaluation
		\item \code{\link{nullRMSEA}} For RMSEA of the null model
	}
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
moreFitIndices(fit)

fit2 <- cfa(HS.model, data=HolzingerSwineford1939, estimator="mlr")
moreFitIndices(fit2)
}
