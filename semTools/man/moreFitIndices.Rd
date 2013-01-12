\name{moreFitIndices}
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
  \item{object}{The lavaan model object provided after running the \code{cfa}, \code{sem}, \code{growth}, or \code{lavaan} functions.}
  \item{nPrior}{The sample size on which prior is based. This argument is used to compute BIC*.}
}
\details{
Normed Fit Index (nfi; West, Taylor, & Wu, 2012) is one of the relative fit indices which can be computed by

\deqn{ nfi =\frac{\chi^{2}_{0} - \chi^{2}_{k}}{\chi^{2}_{0}},}

where \eqn{\chi^{2}_{k}} is the chi-square test statistic value of the target model, \eqn{\chi^{2}_{0}} is the chi-square test statistic value of the null model.  

Incremental Fit Index (ifi; Bollen, 1989), also called Bollen's fit index (BL89; Hu & Bentler, 1998) is one of the relative fit indices which can be computed by

\deqn{ ifi =\frac{\chi^{2}_{0} - \chi^{2}_{k}}{\chi^{2}_{0} - df_{k}},}

where \eqn{df_{k}} is the degree of freedom when fitting the target model.

Relative Noncentrality Index (rni; Hu & Bentler, 1998) is a relative fit index almost identical to CFI, but without the upper limit of 1.  Thus, CFI and RNI are identical except when the fit of the target model is better than its expected value (i.e., the chi-square is less than the degrees of freedom).  It is computed by

\deqn{ rni =\frac{(\chi^{2}_{0} - df_{0}) - (\chi^{2}_{k} - df_{k})}{\chi^{2}_{0} - df_{0}},}

where \eqn{df_{0}} is the degree of freedom when fitting the null model.

McDonald's Centrality Index (Mc; McDonald, 1989) is a global fit index which can be computed by

\deqn{ Mc = exp(-0.5 * ((\chi^{2}_{k} - df_{k}) / N)),}

where \eqn{N} is the sample size.

Gamma Hat (gfi*; West, Taylor, & Wu, 2012) is a global fit index which can be computed by

\deqn{ gfi* =\frac{p}{p + 2 \times \frac{\chi^{2}_{k} - df_{k}}{N - 1}},}

where \eqn{p} is the number of variables in the model. This formula assumes equal number of indicators across groups.

Adjusted Gamma Hat (agfi*; West, Taylor, & Wu, 2012) is a global fit index which can be computed by

\deqn{ agfi* = \left(1 - \frac{K \times p \times (p + 1)}{2 \times df_{k}} \right) \times \left( 1 - gfi* \right) ,}

where \eqn{K} is the number of group (please refer to Dudgeon, 2004 for the multiple-group adjustment for agfi*).

Corrected Akaike Information Criterion (AICc; Burnham & Anderson, 2003) is the corrected version of aic for small sample size:

\deqn{ aicc = f + \frac{2k(k + 1)}{N - k - 1},}

where \eqn{f} is the minimized discrepancy function, which is the product of the log likelihood and -2, and \eqn{k} is the number of parameters in the target model.

Expected Value of Cross-Validation Index (ECVI; West, Taylor, & Wu, 2012) is the average discrepancy in the fitted covariance matrices between two samples of equal sample size across all possible combinations of two samples from the same population:

\deqn{ ecvi = f + \frac{2 \times k}{N},}

Stochastic information criterion (sic; Preacher, 2006) is similar to aic or bic. This index will account for model complexity in the model's function form, in addition to the number of free parameters. This index will be provided only when the chi-squared value is not scaled. sic can be computed by

\deqn{ sic = \frac{1}{2}\left(f - \log{\det{I(\hat{\theta})}}\right),}

where \eqn{I(\hat{\theta})} is the information matrix of the parameters.

Corrected Bayesian Information Criterion (BIC*; Kuha, 2004) is similar to bic but explicitly specifying the sample size on which the prior is based (\eqn{N_{prior}}).

\deqn{ bicc = f + k\log{(1 + N/N_{prior})},}

Hannan-Quinn Information Criterion (hqc; Hannan & Quinn, 1979) is used for model selection similar to aic or bic.

\deqn{ hqc = f + 2k\log{(\log{N})},}

Note that if Satorra-Bentler or Yuan-Bentler's method is used, the fit indices using the scaled chi-square values are also provided. 

See \code{\link{nullRMSEA}} for the further details of the computation of RMSEA of the null model.
}

\value{
\enumerate{
  \item{nfi} Normed Fit Index
  \item{ifi} Incremental Fit Index
  \item{rni} Relative Noncentrality Index
  \item{Mc} McDonald's Centrality Index
  \item{gfiStar} Gamma Hat
  \item{agfiStar} Adjusted Gamma Hat
  \item{baseline.rmsea} RMSEA of the Baseline (Null) Model
  \item{aicc} Corrected Akaike Information Criterion
  \item{ecvi} Expected Value of Cross-Validation Index 
  \item{sic} Stochastic Information Criterion
  \item{bicStar} Bayesian Information Criterion with specifying the prior sample size
  \item{hqc} Hannan-Quinn Information Criterion
  \item{nfi.scaled} Normed Fit Index using Scaled Chi-square
  \item{ifi.scaled} Incremental Fit Index using Scaled Chi-square
  \item{gfiStar.scaled} Gamma Hat using Scaled Chi-square
  \item{agfiStar.scaled} Adjusted Gamma Hat using Scaled Chi-square
  \item{baseline.rmsea.scaled} RMSEA of the Baseline (Null) Model using Scaled Chi-square
}
}
\references{
Bollen, K. A. (1989). A new incremental fit index for general structural equation models. \emph{Sociological Methods and Research, 17}, 303-316.

Burnham, K., & Anderson, D. (2003). \emph{Model selection and multimodel inference: A practical-theoretic approach.} New York, NY: Springer-Verlag.

Dudgeon, P. (2004). A note on extending Steiger's (1998) multiple sample RMSEA adjustment to other noncentrality parameter-based statistic. \emph{Structural Equation Modeling, 11}, 305-319.

Hu, L. & Bentler, P. M. (1998). Fit indices in covariance structure modeling: Sensitivity to underparameterized model misspecification. \emph{Psychological Methods, 3}, 424-453.

Kuha, J. (2004). AIC and BIC: Comparisons of assumptions and performance. \emph{Sociological Methods Research, 33}, 188-229.

McDonald, R. P. (1989). An index of goodnes-of-fit based on noncentrality. \emph{Journal of Classification, 6}, 97-103.

Preacher, K. J. (2006). Quantifying parsimony in structural equation modeling. \emph{Multivariate Behavioral Research, 43}, 227-259.

West, S. G., Taylor, A. B., & Wu, W. (2012). Model fit and model selection in structural equation modeling. In R. H. Hoyle (Ed.), \emph{Handbook of Structural Equation Modeling.} New York: Guilford.
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
	Aaron Boulton (University of Kansas; \email{aboulton@ku.edu})
	Ruben Arslan (Humboldt-University of Berlin, \email{rubenarslan@gmail.com})
  Terrence Jorgensen (University of Kansas; \email{tdj@ku.edu})
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
