\name{miPowerFit}
\alias{miPowerFit}
\alias{miPowerFit}
\title{
Modification indices and their power approach for model fit evaluation
}
\description{
The model fit evaluation approach using modification indices and expected parameter changes.
}
\usage{
miPowerFit(lavaanObj, stdLoad=0.4, cor=0.1, stdBeta=0.1, intcept=0.2, stdDelta=NULL, 
	delta=NULL, cilevel = 0.90)
}
\arguments{
  \item{lavaanObj}{The lavaan model object used to evaluate model fit}
  \item{stdLoad}{The amount of standardized factor loading that one would like to be detected (rejected). The default value is 0.4, which is suggested by Saris and colleagues (2009, p. 571).}
  \item{cor}{The amount of factor or error correlations that one would like to be detected (rejected). The default value is 0.1, which is suggested by Saris and colleagues (2009, p. 571).}
  \item{stdBeta}{The amount of standardized regression coefficients that one would like to be detected (rejected). The default value is 0.1, which is suggested by Saris and colleagues (2009, p. 571).}
  \item{intcept}{The amount of standardized intercept (similar to Cohen's \emph{d} that one would like to be detected (rejected). The default value is 0.2, which is equivalent to a low effect size proposed by Cohen (1988, 1992).}
  \item{stdDelta}{The vector of the standardized parameters that one would like to be detected (rejected). If this argument is specified, the value here will overwrite the other arguments above. The order of the vector must be the same as the row order from modification indices from the \code{lavaan} object. If a single value is specified, the value will be applied to all parameters.}
  \item{delta}{The vector of the unstandardized parameters that one would like to be detected (rejected). If this argument is specified, the value here will overwrite the other arguments above. The order of the vector must be the same as the row order from modification indices from the \code{lavaan} object. If a single value is specified, the value will be applied to all parameters.}
  \item{cilevel}{The confidence level of the confidence interval of expected parameter changes. The confidence intervals are used in the equivalence testing.}
}
\details{
In the lavaan object, one can inspect the modification indices and expected parameter changes. Those values can be used to evaluate model fit by two methods. 

First, Saris, Satorra, and van der Veld (2009, pp. 570-573) used the power to detect modification indices and expected parameter changes to evaluate model fit. First, one should evaluate whether the modification index of each parameter is significant. Second, one should evaluate whether the power to detect a target expected parameter change is high enough. If the modification index is not significant and the power is high, there is no misspecification. If the modification index is significant and the power is low, the fixed parameter is misspecified. If the modification index is significant and the power is high, the expected parameter change is investigated. If the expected parameter change is large (greater than the the target expected parameter change), the parameter is misspecified. If the expected parameter change is low (lower than the target expected parameter change), the parameter is not misspecificied. If the modification index is not significant and the power is low, the decision is inconclusive.

Second, the confidence intervals of the expected parameter changes are formed. These confidence intervals are compared with the range of trivial misspecification, which could be (-\code{delta}, \code{delta}) or (0, \code{delta}) for nonnegative parameters. If the confidence intervals are outside of the range of trivial misspecification, the fixed parameters are severely misspecified. If the confidence intervals are inside the range of trivial misspecification, the fixed parameters are trivially misspecified. If confidence intervals are overlapped the range of trivial misspecification, the decision is inconclusive.
}
\value{
A data frame with these variables:
\enumerate{
  \item{lhs} The left-hand side variable (with respect to the lavaan operator)
  \item{op} The lavaan syntax operator: "~~" represents covariance, "=~" represents factor loading, "~" represents regression, and "~1" represents intercept. 
  \item{rhs} The right-hand side variable (with respect to the lavaan operator)
  \item{group} The group of the parameter
  \item{mi} The modification index of the fixed parameter
  \item{epc} The expected parameter change if the parameter is freely estimated
  \item{target.epc} The target expected parameter change that represents the minimum size of misspecification that one would like to be detected by the test with a high power
  \item{std.epc} The standardized expected parameter change if the parameter is freely estimated
  \item{std.target.epc} The standardized target expected parameter change
  \item{significant.mi} Represents whether the modification index value is significant
  \item{high.power} Represents whether the power is enough to detect the target expected parameter change
  \item{decision.pow} The decision whether the parameter is misspecified or not based on Saris et al's method: \code{"M"} represents the parameter is misspecified, \code{"NM"} represents the parameter is not misspecified, \code{"EPC:M"} represents the parameter is misspecified decided by checking the expected parameter change value, \code{"EPC:NM"} represents the parameter is not misspecified decided by checking the expected parameter change value, and \code{"I"} represents the decision is inconclusive.
  \item{se.epc} The standard errors of the expected parameter changes.
  \item{lower.epc} The lower bound of the confidence interval of expected parameter changes.
  \item{upper.epc} The upper bound of the confidence interval of expected parameter changes.
  \item{lower.std.epc} The lower bound of the confidence interval of standardized expected parameter changes.
  \item{upper.std.epc} The upper bound of the confidence interval of standardized expected parameter changes.
  \item{decision.ci} The decision whether the parameter is misspecified or not based on the confidence interval method: \code{"M"} represents the parameter is misspecified, \code{"NM"} represents the parameter is not misspecified, and \code{"I"} represents the decision is inconclusive.
}
The row numbers matches with the results obtained from the \code{inspect(object, "mi")} function.
}
\references{
Cohen, J. (1988). \emph{Statistical power analysis for the behavioral sciences} (2nd ed.). Hillsdale, NJ: Erlbaum.

Cohen, J. (1992). A power primer. \emph{Psychological Bulletin, 112}, 155-159.

Saris, W. E., Satorra, A., & van der Veld, W. M. (2009). Testing structural equation models or detection of misspecifications? \emph{Structural Equation Modeling, 16}, 561-582.
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\seealso{
	\code{\link{moreFitIndices}} For the additional fit indices information
}
\examples{
library(lavaan)

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939, group="sex", meanstructure=TRUE)
miPowerFit(fit)

model <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
fit2 <- sem(model, data=PoliticalDemocracy, meanstructure=TRUE)
miPowerFit(fit2, stdLoad=0.3, cor=0.2, stdBeta=0.2, intcept=0.5)
}
