\name{reliabilityL2}
\alias{reliabilityL2}
\title{
Calculate the reliability values of a second-order factor
}
\description{
Calculate the reliability values (coefficient omega) of a second-order factor
}
\usage{
reliabilityL2(object, secondFactor)
}
\arguments{
  \item{object}{The lavaan model object provided after running the \code{cfa}, \code{sem}, \code{growth}, or \code{lavaan} functions that has a second-order factor}
  \item{secondFactor}{The name of the second-order factor}
}
\details{
The first formula of the coefficient omega (in the \code{\link{reliability}}) will be mainly used in the calculation. The model-implied covariance matrix of a second-order factor model can be separated into three sources: the second-order factor, the uniqueness of the first-order factor, and the measurement error of indicators: 

\deqn{ \hat{\Sigma} = \Lambda \bold{B} \Phi_2 \bold{B}^{\prime} \Lambda^{\prime} + \Lambda \Psi_{u} \Lambda^{\prime} + \Theta, }

where \eqn{\hat{\Sigma}} is the model-implied covariance matrix, \eqn{\Lambda} is the first-order factor loading, \eqn{\bold{B}} is the second-order factor loading, \eqn{\Phi_2} is the covariance matrix of the second-order factors, \eqn{\Psi_{u}} is the covariance matrix of the unique scores from first-order factors, and \eqn{\Theta} is the covariance matrix of the measurement errors from indicators. Thus, the proportion of the second-order factor explaining the total score, or the coefficient omega at Level 1, can be calculated:

\deqn{ \omega_{L1} = \frac{\bold{1}^{\prime} \Lambda \bold{B} \Phi_2 \bold{B}^{\prime} \Lambda^{\prime} \bold{1}}{\bold{1}^{\prime} \Lambda \bold{B} \Phi_2 \bold{B} ^{\prime} \Lambda^{\prime} \bold{1} + \bold{1}^{\prime} \Lambda \Psi_{u} \Lambda^{\prime} \bold{1} + \bold{1}^{\prime} \Theta \bold{1}}, }

where \eqn{\bold{1}} is the \emph{k}-dimensional vector of 1 and \emph{k} is the number of observed variables. When model-implied covariance matrix among first-order factors (\eqn{\Phi_1}) can be calculated:

\deqn{ \Phi_1 = \bold{B} \Phi_2 \bold{B}^{\prime} + \Psi_{u}, }

Thus, the proportion of the second-order factor explaining the varaince at first-order factor level, or the coefficient omega at Level 2, can be calculated:

\deqn{ \omega_{L2} = \frac{\bold{1_F}^{\prime} \bold{B} \Phi_2 \bold{B}^{\prime} \bold{1_F}}{\bold{1_F}^{\prime} \bold{B} \Phi_2 \bold{B}^{\prime} \bold{1_F} + \bold{1_F}^{\prime} \Psi_{u} \bold{1_F}}, }

where \eqn{\bold{1_F}} is the \emph{F}-dimensional vector of 1 and \emph{F} is the number of first-order factors. 

The partial coefficient omega at Level 1, or the proportion of observed variance explained by the second-order factor after partialling the uniqueness from the first-order factor, can be calculated:

\deqn{ \omega_{L1} = \frac{\bold{1}^{\prime} \Lambda \bold{B} \Phi_2 \bold{B}^{\prime} \Lambda^{\prime} \bold{1}}{\bold{1}^{\prime} \Lambda \bold{B} \Phi_2 \bold{B}^{\prime} \Lambda^{\prime} \bold{1} + \bold{1}^{\prime} \Theta \bold{1}}, }

Note that if the second-order factor has a direct factor loading on some observed variables, the observed variables will be counted as first-order factors.
}
\value{
	Reliability values at Levels 1 and 2 of the second-order factor, as well as the partial reliability value at Level 1
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\seealso{
	\code{\link{reliability}} for the reliability of the first-order factors.
}
\examples{
library(lavaan)

HS.model3 <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 
			  higher =~ visual + textual + speed'

fit6 <- cfa(HS.model3, data=HolzingerSwineford1939)
reliability(fit6) # Should provide a warning for the endogenous variable
reliabilityL2(fit6, "higher")
}
