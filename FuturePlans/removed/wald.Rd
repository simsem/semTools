\name{wald}
\alias{wald}
\title{
Calculate multivariate Wald statistics
}
\description{
Calculate multivariate Wald statistics based on linear combinations of model parameters
}
\usage{
wald(object, syntax)
}
\arguments{
  \item{object}{An output from \code{lavaan}}
  \item{syntax}{Syntax that each line represents one linear constraint. A plus or minus sign is used to separate between each coefficient. An asterisk is used to separate between coefficients and parameters. The coefficient can have a forward slash to represent a division. The parameter names must be matched with the names of lavaan parameters investigated by running the \code{coef} function on a lavaan output. Lines can be separated by semi-colon. A pound sign is allowed for comments. Note that the defined parameters (created by ":=") do not work with this function.}
}
\details{
The formula for multivariate Wald test is

\deqn{ \chi^2 = \left(C\hat{b}\right)^\prime\left[C\hat{V}C^\prime\right]^{-1}\left(C\hat{b}\right),}

where \eqn{C} is the contrast matrix, \eqn{\hat{b}} is the estimated fixed effect, \eqn{\hat{V}} is the asymptotic covariance matrix among fixed effects. 
}
\value{
Chi-square value with \emph{p} value.
}
\author{
Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
# Test the difference in factor loadings
library(lavaan)
HS.model <- ' visual  =~ x1 + con1*x2 + con1*x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + con2*x8 + con2*x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
wald(fit, "con2 - con1")

# Simultaneously test the difference in the influences 
# of x1 and x2 on intercept and slope
model.syntax <- '
    i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
    s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4
    i ~ x1 + x2
    s ~ x1 + x2
    t1 ~ c1
    t2 ~ c2
    t3 ~ c3
    t4 ~ c4
'

fit2 <- growth(model.syntax, data=Demo.growth)
wald.syntax <- '
	i~x1 - i~x2
	1/2*s~x1 - 1/2*s~x2
'
wald(fit2, wald.syntax)

# Mplus example of MODEL TEST
model3 <- ' f1  =~ x1 + p2*x2 + p3*x3 + p4*x4 + p5*x5 + p6*x6
			p4 == 2*p2'

fit3 <- cfa(model3, data=HolzingerSwineford1939)
wald(fit3, "p3; p6 - 0.5*p5")
}
