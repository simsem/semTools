\name{impliedFactorStat}
\alias{impliedFactorMean}
\alias{impliedFactorCov}
\alias{impliedFactorStat}
\title{
Calculate the model-implied factor means and covariance matrix. 
}
\description{
Calculate reliability values of factors by coefficient omega
}
\usage{
impliedFactorStat(object)
impliedFactorMean(object)
impliedFactorCov(object)
}
\arguments{
  \item{object}{The lavaan model object provided after running the \code{cfa}, \code{sem}, \code{growth}, or \code{lavaan} functions.}
}
\details{
The \code{impliedFactorMean} function is used to calculated model-implied factor means:

\deqn{ 
	\mu = \left( \bold{I} - \bold{B} \right)^{-1} \alpha ,
}

where \eqn{\mu} is the model-implied factor mean, \eqn{\bold{I}} is an identity matrix, \eqn{\bold{B}} is an regression coefficient matrix, and \eqn{\alpha} is a vector of factor intercepts.

The \code{impliedFactorCov} function is used to calculated model-implied covariance matrix:

\deqn{ \Phi = \left( \bold{I} - \bold{B} \right)^{-1} \Psi \left(\bold{I} - \bold{B} \right)^{-1\prime} ,}

where \eqn{\Phi} is the model-implied factor covariance matrix, \eqn{\Psi} is the residual factor covariance matrix.

The \code{impliedFactorStat} function is used to provide both model-implied means (if the mean structure is estimated) and covariance matrix.
}
\value{
	Model-implied factor means or model-implied factor covariance matrix, or both
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939, group="school")
impliedFactorStat(fit)
}
