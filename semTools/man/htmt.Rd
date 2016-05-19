\name{htmt}
\alias{htmt}
\title{
	Assessing Discriminant Validity using Heterotrait-Monotrait Ratio
}
\description{
	This function assesses discriminant validity through the heterotrait-monotrait ratio (HTMT) of the correlations (Henseler, Ringlet & Sarstedt, 2015). Specifically, it assesses the average correlation among indicators across constructs (i.e. heterotrait-heteromethod correlations), relative to the average correlation among indicators within the same construct (i.e. monotrait-heteromethod correlations). The resulting HTMT values are interpreted as estimates of inter-construct correlations. Absolute values of the correlations are used to calculate the HTMT matrix.
}
\usage{
htmt(data, model, ...)
}
\arguments{
  \item{data}{
	A desired data set
}
  \item{model}{
	lavaan syntax of a confirmatory factor analysis model where at least two factors are required to indicate indicators measuring the same construct.
  }
  \item{\dots}{
	Other arguments shown in \link[lavaan]{lavCor}
  }
}
\value{
	A matrix showing HTMT values (i.e., discriminant validity) between each pair of factors.
}
\references{
Henseler, J., Ringle, C. M., & Sarstedt, M. (2015). A new criterion for assessing discriminant validity in variance-based structural equation modeling. \emph{Journal of the Academy of Marketing Science, 43}, 115-135.
}
\author{
	Ylenio Longo (University of Nottingham; \email{yleniolongo@gmail.com})
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

dat <- HolzingerSwineford1939[, paste0("x", 1:9)]
htmt(dat, HS.model)
}
