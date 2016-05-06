\name{htmt}
\alias{htmt}
\title{
	Assessing Discriminant Validity using Heterotrait-Monotrait Ratio
}
\description{
	This function assesses the heterotrait-monotrait ratio (HTMT) of the correlations (i.e., the correlations of indicators across constructs measuring different phenomena), relative to the average of the monotrait-heteromethod correlations (i.e., the correlations of indicators within the same construct).
}
\usage{
htmt(data, model)
}
\arguments{
  \item{data}{
	A desired data set
}
  \item{model}{
	lavaan syntax of a confirmatory factor analysis model where at least two factors are required to indicate indicators measuring the same construct.
  }
}
\value{
	A matrix shown htmt values (i.e., discriminant validity) between each pair of factors.
}
\references{
Henseler, J., Ringle, C. M., & Sarstedt, M. (2015). A new criterion for assessing discriminant validity in variance-based structural equation modeling. \emph{Journal of the Academy of Marketing Science, 43}, 115-135.
}
\author{
	Ylenio Longo (University of Nottingham; \email{ylenio.longo@nottingham.ac.uk})
}
\examples{
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

dat <- HolzingerSwineford1939[, paste0("x", 1:9)]
htmt(dat, HS.model)
}
