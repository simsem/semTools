\name{net}
\alias{net}
\title{
	Nesting and Equivalence Testing
}
\description{
This test examines whether models are nested or equivalent based on Bentler and Satorra's (2010) procedure.
}
\usage{
net(..., crit = .0001)
}
\arguments{
  \item{\dots}{
	The \code{lavaan} objects used for test of nesting and equivalence
}
  \item{crit}{
	The upper-bound criterion for testing the equivalence of models. Models are considered nested (or equivalent) if the difference between their chi-squared fit statistics is less than this criterion.
}
}
\details{
	The concept of nesting/equivalence should be the same regardless of estimation method. However, the particular method of testing nesting/equivalence (as described in Bentler & Satorra, 2010) employed by the net function is based on a limited-information estimator (analyzing model-implied means and covariance matrices, not raw data). In the case of robust methods like MLR, the raw data is only utilized for the robust adjustment to SE and chi-sq, and the net function only checks the unadjusted chi-sq for the purposes of testing nesting/equivalence. 
	This method does not apply to models that estimate thresholds for categorical data, so an error message will be issued if such a model is provided.
}
\value{
	The \linkS4class{Net} object representing the outputs for nesting and equivalent testing, including a logical matrix of test results and a vector of degrees of freedom for each model.
}
\author{
    Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\references{
Bentler, P. M., & Satorra, A. (2010). Testing model nesting and equivalence. \emph{Psychological Methods, 15}, 111-123. doi:10.1037/a0019625
}
\examples{
\dontrun{
m1 <- ' visual  =~ x1 + x2 + x3
	textual =~ x4 + x5 + x6
	speed   =~ x7 + x8 + x9 '

	
m2 <- ' f1  =~ x1 + x2 + x3 + x4
	f2 =~ x5 + x6 + x7 + x8 + x9 '

m3 <- ' visual  =~ x1 + x2 + x3
	textual =~ eq*x4 + eq*x5 + eq*x6
	speed   =~ x7 + x8 + x9 '
	
fit1 <- cfa(m1, data = HolzingerSwineford1939)
fit1a <- cfa(m1, data = HolzingerSwineford1939, std.lv = TRUE) # Equivalent to fit1
fit2 <- cfa(m2, data = HolzingerSwineford1939) # Not equivalent to or nested in fit1
fit3 <- cfa(m3, data = HolzingerSwineford1939) # Nested in fit1 and fit1a

tests <- net(fit1, fit1a, fit2, fit3)
tests
summary(tests)
}
}

