\name{exLong}
\alias{exLong}
\title{
Simulated Data set to Demonstrate Longitudinal Measurement Invariance
}
\description{
A simulated data set with 1 factors with 3 indicators in three timepoints
}
\usage{
data(exLong)
}
\format{
  A data frame with 200 observations of 10 variables.
	\describe{
   \item{sex}{Sex of respondents}
   \item{y1t1}{Indicator 1 in Time 1}
   \item{y2t1}{Indicator 2 in Time 1}
   \item{y3t1}{Indicator 3 in Time 1}
   \item{y1t2}{Indicator 1 in Time 2}
   \item{y2t2}{Indicator 2 in Time 2}
   \item{y3t2}{Indicator 3 in Time 2}
   \item{y1t3}{Indicator 1 in Time 3}
   \item{y2t3}{Indicator 2 in Time 3}
   \item{y3t3}{Indicator 3 in Time 3}
   }
}
\source{
Data was generated using the \code{simsem} package.
}
\examples{
head(exLong)
}