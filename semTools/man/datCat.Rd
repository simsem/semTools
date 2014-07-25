\name{datCat}
\alias{datCat}
\title{
Simulated Data set to Demonstrate Categorical Measurement Invariance
}
\description{
A simulated data set with 2 factors with 4 indicators each separated into two groups
}
\usage{
data(datCat)
}
\format{
  A data frame with 200 observations of 9 variables.
	\describe{
   \item{g}{Sex of respondents}
   \item{u1}{Indicator 1}
   \item{u2}{Indicator 2}
   \item{u3}{Indicator 3}
   \item{u4}{Indicator 4}
   \item{u5}{Indicator 5}
   \item{u6}{Indicator 6}
   \item{u7}{Indicator 7}
   \item{u8}{Indicator 8}
   }
}
\source{
Data was generated using the \code{lavaan} package.
}
\examples{
head(datCat)
}