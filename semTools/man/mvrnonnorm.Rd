\name{mvrnonnorm}
\alias{mvrnonnorm}
\title{
Generate Non-normal Data using Vale and Maurelli (1983) method
}
\description{
Generate Non-normal Data using Vale and Maurelli (1983) method. The function is designed to be as similar as the popular \code{mvrnorm} function in the \code{MASS} package. The codes are copied from \code{mvrnorm} function in the \code{MASS} package for argument checking and \code{lavaan} package for data generation using Vale and Maurelli (1983) method.
}
\usage{
mvrnonnorm(n, mu, Sigma, skewness = NULL, kurtosis = NULL, empirical = FALSE)
}
\arguments{
  \item{n}{Sample size}
  \item{mu}{A mean vector}
  \item{Sigma}{A positive-definite symmetric matrix specifying the covariance matrix of the variables}
  \item{skewness}{A vector of skewness of the variables}
  \item{kurtosis}{A vector of excessive kurtosis of the variables}
  \item{empirical}{If \code{TRUE}, \code{mu} and \code{Sigma} specify the empirical not population mean and covariance matrix}
}
\value{
	A data matrix
}
\references{
Vale, C. D. & Maurelli, V. A. (1983) Simulating multivariate nonormal distributions. \emph{Psychometrika, 48}, 465-471.
}
\author{
    The original function is the \code{simulateData} function written by Yves Rosseel in the \code{lavaan} package. The function is adjusted for a convenient usage by Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
mvrnonnorm(100, c(1, 2), matrix(c(10, 2, 2, 5), 2, 2), 
	skewness = c(5, 2), kurtosis = c(3, 3))
}
