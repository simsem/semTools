% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EmpKaiser.R
\name{efa.ekc}
\alias{efa.ekc}
\title{Empirical Kaiser criterion}
\usage{
efa.ekc(data = NULL, sample.cov = NULL, sample.nobs = NULL,
  missing = "default", ordered = NULL, plot = TRUE)
}
\arguments{
\item{data}{A \code{data.frame} or data \code{matrix} containing columns of
variables to be factor-analyzed.}

\item{sample.cov}{A covariance or correlation matrix can be used, instead of
\code{data}, to estimate the eigenvalues.}

\item{sample.nobs}{Number of observations (i.e. sample size) if
\code{is.null(data)} and \verb{sample.cov=} is used.}

\item{missing}{If \code{"listwise"}, incomplete cases are removed listwise from
the \code{data.frame}. If \code{"direct"} or \code{"ml"} or \code{"fiml"} and the \verb{estimator=}
is maximum likelihood, an EM algorithm is used to estimate an unrestricted
covariance matrix (and mean vector). If \code{"pairwise"}, pairwise deletion is
used. If `"default"``, the value is set depending on the estimator and the
mimic option (see \code{\link[lavaan:lavCor]{lavaan::lavCor()}} for details).}

\item{ordered}{\code{character} vector. Only used if object is a \code{data.frame}.
Treat these variables as \verb{ordered=} (ordinal) variables. Importantly, all
other variables will be treated as numeric (unless \code{is.ordered == TRUE} in
\code{data}). (see also \link[lavaan:lavCor]{lavCor})}

\item{plot}{logical. Whether to print a scree plot comparing the sample
eigenvalues with the reference eigenvalues.}
}
\value{
A \code{data.frame} showing the sample and reference eigenvalues.

The number of factors suggested by the Empirical Kaiser Criterion (i.e. the
sample eigenvalues greater than the reference eigenvalues) is returned as an
attribute (see \strong{Examples}).

The number of factors suggested by the original Kaiser Criterion (i.e.
sample eigenvalues > 1) is also printed as a header to the \code{data.frame}
}
\description{
Identify the number of factors to extract based on the Empirical Kaiser
Criterion (EKC). The analysis can be run on a \code{data.frame} or data
\code{matrix} (\code{data}), or on a correlation or covariance matrix
(\code{sample.cov}) and the sample size (\code{sample.nobs}). A
\code{data.frame} is returned with two columns: the eigenvalues from your
data or covariance matrix and the reference eigenvalues. The number of
factors suggested by the Empirical Kaiser Criterion (i.e. the sample
eigenvalues greater than the reference eigenvalues), and the number of
factors suggested by the original Kaiser Criterion
(i.e. sample eigenvalues > 1) is printed above the output.
}
\examples{

## Simulate data with 3 factors
model <- '
  f1 =~ .3*x1 + .5*x2 + .4*x3
  f2 =~ .3*x4 + .5*x5 + .4*x6
  f3 =~ .3*x7 + .5*x8 + .4*x9
'
dat <- simulateData(model, seed = 123)
## save summary statistics
myCovMat <- cov(dat)
myCorMat <- cor(dat)
N <- nrow(dat)

## Run the EKC function
(out <- efa.ekc(dat))

## To extract the recommended number of factors using the EKC:
attr(out, "nfactors")

## If you do not have raw data, you can use summary statistics
(x1 <- efa.ekc(sample.cov = myCovMat, sample.nobs = N, plot = FALSE))
(x2 <- efa.ekc(sample.cov = myCorMat, sample.nobs = N, plot = FALSE))

}
\references{
Braeken, J., & van Assen, M. A. L. M. (2017). An empirical
Kaiser criterion. \emph{Psychological Methods, 22}(3), 450--466.
\doi{10.1037/met0000074}
}
\author{
Ylenio Longo (University of Nottingham;
\email{yleniolongo@gmail.com})

Terrence D. Jorgensen (University of Amsterdam;
\email{TJorgensen314@gmail.com})
}
