% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quark.R
\name{combinequark}
\alias{combinequark}
\title{Combine the results from the quark function}
\usage{
combinequark(quark, percent)
}
\arguments{
\item{quark}{Provide the \code{\link[=quark]{quark()}} object that was returned.  It
should be a list of objects.  Make sure to include it in its entirety.}

\item{percent}{Provide a percentage of variance that you would like to have
explained.  That many components (columns) will be extracted and kept with
the output dataset.  Enter this variable as a number WITHOUT a percentage
sign.}
}
\value{
The output of this function is the original dataset used in quark
combined with enough principal component scores to be able to account for
the amount of variance that was requested.
}
\description{
This function builds upon the \code{\link[=quark]{quark()}} function to provide a
final dataset comprised of the original dataset provided to
\code{\link[=quark]{quark()}} and enough principal components to be able to account
for a certain level of variance in the data.
}
\examples{

set.seed(123321)
dat <- HolzingerSwineford1939[,7:15]
misspat <- matrix(runif(nrow(dat) * 9) < 0.3, nrow(dat))
dat[misspat] <- NA
dat <- cbind(HolzingerSwineford1939[,1:3], dat)

quark.list <- quark(data = dat, id = c(1, 2))

final.data <- combinequark(quark = quark.list, percent = 80)

}
\seealso{
\code{\link[=quark]{quark()}}
}
\author{
Steven R. Chesnut (University of Southern Mississippi
\email{Steven.Chesnut@usm.edu})
}
