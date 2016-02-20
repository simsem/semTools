\name{loadingFromAlpha}
\alias{loadingFromAlpha}
\title{Find standardized factor loading from coefficient alpha}
\description{
    Find standardized factor loading from coefficient alpha assuming that all items have equal loadings.
}
\usage{
loadingFromAlpha(alpha, ni)
}
\arguments{
  \item{alpha}{A desired coefficient alpha value.}
  \item{ni}{A desired number of items.}
}
\value{
    \item{result}{The standardized factor loadings that make desired coefficient alpha with specified number of items.}
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
    loadingFromAlpha(0.8, 4)
}
