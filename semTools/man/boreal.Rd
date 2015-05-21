\name{boreal}
\alias{boreal}
\title{
The Boreal Vegetation Dataset 
}
\description{
A part of the Boreal Vegetation dataset from Zuur et al. (2009)
}
\usage{
data(boreal)
}
\format{
  A data frame with 533 observations of 7 variables.
	\describe{
   \item{point}{Point ID}
   \item{x}{X-coordinate}
   \item{y}{Y-coordinate}
   \item{nTot}{species richness}
   \item{NDVI}{The normalized difference vegetation index from satellite data}
   \item{T61}{climate}
   \item{Wet}{Wetness}
   }
}
\references{
Zuur, A., Ieno, E. N., Walker, N., Saveliev, A. A., & Smith, G. M. (2009). \emph{Mixed effects models and extensions in ecology with R}. New York: Springer.
}
\examples{
head(boreal)
}