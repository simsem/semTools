\name{miPoolChi}
\alias{miPoolChi}
\title{
Function to pool chi-square statistics from the result from multiple imputation
}
\description{
The function combines likelihood ratio chi-square statistics from an analysis of multiply imputed data sets using the method proposed by Li, Meng, Raghunathan, and Rubin (1991, p. 74).
}
\usage{
miPoolChi(chis, df)
}
\arguments{
  \item{chis}{
A vector of chi-square statistics
}
  \item{df}{
Degree of freedom that the chi-square statistics are based on
}
}
\details{
The chi-square statistics are pooled by the procedure proposed by Li, Meng, Raghunathan, and Rubin (1991; Equations 2.1, 2.2, 2.16, and 2.17). 
}
\value{
The resulting value from the pooled chi-square is F-statistic. If the denominator degree of freedom is large, the F value multiplied by the numerator degree of freedom will approximate the chi-square statistics.
}
\references{
Li, K. H., Meng, X. L., Raghunathan, T. E., & Rubin, D. B. (1991). Significance levels from repeated p-values with multiply-imputed data. \emph{Statistica Sinica, 1,} 65-92.
}
\author{
Craig Enders originally wrote this function in SAS, \url{http://psychology.clas.asu.edu/files/CombiningLikelihoodRatioChi-SquareStatisticsFromaMIAnalysis.sas}. Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu}) modified the function to run in R.
}
\examples{
miPoolChi(c(89.864, 81.116, 71.500, 49.022, 61.986, 64.422, 55.256, 57.890, 79.416, 63.944), 2)
}
