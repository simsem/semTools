\name{plotRMSEApower}
\alias{plotRMSEApower}
\title{
Plot power curves for RMSEA
}
\description{
Plots power of RMSEA over a range of sample sizes
}
\usage{
plotRMSEApower(rmsea0, rmseaA, df, nlow, nhigh, steps=1, alpha=.05, group=1, ...)
}
\arguments{
  \item{rmsea0}{Null RMSEA}
  \item{rmseaA}{Alternative RMSEA}
  \item{df}{Model degrees of freedom}
  \item{nlow}{Lower sample size}
  \item{nhigh}{Upper sample size}
  \item{steps}{Increase in sample size for each iteration. Smaller values of steps will lead to more precise plots. However, smaller step sizes means a longer run time.}
  \item{alpha}{Alpha level used in power calculations}
  \item{group}{The number of group that is used to calculate RMSEA.}
  \item{\dots}{The additional arguments for the plot function.}
  }
\details{
This function creates plot of power for RMSEA against a range of sample sizes. The plot places sample size on the horizontal axis and power on the vertical axis. The user should indicate the lower and upper values for sample size and the sample size between each estimate ("step size") We strongly urge the user to read the sources below (see References) before proceeding.  A web version of this function is available at: \url{http://quantpsy.org/rmsea/rmseaplot.htm}.
}

\value{
	\enumerate{
		\item{plot} Plot of power for RMSEA against a range of sample sizes
	}
}
\references{
MacCallum, R. C., Browne, M. W., & Cai, L. (2006). Testing differences between nested covariance structure models: Power analysis and null hypotheses. \emph{Psychological Methods, 11,} 19-35.

MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis and determination of sample size for covariance structure modeling. \emph{Psychological Methods, 1,} 130-149.

MacCallum, R. C., Lee, T., & Browne, M. W. (2010). The issue of isopower in power analysis for tests of structural equation models. \emph{Structural Equation Modeling, 17,} 23-41.

Preacher, K. J., Cai, L., & MacCallum, R. C. (2007). Alternatives to traditional model comparison strategies for covariance structure models. In T. D. Little, J. A. Bovaird, & N. A. Card (Eds.), \emph{Modeling contextual effects in longitudinal studies} (pp. 33-62). Mahwah, NJ: Lawrence Erlbaum Associates.

Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA fit index. \emph{Structural Equation Modeling, 5,} 411-419.

Steiger, J. H., & Lind, J. C. (1980, June). \emph{Statistically based tests for the number of factors.} Paper presented at the annual meeting of the Psychometric Society, Iowa City, IA.
}
\author{
    Alexander M. Schoemann (East Carolina University; \email{schoemanna@ecu.edu})
	Kristopher J. Preacher (Vanderbilt University; \email{kris.preacher@vanderbilt.edu})
	Donna L. Coffman (Pennsylvania State University; \email{dlc30@psu.edu.})
}
\seealso{
	\itemize{
		\item \code{\link{plotRMSEAdist}} to visualize the RMSEA distributions
		\item \code{\link{findRMSEApower}} to find the statistical power based on population RMSEA given a sample size
		\item \code{\link{findRMSEAsamplesize}} to find the minium sample size for a given statistical power based on population RMSEA
	}
}
\examples{
plotRMSEApower(.025, .075, 23, 100, 500, 10)
}