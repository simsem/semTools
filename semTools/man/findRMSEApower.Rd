% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerAnalysisRMSEA.R
\name{findRMSEApower}
\alias{findRMSEApower}
\title{Find the statistical power based on population RMSEA}
\usage{
findRMSEApower(rmsea0, rmseaA, df, n, alpha = 0.05, group = 1)
}
\arguments{
\item{rmsea0}{Null RMSEA}

\item{rmseaA}{Alternative RMSEA}

\item{df}{Model degrees of freedom}

\item{n}{Sample size of a dataset}

\item{alpha}{Alpha level used in power calculations}

\item{group}{The number of group that is used to calculate RMSEA.}
}
\description{
Find the proportion of the samples from the sampling distribution of RMSEA
in the alternative hypothesis rejected by the cutoff dervied from the
sampling distribution of RMSEA in the null hypothesis. This function can be
applied for both test of close fit and test of not-close fit (MacCallum,
Browne, & Suguwara, 1996)
}
\details{
This function find the proportion of sampling distribution derived from the
alternative RMSEA that is in the critical region derived from the sampling
distribution of the null RMSEA. If \code{rmseaA} is greater than
\code{rmsea0}, the test of close fit is used and the critical region is in
the right hand side of the null sampling distribution. On the other hand, if
\code{rmseaA} is less than \code{rmsea0}, the test of not-close fit is used
and the critical region is in the left hand side of the null sampling
distribution (MacCallum, Browne, & Suguwara, 1996).

There is also a Shiny app called "power4SEM" that provides a graphical user
interface for this functionality (Jak et al., in press).  It can be accessed
at \url{https://sjak.shinyapps.io/power4SEM/}.
}
\examples{

findRMSEApower(rmsea0 = .05, rmseaA = .08, df = 20, n = 200)

}
\references{
MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis
and determination of sample size for covariance structure modeling.
\emph{Psychological Methods, 1}(2), 130--149. \doi{10.1037/1082-989X.1.2.130}

Jak, S., Jorgensen, T. D., Verdam, M. G., Oort, F. J., & Elffers, L.
(2021). Analytical power calculations for structural equation modeling:
A tutorial and Shiny app. \emph{Behavior Research Methods, 53}, 1385--1406.
\doi{10.3758/s13428-020-01479-0}
}
\seealso{
\itemize{
\item \code{\link[=plotRMSEApower]{plotRMSEApower()}} to plot the statistical power based on
population RMSEA given the sample size
\item \code{\link[=plotRMSEAdist]{plotRMSEAdist()}} to visualize the RMSEA distributions
\item \code{\link[=findRMSEAsamplesize]{findRMSEAsamplesize()}} to find the minium sample size for
a given statistical power based on population RMSEA
}
}
\author{
Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
