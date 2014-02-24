\name{fmi}
\alias{fmi}
\title{
Fraction of Missing Information.
}
\description{
	This function takes a list of imputed data sets and estimates the Fraction of Missing Information of the Variances and Means for each variable.
}
\usage{
fmi(dat.imp, method="saturated", varnames=NULL, group=NULL, exclude=NULL, 
digits=3)
}
\arguments{
  \item{dat.imp}{
	List of imputed data sets, the function only accept a list of data frames.
}
  \item{method}{
	Specified the model used to estimated the variances and means. Can be one of the following: \code{"saturated"} (\code{"sat"}) or \code{"null"}, the default is \code{"saturated"}. See Details for more information.
}
 \item{varnames}{
	A vector of variables names. This argument allow the user to get the fmi of a subset of variables. The function by default will estimate the fmi for all the variables.
}
 \item{group}{
  A variable name defining the groups. This will give the fmi for each group. 
}
 \item{exclude}{
  A vector of variables names. These variables will be excluded from the analysis.
}
 \item{digits}{
	Number of decimals to print in the results.
}
}
\details{
The function estimates a variance/covariance model for each data set using lavaan. If method = \code{"saturated"} the function will estimate all the variances and covariances, 
if method = \code{"null"} the function will only estimate the variances. The saturated model gives more reliable estimates. 
With big data sets using the saturated model could take a lot of time.
In the case of having problems with big data sets it is helpful to select a subset of variables with \code{varnames} and/or use the \code{"null"} model.
The function does not accept character variables.
}
\value{
fmi returns a list with the Fraction of Missing Information of the Variances and Means for each variable in the data set.
\item{Variances}{The estimated variance for each variable, and the respective standard error.
Two estimates Fraction of Missing Information of the Variances. The first estimate of fmi (fmi.1) is asymptotic fmi and the second estimate of fmi (fmi.2) is corrected for small numbers of imputations}
\item{Means}{The estimated mean for each variable, and the respective standard error.
Two estimates Fraction of Missing Information of the Means. The first estimate of fmi (fmi.1) is asymptotic fmi and the second estimate of fmi (fmi.2) is corrected for small numbers of imputations}
}
\references{
Rubin, D.B. (1987) \emph{Multiple Imputation for Nonresponse in Surveys.} J. Wiley & Sons, New York.

Savalei, V. & Rhemtulla, M. (2012) On Obtaining Estimates of the Fraction
of Missing Information From Full Information Maximum Likelihood, \emph{Structural Equation Modeling: A Multidisciplinary Journal, 19:3}, 477-494.

Wagner, J. (2010) The Fraction of Missing Information as a Tool for Monitoring the Quality of Survey Data, \emph{Public Opinion Quarterly, 74:2}, 223-243.
}
\author{Mauricio Garnier Villarreal (University of Kansas; \email{mgv@ku.edu})
}
\examples{
library(Amelia)
library(lavaan)

modsim <- '
f1 =~ 0.7*y1+0.7*y2+0.7*y3
f2 =~ 0.7*y4+0.7*y5+0.7*y6
f3 =~ 0.7*y7+0.7*y8+0.7*y9'

datsim <- simulateData(modsim,model.type="cfa", meanstructure=TRUE, 
                       std.lv=TRUE, sample.nobs=c(200,200))
randomMiss2 <- rbinom(prod(dim(datsim)), 1, 0.1)
randomMiss2 <- matrix(as.logical(randomMiss2), nrow=nrow(datsim))
randomMiss2[,10] <- FALSE
datsim[randomMiss2] <- NA
datsimMI <- amelia(datsim,m=3,idvars="group")

out1 <- fmi(datsimMI$imputations, exclude="group")
out1
                       
out2 <- fmi(datsimMI$imputations, exclude="group", method="null")
out2
                       
out3 <- fmi(datsimMI$imputations, varnames=c("y1","y2","y3","y4"))
out3

out4 <- fmi(datsimMI$imputations, group="group")
out4

}
