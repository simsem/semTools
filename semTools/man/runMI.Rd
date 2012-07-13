\name{runMI}
\alias{runMI}
\title{
Multiply impute and analyze data using lavaan
}
\description{
This function takes data with missing observations, multiple imputes the data, runs a SEM using lavaan and combines the results using Rubin's rules.
}
\usage{
runMI(data.mat,data.model, m, miPackage="Amelia", digits=3, seed=12345, 
    std.lv = FALSE, estimator = "ML", group = NULL, group.equal = "", ...)
}
\arguments{
  \item{data.mat}{
Data frame with missing observations or a list of data frames where each data frame is one imputed data set (for imputed data generated outside of the function). If a list of data frames is supplied, then other options can be left at the default.
}
  \item{data.model}{
lavaan syntax for the the model to be analyzed. 
}
  \item{m}{
Number of imputations wanted.
}
 \item{miPackage}{
Package to be used for imputation. Currently runMI only uses Amelia or mice for imputation. 
}
 \item{digits}{
	Number of digits to print in the results.
}
 \item{seed}{
	Random number seed to be used in imputations.
}
 \item{std.lv}{
	lavaan option. If TRUE, the metric of each latent variable is determined by fixing their variances to 1.0. If FALSE, the metric of each latent variable is determined by fixing the factor loading of the first indicator to 1.0.
}
 \item{estimator}{
	lavaan option. The estimator to be used. Can be one of the following: "ML" for maximum likelihood, "GLS" for generalized least squares, "WLS" for weighted least squares (sometimes called ADF estimation), "MLM" for maximum likelihood estimation with robust standard errors and a Satorra-Bentler scaled test statistic, "MLF" for maximum likelihood estimation with standard errors based on first-order derivatives and a conventional test statistic, "MLR" for maximum likelihood estimation with robust 'Huber-White' standard errors and a scaled test statistic which is asymptotically equivalent to the Yuan-Bentler T2-star test statistic. Note that the "MLM", "MLF" and "MLR" choices only affect the standard errors and the test statistic. 
}
 \item{group}{
	lavaan option. A variable name in the data frame defining the groups in a multiple group analysis.
}
 \item{group.equal}{
	lavaan option. A vector of character strings. Only used in a multiple group analysis. Can be one or more of the following: "loadings", "intercepts", "means", "regressions", "residuals", "residual.covariances", "lv.variances" or "lv.covariances", specifying the pattern of equality constraints across multiple groups.
}
 \item{...}{
	Other arguments to be passed to the imputation package
}
}
\value{
runMI returns a list with pooled fit indices, estimates, standard errors and fraction missing information.
\item{fit}{Pooled fit information. The first set of fit information are simply averaged across imputations and are not trustworthy. The second set of fit information, is a pooled Chi-square statistic based on Li, Meng, Raghunathan, & Rubin (1991) }
\item{parameters}{Pooled parameter estimates and standard errors. Wald statistics and p values are computed from the pooled estimates and standard errors. Also contains two estimates of Fraction of Missing Information (FMI). The first estimate of FMI (FMI.1) is asymptotic FMI and the second estimate of FMI (FMI.2) is corrected for small numbers of imputation}
}
\references{
Li, K.H., Meng, X.-L., Raghunathan, T.E. and Rubin, D.B. (1991). Significance Levels From Repeated p-values with Multiply-Imputed Data. \emph{Statistica Sinica, 1}, 65-92.

Rubin, D.B. (1987) \emph{Multiple Imputation for Nonresponse in Surveys.} J. Wiley & Sons, New York.
}
\author{Alexander M. Schoemann (University of Kansas;  \email{schoemann@ku.edu})
Patrick Miller (University of Kansas; \email{patr1ckm@ku.edu})
Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
Mijke Rhemtulla (University of Kansas; \email{mijke@ku.edu})
Alexander Robitzsch (Federal Institute for Education Research, Innovation, and Development of the Austrian School System, Salzburg, Austria; \email{a.robitzsch@bifie.at})
Craig Enders (Arizona State University; \email{Craig.Enders@asu.edu})
Mauricio Garnier Villarreal (University of Kansas; \email{mgv@ku.edu})
}
\examples{

library(lavaan)

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

HSMiss <- HolzingerSwineford1939[,paste("x", 1:9, sep="")]
randomMiss <- rbinom(prod(dim(HSMiss)), 1, 0.1)
randomMiss <- matrix(as.logical(randomMiss), nrow=nrow(HSMiss))
HSMiss[randomMiss] <- NA

out <- runMI(HSMiss, HS.model, m = 3)

HSMiss2 <- cbind(HSMiss, school = HolzingerSwineford1939[,"school"])
out2 <- runMI(HSMiss2, HS.model, m = 3, group="school", noms="school")

library(Amelia)

modsim <- '
f1 =~ 0.7*y1+0.7*y2+0.7*y3
f2 =~ 0.7*y4+0.7*y5+0.7*y6
f3 =~ 0.7*y7+0.7*y8+0.7*y9'

mod <- '
f1 =~ y1+y2+y3
f2 =~ y4+y5+y6
f3 =~ y7+y8+y9'

datsim <- simulateData(modsim,model.type="cfa", meanstructure=TRUE, 
	std.lv=TRUE, sample.nobs=c(200,200))
randomMiss2 <- rbinom(prod(dim(datsim)), 1, 0.1)
randomMiss2 <- matrix(as.logical(randomMiss2), nrow=nrow(datsim))
datsim[randomMiss2] <- NA
datsimMI <- amelia(datsim,m=3, noms="group")

out3 <- runMI(datsimMI$imputations, mod, group="group")

}
