\name{runMI}
\alias{runMI}
\title{
Multiply impute and analyze data using lavaan
}
\description{
This function takes data with missing observations, multiple imputes the data, runs a SEM using lavaan and combines the results using Rubin's rules.
}
\usage{
runMI(data.mat, data.model, m, miPackage="Amelia", digits = 3, ...)
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
	lavaan option. The estimator to be used. Can be one of the following: "ML" for maximum likelihood, "GLS" for generalized least squares, "WLS" for weighted least squares (sometimes called ADF estimation), "MLM" for maximum likelihood estimation with robust standard errors and a Satorra-Bentler scaled test statistic, "MLF" for maximum likelihood estimation with standard errors based on first-order derivatives and a conventional test statistic, "MLR" for maximum likelihood estimation with robust ‘Huber-White’ standard errors and a scaled test statistic which is asymptotically equivalent to the Yuan-Bentler T2-star test statistic. Note that the "MLM", "MLF" and "MLR" choices only affect the standard errors and the test statistic. 
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
Li, K.H., Meng, X.-L., Raghunathan, T.E. and Rubin, D.B. (1991). Significance Levels From Repeated p-values with Multiply-Imputed Data. Statistica Sinica, 1, 65-92.
Rubin, D.B. (1987) Multiple Imputation for Nonresponse in Surveys. J. Wiley & Sons, New York.
}
\author{Alexander M. Schoemann (University of Kansas;  \email{schoemann@ku.edu})
Patrick Miller(University of Kansas; \email{patr1ckm@ku.edu})
Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
Mijke Rhemtulla (University of Kansas; \email{mijke@ku.edu})
Alexander Robitzsch (Federal Institute for Education Research, Innovation, and
Development of the Austrian School System, Salzburg, Austria; \email{a.robitzsch@bifie.at})
}
\examples{
## No Example

}
