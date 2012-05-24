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
 \item{...}{
	Other arguments to be passed to the imputation package
}
}
\value{
runMI returns a list with pooled fit indices, estimates, standard errors and fraction missing information.
\item{fit}{Pooled fit information. The first set of fit information are simply averaged across imputations and are not trustworthy. The second set of fit information, is a pooled Chi-square statistic based on Li, Meng, Raghunathan, & Rubin (1991) }
\item{parameters}{Pooled parameter estimates and standard errors. Wald statistics and p values are computed from the pooled estimates and standard errors. Also contains two estimates of Fraction of Missing Information (FMI). Includes asymptotic FMI (FMI.1) and FMI that is corrected for small numbers of imputation (FMI.2)}
}
\references{
Li, K.H., Meng, X.-L., Raghunathan, T.E. and Rubin, D.B. (1991). Significance Levels From Repeated p-values with Multiply-Imputed Data. Statistica Sinica, 1, 65-92.
Rubin, D.B. (1987) Multiple Imputation for Nonresponse in Surveys. J. Wiley & Sons, New York.
}
\author{Alexander M. Schoemann (University of Kansas;  \email{schoemann@ku.edu})
Patrick Miller(University of Kansas; \email{patr1ckm@ku.edu})
Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
Mijke Rhemtulla (University of Kansas; \email{mijke@ku.edu})
Alexander Robitzsch (UNIVERSITY?; \email{a.robitzsch@bifie.at})
}
\examples{
## No Example

}
