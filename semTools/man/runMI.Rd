\name{runMI}
\alias{runMI}
\alias{cfa.mi}
\alias{sem.mi}
\alias{growth.mi}
\alias{lavaan.mi}
\title{
Multiply impute and analyze data using lavaan
}
\description{
This function takes data with missing observations, multiple imputes the data, runs a SEM using lavaan and combines the results using Rubin's rules. Note that parmeter estimates and standard errors are pooled by the Rubin's (1987) rule. The chi-square statistics and the related fit indices are pooled by the method described in \code{"chi"} argument. SRMR is calculated based on the average model-implied means and covariance matrices across imputations.
}
\usage{
runMI(model, data, m, miArgs=list(), chi="all", miPackage="Amelia", 
	seed=12345, fun, ...) 
cfa.mi(model, data, m, miArgs=list(), miPackage="Amelia", chi="all", 
	seed=12345, ...)
sem.mi(model, data, m, miArgs=list(), miPackage="Amelia", chi="all",  
	seed=12345, ...)
growth.mi(model, data, m, miArgs=list(), miPackage="Amelia", chi="all", 
	seed=12345, ...)
lavaan.mi(model, data, m, miArgs=list(), miPackage="Amelia", chi="all", 
	seed=12345, ...)
}
\arguments{
  \item{model}{
lavaan syntax for the the model to be analyzed. 
}
  \item{data}{
Data frame with missing observations or a list of data frames where each data frame is one imputed data set (for imputed data generated outside of the function). If a list of data frames is supplied, then other options can be left at the default.
}
  \item{m}{
Number of imputations wanted.
}
  \item{miArgs}{
Addition arguments for the multiple-imputation function. The arguments should be put in a list (see example below).
}
 \item{miPackage}{
Package to be used for imputation. Currently these functions only support \code{"Amelia"} or \code{"mice"} for imputation. 
}
 \item{chi}{
The method to combine the chi-square. Can be one of the following: \code{"MR"} for the method proposed for Meng & Rubin (1992), \code{"Mplus"} for the method used in Mplus (Asparouhov & Muthen, 2010), \code{"LMRR"} for the method proposed by Li, Meng, Raghunathan, & Rubin (1991), and \code{"all"} to show the three methods in the output. The default is \code{"all"}. 
}
 \item{seed}{
	Random number seed to be used in imputations.
}
 \item{fun}{
	The character of the function name used in running lavaan model (\code{"cfa"}, \code{"sem"}, \code{"growth"}, \code{"lavaan"}). 
}
 \item{...}{
	Other arguments to be passed to the specified lavaan function (\code{"cfa"}, \code{"sem"}, \code{"growth"}, \code{"lavaan"}).
}
}
\value{
	The \code{\linkS4class{lavaanStar}} object which contains the original \code{lavaan} object (where the appropriate parameter estimates, appropriate standard errors, and chi-squares are filled), the additional fit-index values of the null model, which need to be adjusted to multiple datasets, and the information from pooling multiple results.
}
\references{
Asparouhov T. & Muthen B. (2010).\emph{Chi-Square Statistics with Multiple Imputation}. Technical Report. www.statmodel.com. 

Li, K.H., Meng, X.-L., Raghunathan, T.E. and Rubin, D.B. (1991). Significance Levels From Repeated p-values with Multiply-Imputed Data. \emph{Statistica Sinica, 1}, 65-92.

Meng, X.L. & Rubin, D.B. (1992). Performing likelihood ratio tests with multiply-imputed data sets. \emph{Biometrika, 79}, 103 - 111.

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

out <- cfa.mi(HS.model, data=HSMiss, m = 3, chi="all")
summary(out)
inspect(out, "fit")
inspect(out, "impute")

\dontrun{
##Multiple group example
HSMiss2 <- cbind(HSMiss, school = HolzingerSwineford1939[,"school"])
out2 <- cfa.mi(HS.model, data=HSMiss2, m = 3, miArgs=list(noms="school"), chi="MR", group="school")
summary(out2)
inspect(out2, "fit")
inspect(out2, "impute")

##Example using previously imputed data with runMI
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

out3 <- runMI(mod, data=datsimMI$imputations, chi="LMRR", group="group", fun="cfa")
summary(out3)
inspect(out3, "fit")
inspect(out3, "impute")
}
}
