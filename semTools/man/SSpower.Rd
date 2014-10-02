\name{SSpower}
\alias{SSpower}
\title{
	Power for model parameters 
}
\description{
Determines power for model parameters using the Satorra & Sarris (1985) method
}
\usage{
SSpower(popModel, n, powerModel, fun = "cfa", nparam = 1, alpha = .05, ...)
}
\arguments{
\item{popModel}{
    lavaan syntax for the population model. This model should specify population  values for all paramters in the model.
}
\item{n}{
    Sample size used in power calculation 
}
\item{powerModel}{
    lavaan syntax for the model to be analyzed. This syntax should have the parameter(s) of interest fixed to 0 (or some other number).
}
  \item{fun}{
	The character of the function name used in running lavaan model (\code{"cfa"}, \code{"sem"}, \code{"growth"}, \code{"lavaan"}). 
}
\item{nparam}{
    The number of parameters one is constrained in \code{powerModel}.
}
\item{alpha}{
    The Type I error rate used to  assess power
}
 \item{...}{
	Other arguments to be passed to the specified lavaan function (\code{"cfa"}, \code{"sem"}, \code{"growth"}, \code{"lavaan"}).
}
}

\author{
    Alexander M. Schoemann (East Carolina University; \email{schoemanna@ecu.edu})
}
\references{
Satorra, A., & Saris, W. E. (1985). Power of the likelihood ratio test in covariance structure analysis. \emph{Psychometrika, 50}, 83-90.}

\examples{
library(lavaan)

#Specify population values. Note every paramter has a fixed value
modelP <- '
         f1 =~ .7*V1 + .7*V2 + .7*V3 + .7*V4
         f2 =~ .7*V5 + .7*V6 + .7*V7 + .7*V8
         
         f1 ~~ .3*f2
         f1 ~~ 1*f1
         f2 ~~ 1*f2
         
         V1 ~~ .51*V1
         V2 ~~ .51*V2
         V3 ~~ .51*V3
         V4 ~~ .51*V4
         V5 ~~ .51*V5
         V6 ~~ .51*V6
         V7 ~~ .51*V7
         V8 ~~ .51*V8  
         '

#Specify model to be analyzed. Note parameter of interest f1~~f2 is fixed to 0.
modelA <- '
         f1 =~ V1 + V2 + V3 + V4
         f2 =~ V5 + V6 + V7 + V8
         
         f1 ~~ 0*f2

  
         '

	
SSpower(modelP, 150, modelA, std.lv=TRUE)

##Get power for a range of values

Ns <- seq(100, 500, 40)
powVals <- rep(NA, length(Ns))
for(i in 1:length(Ns)){
powVals[i] <- SSpower(modelP, Ns[i], modelA)
}
plot(Ns, powVals, type = 'l')

}

