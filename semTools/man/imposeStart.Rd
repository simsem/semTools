\name{imposeStart}
\alias{imposeStart}
\title{
	Specify starting values from a lavaan output
}
\description{
This function will save the parameter estimates of a lavaan output and impose those parameter estimates as starting values for another analysis model. The free parameters with the same names or the same labels across two models will be imposed the new starting values. This function may help to increase the chance of convergence in a complex model (e.g., multitrait-multimethod model or complex longitudinal invariance model).
}
\usage{
imposeStart(out, expr, silent = TRUE)
}
\arguments{
  \item{out}{
	The \code{lavaan} output that users wish to use the parameter estimates as staring values for an analysis model
}
  \item{expr}{
	The original code that users use to run a lavaan model
}
  \item{silent}{
	Logical to print the parameter table with new starting values
}
}
\value{
	A fitted lavaan model
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
# The following example show that the longitudinal weak invariance model
# using effect coding was not convergent with three time points but convergent
# with two time points. Thus, the parameter estimates from the model with
# two time points are used as starting values of the three time points.
# The model with new starting values is convergent properly.

weak2time <- ' 
	# Loadings
	f1t1 =~ LOAD1*y1t1 + LOAD2*y2t1 + LOAD3*y3t1
    f1t2 =~ LOAD1*y1t2 + LOAD2*y2t2 + LOAD3*y3t2
	
	# Factor Variances
	f1t1 ~~ f1t1
	f1t2 ~~ f1t2
	
	# Factor Covariances
	f1t1 ~~ f1t2 
	
	# Error Variances
	y1t1 ~~ y1t1
	y2t1 ~~ y2t1
	y3t1 ~~ y3t1
	y1t2 ~~ y1t2
	y2t2 ~~ y2t2
	y3t2 ~~ y3t2
	
	# Error Covariances
	y1t1 ~~ y1t2 
	y2t1 ~~ y2t2 
	y3t1 ~~ y3t2
	
	# Factor Means
	f1t1 ~ NA*1
	f1t2 ~ NA*1
	
	# Measurement Intercepts
	y1t1 ~ INT1*1
	y2t1 ~ INT2*1
	y3t1 ~ INT3*1
	y1t2 ~ INT4*1
	y2t2 ~ INT5*1
	y3t2 ~ INT6*1
	
	# Constraints for Effect-coding Identification
	LOAD1 == 3 - LOAD2 - LOAD3
	INT1 == 0 - INT2 - INT3
	INT4 == 0 - INT5 - INT6
'
model2time <- lavaan(weak2time, data = exLong)

weak3time <- ' 
	# Loadings
	f1t1 =~ LOAD1*y1t1 + LOAD2*y2t1 + LOAD3*y3t1
    f1t2 =~ LOAD1*y1t2 + LOAD2*y2t2 + LOAD3*y3t2
    f1t3 =~ LOAD1*y1t3 + LOAD2*y2t3 + LOAD3*y3t3
	
	# Factor Variances
	f1t1 ~~ f1t1
	f1t2 ~~ f1t2
	f1t3 ~~ f1t3
	
	# Factor Covariances
	f1t1 ~~ f1t2 + f1t3
	f1t2 ~~ f1t3 
	
	# Error Variances
	y1t1 ~~ y1t1
	y2t1 ~~ y2t1
	y3t1 ~~ y3t1
	y1t2 ~~ y1t2
	y2t2 ~~ y2t2
	y3t2 ~~ y3t2
	y1t3 ~~ y1t3
	y2t3 ~~ y2t3
	y3t3 ~~ y3t3
	
	# Error Covariances
	y1t1 ~~ y1t2 
	y2t1 ~~ y2t2 
	y3t1 ~~ y3t2
	y1t1 ~~ y1t3 
	y2t1 ~~ y2t3 
	y3t1 ~~ y3t3
	y1t2 ~~ y1t3
	y2t2 ~~ y2t3 
	y3t2 ~~ y3t3
	
	# Factor Means
	f1t1 ~ NA*1
	f1t2 ~ NA*1
	f1t3 ~ NA*1
	
	# Measurement Intercepts
	y1t1 ~ INT1*1
	y2t1 ~ INT2*1
	y3t1 ~ INT3*1
	y1t2 ~ INT4*1
	y2t2 ~ INT5*1
	y3t2 ~ INT6*1
	y1t3 ~ INT7*1
	y2t3 ~ INT8*1
	y3t3 ~ INT9*1
	
	# Constraints for Effect-coding Identification
	LOAD1 == 3 - LOAD2 - LOAD3
	INT1 == 0 - INT2 - INT3
	INT4 == 0 - INT5 - INT6
	INT7 == 0 - INT8 - INT9
'
### The following command does not provide convergent result
# model3time <- lavaan(weak3time, data = exLong)

### Use starting values from the model with two time points
model3time <- imposeStart(model2time, lavaan(weak3time, data = exLong))
summary(model3time)
}
