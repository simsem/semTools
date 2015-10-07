###Function to do power analysis for parameters with Satorra & Sarris method
###Alexander M. Schoemann
###11/4/2014


##Steps:
##1. Specify model (use lavaan syntax based on simulateData)
##2. get model implied covariance matrix
##3. Fit model with parameter constrained to 0 (or take a model specification for multiparameter tests?)
##4. Use chi square from step 3 as non-centrality parameter to get power.


	
##Function to return power for a given model parameter
#inputs: popModel = lavaan syntax specifying data generating model (be sure to provide a value for all parameters), n =sample size (either scalar or vector), powerModel=Model to be fit, with parameter of interest fixed to 0, fun = lavaan function to use, nparam = number of parameters fixed in the power Model, ... additional arguments to pass to lavaan
SSpower <- function(popModel, n, powerModel, fun = "cfa", nparam = 1, alpha = .05, ...) {
##Two item list, first item is covariance matrix, second item is mean vector
popCov <- lavaan::fitted(do.call(fun, list(model=popModel)))

##Fit model with parameter(s) fixed to 0
out <- list(model=powerModel, sample.cov=popCov[[1]], sample.mean = popCov[[2]], sample.nobs=n)
out <- c(out, list(...))
mod <- do.call(fun, out)

##get NCP from chi square
ncp <- lavaan::fitmeasures(mod)["chisq"]
critVal <- qchisq(1-alpha, nparam)

1-pchisq(critVal, nparam, ncp)
 
}

#Test the function
# model <- '
         # f1 =~ .7?V1 + .7?V2 + .7?V3 + .7?V4
         # f2 =~ .7?V5 + .7?V6 + .7?V7 + .7?V8
         
         # f1 ~~ .3?f2
         # f1 ~~ 1*f1
         # f2 ~~ 1*f2
         
         # V1 ~~ .51?V1
         # V2 ~~ .51?V2
         # V3 ~~ .51?V3
         # V4 ~~ .51?V4
         # V5 ~~ .51?V5
         # V6 ~~ .51?V6
         # V7 ~~ .51?V7
         # V8 ~~ .51?V8  
         # '


# model2 <- '
         # f1 =~ .7?V1 + .7?V2 + .7?V3 + .7?V4
         # f2 =~ .7?V5 + .7?V6 + .7?V7 + .7?V8
         
         # f1 ~~ 0*f2
         # f1 ~~ 1*f1
         # f2 ~~ 1*f2
         
         # V1 ~~ .51?V1
         # V2 ~~ .51?V2
         # V3 ~~ .51?V3
         # V4 ~~ .51?V4
         # V5 ~~ .51?V5
         # V6 ~~ .51?V6
         # V7 ~~ .51?V7
         # V8 ~~ .51?V8  
         # '

	
# SSpower(model, 150, model2)

#Get power for a range of values

# powVals <- NULL
# Ns <- seq(120, 500, 10)
# for(i in Ns){
# powVals <- c(powVals, SSpower(model, i, model2))
# }
# plot(Ns, powVals, type = 'l')

#Test with multiple params
# model3 <- '
         # f1 =~ 1*V1 + 1*V2 + 1*V3 + 1*?V4
         # f2 =~ .7?V5 + .7?V6 + .7?V7 + .7?V8
         
         # f1 ~~ f2
         # f1 ~~ 1*f1
         # f2 ~~ 1*f2
         
         # V1 ~~ .51?V1
         # V2 ~~ .51?V2
         # V3 ~~ .51?V3
         # V4 ~~ .51?V4
         # V5 ~~ .51?V5
         # V6 ~~ .51?V6
         # V7 ~~ .51?V7
         # V8 ~~ .51?V8  
         # '
# SSpower(model, 150, model3, nparam=4)
