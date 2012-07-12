## Title: Compute more fit indices
## Author: Sunthud Pornprasertmanit <psunthud@ku.edu>, Aaron Boulton <aboulton@ku.edu>
## Description: Calculations for promising alternative fit indices
##----------------------------------------------------------------------------##

moreFitIndices <- function(object, nPrior = 1) {
	# Extract fit indices information from lavaan object
	fit <- inspect(object, "fit")
	# Get the number of variable
	p <- length(fitted(object)$mean)
	
	# Get the number of parameters
	nParam <- fit["baseline.df"] - fit["df"]
	
	# Get number of observations
	n <- fit["ntotal"]
	
	# Calculate the minimized discrepancy function
	f <- -2 * fit["logl"]
	
	# Find the number of groups
	ngroup <- object@Data@ngroups
	
	# Compute fit indices
	gfiStarValue <- p / (p + 2 * ((fit["chisq"] - fit["df"]) / (n - 1)))
	agfiStarValue <- 1 - (((ngroup * p * (p + 1)) / 2) / fit["df"]) * (1 - gfiStarValue)
	sicValue <- NA
	estSpec <- as.list(object@call)$estimator
	if(!(!is.null(estSpec) && (estSpec %in% c("mlr", "mlm", "mlf")))) try(sicValue <- sic(f, object), silent=TRUE)
	nfiValue <- (fit["baseline.chisq"] - fit["chisq"])/fit["baseline.chisq"]
	ifiValue <- (fit["baseline.chisq"] - fit["chisq"])/(fit["baseline.chisq"] - fit["df"])
	ciacValue <- f + (2 * nParam * (nParam + 1)) / (n - nParam - 1)
	ecviValue <- f + (2 * nParam / n)
	bicStarValue <- f + log(1 + n/nPrior) * nParam
	hqcValue <- f + 2 * log(log(n)) * nParam
	
	# Vector of result
	result <- c(nfiValue, ifiValue, gfiStarValue, agfiStarValue, ciacValue, ecviValue, sicValue, bicStarValue, hqcValue)
	names(result) <- c("nfi", "ifi", "gfi*", "agfi*", "ciac", "ecvi", "sic", "bic*", "hqc")
	return(result)
}

## Stochastic Information Criterion
## f = minimized discrepancy function
## lresults = lavaan sem output object

sic <- function(f, lresults = NULL) {
	expinf <- NA
	v <- NA
	try(v <- vcov(lresults), silent=TRUE)
	ifelse(is.na(v), return(NA), try(expinf <- solve(v) / lresults@SampleStats@ntotal, silent=TRUE))
	sic <- as.numeric(f + log(det(lresults@SampleStats@ntotal * (expinf))))/2
	return(sic)
}
