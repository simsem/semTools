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
	nfiValue <- (fit["baseline.chisq"] - fit["chisq"])/fit["baseline.chisq"]
	ifiValue <- (fit["baseline.chisq"] - fit["chisq"])/(fit["baseline.chisq"] - fit["df"])
	ciacValue <- f + (2 * nParam * (nParam + 1)) / (n - nParam - 1)
	ecviValue <- f + (2 * nParam / n)
	bicStarValue <- f + log(1 + n/nPrior) * nParam
	hqcValue <- f + 2 * log(log(n)) * nParam
	
	# Vector of result
	result <- c(nfiValue, ifiValue, gfiStarValue, agfiStarValue, ciacValue, ecviValue, bicStarValue, hqcValue)
	names(result) <- c("nfi", "ifi", "gfi*", "agfi*", "ciac", "ecvi", "bic*", "hqc")
	if(object@Options$test %in% c("satorra.bentler", "yuan.bentler")) {
		gfiStarScaledValue <- p / (p + 2 * ((fit["chisq.scaled"] - fit["df.scaled"]) / (n - 1)))
		agfiStarScaledValue <- 1 - (((ngroup * p * (p + 1)) / 2) / fit["df.scaled"]) * (1 - gfiStarValue)
		nfiScaledValue <- (fit["baseline.chisq.scaled"] - fit["chisq.scaled"])/fit["baseline.chisq.scaled"]
		ifiScaledValue <- (fit["baseline.chisq.scaled"] - fit["chisq.scaled"])/(fit["baseline.chisq.scaled"] - fit["df.scaled"])
		resultScaled <- c(gfiStarScaledValue, agfiStarScaledValue, nfiScaledValue, ifiScaledValue)
		names(resultScaled) <- c("nfi.scaled", "ifi.scaled", "gfi*.scaled", "agfi*.scaled")
		result <- c(result, resultScaled)
    } else {
		sicValue <- sic(f, object)
		result <- c(result, "sic" = sicValue)
	}
	
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
