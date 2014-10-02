## Title: Compute more fit indices
## Author: Sunthud Pornprasertmanit <psunthud@ku.edu>, Aaron Boulton <aboulton@ku.edu>, Ruben Arslan <rubenarslan@gmail.com>, Terrence Jorgensen <tdj@ku.edu>
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
	n <- object@SampleStats@ntotal
	
	# Calculate the minimized discrepancy function
	f <- -2 * fit["logl"]
	# Find the number of groups
	ngroup <- object@Data@ngroups
	
	# Compute fit indices
	gammaHatValue <- p / (p + 2 * ((fit["chisq"] - fit["df"]) / (n - 1)))
	adjGammaHatValue <- 1 - (((ngroup * p * (p + 1)) / 2) / fit["df"]) * (1 - gammaHatValue)
	nullRmseaValue <- nullRMSEA(object, silent = TRUE)
	result <- c(gammaHatValue, adjGammaHatValue, nullRmseaValue)
	names(result) <- c("gammaHat", "adjGammaHat", "baseline.rmsea")

	if(!is.na(f)) {
		aiccValue <- f + (2 * nParam * (nParam + 1)) / (n - nParam - 1)
		bicStarValue <- f + log(1 + n/nPrior) * nParam
		hqcValue <- f + 2 * log(log(n)) * nParam
		temp <- c(aiccValue, bicStarValue, hqcValue)
		names(temp) <- c("aic.smallN", "bic.priorN", "hqc")
		result <- c(result, temp)
	}
	
	# Vector of result
	if(object@Options$test %in% c("satorra.bentler", "yuan.bentler")) {
		gammaHatScaledValue <- p / (p + 2 * ((fit["chisq.scaled"] - fit["df.scaled"]) / (n - 1)))
		adjGammaHatScaledValue <- 1 - (((ngroup * p * (p + 1)) / 2) / fit["df.scaled"]) * (1 - gammaHatScaledValue)
		nullRmseaScaledValue <- nullRMSEA(object, scaled = TRUE, silent = TRUE)
		resultScaled <- c(gammaHatScaledValue, adjGammaHatScaledValue, nullRmseaScaledValue)
		names(resultScaled) <- c("gammaHat.scaled", "adjGammaHat.scaled", "baseline.rmsea.scaled")
		result <- c(result, resultScaled)
    } else {
		if(!is.na(f)) {
			sicValue <- sic(f, object)
			result <- c(result, "sic" = sicValue)
		}
	}
	
	return(result)
}

## Stochastic Information Criterion
## f = minimized discrepancy function
## lresults = lavaan sem output object

sic <- function(f, lresults = NULL) {
	expinf <- NA
	v <- NA
	try(v <- vcov(lresults), silent = TRUE)
	ifelse(is.na(v) || det(v) == 0, return(NA), try(expinf <- solve(v) / lresults@SampleStats@ntotal, silent = TRUE))
	sic <- as.numeric(f + log(det(lresults@SampleStats@ntotal * (expinf))))/2
	return(sic)
}


nullRMSEA <- function (object, scaled = FALSE, silent = FALSE) { 
	# return RMSEA of the null model, warn if it is lower than 0.158, because it makes the TLI/CLI hard to interpret
	test <- object@Options$test 
	
	fits <- lavaan::fitMeasures(object)
	N <- object@SampleStats@ntotal # sample size
	
	X2 <- as.numeric ( fits['baseline.chisq'] ) # get baseline chisq
	df <- as.numeric ( fits['baseline.df'] ) # get baseline df 
	G <- object@Data@ngroups # number of groups
	
	### a simple rip from fit.measures.R in lavaan's codebase.
	N.RMSEA <- max(N, X2*4) # Check with lavaan
        # RMSEA
	if(df > 0) {
		if(scaled) {
			d <- sum(object@Fit@test[[2]]$trace.UGamma)
		} 
		if(object@Options$mimic %in% c("Mplus", "lavaan")) {
			GG <- 0
			RMSEA <- sqrt( max( c((X2/N)/df - 1/(N-GG), 0) ) ) * sqrt(G)
			if(scaled && test != "scaled.shifted") {
				RMSEA.scaled <- 
					 sqrt( max( c((X2/N)/d - 1/(N-GG), 0) ) ) * sqrt(G)
			} else if(test == "scaled.shifted") {
				RMSEA.scaled <-
					 sqrt( max(c((as.numeric(fits["baseline.chisq.scaled"])/N)/df - 1/(N-GG), 0))) * sqrt(G)
			}
		} else {
			RMSEA <- sqrt( max( c((X2/N)/df - 1/N, 0) ) )
			if(scaled) {
				RMSEA.scaled <- sqrt( max( c((X2/N)/d - 1/N, 0) ) )
			}
		}
	} else {
		RMSEA <- RMSEA.scaled <- 0
	}
	if(scaled) {
		RMSEA <- RMSEA.scaled
	}
	if(!silent) {
		if(RMSEA < 0.158 ) { 
			cat(paste0("TLI and other incremental fit indices may not be that informative, because the RMSEA of the baseline model is lower than 0.158 (Kenny, Kaniskan, & McCoach, 2011). The baseline RMSEA is ",round(RMSEA,3), "\n"))
		} else {
			cat(paste0("Baseline RMSEA: ",round(RMSEA,3), "\n"))
		}
	}
	invisible(RMSEA)
}

