## Title: Compute more fit indices
## Author: Sunthud Pornprasertmanit <psunthud@ku.edu>, Aaron Boulton <aboulton@ku.edu>, Ruben Arslan <rubenarslan@gmail.com>, Terrence Jorgensen <TJorgensen314@gmail.com>
## Description: Calculations for promising alternative fit indices
##----------------------------------------------------------------------------##

moreFitIndices <- function(object, fit.measures = "all", nPrior = 1) {
  ## check for validity of user-specified "fit.measures" argument
  fit.choices <- c("gammaHat","adjGammaHat","baseline.rmsea","aic.smallN","bic.priorN","hqc",
                   "sic","gammaHat.scaled","adjGammaHat.scaled","baseline.rmsea.scaled")
  flags <- setdiff(fit.measures, c("all", fit.choices))
  if (length(flags)) stop(paste("Argument 'fit.measures' includes invalid options:",
                                paste(flags, collapse = ", "),
                                "Please choose 'all' or among the following:",
                                paste(fit.choices, collapse = ", "), sep = "\n"))
  if("all" %in% fit.measures) fit.measures <- fit.choices

  # Extract fit indices information from lavaan object
  fit <- inspect(object, "fit")
  # Get the number of variable
  p <- length(object@Data@ov.names[[1]])
  # Get the number of parameters
  nParam <- fit["npar"]
	
  # Get number of observations
  n <- object@SampleStats@ntotal
  # Find the number of groups
  ngroup <- object@Data@ngroups
  
  # Calculate -2*log(likelihood)
  f <- -2 * fit["logl"]
  
  # Compute fit indices
  result <- list()
  if (length(grep("gamma", fit.measures, ignore.case = TRUE))) {
    gammaHatValue <- p / (p + 2 * ((fit["chisq"] - fit["df"]) / (n - 1)))
    adjGammaHatValue <- 1 - (((ngroup * p * (p + 1)) / 2) / fit["df"]) * (1 - gammaHatValue)
    result["gammaHat"] <- gammaHatValue
    result["adjGammaHat"] <- adjGammaHatValue
    if(object@Options$test %in% c("satorra.bentler", "yuan.bentler")) {
      gammaHatScaledValue <- p / (p + 2 * ((fit["chisq.scaled"] - fit["df.scaled"]) / (n - 1)))
      adjGammaHatScaledValue <- 1 - (((ngroup * p * (p + 1)) / 2) / fit["df.scaled"]) * (1 - gammaHatScaledValue)
      result["gammaHat.scaled"] <- gammaHatScaledValue
      result["adjGammaHat.scaled"] <- adjGammaHatScaledValue
    }
  }
  if (length(grep("rmsea", fit.measures))) {
    result["baseline.rmsea"] <- nullRMSEA(object, silent = TRUE)
    if(object@Options$test %in% c("satorra.bentler", "yuan.bentler")) {
      result["baseline.rmsea.scaled"] <- nullRMSEA(object, scaled = TRUE, silent = TRUE)
    }
  }
  if(!is.na(f)) {
    if("aic.smallN" %in% fit.measures) result["aic.smallN"] <- f + (2 * nParam * (nParam + 1)) / (n - nParam - 1)
    if("bic.priorN" %in% fit.measures) result["bic.priorN"] <- f + log(1 + n/nPrior) * nParam
    if("hqc" %in% fit.measures) result["hqc"] <- f + 2 * log(log(n)) * nParam
    if("sic" %in% fit.measures) result["sic"] <- sic(f, object)
  }
  unlist(result[fit.measures])
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

## Stochastic Information Criterion
## f = minimized discrepancy function
## lresults = lavaan sem output object

sic <- function(f, lresults = NULL) {

  E.inv <- lavaan::lavTech(lresults, "inverted.information")
  if(inherits(E.inv, "try-error")) {
    return(as.numeric(NA))
  }
  E <- MASS::ginv(E.inv) * nobs(lresults)

  eigvals <- eigen(E, symmetric = TRUE, only.values = TRUE)$values
  # only positive ones
  eigvals <- eigvals[ eigvals > sqrt(.Machine$double.eps)]

  DET <- prod(eigvals)

  ## check singular
  if (DET <= 0) return(NA)

  ## return SIC
  f + log(DET)
}
