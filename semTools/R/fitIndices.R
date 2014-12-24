## Title: Compute more fit indices
## Author: Sunthud Pornprasertmanit <psunthud@ku.edu>, Aaron Boulton <aboulton@ku.edu>, Ruben Arslan <rubenarslan@gmail.com>, Terrence Jorgensen <tdj@ku.edu>
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
  ## calculate Jacobian
  R <- lavaan::lav_func_jacobian_complex(func = lresults@Model@ceq.function, x = lresults@Fit@x)
  
  ## FIXME: CRAN does not allow hidden functions. Yves can include a new method
  ## in lavaan 0.5-18, inspect(lresults, "information")
  ## calculate Fisher Information Matrix
  if(lresults@Options$estimator == "ML" && lresults@Options$mimic == "Mplus") {
    E <- computeExpectedInformationMLM(lresults@Model, lavsamplestats = lresults@SampleStats)
  } else {
    E <- computeExpectedInformation(lresults@Model, lavsamplestats = lresults@SampleStats)
  }
  
  ## calculate Fisher Information Matrix of only the non-redundant parameters
  QR <- qr(t(R))
  ranK <- QR$rank
  Q <- qr.Q(QR, complete = TRUE)
  Q1 <- Q[ , 1:ranK, drop = FALSE]         # range space
  Q2 <- Q[ , -seq_len(ranK), drop = FALSE] # null space
  FIM <- t(Q2) %*% E %*% Q2
  
  ## check that the non-redundant Fisher Information Matrix is not singular
  if (det(FIM) <= 0) return(NA)
  
  ## return SIC
  f + log(det(lresults@SampleStats@ntotal * FIM))
}

## functions to compute model information, borrowed from lavaan source code (lav_information.R)
computeExpectedInformation <- function(lavmodel       = NULL, 
                                       lavsamplestats = NULL, 
                                       lavdata        = NULL,
                                       estimator      = "ML",
                                       # is no Delta is provided, we compute 
                                       # Delta for the free parameters only
                                       Delta = computeDelta(lavmodel=lavmodel), 
                                       lavcache       = NULL,
                                       extra          = FALSE) {
  
  # compute/get WLS.V
  # if DWLS or ULS, this is the diagonal only! (since 0.5-17)
  WLS.V <- lav_model_wls_v(lavmodel       = lavmodel, 
                           lavsamplestats = lavsamplestats,
                           estimator      = estimator,
                           lavdata        = lavdata)
  
  # compute Information per group
  Info.group  <- vector("list", length=lavsamplestats@ngroups)
  for(g in 1:lavsamplestats@ngroups) {
    # note LISREL documentation suggest (Ng - 1) instead of Ng...
    fg <- lavsamplestats@nobs[[g]]/lavsamplestats@ntotal
    # compute information for this group
    if(estimator %in% c("DWLS", "ULS")) {
      # diagonal weight matrix
      Delta2 <- sqrt(WLS.V[[g]]) * Delta[[g]]
      Info.group[[g]] <- fg * crossprod(Delta2)
    } else {
      # full weight matrix
      # Info.group[[g]] <- 
      # fg * (t(Delta[[g]]) %*% WLS.V[[g]] %*% Delta[[g]])
      Info.group[[g]] <- 
        fg * ( crossprod(Delta[[g]], WLS.V[[g]]) %*% Delta[[g]] )
    }
  }
  
  # assemble over groups
  Information <- Info.group[[1]]
  if(lavsamplestats@ngroups > 1) {
    for(g in 2:lavsamplestats@ngroups) {
      Information <- Information + Info.group[[g]]
    }
  }
  
  if(extra) {
    attr(Information, "Delta") <- Delta
    attr(Information, "WLS.V") <- WLS.V # unweighted
  }
  
  Information
}

# only for Mplus MLM
computeExpectedInformationMLM <- function(lavmodel = NULL, 
                                          lavsamplestats = NULL, 
                                          Delta = computeDelta(lavmodel = 
                                                                 lavmodel)) {
  
  # compute WLS.V 
  WLS.V <- vector("list", length=lavsamplestats@ngroups)
  if(lavmodel@group.w.free) {
    GW <- unlist(computeGW(lavmodel = lavmodel))
  }
  for(g in 1:lavsamplestats@ngroups) {
    WLS.V[[g]] <- compute.A1.sample(lavsamplestats=lavsamplestats, group=g,
                                    meanstructure=TRUE, 
                                    information="expected")
    # the same as GLS... (except for the N/N-1 scaling)
    if(lavmodel@group.w.free) {
      # unweight!!
      a <- exp(GW[g]) / lavsamplestats@nobs[[g]]
      # a <- exp(GW[g]) * lavsamplestats@ntotal / lavsamplestats@nobs[[g]]
      WLS.V[[g]] <- bdiag( matrix(a,1,1), WLS.V[[g]])
    }
  }
  
  # compute Information per group
  Info.group  <- vector("list", length=lavsamplestats@ngroups)
  for(g in 1:lavsamplestats@ngroups) {
    fg <- lavsamplestats@nobs[[g]]/lavsamplestats@ntotal
    # compute information for this group
    Info.group[[g]] <- fg * (t(Delta[[g]]) %*% WLS.V[[g]] %*% Delta[[g]])
  }
  
  # assemble over groups
  Information <- Info.group[[1]]
  if(lavsamplestats@ngroups > 1) {
    for(g in 2:lavsamplestats@ngroups) {
      Information <- Information + Info.group[[g]]
    }
  }
  
  # always
  attr(Information, "Delta") <- Delta
  attr(Information, "WLS.V") <- WLS.V # unweighted
  
  Information
}
