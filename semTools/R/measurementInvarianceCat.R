measurementInvarianceCat <- function(..., std.lv = FALSE, strict=FALSE, quiet=FALSE) {
	List <- list(...)
	if(!is.null(List$parameterization) && tolower(List$parameterization) != "theta") warnings("The parameterization is set to 'theta' by default.")
	List$parameterization <- "theta"

	# Find the number of groups
	ngroups <- 1
	if(!is.null(List$group)) {
		group <- List$group
		ngroups <- length(setdiff(unique(List$data[,group]), NA))
	} else {
		stop("Please specify the group variable")
	}
	
	# Get the lavaan parameter table
	template <- do.call("cfa", c(List, do.fit=FALSE))
    lavaanParTable <- parTable(template)
	
	# Check whether all variables are categorical
	sampstat <- inspect(template, "samp")[[1]]
	meanname <- names(sampstat$mean)
	thname <- names(sampstat$th)
	if(any(is.na(charmatch(meanname, thname)))) stop("Some variables in your model are not identified as categorical.")
	
	varList <- lavaanParTable$rhs[lavaanParTable$op == "=~"]
	facName <- lavaanParTable$lhs[(lavaanParTable$op == "=~") & (lavaanParTable$rhs %in% varList)]
	if(length(unique(sapply(facName, function(x) length(x)))) > 1) stop("The numbers of variables in each element are not equal.")
	varList <- unique(varList)
	facName <- unique(facName)
	
	# Check whether the factor configuration is the same across gorups
	groupParTable <- split(lavaanParTable, lavaanParTable$group)
	group1pt <- groupParTable[[1]]
	groupParTable <- lapply(groupParTable, "[", c("lhs", "op", "rhs"))
	if(!multipleAllEqualList(lapply(groupParTable, function(x) sapply(x, "[", x$op == "=~")))) stop("Factor configuration is not the same across groups")
	
	# Extract the number of thresholds
	numThreshold <- table(sapply(group1pt, "[", group1pt$op == "|")[,"lhs"])
	
	# Find the indicators of each factor
	group1facload <- sapply(group1pt, "[", group1pt$op == "=~")
	factorRep <- split(group1facload[,"rhs"], group1facload[,"lhs"])
	
	# Find marker variables
	marker <- rep(NA, length(factorRep))
	numThresholdMarker <- rep(NA, length(factorRep))
	for(i in seq_along(factorRep)) {
		temp <- sapply(group1pt, "[", group1pt$rhs %in% factorRep[[i]] & group1pt$op == "=~" & group1pt$lhs == names(factorRep)[i])
		marker[i] <- temp[!is.na(temp[,"ustart"]), "rhs"]
		numThresholdMarker[i] <- numThreshold[marker[i]]
	}
	
	numThresholdFactorRep <- lapply(factorRep, function(x) numThreshold[x])
	constraintSecondThreshold <- unlist(lapply(numThresholdFactorRep, function(x) names(which(x > 1)[1])))
	constraintSecondThreshold <- constraintSecondThreshold[!is.na(constraintSecondThreshold)]
			# Find the marker variable of each facto

	for(i in names(numThreshold)) {
		lavaanParTable <- constrainParTable(lavaanParTable, i, "|", "t1", 1:ngroups)
	}
	
	if(length(constraintSecondThreshold) > 0) {
		for(i in constraintSecondThreshold) {
			lavaanParTable <- constrainParTable(lavaanParTable, i, "|", "t2", 1:ngroups)
		}
	}
	
	# Group 1
	for(i in facName) {
		lavaanParTable <- fixParTable(lavaanParTable, i, "~1", "", 1, 0) # Fix factor means as 0
		if(std.lv) {
			lavaanParTable <- fixParTable(lavaanParTable, i, "~~", i, 1, 1)
		} else {
			lavaanParTable <- freeParTable(lavaanParTable, i, "~~", i, 1, NA) # Free factor variances
		}
		# Assuming that all factor covariances are freeParTable
	}
	for(i in varList) {
		lavaanParTable <- fixParTable(lavaanParTable, i, "~~", i, 1, 1)
	}
	
	# Other groups
	for(k in 2:ngroups) {
		for(i in facName) {
			lavaanParTable <- freeParTable(lavaanParTable, i, "~1", "", k, NA)
			if(std.lv) {
				lavaanParTable <- fixParTable(lavaanParTable, i, "~~", i, k, 1)
			} else {
				lavaanParTable <- freeParTable(lavaanParTable, i, "~~", i, k, NA)
			}
		}	
		for(i in varList) {
			lavaanParTable <- freeParTable(lavaanParTable, i, "~~", i, k, NA)
		}
		# Fix the indicator variances of marker variables with two categories as 1
		for(i in seq_along(marker)) {
			if(numThresholdMarker[i] == 1)  lavaanParTable <- fixParTable(lavaanParTable, marker[i], "~~", marker[i], k, 1)
		}
	}
	
	if(std.lv) {
		for(i in seq_along(factorRep)) {
			lavaanParTable <- freeParTable(lavaanParTable, names(factorRep)[i], "=~", marker[i], 1:ngroups, NA)
		}
	}
	# Fit configural invariance
	ListConfigural <- List
	ListConfigural$model <- lavaanParTable
	fitConfigural <- do.call("lavaan", ListConfigural)
	
	# Create the parameter table for metric invariance
	ptMetric <- lavaanParTable
	
	for(i in seq_along(factorRep)) {
		varwithin <- factorRep[[i]]
		if(!std.lv) {
			varwithin <- setdiff(varwithin, marker[i])
		}
		for(j in seq_along(varwithin)) {
			ptMetric <- constrainParTable(ptMetric, names(factorRep)[i], "=~", varwithin[j], 1:ngroups)
		}
	}
	if(std.lv) {
		for(k in 2:ngroups) {
			for(i in facName) {
				ptMetric <- freeParTable(ptMetric, i, "~~", i, k, NA)
			}
		}
	}
	
	ListMetric <- List
	ListMetric$model <- ptMetric
	fitMetric <- do.call("lavaan", ListMetric)

	ptMeans <- ptStrict <- ptMetric
	
	nonMarker <- setdiff(names(numThreshold), marker)
	nonDichoMarker <- numThreshold[which(numThreshold[nonMarker] > 1)]
	scalar <- length(nonDichoMarker) > 0
	if(scalar) {
		ptScalar <- ptMetric
		for(i in seq_along(numThreshold)) {
			thresholdName <- paste0("t", 1:numThreshold[i])
			for(j in seq_along(thresholdName)) {
				ptScalar <- constrainParTable(ptScalar, names(numThreshold)[i], "|", thresholdName[j], 1:ngroups)
			}
		}
		ListScalar <- List
		ListScalar$model <- ptScalar
		fitScalar <- do.call("lavaan", ListScalar)		
		ptMeans <- ptStrict <- ptScalar
	}
	

	fitStrict <- NULL
	# Create the parameter table for strict invariance if specified
	if(strict) {
		ptStrict <- ptScalar
		for(k in 2:ngroups) {
			# Constrain measurement error variances
			for(i in varList) {
				ptStrict <- fixParTable(ptStrict, i, "~~", i, k, 1)
			}
		}
		ListStrict <- List
		ListStrict$model <- ptStrict
		fitStrict <- do.call("lavaan", ListStrict)		
		ptMeans <- ptStrict
	} 
	
	# Create the parameter table for mean equality
	
	# Constrain factor means to be equal
	for(k in 2:ngroups) {
		ptMeans <- fixParTable(ptMeans, facName, "~1", "", k, ustart = 0)
	}
	ListMeans <- List
	ListMeans$model <- ptMeans
	fitMeans <- do.call("lavaan", ListMeans)		

	# Modify these functions from measurementInvariance function
	if(!quiet) {
		cat("\n#################### Measurement invariance tests ####################\n")
		if(!scalar) {
			cat("\nNOTE: All thresholds are equal across groups for scale identification.\n")
		}
		cat("\n#################### Model 1: configural invariance:\n")
		printFitLine(fitConfigural)

		cat("\n#################### Model 2: weak invariance (equal loadings):\n")
		printFitLine(fitMetric)

		cat("\n[Model 1 versus model 2]\n")
		difftest(fitConfigural, fitMetric)
		
		if(scalar & strict) {
			cat("\n#################### Model 3: strong invariance (equal loadings + thresholds):\n")
			printFitLine(fitScalar)
			cat("\n[Model 1 versus model 3]\n")
			difftest(fitConfigural, fitScalar)
			cat("\n[Model 2 versus model 3]\n")
			difftest(fitMetric, fitScalar)
            cat("\n#################### Model 4: strict invariance (equal loadings + thresholds + residuals):\n")
            printFitLine(fitStrict)
            cat("\n[Model 1 versus model 4]\n")
            difftest(fitConfigural, fitStrict)
            cat("\n[Model 2 versus model 4]\n")
            difftest(fitMetric, fitStrict)
            cat("\n[Model 3 versus model 4]\n")
            difftest(fitScalar, fitStrict)
  
            cat("\n#################### Model 5: equal loadings + thresholds + residuals + means:\n")
            printFitLine(fitMeans, horizontal=TRUE)
            cat("\n[Model 1 versus model 5]\n")
            difftest(fitConfigural, fitMeans)
            cat("\n[Model 2 versus model 5]\n")
            difftest(fitMetric, fitMeans)
            cat("\n[Model 3 versus model 5]\n")
            difftest(fitScalar, fitMeans)
            cat("\n[Model 4 versus model 5]\n")
            difftest(fitStrict, fitMeans)		
		} else if (!scalar & strict) {
            cat("\n#################### Model 3: strict invariance (equal loadings + residuals):\n")
            printFitLine(fitStrict)
            cat("\n[Model 1 versus model 3]\n")
            difftest(fitConfigural, fitStrict)
            cat("\n[Model 2 versus model 3]\n")
            difftest(fitMetric, fitStrict)
            cat("\n#################### Model 4: equal loadings + residuals + means:\n")
            printFitLine(fitMeans, horizontal=TRUE)
            cat("\n[Model 1 versus model 4]\n")
            difftest(fitConfigural, fitMeans)
            cat("\n[Model 2 versus model 4]\n")
            difftest(fitMetric, fitMeans)
            cat("\n[Model 3 versus model 4]\n")
            difftest(fitStrict, fitMeans)		
		} else if (scalar & !strict) {
			cat("\n#################### Model 3: strong invariance (equal loadings + thresholds):\n")
			printFitLine(fitScalar)
			cat("\n[Model 1 versus model 3]\n")
			difftest(fitConfigural, fitScalar)
			cat("\n[Model 2 versus model 3]\n")
			difftest(fitMetric, fitScalar)
            cat("\n#################### Model 4: equal loadings + thresholds + means:\n")
            printFitLine(fitMeans)
            cat("\n[Model 1 versus model 4]\n")
            difftest(fitConfigural, fitMeans)
            cat("\n[Model 2 versus model 4]\n")
            difftest(fitMetric, fitMeans)
            cat("\n[Model 3 versus model 4]\n")
            difftest(fitScalar, fitMeans)
		} else {
            cat("\n#################### Model 3: equal loadings + means:\n")
            printFitLine(fitMeans)
            cat("\n[Model 1 versus model 3]\n")
            difftest(fitConfigural, fitMeans)
            cat("\n[Model 2 versus model 3]\n")
            difftest(fitMetric, fitMeans)
		}
	}
	return(invisible(list(configural = fitConfigural, metric = fitMetric, scalar = fitScalar, strict = fitStrict, means = fitMeans)))
}

multipleAllEqual <- function(...) {
    obj <- list(...)
    multipleAllEqualList(obj)
}

multipleAllEqualList <- function(obj) {
    for (i in 2:length(obj)) {
        for (j in 1:(i - 1)) {
            temp <- isTRUE(all.equal(obj[[i]], obj[[j]]))
            if (!temp) 
                return(FALSE)
        }
    }
    return(TRUE)
} 

multipleAnyEqual <- function(...) {
    obj <- list(...)
    multipleAnyEqualList(obj)
}

multipleAnyEqualList <- function(obj) {
    for (i in 2:length(obj)) {
        for (j in 1:(i - 1)) {
            temp <- isTRUE(all.equal(obj[[i]], obj[[j]]))
            if (temp) 
                return(TRUE)
        }
    }
    return(FALSE)
} 


