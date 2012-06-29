## Title: Longitudinal (or within-group, such as dyadic data) measurement invariance
## Author: Sunthud Pornprasertmanit <psunthud@ku.edu>
## Description: Test measurement invariance and save the fitted objects
##----------------------------------------------------------------------------##

longInvariance <- function(model, varList, auto = "all", constrainAuto = FALSE, fixed.x = TRUE, std.lv = FALSE, group=NULL, group.equal="", group.partial="", warn=TRUE, debug=FALSE, strict = FALSE, quiet = FALSE, ...) {

	List <- list(...)

	# Find the number of groups
	ngroups <- 1
	if(!is.null(group)) {
		if(!is.null(List$data)) {
			ngroups <- length(unique(List$data[,group]))
		} else if (!is.null(List$sample.cov)) {
			ngroups <- length(List$sample.cov)
		} else {
			stop("Cannot find the specifying variable name in the 'group' argument.")
		}
	}
	
	# Get the lavaan parameter table
    if(is.character(model)) {
        lavaanParTable <- 
            lavaanify(model           = model,
                      meanstructure   = TRUE, 
                      int.ov.free     = TRUE,
                      int.lv.free     = FALSE,
                      orthogonal      = FALSE, 
                      fixed.x         = fixed.x,
                      std.lv          = std.lv,

                      auto.fix.first  = ifelse(std.lv, FALSE, TRUE),
                      auto.fix.single = TRUE,
                      auto.var        = TRUE,
                      auto.cov.lv.x   = TRUE,
                      auto.cov.y      = TRUE,

                      ngroups         = ngroups,
                      group.equal     = group.equal, 
                      group.partial   = group.partial,
                      debug           = debug,
                      warn            = warn,
                      as.data.frame.  = TRUE)
    } else if(is.list(model)) {
        if(!is.null(model$lhs) && !is.null(model$op)  &&
           !is.null(model$rhs) && !is.null(model$free)) {
            lavaanParTable <- model
        } else if(is.character(model[[1]])) {
            stop("lavaan ERROR: model is a list, but not a parameterTable?")
        }
    } else {
        cat("model type: ", class(model), "\n")
        stop("lavaan ERROR: model is not of type character or list")
    }
	
	# Error checking on the varList argument and get the factor name corresponding to each elements of the list
	facName <- lapply(varList, function(vec, pt) pt$lhs[(pt$op == "=~") & (pt$rhs %in% vec)], pt=lavaanParTable)
	if(any(sapply(facName, function(x) length(unique(x)) > 1))) stop("The factor names of the same element of the 'varList' are not the same.")
	if(length(unique(sapply(facName, function(x) length(x)))) > 1) stop("The numbers of variables in each element are not equal.")
	facName <- unlist(lapply(facName, unique))
	
	# Impose the autocorrelation in the parameter table
	if(auto != 0) {
		if(is.numeric(auto) && auto >= length(varList)) stop("The number of lag in auto-correlation is not possible in the current number of timepoints.")
		if(auto == "all") auto <- length(varList) - 1
		for(k in 1:ngroups) {
			for(i in 1:length(varList[[1]])) {
				name <- sapply(varList, function(x, element) x[element], element = i)
				for(j in 1:auto) {
					vec <- 1:(length(varList) - j)
					lavaanParTable <- freeParTable(lavaanParTable, name[vec], "~~", name[vec + j], k, ustart = NA)
					if(constrainAuto & (length(vec) > 1)) lavaanParTable <- constrainParTable(lavaanParTable, name[vec], "~~", name[vec + j], k)
				}
			}
		}
	}
	
	# Fit configural invariance
	fitConfigural <- lavaan(lavaanParTable, ..., group=group, group.equal=group.equal, group.partial=group.partial, warn=TRUE, debug=FALSE)
	
	# Create the parameter table for metric invariance
	ptMetric <- lavaanParTable
	if(std.lv) {
		for(k in 1:ngroups) {
			# Free variances of factor 2, 3, ...
			ptMetric <- freeParTable(ptMetric, facName[-1], "~~", facName[-1], k, ustart = NA)
			
			# Constrain factor loadings
			for(i in 1:length(varList[[1]])) {
				ptMetric <- constrainParTable(ptMetric, facName, "=~", sapply(varList, function(x, element) x[element], element = i), k)
			}
		}
		ptMetric$ustart[(ptMetric$op == "=~") & (ptMetric$rhs %in% sapply(varList, function(x, element) x[element], element = 1))] <- 1
		
	} else {
		for(k in 1:ngroups) {
			# Constrain factor loadings but keep marker variables
			for(i in 2:length(varList[[1]])) {
				ptMetric <- constrainParTable(ptMetric, facName, "=~", sapply(varList, function(x, element) x[element], element = i), k)
			}
		}
	}
	fitMetric <- lavaan(ptMetric, ..., group=group, group.equal=group.equal, group.partial=group.partial, warn=TRUE, debug=FALSE)
	
	# Create the parameter table for scalar invariance
	ptScalar <- ptMetric
	for(k in 1:ngroups) {
		# Free means of factors 2, 3, ...
		ptScalar <- freeParTable(ptScalar, facName[-1], "~1", "", k, ustart = NA)
		
		# Constrain measurement intercepts
		for(i in 1:length(varList[[1]])) {
			ptScalar <- constrainParTable(ptScalar, sapply(varList, function(x, element) x[element], element = i), "~1", "", k)
		}
	}
	ptScalar$ustart[(ptMetric$op == "~1") & (ptMetric$rhs %in% facName)] <- 0
	fitScalar <- lavaan(ptScalar, ..., group=group, group.equal=group.equal, group.partial=group.partial, warn=TRUE, debug=FALSE)
	
	ptMeans <- ptScalar
	
	# Create the parameter table for strict invariance if specified
	ptStrict <- ptScalar
	fitStrict <- NULL
	if(strict) {
		ptStrict <- ptScalar
		for(k in 1:ngroups) {
			# Constrain measurement error variances
			for(i in 1:length(varList[[1]])) {
				name <- sapply(varList, function(x, element) x[element], element = i)
				ptStrict <- constrainParTable(ptStrict, name, "~~", name, k)
			}
		}
		fitStrict <- lavaan(ptStrict, ..., group=group, group.equal=group.equal, group.partial=group.partial, warn=TRUE, debug=FALSE)
		ptMeans <- ptStrict
	} 
	
	# Create the parameter table for mean equality
	
	# Constrain factor means to be equal
	for(k in 1:ngroups) {
		ptMeans <- fixParTable(ptMeans, facName[-1], "~1", "", k, ustart = 0)
	}
	fitMeans <- lavaan(ptMeans, ..., group=group, group.equal=group.equal, group.partial=group.partial, warn=TRUE, debug=FALSE)
	
	# Modify these functions from measurementInvariance function
	if(!quiet) {
		cat("\n#################### Measurement invariance tests ####################\n")
		cat("\nThe order of autocorrelation: ", auto, "\n")
		cat("\n#################### Model 1: configural invariance:\n")
		printFitLine(fitConfigural)

		cat("\n#################### Model 2: weak invariance (equal loadings):\n")
		printFitLine(fitMetric)

		cat("\n[Model 1 versus model 2]\n")
		difftest(fitConfigural, fitMetric)

		cat("\n#################### Model 3: strong invariance (equal loadings + intercepts):\n")
		printFitLine(fitScalar)
		cat("\n[Model 1 versus model 3]\n")
		difftest(fitConfigural, fitScalar)
		cat("\n[Model 2 versus model 3]\n")
		difftest(fitMetric, fitScalar)
		if(strict) {
            cat("\n#################### Model 4: strict invariance (equal loadings + intercepts + residuals):\n")
            printFitLine(fitStrict)
            cat("\n[Model 1 versus model 4]\n")
            difftest(fitConfigural, fitStrict)
            cat("\n[Model 2 versus model 4]\n")
            difftest(fitMetric, fitStrict)
            cat("\n[Model 3 versus model 4]\n")
            difftest(fitScalar, fitStrict)
  
            cat("\n#################### Model 5: equal loadings + intercepts + residuals + means:\n")
            printFitLine(fitMeans, horizontal=TRUE)
            cat("\n[Model 1 versus model 5]\n")
            difftest(fitConfigural, fitMeans)
            cat("\n[Model 2 versus model 5]\n")
            difftest(fitMetric, fitMeans)
            cat("\n[Model 3 versus model 5]\n")
            difftest(fitScalar, fitMeans)
            cat("\n[Model 4 versus model 5]\n")
            difftest(fitStrict, fitMeans)
        } else {
            cat("\n#################### Model 4: equal loadings + intercepts + means:\n")
            printFitLine(fitMeans)
            cat("\n[Model 1 versus model 4]\n")
            difftest(fitConfigural, fitMeans)
            cat("\n[Model 2 versus model 4]\n")
            difftest(fitMetric, fitMeans)
            cat("\n[Model 3 versus model 4]\n")
            difftest(fitScalar, fitMeans)
        }
	}
	return(invisible(list(configural = fitConfigural, metric = fitMetric, scalar = fitScalar, strict = fitStrict, means = fitMeans)))
}

# freeParTable: Free elements in parameter table

freeParTable <- function(parTable, lhs, op, rhs, group, ustart = NA) {
	target <- cbind(lhs, op, rhs, group)
	element <- apply(target, 1, matchElement, parTable=parTable)
	for(i in 1:nrow(target)) {
		if((length(element[i]) == 0) | is.na(element[i])) {
			parTable <- as.list(parTable)
			parTable$id <- c(parTable$id, as.integer(max(parTable$id) + 1))
			parTable$lhs <- c(parTable$lhs, as.character(target[i, 1]))
			parTable$op <- c(parTable$op, as.character(target[i, 2]))
			parTable$rhs <- c(parTable$rhs, as.character(target[i, 3]))
			parTable$user <- c(parTable$user, as.integer(1))
			parTable$group <- c(parTable$group, as.integer(target[i, 4]))
			parTable$free <- c(parTable$free, as.integer(max(parTable$free) + 1))
			parTable$ustart <- c(parTable$ustart, as.numeric(NA))
			parTable$exo <- c(parTable$exo, as.integer(0))
			parTable$eq.id <- c(parTable$eq.id, as.integer(0))
			parTable$label <- c(parTable$label, as.character(""))
			parTable$unco <- c(parTable$unco, as.integer(max(parTable$unco) + 1))
		} else {
			if(parTable$free[element[i]] != 0) warnings(paste("The", lhs, op, rhs, group, "is free already."))
			parTable$unco[element[i]] <- max(parTable$unco) + 1
			parTable$ustart[element[i]] <- ustart
			parTable$user[element[i]] <- 1
			parTable$free[element[i]] <- max(parTable$free) + 1
		}
	}
	parTable$unco <- rearrangeFreeElement(parTable$unco)
	parTable$free <- rearrangeFreeElement(parTable$free)	
	return(parTable)
}

# fixParTable: Fix elements in parameter table

fixParTable <- function(parTable, lhs, op, rhs, group, ustart = NA) {
	target <- cbind(lhs, op, rhs, group)
	element <- apply(target, 1, matchElement, parTable=parTable)
	for(i in 1:nrow(target)) {
		if(parTable$free[element[i]] == 0) warnings(paste("The", lhs, op, rhs, group, "is fixed already."))
		parTable$unco[element[i]] <- 0
		parTable$ustart[element[i]] <- ustart
		parTable$user[element[i]] <- 1
		parTable$free[element[i]] <- 0
	}
	parTable$unco <- rearrangeFreeElement(parTable$unco)
	parTable$free <- rearrangeFreeElement(parTable$free)	
	return(parTable)
}

# constrainParTable: Impose equality constraints in any set of elements in the parameter table

constrainParTable <- function(parTable, lhs, op, rhs, group) {
	target <- cbind(lhs, op, rhs, group)
	element <- apply(target, 1, matchElement, parTable=parTable)
	parTable$user[element[1]] <- 1

	free <- parTable$free[element[1]]
	if(parTable$eq.id[element[1]] == 0) parTable$eq.id[element[1]] <- max(parTable$eq.id) + 1
	eq.id <- parTable$eq.id[element[1]]
	for(i in 2:nrow(target)) {
		parTable$user[element[i]] <- 1
		parTable$free[element[i]] <- free
		parTable$eq.id[element[i]] <- eq.id
	}
	parTable$unco <- rearrangeFreeElement(parTable$unco)
	parTable$free <- rearrangeFreeElement(parTable$free)	
	return(parTable)	
}

# matchElement: Find the number of row that have the specification in vec (lhs, op, rhs, group)

matchElement <- function(parTable, vec) {
	which((parTable$lhs == vec[1]) & (parTable$op == vec[2]) & (parTable$rhs == vec[3]) & (parTable$group == vec[4]))
}

# rearrangeFreeElement: Rearrange the number listed in 'free' or 'unco' in parameter tables 

rearrangeFreeElement <- function(vec) {
	vec2 <- vec
	vec <- vec[vec != 0]
	uvec <- unique(vec)
	newvec <- 1:length(unique(vec))
	vec2[vec2 != 0] <- newvec[match(vec, uvec)]
	class(vec2) <- "integer"
	return(vec2)
}
