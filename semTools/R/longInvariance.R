## Title: Longitudinal (or within-group, such as dyadic data) measurement invariance
## Author: Sunthud Pornprasertmanit <psunthud@ku.edu>
## Description: Test measurement invariance and save the fitted objects
##----------------------------------------------------------------------------##

longInvariance <- function(model, varList, auto = "all", constrainAuto = FALSE, fixed.x = TRUE, std.lv = FALSE, group=NULL, group.equal="", group.partial="", warn=TRUE, debug=FALSE, strict = FALSE, quiet = FALSE, fit.measures = "default", method = "satorra.bentler.2001", ...) {

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
            lavaan::lavaanify(model           = model,
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
	fitConfigural <- lavaan::lavaan(lavaanParTable, ..., group=group, group.equal=group.equal, group.partial=group.partial, warn=TRUE, debug=FALSE)
	
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
	fitMetric <- lavaan::lavaan(ptMetric, ..., group=group, group.equal=group.equal, group.partial=group.partial, warn=TRUE, debug=FALSE)
	
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
	fitScalar <- lavaan::lavaan(ptScalar, ..., group=group, group.equal=group.equal, group.partial=group.partial, warn=TRUE, debug=FALSE)
	
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
		fitStrict <- lavaan::lavaan(ptStrict, ..., group=group, group.equal=group.equal, group.partial=group.partial, warn=TRUE, debug=FALSE)
		ptMeans <- ptStrict
	} 
	
	# Create the parameter table for mean equality
	
	# Constrain factor means to be equal
	for(k in 1:ngroups) {
		ptMeans <- fixParTable(ptMeans, facName[-1], "~1", "", k, ustart = 0)
	}
	fitMeans <- lavaan::lavaan(ptMeans, ..., group=group, group.equal=group.equal, group.partial=group.partial, warn=TRUE, debug=FALSE)

	FIT <- invisible(list(fit.configural = fitConfigural, fit.loadings = fitMetric, fit.thresholds = fitScalar, fit.residuals = fitStrict, fit.means = fitMeans))
	FIT <- FIT[!sapply(FIT, is.null)]

	if(!quiet) {
        printInvarianceResult(FIT, fit.measures, method)
    }

    invisible(FIT)
	
	# Modify these functions from measurementInvariance function
	# if(!quiet) {
		# cat("\n#################### Measurement invariance tests ####################\n")
		# cat("\nThe order of autocorrelation: ", auto, "\n")
		# cat("\n#################### Model 1: configural invariance:\n")
		# printFitLine(fitConfigural)

		# cat("\n#################### Model 2: weak invariance (equal loadings):\n")
		# printFitLine(fitMetric)

		# cat("\n[Model 1 versus model 2]\n")
		# difftest(fitConfigural, fitMetric)

		# cat("\n#################### Model 3: strong invariance (equal loadings + intercepts):\n")
		# printFitLine(fitScalar)
		# cat("\n[Model 1 versus model 3]\n")
		# difftest(fitConfigural, fitScalar)
		# cat("\n[Model 2 versus model 3]\n")
		# difftest(fitMetric, fitScalar)
		# if(strict) {
            # cat("\n#################### Model 4: strict invariance (equal loadings + intercepts + residuals):\n")
            # printFitLine(fitStrict)
            # cat("\n[Model 1 versus model 4]\n")
            # difftest(fitConfigural, fitStrict)
            # cat("\n[Model 2 versus model 4]\n")
            # difftest(fitMetric, fitStrict)
            # cat("\n[Model 3 versus model 4]\n")
            # difftest(fitScalar, fitStrict)
  
            # cat("\n#################### Model 5: equal loadings + intercepts + residuals + means:\n")
            # printFitLine(fitMeans, horizontal=TRUE)
            # cat("\n[Model 1 versus model 5]\n")
            # difftest(fitConfigural, fitMeans)
            # cat("\n[Model 2 versus model 5]\n")
            # difftest(fitMetric, fitMeans)
            # cat("\n[Model 3 versus model 5]\n")
            # difftest(fitScalar, fitMeans)
            # cat("\n[Model 4 versus model 5]\n")
            # difftest(fitStrict, fitMeans)
        # } else {
            # cat("\n#################### Model 4: equal loadings + intercepts + means:\n")
            # printFitLine(fitMeans)
            # cat("\n[Model 1 versus model 4]\n")
            # difftest(fitConfigural, fitMeans)
            # cat("\n[Model 2 versus model 4]\n")
            # difftest(fitMetric, fitMeans)
            # cat("\n[Model 3 versus model 4]\n")
            # difftest(fitScalar, fitMeans)
        # }
	# }
	# return(invisible(list(fit.configural = fitConfigural, fit.loadings = fitMetric, fit.intercepts = fitScalar, fit.residuals = fitStrict, fit.means = fitMeans)))
}

# freeParTable: Free elements in parameter table

freeParTable <- function(parTable, lhs, op, rhs, group, ustart = NA) {
	parTable$start <- parTable$est <- parTable$se <- NULL
	target <- cbind(lhs, op, rhs, group)
	for(i in 1:nrow(target)) {
		targetElem <- matchElement(parTable = parTable, vec = target[i,])
		ptargetElem <- parTable$plabel[targetElem]
		if((length(targetElem) == 0) || is.na(targetElem)) {
			newline <- list(lhs = as.character(target[i, 1]),
				op = as.character(target[i, 2]),
				rhs = as.character(target[i, 3]),
				group = as.integer(target[i, 4]),
				free = as.integer(max(parTable$free) + 1),
				ustart = as.numeric(NA))
			parTable <- patMerge(pt1 = parTable, pt2 = newline)
		} else {
			if(parTable$free[targetElem] == 0) {
				parTable$ustart[targetElem] <- ustart
				parTable$user[targetElem] <- 1
				parTable$free[targetElem] <- max(parTable$free) + 1
			}
			equalelement <- which(parTable$op == "==")
			rmelem <- intersect(union(match(ptargetElem, parTable$lhs), match(ptargetElem, parTable$rhs)), equalelement)
			if(length(rmelem) > 0) parTable <- removeEqCon(parTable, rmelem)
		}
	}
	parTable <- rearrangept(parTable)	
	return(parTable)
}

# removeEqCon: Remove equality constraints

removeEqCon <- function(pt, element) {
	pt <- lapply(pt, "[", -element)
	pt$id <- seq_along(pt$id)
	pt
}

# fixParTable: Fix elements in parameter table

fixParTable <- function(parTable, lhs, op, rhs, group, ustart = NA) {
	parTable$start <- parTable$est <- parTable$se <- NULL
	target <- cbind(lhs, op, rhs, group)
	element <- apply(target, 1, matchElement, parTable=parTable)
	for(i in 1:nrow(target)) {
		if(parTable$free[element[i]] == 0) warnings(paste("The", lhs, op, rhs, group, "is fixed already."))
		
		# equalelement <- which(parTable$op == "==")
		# targetElem <- matchElement(parTable = parTable, vec = target[i,])
		# ptargetElem <- parTable$plabel[targetElem]
		# rmelem <- intersect(union(match(ptargetElem, parTable$lhs), match(ptargetElem, parTable$rhs)), equalelement)
		# if(length(rmelem) > 0) parTable <- removeEqCon(parTable, rmelem)
		
		parTable$ustart[element[i]] <- ustart
		parTable$user[element[i]] <- 1
		parTable$free[element[i]] <- 0
	}
	parTable <- rearrangept(parTable)
	# rearrangePlabel with change all equality constraints
	return(parTable)
}

# constrainParTable: Impose equality constraints in any set of elements in the parameter table

constrainParTable <- function(parTable, lhs, op, rhs, group) {
	parTable$start <- parTable$est <- parTable$se <- NULL
	target <- cbind(lhs, op, rhs, group)
	element <- apply(target, 1, matchElement, parTable=parTable)

	#     id lhs  op rhs user group free ustart exo label plabel  start	
	for(i in 2:length(element)) {
		len <- length(parTable$id)
		newline <- list(lhs = parTable$plabel[element[1]],
			op = "==",
			rhs = parTable$plabel[element[i]])
		if(!any(parTable$lhs == newline$lhs & parTable$op == newline$op & parTable$rhs == newline$rhs)) parTable <- patMerge(pt1 = parTable, pt2 = newline)
	}
	return(parTable)	
}

# matchElement: Find the number of row that have the specification in vec (lhs, op, rhs, group)

matchElement <- function(parTable, vec) {
	if(is.null(parTable$group)) {
		return(which((parTable$lhs == vec[1]) & (parTable$op == vec[2]) & (parTable$rhs == vec[3])))
	} else {
		return(which((parTable$lhs == vec[1]) & (parTable$op == vec[2]) & (parTable$rhs == vec[3]) & (parTable$group == vec[4])))
	}
}

# rearrangeFreeElement: Rearrange the number listed in 'free' in parameter tables 

rearrangeFreeElement <- function(vec) {
	vec2 <- vec
	vec <- vec[vec != 0]
	uvec <- unique(vec)
	newvec <- 1:length(unique(vec))
	vec2[vec2 != 0] <- newvec[match(vec, uvec)]
	class(vec2) <- "integer"
	return(vec2)
}

createplabel <- function(num) {
	result <- paste0(".p", num, ".")
	result[num == 0] <- ""
	result
}

# rearrangept: Rearrange parameter table and plabel

rearrangept <- function(pt) {
	oldfree <- pt$free
	newfree <- rearrangeFreeElement(oldfree)
	oldplabel <- pt$plabel
	newplabel <- createplabel(seq_along(pt$op))
	eqpos <- which(pt$op == "==")
	newplabel[eqpos] <- ""
	if(length(eqpos) > 0) {
		eqlhs <- pt$lhs[eqpos]
		eqrhs <- pt$rhs[eqpos]
		matchlhs <- match(eqlhs, oldplabel)
		matchrhs <- match(eqrhs, oldplabel)
		neweqlhs <- newplabel[matchlhs]
		neweqrhs <- newplabel[matchrhs]
		neweqlhs[is.na(matchlhs)] <- eqlhs[is.na(matchlhs)]
		neweqrhs[is.na(matchrhs)] <- eqrhs[is.na(matchrhs)]
		pt$lhs[eqpos] <- neweqlhs
		pt$rhs[eqpos] <- neweqrhs
	}
	pt$free <- newfree
	pt$plabel <- newplabel
	pt
}

getValue <- function(parTable, est, lhs, op, rhs, group) {
	target <- cbind(lhs, op, rhs, group)
	element <- apply(target, 1, matchElement, parTable=parTable)
	free <- parTable$free[element]
	out <- parTable$ustart[element]
	out[free != 0] <- est[free[free != 0]]
	out
}

patMerge <- function (pt1 = NULL, pt2 = NULL, remove.duplicated = FALSE, 
    fromLast = FALSE, warn = TRUE) 
{
    pt1 <- as.data.frame(pt1, stringsAsFactors = FALSE)
    pt2 <- as.data.frame(pt2, stringsAsFactors = FALSE)
    stopifnot(!is.null(pt1$lhs), !is.null(pt1$op), !is.null(pt1$rhs), 
        !is.null(pt2$lhs), !is.null(pt2$op), !is.null(pt2$rhs))
    if (is.null(pt1$group) && is.null(pt2$group)) {
        TMP <- rbind(pt1[, c("lhs", "op", "rhs", "group")], pt2[, 
            c("lhs", "op", "rhs", "group")])
    }
    else {
        if (is.null(pt1$group) && !is.null(pt2$group)) {
            pt1$group <- rep(1L, length(pt1$lhs))
        }
        else if (is.null(pt2$group) && !is.null(pt1$group)) {
            pt2$group <- rep(1L, length(pt2$lhs))
        }
        TMP <- rbind(pt1[, c("lhs", "op", "rhs", "group")], pt2[, 
            c("lhs", "op", "rhs", "group")])
    }
    if (is.null(pt1$user) && !is.null(pt2$user)) {
        pt1$user <- rep(0L, length(pt1$lhs))
    }
    else if (is.null(pt2$user) && !is.null(pt1$user)) {
        pt2$user <- rep(0L, length(pt2$lhs))
    }
    if (is.null(pt1$free) && !is.null(pt2$free)) {
        pt1$free <- rep(0L, length(pt1$lhs))
    }
    else if (is.null(pt2$free) && !is.null(pt1$free)) {
        pt2$free <- rep(0L, length(pt2$lhs))
    }
    if (is.null(pt1$ustart) && !is.null(pt2$ustart)) {
        pt1$ustart <- rep(0, length(pt1$lhs))
    }
    else if (is.null(pt2$ustart) && !is.null(pt1$ustart)) {
        pt2$ustart <- rep(0, length(pt2$lhs))
    }
    if (is.null(pt1$exo) && !is.null(pt2$exo)) {
        pt1$exo <- rep(0L, length(pt1$lhs))
    }
    else if (is.null(pt2$exo) && !is.null(pt1$exo)) {
        pt2$exo <- rep(0L, length(pt2$lhs))
    }
    if (is.null(pt1$label) && !is.null(pt2$label)) {
        pt1$label <- rep("", length(pt1$lhs))
    }
    else if (is.null(pt2$label) && !is.null(pt1$label)) {
        pt2$label <- rep("", length(pt2$lhs))
    }
    if (is.null(pt1$plabel) && !is.null(pt2$plabel)) {
        pt1$plabel <- rep("", length(pt1$lhs))
    }
    else if (is.null(pt2$plabel) && !is.null(pt1$plabel)) {
        pt2$plabel <- rep("", length(pt2$lhs))
    }
    if (is.null(pt1$start) && !is.null(pt2$start)) {
        pt1$start <- rep(as.numeric(NA), length(pt1$lhs))
    }
    else if (is.null(pt2$start) && !is.null(pt1$start)) {
        pt2$start <- rep(as.numeric(NA), length(pt2$lhs))
    }
	if(!is.null(pt1$est)) pt1$est <- NULL
	if(!is.null(pt2$est)) pt2$est <- NULL
	if(!is.null(pt1$se)) pt1$se <- NULL
	if(!is.null(pt2$se)) pt2$se <- NULL
    if (remove.duplicated) {
        idx <- which(duplicated(TMP, fromLast = fromLast))
        if (length(idx)) {
            if (warn) {
                warning("lavaan WARNING: duplicated parameters are ignored:\n", 
                  paste(apply(pt1[idx, c("lhs", "op", "rhs")], 
                    1, paste, collapse = " "), collapse = "\n"))
            }
            if (fromLast) {
                pt1 <- pt1[-idx, ]
            }
            else {
                idx <- idx - nrow(pt1)
                pt2 <- pt2[-idx, ]
            }
        }
    } else if (!is.null(pt1$start) && !is.null(pt2$start)) {
        for (i in 1:length(pt1$lhs)) {
            idx <- which(pt2$lhs == pt1$lhs[i] & pt2$op == pt1$op[i] & 
                pt2$rhs == pt1$rhs[i] & pt2$group == pt1$group[i])
            pt2$start[idx] <- pt1$start[i]
        }
    }
    if (is.null(pt1$id) && !is.null(pt2$id)) {
        nid <- max(pt2$id)
        pt1$id <- (nid + 1L):(nid + nrow(pt1))
    }
    else if (is.null(pt2$id) && !is.null(pt1$id)) {
        nid <- max(pt1$id)
        pt2$id <- (nid + 1L):(nid + nrow(pt2))
    }
    NEW <- base::merge(pt1, pt2, all = TRUE, sort = FALSE)
    NEW
}
