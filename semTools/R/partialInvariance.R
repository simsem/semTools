# Work with only with congeneric models

partialInvariance <- function(fit, type, free = NULL, fix = NULL, return.fit = FALSE) { 
	type <- tolower(type)
	numType <- 0
	fit1 <- fit0 <- NULL
	# fit0 = Nested model, fit1 = Parent model
	if(type %in% c("metric", "weak", "loading", "loadings")) {
		numType <- 1
		if(all(c("configural", "metric") %in% names(fit))) {
			fit1 <- fit$configural
			fit0 <- fit$metric
		} else {
			stop("Both configural and metric invariance models are needed in the 'fit' argument")
		}
	} else if (type %in% c("scalar", "strong", "intercept", "intercepts", "threshold", "thresholds")) {
		numType <- 2
		if(all(c("metric", "scalar") %in% names(fit))) {
			fit1 <- fit$metric
			fit0 <- fit$scalar
		} else {
			stop("Both metric and scalar invariance models are needed in the 'fit' argument")
		}
	} else if (type %in% c("strict", "residual", "residuals", "error", "errors")) {
		numType <- 3
		if(all(c("scalar", "strict") %in% names(fit))) {
			fit1 <- fit$scalar
			fit0 <- fit$strict
		} else {
			stop("Both scalar and strict invariance models are needed in the 'fit' argument")
		}
	} else if (type %in% c("means", "mean")) {
		numType <- 4
		if("means" %in% names(fit)) {
			fit0 <- fit$means
			if("strict" %in% names(fit)) {
				fit1 <- fit$strict			
			} else if ("scalar" %in% names(fit)) {
				fit1 <- fit$scalar			
			} else {
				stop("Either scalar or strict invariance models is needed in the 'fit' argument")
			}
		} else {
			stop("Mean invariance models is needed in the 'fit' argument")
		}
	} else {
		stop("Please specify the correct type of measurement invariance. See the help page.")
	}
	pt1 <- partable(fit1)
	pt0 <- partable(fit0)
	namept1 <- paramNameFromPt(pt1)
	namept0 <- paramNameFromPt(pt0)
	if(length(table(table(pt0$rhs[pt0$op == "=~"]))) != 1) stop("The model is not congeneric. This function does not support non-congeneric model.")
	varfree <- varnames <- unique(pt0$rhs[pt0$op == "=~"])
	facnames <- unique(pt0$lhs[(pt0$op == "=~") & (pt0$rhs %in% varnames)])
	facrepresent <- table(pt0$lhs[(pt0$op == "=~") & (pt0$rhs %in% varnames)], pt0$rhs[(pt0$op == "=~") & (pt0$rhs %in% varnames)])
	if(any(apply(facrepresent, 2, function(x) sum(x != 0)) > 1)) stop("The model is not congeneric. This function does not support non-congeneric model.")
	facList <- list()
	for(i in 1:nrow(facrepresent)) {
		facList[[i]] <- colnames(facrepresent)[facrepresent[i,] > 0]
	}
	names(facList) <- rownames(facrepresent)
	facList <- facList[match(names(facList), facnames)]
	fixLoadingFac <- list()
	for(i in seq_along(facList)) {
		select <- pt1$lhs == names(facList)[i] & pt1$op == "=~" & pt1$rhs %in% facList[[i]] & pt1$group == 1 & pt1$free == 0 & (!is.na(pt1$ustart) & pt1$ustart > 0)
		fixLoadingFac[[i]] <- pt1$rhs[select]
	}
	names(fixLoadingFac) <- names(facList)
	fixIntceptFac <- list()
	for(i in seq_along(facList)) {
		select <- pt1$op == "~1" & pt1$rhs %in% facList[[i]] & pt1$group == 1 & pt1$free == 0
		fixIntceptFac[[i]] <- pt1$rhs[select]
	}
	names(fixIntceptFac) <- names(facList)
	
	ngroups <- max(pt0$group)
	if(ngroups <= 1) stop("Well, the number of groups is 1. Measurement invariance across 'groups' cannot be done.")

	if(numType == 4) {
		if(!all(c(free, fix) %in% facnames)) stop("'free' and 'fix' arguments should consist of factor names because mean invariance is tested.")
	} else {
		if(!all(c(free, fix) %in% varnames)) stop("'free' and 'fix' arguments should consist of variable names.")
	}
	result <- fixCon <- freeCon <- NULL
	listFreeCon <- listFixCon <- list()		
	beta <- coef(fit1)
	waldMat <- matrix(0, ngroups - 1, length(beta))
	if(numType == 1) {
		if(!is.null(free) | !is.null(fix)) {
			if(!is.null(fix)) {
				facinfix <- findFactor(fix, facList)
				dup <- duplicated(facinfix)
				for(i in seq_along(fix)) {
					if(dup[i]) {
						pt0 <- constrainParTable(pt0, facinfix[i], "=~", fix[i], 1:ngroups)
						pt1 <- constrainParTable(pt1, facinfix[i], "=~", fix[i], 1:ngroups)					
					} else {
						oldmarker <- fixLoadingFac[[facinfix[i]]]
						if(length(oldmarker) > 0) {
							oldmarkerval <- pt1$ustart[pt1$lhs == facinfix[i] & pt1$op == "=~" & pt1$rhs == oldmarker & pt1$group == 1]
							if(oldmarker == fix[i]) {
								pt0 <- fixParTable(pt0, facinfix[i], "=~", fix[i], 1:ngroups, oldmarkerval)
								pt1 <- fixParTable(pt1, facinfix[i], "=~", fix[i], 1:ngroups, oldmarkerval)
							} else {
								pt0 <- freeParTable(pt0, facinfix[i], "=~", oldmarker, 1:ngroups)
								pt0 <- constrainParTable(pt0, facinfix[i], "=~", oldmarker, 1:ngroups)
								pt1 <- freeParTable(pt1, facinfix[i], "=~", oldmarker, 1:ngroups)
								pt0 <- fixParTable(pt0, facinfix[i], "=~", fix[i], 1:ngroups, oldmarkerval)
								pt1 <- fixParTable(pt1, facinfix[i], "=~", fix[i], 1:ngroups, oldmarkerval)
								fixLoadingFac[[facinfix[i]]] <- fix[i]
							}
						} else {
							pt0 <- constrainParTable(pt0, facinfix[i], "=~", fix[i], 1:ngroups)
							pt1 <- constrainParTable(pt1, facinfix[i], "=~", fix[i], 1:ngroups)					
						}
					}
				}
			}
			if(!is.null(free)) {
				facinfree <- findFactor(free, facList)
				for(i in seq_along(free)) {
					# Need to change marker variable if fixed
					oldmarker <- fixLoadingFac[[facinfree[i]]]
					if(length(oldmarker) > 0 && oldmarker == free[i]) {
						oldmarkerval <- pt1$ustart[pt1$lhs == facinfix[i] & pt1$op == "=~" & pt1$rhs == oldmarker & pt1$group == 1]
						candidatemarker <- setdiff(facList[[facinfree[i]]], free[i])[1]
						pt0 <- freeParTable(pt0, facinfree[i], "=~", free[i], 1:ngroups)
						pt1 <- freeParTable(pt1, facinfree[i], "=~", free[i], 1:ngroups)
						pt0 <- fixParTable(pt0, facinfix[i], "=~", candidatemarker, 1:ngroups, oldmarkerval)
						pt1 <- fixParTable(pt1, facinfix[i], "=~", candidatemarker, 1:ngroups, oldmarkerval)
						fixLoadingFac[[facinfix[i]]] <- candidatemarker
					} else {
						pt0 <- freeParTable(pt0, facinfree[i], "=~", free[i], 1:ngroups)
						pt1 <- freeParTable(pt1, facinfree[i], "=~", free[i], 1:ngroups)
					}
				}
			}
			namept1 <- paramNameFromPt(pt1)
			namept0 <- paramNameFromPt(pt0)
			fit0 <- refit(pt0, fit0)
			fit1 <- refit(pt1, fit1)
			beta <- coef(fit1)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		fixCon <- freeCon <- waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p")
		colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")
		index <- which((pt1$rhs %in% varfree) & (pt1$op == "=~") & (pt1$group == 1))
		facinfix <- findFactor(fix, facList)
		varinfixvar <- unlist(facList[facinfix])
		varinfixvar <- setdiff(varinfixvar, setdiff(varinfixvar, varfree))
		indexfixvar <- which((pt1$rhs %in% varinfixvar) & (pt1$op == "=~") & (pt1$group == 1))
		varnonfixvar <- setdiff(varfree, varinfixvar)
		indexnonfixvar <- setdiff(index, indexfixvar)
		
		pos <- 1
		for(i in seq_along(indexfixvar)) {
			runnum <- indexfixvar[i]
			temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
			tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
			if(!is(tryresult, "try-error")) {
				compresult <- try(modelcomp <- lavTestLRT(tempfit, fit1), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[pos,] <- unlist(modelcomp[2,5:7])
			}
			listFixCon <- c(listFixCon, tryresult)
			temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavTestLRT(tempfit0, fit0), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[pos,] <- unlist(modelcomp0[2,5:7])
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			waldCon[pos,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
			pos <- pos + 1
		}

		
		facinvarfree <- findFactor(varnonfixvar, facList)
		for(i in seq_along(indexnonfixvar)) {
			runnum <- indexnonfixvar[i]
			# Need to change marker variable if fixed
			oldmarker <- fixLoadingFac[[facinvarfree[i]]]
			if(length(oldmarker) > 0 && oldmarker == varnonfixvar[i]) {
				candidatemarker <- setdiff(facList[[facinvarfree[i]]], varnonfixvar[i])[1]
				temp <- freeParTable(pt1, facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups)
				temp <- constrainParTable(temp, facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups)
				temp <- fixParTable(temp, facinvarfree[i], "=~", candidatemarker, 1:ngroups)
				newparent <- freeParTable(pt1, facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups)
				newparent <- fixParTable(newparent, facinvarfree[i], "=~", candidatemarker, 1:ngroups)
				newparentfit <- refit(newparent, fit1)
				
				tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
				if(!is(tryresult, "try-error")) {
					compresult <- try(modelcomp <- lavTestLRT(tempfit, newparentfit), silent = TRUE)
					if(!is(compresult, "try-error"))  fixCon[pos,] <- unlist(modelcomp[2,5:7])
				}
				waldCon[pos,] <- waldConstraint(newparentfit, newparent, waldMat, cbind(facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups))
			} else {
				temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
				tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
				if(!is(tryresult, "try-error")) {
					compresult <- try(modelcomp <- lavTestLRT(tempfit, fit1), silent = TRUE)
					if(!is(compresult, "try-error"))  fixCon[pos,] <- unlist(modelcomp[2,5:7])
				}
				waldCon[pos,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
			}
			listFixCon <- c(listFixCon, tryresult)
			if(length(oldmarker) > 0 && oldmarker == varnonfixvar[i]) {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
			} else {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			}
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavTestLRT(tempfit0, fit0), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[pos,] <- unlist(modelcomp0[2,5:7])
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			pos <- pos + 1
		}

		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- namept1[c(indexfixvar, indexnonfixvar)]
		result <- cbind(freeCon, fixCon, waldCon)		
	} else if (numType == 2) {
		if(!is.null(free) | !is.null(fix)) {
			if(!is.null(fix)) {
				facinfix <- findFactor(fix, facList)
				dup <- duplicated(facinfix)
				for(i in seq_along(fix)) {
					if(dup[i]) {
						pt0 <- constrainParTable(pt0, fix[i], "~1", "", 1:ngroups)
						pt1 <- constrainParTable(pt1, fix[i], "~1", "", 1:ngroups)					
					} else {
						oldmarker <- fixIntceptFac[[facinfix[i]]]
						if(length(oldmarker) > 0) {
							oldmarkerval <- pt1$ustart[pt1$lhs == fix[i] & pt1$op == "~1" & pt1$rhs == "" & pt1$group == 1]
							if(oldmarker == fix[i]) {
								pt0 <- fixParTable(pt0, fix[i], "~1", "", 1:ngroups, oldmarkerval)
								pt1 <- fixParTable(pt1, fix[i], "~1", "", 1:ngroups, oldmarkerval)
							} else {
								pt0 <- freeParTable(pt0, oldmarker, "~1", "", 1:ngroups)
								pt0 <- constrainParTable(pt0, oldmarker, "~1", "", 1:ngroups)
								pt1 <- freeParTable(pt1, oldmarker, "~1", "", 1:ngroups)
								pt0 <- fixParTable(pt0, fix[i], "~1", "", 1:ngroups, oldmarkerval)
								pt1 <- fixParTable(pt1, fix[i], "~1", "", 1:ngroups, oldmarkerval)
								fixIntceptFac[[facinfix[i]]] <- fix[i]
							}
						} else {
							pt0 <- constrainParTable(pt0, fix[i], "~1", "", 1:ngroups)
							pt1 <- constrainParTable(pt1, fix[i], "~1", "", 1:ngroups)						
						}
					}
				}
			}
			if(!is.null(free)) {
				facinfree <- findFactor(free, facList)
				for(i in seq_along(free)) {
					# Need to change marker variable if fixed
					oldmarker <- fixIntceptFac[[facinfree[i]]]
					if(length(oldmarker) > 0 && oldmarker == free[i]) {
						oldmarkerval <- pt1$ustart[pt1$lhs == oldmarker & pt1$op == "~1" & pt1$rhs == "" & pt1$group == 1]
						candidatemarker <- setdiff(facList[[facinfree[i]]], free[i])[1]
						pt0 <- freeParTable(pt0, free[i], "~1", "", 1:ngroups)
						pt1 <- freeParTable(pt1, free[i], "~1", "", 1:ngroups)
						pt0 <- fixParTable(pt0, candidatemarker, "~1", "", 1:ngroups, oldmarkerval)
						pt1 <- fixParTable(pt1, candidatemarker, "~1", "", 1:ngroups, oldmarkerval)
						fixIntceptFac[[facinfix[i]]] <- candidatemarker
					} else {
						pt0 <- freeParTable(pt0, free[i], "~1", "", 1:ngroups)
						pt1 <- freeParTable(pt1, free[i], "~1", "", 1:ngroups)
					}
				}
			}
			namept1 <- paramNameFromPt(pt1)
			namept0 <- paramNameFromPt(pt0)
			fit0 <- refit(pt0, fit0)
			fit1 <- refit(pt1, fit1)
			beta <- coef(fit1)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		fixCon <- freeCon <- waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p")
		colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")
		index <- which((pt1$lhs %in% varfree) & (pt1$op == "~1") & (pt1$group == 1))
		facinfix <- findFactor(fix, facList)
		varinfixvar <- unlist(facList[facinfix])
		varinfixvar <- setdiff(varinfixvar, setdiff(varinfixvar, varfree))
		indexfixvar <- which((pt1$lhs %in% varinfixvar) & (pt1$op == "~1") & (pt1$group == 1))
		varnonfixvar <- setdiff(varfree, varinfixvar)
		indexnonfixvar <- setdiff(index, indexfixvar)
	
		pos <- 1
		for(i in seq_along(varinfixvar)) {
			runnum <- index[i]
			temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
			tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
			if(!is(tryresult, "try-error")) {
				compresult <- try(modelcomp <- lavTestLRT(tempfit, fit1), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[pos,] <- unlist(modelcomp[2,5:7])
			}
			listFixCon <- c(listFixCon, tryresult)
			temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavTestLRT(tempfit0, fit0), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[pos,] <- unlist(modelcomp0[2,5:7])
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			waldCon[pos,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
			pos <- pos + 1
		}

		facinvarfree <- findFactor(varfree, facList)
		for(i in seq_along(varnonfixvar)) {
			runnum <- index[i]
			# Need to change marker variable if fixed
			oldmarker <- fixIntceptFac[[facinvarfree[i]]]
			if(length(oldmarker) > 0 && oldmarker == varfree[i]) {
				candidatemarker <- setdiff(facList[[facinvarfree[i]]], varfree[i])[1]
				temp <- freeParTable(pt1, varfree[i], "~1", "", 1:ngroups)
				temp <- constrainParTable(temp, varfree[i], "~1", "", 1:ngroups)
				temp <- fixParTable(temp, candidatemarker, "~1", "", 1:ngroups)
				newparent <- freeParTable(pt1, varfree[i], "~1", "", 1:ngroups)
				newparent <- fixParTable(newparent, candidatemarker, "~1", "", 1:ngroups)
				newparentfit <- refit(newparent, fit1)
				
				tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
				if(!is(tryresult, "try-error")) {
					compresult <- try(modelcomp <- lavTestLRT(tempfit, newparentfit), silent = TRUE)
					if(!is(compresult, "try-error"))  fixCon[pos,] <- unlist(modelcomp[2,5:7])
				}
				waldCon[pos,] <- waldConstraint(newparentfit, newparent, waldMat, cbind(varfree[i], "~1", "", 1:ngroups))
			} else {
				temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
				tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
				if(!is(tryresult, "try-error")) {
					compresult <- try(modelcomp <- lavTestLRT(tempfit, fit1), silent = TRUE)
					if(!is(compresult, "try-error"))  fixCon[pos,] <- unlist(modelcomp[2,5:7])
				}
				waldCon[pos,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
			}
			listFixCon <- c(listFixCon, tryresult)
			if(length(oldmarker) > 0 && oldmarker == varfree[i]) {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
			} else {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			}
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavTestLRT(tempfit0, fit0), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[pos,] <- unlist(modelcomp0[2,5:7])
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			pos <- pos + 1
		}

		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- namept1[c(indexfixvar, indexnonfixvar)]
		result <- cbind(freeCon, fixCon, waldCon)		
	} else if (numType == 3) {
		if(!is.null(free) | !is.null(fix)) {
			if(!is.null(fix)) {
				for(i in seq_along(fix)) {
					pt0 <- constrainParTable(pt0, fix[i], "~~", fix[i], 1:ngroups)
					pt1 <- constrainParTable(pt1, fix[i], "~~", fix[i], 1:ngroups)					
				}
			}
			if(!is.null(free)) {
				for(i in seq_along(free)) {
					pt0 <- freeParTable(pt0, free[i], "~~", free[i], 1:ngroups)
					pt1 <- freeParTable(pt1, free[i], "~~", free[i], 1:ngroups)
				}
			}
			namept1 <- paramNameFromPt(pt1)
			namept0 <- paramNameFromPt(pt0)
			fit0 <- refit(pt0, fit0)
			fit1 <- refit(pt1, fit1)
			beta <- coef(fit1)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		fixCon <- freeCon <- waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p")
		colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")
		index <- which((pt1$lhs %in% varfree) & (pt1$op == "~~") & (pt1$lhs == pt1$rhs) & (pt1$group == 1))
		for(i in seq_along(index)) {
			runnum <- index[i]
			temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
			tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
			if(!is(tryresult, "try-error")) {
				compresult <- try(modelcomp <- lavTestLRT(tempfit, fit1), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[i,] <- unlist(modelcomp[2,5:7])
			}
			listFixCon <- c(listFixCon, tryresult)
			temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavTestLRT(tempfit0, fit0), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[i,] <- unlist(modelcomp0[2,5:7])
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			waldCon[i,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
		}
		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- namept1[index]
		result <- cbind(freeCon, fixCon, waldCon)		
	} else if (numType == 4) {
		varfree <- facnames
		if(!is.null(free) | !is.null(fix)) {
			if(!is.null(fix)) {
				for(i in seq_along(fix)) {
					pt0 <- constrainParTable(pt0, fix[i], "~1", "", 1:ngroups)
					pt1 <- constrainParTable(pt1, fix[i], "~1", "", 1:ngroups)					
				}
			}
			if(!is.null(free)) {
				for(i in seq_along(free)) {
					pt0 <- freeParTable(pt0, free[i], "~1", "", 1:ngroups)
					pt1 <- freeParTable(pt1, free[i], "~1", "", 1:ngroups)
				}
			}
			namept1 <- paramNameFromPt(pt1)
			namept0 <- paramNameFromPt(pt0)
			fit0 <- refit(pt0, fit0)
			fit1 <- refit(pt1, fit1)
			beta <- coef(fit1)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		fixCon <- freeCon <- waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p")
		colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")
		index <- which((pt1$lhs %in% varfree) & (pt1$op == "~1") & (pt1$group == 1))
		for(i in seq_along(index)) {
			runnum <- index[i]
			isfree <- pt1$free[runnum] != 0
			if(isfree) {
				temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
			} else {
				temp <- fixParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 2:ngroups, ustart = pt1$ustart[runnum])
			}
			tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
			if(!is(tryresult, "try-error")) {
				compresult <- try(modelcomp <- lavTestLRT(tempfit, fit1), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[i,] <- unlist(modelcomp[2,5:7])
			}
			listFixCon <- c(listFixCon, tryresult)
			isfree0 <- pt0$free[runnum] != 0
			if(isfree0) {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			} else {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
			}
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavTestLRT(tempfit0, fit0), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[i,] <- unlist(modelcomp0[2,5:7])
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			waldCon[i,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
		}
		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- namept1[index]
		result <- cbind(freeCon, fixCon, waldCon)		
	}
	if(return.fit) {
		return(invisible(list(result = result, models = list(free = listFreeCon, fix = listFixCon))))
	} else {
		return(result)
	}
}

findFactor <- function(var, facList) {
	tempfac <- lapply(facList, intersect, var)
	facinvar <- rep(names(tempfac), sapply(tempfac, length))
	facinvar[match(unlist(tempfac), var)]
}

waldConstraint <- function(fit, pt, mat, ...) {
	dotdotdot <- list(...)
	overallMat <- NULL
	for(i in seq_along(dotdotdot)) {
		target <- dotdotdot[[i]]
		tempMat <- mat
		element <- apply(target, 1, matchElement, parTable=pt)
		freeIndex  <- pt$free[element]
		tempMat[,freeIndex[1]] <- -1
		for(m in 2:length(freeIndex)) {
			tempMat[m - 1, freeIndex[m]] <- 1
		}
		overallMat <- rbind(overallMat, tempMat)
	}
	if(any(apply(overallMat, 1, sum) != 0)) {
		return(rep(NA, 3))
	} else {
		return(waldContrast(fit, overallMat))
	}
}
