# Wald stat did not show up

partialInvarianceCat <- function(fit, type, free = NULL, fix = NULL, refgroup = 1, poolvar = TRUE, p.adjust = "none", return.fit = FALSE, method = "satorra.bentler.2001") { 
	# model <- ' f1 =~ u1 + u2 + u3 + u4
			   # f2 =~ u5 + u6 + u7 + u8'

	# modelsCat2 <- measurementInvarianceCat(model, data = datCat, group = "g", parameterization="theta", 
		# estimator="wlsmv", strict = TRUE)
	# fit <- modelsCat2
	# type <- "weak"
	# free <- NULL
	# fix <- NULL
	# refgroup <- 1
	# poolvar <- TRUE
	# p.adjust <- "none"
	# return.fit <- FALSE
	# method = "satorra.bentler.2001"

	type <- tolower(type)
	numType <- 1
	fit1 <- fit0 <- NULL
	# fit0 = Nested model, fit1 = Parent model
	if(type %in% c("metric", "weak", "loading", "loadings")) {
		numType <- 1
		if(all(c("fit.configural", "fit.loadings") %in% names(fit))) {
			fit1 <- fit$fit.configural
			fit0 <- fit$fit.loadings
		} else {
			stop("The elements named 'fit.configural' and 'fit.loadings' are needed in the 'fit' argument")
		}
	} else if (type %in% c("scalar", "strong", "intercept", "intercepts", "threshold", "thresholds")) {
		numType <- 2
		if(all(c("fit.loadings", "fit.thresholds") %in% names(fit))) {
			fit1 <- fit$fit.loadings
			fit0 <- fit$fit.thresholds
		} else {
			stop("The elements named 'fit.loadings' and 'fit.thresholds' are needed in the 'fit' argument")
		}
	} else if (type %in% c("strict", "residual", "residuals", "error", "errors")) {
		numType <- 3
		if("fit.residuals" %in% names(fit)) {
			fit0 <- fit$fit.residuals
			if("fit.thresholds" %in% names(fit)) {
				fit1 <- fit$fit.thresholds			
			} else if ("fit.loadings" %in% names(fit)) {
				fit1 <- fit$fit.loadings			
			} else {
				stop("The element named either 'fit.thresholds' or 'fit.loadings' is needed in the 'fit' argument")
			}			
		} else {
			stop("The element named 'fit.residuals' is needed in the 'fit' argument")
		}
	} else if (type %in% c("means", "mean")) {
		numType <- 4
		if("fit.means" %in% names(fit)) {
			fit0 <- fit$fit.means
			if("fit.residuals" %in% names(fit)) {
				fit1 <- fit$fit.residuals			
			} else if ("fit.thresholds" %in% names(fit)) {
				fit1 <- fit$fit.thresholds			
			} else if ("fit.loadings" %in% names(fit)) {
				fit1 <- fit$fit.loadings			
			} else {
				stop("The element named either 'fit.residuals', 'fit.thresholds', or 'fit.loadings' is needed in the 'fit' argument")
			}
		} else {
			stop("The element named 'fit.means' is needed in the 'fit' argument")
		}
	} else {
		stop("Please specify the correct type of measurement invariance. See the help page.")
	}
	pt1 <- lavaan::partable(fit1)
	pt0 <- lavaan::partable(fit0)
	pt0$start <- pt0$est <- pt0$se <- NULL
	pt1$start <- pt1$est <- pt1$se <- NULL
	
	pt1$label[substr(pt1$label, 1, 1) == "." & substr(pt1$label, nchar(pt1$label), nchar(pt1$label)) == "."] <- ""
	pt0$label[substr(pt0$label, 1, 1) == "." & substr(pt0$label, nchar(pt0$label), nchar(pt0$label)) == "."] <- ""
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
	
	# Find the number of thresholds
	# Check whether the factor configuration is the same across gorups
	
	conParTable <- lapply(pt1, "[", pt1$op == "==")
	group1pt <- lapply(pt1, "[", pt1$group != 1)
	
	numThreshold <- table(sapply(group1pt, "[", group1pt$op == "|")[,"lhs"])
	plabelthres <- split(group1pt$plabel[group1pt$op == "|"], group1pt$lhs[group1pt$op == "|"])
	numFixedThreshold <- sapply(lapply(plabelthres, function(vec) !is.na(match(vec, conParTable$lhs)) | !is.na(match(vec, conParTable$rhs))), sum)[names(numThreshold)]  
	
	#numFixedThreshold <- table(sapply(group1pt, "[", group1pt$op == "|" & group1pt$eq.id != 0)[,"lhs"])
	fixIntceptFac <- list()
	for(i in seq_along(facList)) {
		tmp <- numFixedThreshold[facList[[i]]]
		if(all(tmp > 1)) {
			fixIntceptFac[[i]] <- integer(0)
		} else {
			fixIntceptFac[[i]] <- names(which.max(tmp))[1]
		}
	}
	names(fixIntceptFac) <- names(facList)
	
	ngroups <- max(pt0$group)
	neach <- lavaan::lavInspect(fit0, "nobs")
	groupvar <- lavaan::lavInspect(fit0, "group")
	grouplab <- lavaan::lavInspect(fit0, "group.label")
	if(!is.numeric(refgroup)) refgroup <- which(refgroup == grouplab)
	grouporder <- 1:ngroups
	grouporder <- c(refgroup, setdiff(grouporder, refgroup))
	grouplaborder <- grouplab[grouporder]
	complab <- paste(grouplaborder[2:ngroups], "vs.", grouplaborder[1])
	if(ngroups <= 1) stop("Well, the number of groups is 1. Measurement invariance across 'groups' cannot be done.")

	if(numType == 4) {
		if(!all(c(free, fix) %in% facnames)) stop("'free' and 'fix' arguments should consist of factor names because mean invariance is tested.")
	} else {
		if(!all(c(free, fix) %in% varnames)) stop("'free' and 'fix' arguments should consist of variable names.")
	}
	result <- fixCon <- freeCon <- NULL
	estimates <- NULL
	listFreeCon <- listFixCon <- list()	
	beta <- lavaan::coef(fit1)
	beta0 <- lavaan::coef(fit0)
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
			beta <- lavaan::coef(fit1)
			beta0 <- lavaan::coef(fit0)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		estimates <- matrix(NA, length(varfree), ngroups + 1)
		stdestimates <- matrix(NA, length(varfree), ngroups)
		colnames(estimates) <- c("poolest", paste0("load:", grouplab))
		colnames(stdestimates) <- paste0("std:", grouplab)
		esstd <- esz <- matrix(NA, length(varfree), ngroups - 1)
		colnames(esstd) <- paste0("diff_std:", complab)
		colnames(esz) <- paste0("q:", complab)
		fixCon <- freeCon <- matrix(NA, length(varfree), 4)
		waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
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
				compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
			}
			listFixCon <- c(listFixCon, tryresult)
			temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			estimates[pos, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				loadVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[pos, 2:ncol(estimates)] <- loadVal
				facVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~~", pt0$lhs[runnum], 1:ngroups)
				totalVal <- sapply(thetaImpliedTotalVar(tempfit0), function(x, v) x[v, v], v = pt0$rhs[runnum])
				names(facVal) <- names(totalVal) <- grouplab
				ifelse(poolvar, refFacVal <- poolVariance(facVal, neach), refFacVal <- facVal[refgroup])
				ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
				stdLoadVal <- loadVal * sqrt(refFacVal) / sqrt(refTotalVal)
				stdestimates[pos,] <- stdLoadVal
				stdLoadVal <- stdLoadVal[grouporder]
				esstd[pos,] <- stdLoadVal[2:ngroups] - stdLoadVal[1]
				if(any(abs(stdLoadVal) > 0.9999)) warning(paste("Standardized Loadings of", pt0$rhs[runnum], "in some groups are less than -1 or over 1. The standardized loadings used in Fisher z transformation are changed to -0.9999 or 0.9999."))
				stdLoadVal[stdLoadVal > 0.9999] <- 0.9999
				stdLoadVal[stdLoadVal < -0.9999] <- -0.9999
				zLoadVal <- atanh(stdLoadVal)
				esz[pos,] <- zLoadVal[2:ngroups] - zLoadVal[1]
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
				newparentresult <- try(newparentfit <- refit(newparent, fit1), silent = TRUE)
				if(!is(newparentresult, "try-error")) {
					tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
					if(!is(tryresult, "try-error")) {
						compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, newparentfit, method = method), silent = TRUE)
						if(!is(compresult, "try-error")) fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(newparentfit, tempfit))
					}
					waldCon[pos,] <- waldConstraint(newparentfit, newparent, waldMat, cbind(facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups))
				}
			} else {
				temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
				tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
				if(!is(tryresult, "try-error")) {
					compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
					if(!is(compresult, "try-error"))  fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
				}
				waldCon[pos,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
			}
			listFixCon <- c(listFixCon, tryresult)
			if(length(oldmarker) > 0 && oldmarker == varnonfixvar[i]) {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
			} else {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			}
			estimates[pos, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error")) freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				loadVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[pos, 2:ncol(estimates)] <- loadVal
				facVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~~", pt0$lhs[runnum], 1:ngroups)
				totalVal <- sapply(thetaImpliedTotalVar(tempfit0), function(x, v) x[v, v], v = pt0$rhs[runnum])
				names(facVal) <- names(totalVal) <- grouplab
				ifelse(poolvar, refFacVal <- poolVariance(facVal, neach), refFacVal <- facVal[refgroup])
				ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
				stdLoadVal <- loadVal * sqrt(refFacVal) / sqrt(refTotalVal)
				stdestimates[pos,] <- stdLoadVal
				stdLoadVal <- stdLoadVal[grouporder]
				esstd[pos,] <- stdLoadVal[2:ngroups] - stdLoadVal[1]
				if(any(abs(stdLoadVal) > 0.9999)) warning(paste("Standardized Loadings of", pt0$rhs[runnum], "in some groups are less than -1 or over 1. The standardized loadings used in Fisher z transformation are changed to -0.9999 or 0.9999."))
				stdLoadVal[stdLoadVal > 0.9999] <- 0.9999
				stdLoadVal[stdLoadVal < -0.9999] <- -0.9999
				zLoadVal <- atanh(stdLoadVal)
				esz[pos,] <- zLoadVal[2:ngroups] - zLoadVal[1]
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			pos <- pos + 1
		}
		freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
		fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
		waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)

		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[c(indexfixvar, indexnonfixvar)]
		estimates <- cbind(estimates, stdestimates, esstd, esz)
		result <- cbind(freeCon, fixCon, waldCon)		
	} else if (numType == 2) {
		if(!is.null(free) | !is.null(fix)) {
			if(!is.null(fix)) {
				facinfix <- findFactor(fix, facList)
				dup <- duplicated(facinfix)
				for(i in seq_along(fix)) {
					numfixthres <- numThreshold[fix[i]]
					if(numfixthres > 1) {
						if(dup[i]) {
							for(s in 2:numfixthres) {
								pt0 <- constrainParTable(pt0, fix[i], "|", paste0("t", s), 1:ngroups)
								pt1 <- constrainParTable(pt1, fix[i], "|", paste0("t", s), 1:ngroups)	
							}
						} else {
							oldmarker <- fixIntceptFac[[facinfix[i]]]
							numoldthres <- numThreshold[oldmarker]
							if(length(oldmarker) > 0) {
								if(oldmarker == fix[i]) {
									for(s in 2:numfixthres) {
										pt0 <- constrainParTable(pt0, fix[i], "|", paste0("t", s), 1:ngroups)
										pt1 <- constrainParTable(pt1, fix[i], "|", paste0("t", s), 1:ngroups)			
									}	
								} else {
									for(r in 2:numoldthres) {
										pt1 <- freeParTable(pt1, oldmarker, "|", paste0("t", r), 1:ngroups)									
									}
									for(s in 2:numfixthres) {
										pt0 <- constrainParTable(pt0, fix[i], "|", paste0("t", s), 1:ngroups)		
										pt1 <- constrainParTable(pt1, fix[i], "|", paste0("t", s), 1:ngroups)	
									}	
									fixIntceptFac[[facinfix[i]]] <- fix[i]
								}
							} else {
								for(s in 2:numfixthres) {
									pt0 <- constrainParTable(pt0, fix[i], "|", paste0("t", s), 1:ngroups)
									pt1 <- constrainParTable(pt1, fix[i], "|", paste0("t", s), 1:ngroups)			
								}				
							}
						}
					}
				}
			}
			if(!is.null(free)) {
				facinfree <- findFactor(free, facList)
				for(i in seq_along(free)) {
					numfreethres <- numThreshold[free[i]]
					# Need to change marker variable if fixed
					oldmarker <- fixIntceptFac[[facinfree[i]]]
					numoldthres <- numThreshold[oldmarker]
					if(length(oldmarker) > 0 && oldmarker == free[i]) {
						candidatemarker <- setdiff(facList[[facinfree[i]]], free[i])
						candidatemarker <- candidatemarker[numThreshold[candidatemarker] > 1][1]
						numcandidatethres <- numThreshold[candidatemarker]
						pt0 <- constrainParTable(pt0, candidatemarker, "|", "t2", 1:ngroups)
						pt1 <- constrainParTable(pt1, candidatemarker, "|", "t2", 1:ngroups)
						for(s in 2:numfixthres) {
							pt0 <- freeParTable(pt0, free[i], "|", paste0("t", s), 1:ngroups)
							pt1 <- freeParTable(pt1, free[i], "|", paste0("t", s), 1:ngroups)			
						}	
						fixIntceptFac[[facinfix[i]]] <- candidatemarker
					} else {
						for(s in 2:numfixthres) {
							pt0 <- freeParTable(pt0, free[i], "|", paste0("t", s), 1:ngroups)
							pt1 <- freeParTable(pt1, free[i], "|", paste0("t", s), 1:ngroups)			
						}	
					}
				}
			}
			namept1 <- paramNameFromPt(pt1)
			namept0 <- paramNameFromPt(pt0)
			fit0 <- refit(pt0, fit0)
			fit1 <- refit(pt1, fit1)
			beta <- lavaan::coef(fit1)
			beta0 <- lavaan::coef(fit0)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		maxcolumns <- max(numThreshold[varfree]) - 1
		tname <- paste0("t", 2:(maxcolumns + 1))
		estimates <- matrix(NA, length(varfree), (ngroups * length(tname)) + length(tname))
		stdestimates <- matrix(NA, length(varfree), ngroups * length(tname))
		tnameandlab <- expand.grid(tname, grouplab)
		colnames(estimates) <- c(paste0("pool:", tname), paste0(tnameandlab[,1], ":", tnameandlab[,2]))
		colnames(stdestimates) <- paste0("std:", tnameandlab[,1], ":", tnameandlab[,2])
		esstd <- matrix(NA, length(varfree), (ngroups - 1)* length(tname))
		tnameandcomplab <- expand.grid(tname, complab)
		colnames(esstd) <- paste0("diff_std:", tnameandcomplab[,1], ":", tnameandcomplab[,2])
		fixCon <- freeCon <- matrix(NA, length(varfree), 4)
		waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
		colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")

		facinfix <- findFactor(fix, facList)
		varinfixvar <- unlist(facList[facinfix])
		varinfixvar <- setdiff(varinfixvar, setdiff(varinfixvar, varfree))
		varnonfixvar <- setdiff(varfree, varinfixvar)
		
		pos <- 1
		for(i in seq_along(varinfixvar)) {
			temp <- pt1
			for(s in 2:numThreshold[varinfixvar[i]]) {
				runnum <- which((pt1$lhs == varfree[i]) & (pt1$op == "|") & (pt1$rhs == paste0("t", s)) & (pt1$group == 1))
				temp <- constrainParTable(temp, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
			}
			tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
			if(!is(tryresult, "try-error")) {
				compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
			}
			listFixCon <- c(listFixCon, tryresult)
			temp0 <- pt0
			for(s in 2:numThreshold[varinfixvar[i]]) {
				runnum <- which((pt0$lhs == varfree[i]) & (pt0$op == "|") & (pt0$rhs == paste0("t", s)) & (pt0$group == 1))
				temp0 <- freeParTable(temp0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[pos, s - 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			}			
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				for(s in 2:numThreshold[varinfixvar[i]]) {
					runnum <- which((pt0$lhs == varfree[i]) & (pt0$op == "|") & (pt0$rhs == paste0("t", s)) & (pt0$group == 1))
					thresVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
					estimates[pos, maxcolumns*(1:ngroups) + (s - 1)] <- thresVal
					totalVal <- sapply(thetaImpliedTotalVar(tempfit0), function(x, v) x[v, v], v = pt0$lhs[runnum])
					ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
					stdIntVal <- thresVal / sqrt(refTotalVal)
					stdestimates[pos, maxcolumns*(1:ngroups - 1) + (s - 1)] <- stdIntVal
					stdIntVal <- stdIntVal[grouporder]
					esstd[pos, maxcolumns*(1:length(complab) - 1) + (s - 1)] <- stdIntVal[2:ngroups] - stdIntVal[1]
				}
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			args <- list(fit1, pt1, waldMat)
			for(s in 2:numThreshold[varinfixvar[i]]) {
				runnum <- which((pt1$lhs == varfree[i]) & (pt1$op == "|") & (pt1$rhs == paste0("t", s)) & (pt1$group == 1))
				args <- c(args, list(cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)))
			}
			waldCon[pos,] <- do.call(waldConstraint, args)
			pos <- pos + 1
		}

		facinvarfree <- findFactor(varnonfixvar, facList)
		for(i in seq_along(varnonfixvar)) {
			# Need to change marker variable if fixed
			oldmarker <- fixIntceptFac[[facinvarfree[i]]]
			if(length(oldmarker) > 0 && oldmarker == varfree[i]) {
				candidatemarker <- setdiff(facList[[facinvarfree[i]]], varnonfixvar[i])
				candidatemarker <- candidatemarker[numThreshold[candidatemarker] > 1][1]
				numcandidatethres <- numThreshold[candidatemarker]
				newparent <- constrainParTable(pt1, candidatemarker, "|", "t2", 1:ngroups)
				for(s in 2:numcandidatethres) {
					newparent <- freeParTable(newparent, varnonfixvar[i], "|", paste0("t", s), 1:ngroups)		
				}	
				temp <- newparent
				for(s in 2:numThreshold[varnonfixvar[i]]) {
					runnum <- which((newparent$lhs == varnonfixvar[i]) & (newparent$op == "|") & (newparent$rhs == paste0("t", s)) & (newparent$group == 1))
					temp <- constrainParTable(temp, newparent$lhs[runnum], newparent$op[runnum], newparent$rhs[runnum], 1:ngroups)
				}
				newparentresult <- try(newparentfit <- refit(newparent, fit1), silent = TRUE)
				if(!is(newparentresult, "try-error")) {
					tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
					if(!is(tryresult, "try-error")) {
						compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, newparentfit, method = method), silent = TRUE)
						if(!is(compresult, "try-error")) fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(newparentfit, tempfit))
					}
					args <- list(newparentfit, newparent, waldMat)
					for(s in 2:numThreshold[varnonfixvar[i]]) {
						runnum <- which((newparent$lhs == varnonfixvar[i]) & (newparent$op == "|") & (newparent$rhs == paste0("t", s)) & (newparent$group == 1))
						args <- c(args, list(cbind(newparent$lhs[runnum], newparent$op[runnum], newparent$rhs[runnum], 1:ngroups)))
					}
					waldCon[pos,] <- do.call(waldConstraint, args)
				}
			} else {
				temp <- pt1
				for(s in 2:numThreshold[varnonfixvar[i]]) {
					runnum <- which((pt1$lhs == varfree[i]) & (pt1$op == "|") & (pt1$rhs == paste0("t", s)) & (pt1$group == 1))
					temp <- constrainParTable(temp, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
				}
				tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
				if(!is(tryresult, "try-error")) {
					compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
					if(!is(compresult, "try-error"))  fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
				}
				args <- list(fit1, pt1, waldMat)
				for(s in 2:numThreshold[varnonfixvar[i]]) {
					runnum <- which((pt1$lhs == varfree[i]) & (pt1$op == "|") & (pt1$rhs == paste0("t", s)) & (pt1$group == 1))
					args <- c(args, list(cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)))
				}
				waldCon[pos,] <- do.call(waldConstraint, args)
			}
			listFixCon <- c(listFixCon, tryresult)
			
			temp0 <- pt0
			for(s in 2:numThreshold[varnonfixvar[i]]) {
				runnum <- which((pt0$lhs == varfree[i]) & (pt0$op == "|") & (pt0$rhs == paste0("t", s)) & (pt0$group == 1))
				temp0 <- freeParTable(temp0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[pos, s - 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			}
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				for(s in 2:numThreshold[varnonfixvar[i]]) {
					runnum <- which((pt0$lhs == varfree[i]) & (pt0$op == "|") & (pt0$rhs == paste0("t", s)) & (pt0$group == 1))
					thresVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
					estimates[pos, maxcolumns*(1:ngroups) + (s - 1)] <- thresVal
					totalVal <- sapply(thetaImpliedTotalVar(tempfit0), function(x, v) x[v, v], v = pt0$lhs[runnum])
					ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
					stdIntVal <- thresVal / sqrt(refTotalVal)
					stdestimates[pos, maxcolumns*(1:ngroups - 1) + (s - 1)] <- stdIntVal
					stdIntVal <- stdIntVal[grouporder]
					esstd[pos, maxcolumns*(1:length(complab) - 1) + (s - 1)] <- stdIntVal[2:ngroups] - stdIntVal[1]
				}
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			pos <- pos + 1
		}
		freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
		fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
		waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)
		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- paste0(c(varinfixvar, varnonfixvar), "|")
		estimates <- cbind(estimates, stdestimates, esstd)
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
			beta <- lavaan::coef(fit1)
			beta0 <- lavaan::coef(fit0)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		estimates <- matrix(NA, length(varfree), ngroups + 1)
		stdestimates <- matrix(NA, length(varfree), ngroups)
		colnames(estimates) <- c("poolest", paste0("errvar:", grouplab))
		colnames(stdestimates) <- paste0("std:", grouplab)
		esstd <- esz <- matrix(NA, length(varfree), ngroups - 1)
		colnames(esstd) <- paste0("diff_std:", complab)
		colnames(esz) <- paste0("h:", complab)
		fixCon <- freeCon <- matrix(NA, length(varfree), 4)
		waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
		colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")
		index <- which((pt1$lhs %in% varfree) & (pt1$op == "~~") & (pt1$lhs == pt1$rhs) & (pt1$group == 1))
		for(i in seq_along(index)) {
			runnum <- index[i]
			ustart <- getValue(pt1, beta, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1)
			temp <- fixParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 2:ngroups, ustart)
			tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
			if(!is(tryresult, "try-error")) {
				compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[i,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
			}
			listFixCon <- c(listFixCon, tryresult)
			temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
			estimates[i, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[i,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				errVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[i, 2:ncol(estimates)] <- errVal
				totalVal <- sapply(thetaImpliedTotalVar(tempfit0), function(x, v) x[v, v], v = pt0$rhs[runnum])
				ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
				stdErrVal <- errVal / sqrt(refTotalVal)
				stdestimates[i,] <- stdErrVal
				stdErrVal <- stdErrVal[grouporder]
				esstd[i,] <- stdErrVal[2:ngroups] - stdErrVal[1]
				if(any(abs(stdErrVal) > 0.9999)) warning(paste("The uniqueness of", pt0$rhs[runnum], "in some groups are over 1. The uniqueness used in arctan transformation are changed to 0.9999."))
				stdErrVal[stdErrVal > 0.9999] <- 0.9999
				zErrVal <- asin(sqrt(stdErrVal))
				esz[i,] <- zErrVal[2:ngroups] - zErrVal[1]
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			waldCon[i,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
		}
		freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
		fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
		waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)
		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[index]
		estimates <- cbind(estimates, stdestimates, esstd, esz)
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
			beta <- lavaan::coef(fit1)
			beta0 <- lavaan::coef(fit0)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		estimates <- matrix(NA, length(varfree), ngroups + 1)
		stdestimates <- matrix(NA, length(varfree), ngroups)
		colnames(estimates) <- c("poolest", paste0("mean:", grouplab))
		colnames(stdestimates) <- paste0("std:", grouplab)
		esstd <- matrix(NA, length(varfree), ngroups - 1)
		colnames(esstd) <- paste0("diff_std:", complab)
		fixCon <- freeCon <- matrix(NA, length(varfree), 4)
		waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
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
				compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[i,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
			}
			listFixCon <- c(listFixCon, tryresult)
			isfree0 <- pt0$free[runnum] != 0
			if(isfree0) {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			} else {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
			}
			estimates[i, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[i,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				meanVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[i, 2:ncol(estimates)] <- meanVal
				facVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~~", pt0$lhs[runnum], 1:ngroups)
				ifelse(poolvar, refFacVal <- poolVariance(facVal, neach), refFacVal <- facVal[refgroup])
				stdMeanVal <- meanVal / sqrt(refFacVal)
				stdestimates[i,] <- stdMeanVal
				stdMeanVal <- stdMeanVal[grouporder]
				esstd[i,] <- stdMeanVal[2:ngroups] - stdMeanVal[1]
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			waldCon[i,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
		}
		freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
		fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
		waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)
		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[index]
		estimates <- cbind(estimates, stdestimates, esstd)
		result <- cbind(freeCon, fixCon, waldCon)			
	}
	if(return.fit) {
		return(invisible(list(estimates = estimates, results = result, models = list(free = listFreeCon, fix = listFixCon, nested = fit0, parent = fit1))))
	} else {
		return(list(estimates = estimates, results = result))
	}
}

thetaImpliedTotalVar <- function(object) {
	param <- lavaan::lavInspect(object, "coef")
	ngroup <- lavaan::lavInspect(object, "ngroups")
	name <- names(param)
	if(ngroup == 1) {
		ly <- param[name == "lambda"]
	} else {
		ly <- lapply(param, "[[", "lambda")
	}
	ps <- lavaan::lavInspect(object, "cov.lv")
	if(ngroup == 1) ps <- list(ps)
	if(ngroup == 1) {
		te <- param[name == "theta"]
	} else {
		te <- lapply(param, "[[", "theta")
	}
	result <- list()
	for(i in 1:ngroup) {
		result[[i]] <- ly[[i]]%*%ps[[i]]%*%t(ly[[i]]) + te[[i]]
	}
	result
}
