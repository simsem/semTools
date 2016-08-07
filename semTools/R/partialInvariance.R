# Work with only with congeneric models

partialInvariance <- function(fit, type, free = NULL, fix = NULL, refgroup = 1, poolvar = TRUE, p.adjust = "none", fbound = 2, return.fit = FALSE, method = "satorra.bentler.2001") { 
	# fit <- measurementInvariance(HW.model, data=HolzingerSwineford1939, group="school", strict = TRUE)
	# type <- "weak"
	# free <- NULL
	# fix <- "x1"
	# refgroup <- 1
	# poolvar <- TRUE
	# p.adjust <- "none"
	# return.fit <- FALSE
	# fbound <- 2
	# method <- "satorra.bentler.2001"


	type <- tolower(type)
	numType <- 0
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
		if(all(c("fit.loadings", "fit.intercepts") %in% names(fit))) {
			fit1 <- fit$fit.loadings
			fit0 <- fit$fit.intercepts
		} else {
			stop("The elements named 'fit.loadings' and 'fit.intercepts' are needed in the 'fit' argument")
		}
	} else if (type %in% c("strict", "residual", "residuals", "error", "errors")) {
		numType <- 3
		if(all(c("fit.intercepts", "fit.residuals") %in% names(fit))) {
			fit1 <- fit$fit.intercepts
			fit0 <- fit$fit.residuals
		} else {
			stop("The elements named 'fit.intercepts' and 'fit.residuals' are needed in the 'fit' argument")
		}
	} else if (type %in% c("means", "mean")) {
		numType <- 4
		if("fit.means" %in% names(fit)) {
			fit0 <- fit$fit.means
			if("fit.residuals" %in% names(fit)) {
				fit1 <- fit$fit.residuals		
			} else if ("fit.intercepts" %in% names(fit)) {
				fit1 <- fit$fit.intercepts			
			} else {
				stop("The elements named either 'fit.residuals' or 'fit.intercepts	' is needed in the 'fit' argument")
			}
		} else {
			stop("The elements named 'fit.means' is needed in the 'fit' argument")
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
	fixIntceptFac <- list()
	for(i in seq_along(facList)) {
		select <- pt1$op == "~1" & pt1$rhs %in% facList[[i]] & pt1$group == 1 & pt1$free == 0
		fixIntceptFac[[i]] <- pt1$rhs[select]
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
		
		obsmean <- sapply(lavaan::lavInspect(fit0, "sampstat"), "[[", "mean")
		obsmean <- obsmean[,grouporder]
		obsdiff <- obsmean[,2:ngroups, drop = FALSE] - matrix(obsmean[,1], nrow(obsmean), ngroups - 1)
		obsdiff <- obsdiff[varfree, , drop = FALSE]
		colnames(obsdiff) <- paste0("diff_mean:", complab)
		
		estimates <- matrix(NA, length(varfree), ngroups + 1)
		stdestimates <- matrix(NA, length(varfree), ngroups)
		colnames(estimates) <- c("poolest", paste0("load:", grouplab))
		colnames(stdestimates) <- paste0("std:", grouplab)
		esstd <- esz <- matrix(NA, length(varfree), ngroups - 1)
		colnames(esstd) <- paste0("diff_std:", complab)
		colnames(esz) <- paste0("q:", complab)
		esdiff <- matrix(NA, length(varfree), ngroups - 1)
		
		# Extract facmean, facsd, load, tau -> lowdiff, highdiff
		lowdiff <- matrix(NA, length(varfree), ngroups - 1)
		highdiff <- matrix(NA, length(varfree), ngroups - 1)
		colnames(lowdiff) <- paste0("low_fscore:", complab)
		colnames(highdiff) <- paste0("high_fscore:", complab)
		
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
				totalVal <- sapply(lavaan::fitted.values(tempfit0), function(x, v) x$cov[v, v], v = pt0$rhs[runnum])
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
				
				facMean <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~1", "", 1:ngroups)
				wlow <- min(facMean - fbound * sqrt(facVal))
				whigh <- max(facMean + fbound * sqrt(facVal))
				intVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$rhs[runnum], "~1", "", 1:ngroups)
				loadVal <- loadVal[grouporder]
				intVal <- intVal[grouporder]
				loaddiff <- loadVal[2:ngroups] - loadVal[1]
				intdiff <- intVal[2:ngroups] - intVal[1]
				lowdiff[pos,] <- intdiff + wlow * loaddiff
				highdiff[pos,] <- intdiff + whigh * loaddiff				
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
				temp <- fixParTable(temp, facinvarfree[i], "=~", candidatemarker, 1:ngroups, ustart = 1)
				temp <- constrainParTable(temp, facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups)
				newparent <- freeParTable(pt1, facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups)
				newparent <- fixParTable(newparent, facinvarfree[i], "=~", candidatemarker, 1:ngroups, ustart = 1)
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
				totalVal <- sapply(lavaan::fitted.values(tempfit0), function(x, v) x$cov[v, v], v = pt0$rhs[runnum])
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
				
				facMean <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~1", "", 1:ngroups)
				wlow <- min(facMean - fbound * sqrt(facVal))
				whigh <- max(facMean + fbound * sqrt(facVal))
				intVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$rhs[runnum], "~1", "", 1:ngroups)
				loadVal <- loadVal[grouporder]
				intVal <- intVal[grouporder]
				loaddiff <- loadVal[2:ngroups] - loadVal[1]
				intdiff <- intVal[2:ngroups] - intVal[1]
				lowdiff[pos,] <- intdiff + wlow * loaddiff
				highdiff[pos,] <- intdiff + whigh * loaddiff
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			pos <- pos + 1
		}
		freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
		fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
		waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)
		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[c(indexfixvar, indexnonfixvar)]
		estimates <- cbind(estimates, stdestimates, esstd, esz, obsdiff, lowdiff, highdiff)
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
			beta <- lavaan::coef(fit1)
			beta0 <- lavaan::coef(fit0)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		obsmean <- sapply(lavaan::lavInspect(fit0, "sampstat"), "[[", "mean")
		obsmean <- obsmean[,grouporder]
		obsdiff <- obsmean[,2:ngroups, drop = FALSE] - matrix(obsmean[,1], nrow(obsmean), ngroups - 1)
		obsdiff <- obsdiff[varfree, , drop = FALSE]
		colnames(obsdiff) <- paste0("diff_mean:", complab)

		# Prop diff
		propdiff <- matrix(NA, length(varfree), ngroups - 1)
		colnames(propdiff) <- paste0("propdiff:", complab)

		estimates <- matrix(NA, length(varfree), ngroups + 1)
		stdestimates <- matrix(NA, length(varfree), ngroups)
		colnames(estimates) <- c("poolest", paste0("int:", grouplab))
		colnames(stdestimates) <- paste0("std:", grouplab)
		esstd <- matrix(NA, length(varfree), ngroups - 1)
		colnames(esstd) <- paste0("diff_std:", complab)
		fixCon <- freeCon <- matrix(NA, length(varfree), 4)
		waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
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
				intVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[pos, 2:ncol(estimates)] <- intVal
				totalVal <- sapply(lavaan::fitted.values(tempfit0), function(x, v) x$cov[v, v], v = pt0$lhs[runnum])
				ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
				stdIntVal <- intVal / sqrt(refTotalVal)
				stdestimates[pos,] <- stdIntVal
				stdIntVal <- stdIntVal[grouporder]
				esstd[pos,] <- stdIntVal[2:ngroups] - stdIntVal[1]
				
				intVal <- intVal[grouporder]
				propdiff[pos,] <- (intVal[2:ngroups] - intVal[1]) / obsdiff[pos,]
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			waldCon[pos,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
			pos <- pos + 1
		}

		facinvarfree <- findFactor(varfree, facList)
		for(i in seq_along(varnonfixvar)) {
			runnum <- indexnonfixvar[i]
			# Need to change marker variable if fixed
			oldmarker <- fixIntceptFac[[facinvarfree[i]]]
			if(length(oldmarker) > 0 && oldmarker == varfree[i]) {
				candidatemarker <- setdiff(facList[[facinvarfree[i]]], varfree[i])[1]
				temp <- freeParTable(pt1, varfree[i], "~1", "", 1:ngroups)
				temp <- constrainParTable(temp, varfree[i], "~1", "", 1:ngroups)
				temp <- fixParTable(temp, candidatemarker, "~1", "", 1:ngroups)
				newparent <- freeParTable(pt1, varfree[i], "~1", "", 1:ngroups)
				newparent <- fixParTable(newparent, candidatemarker, "~1", "", 1:ngroups)
				newparentresult <- try(newparentfit <- refit(newparent, fit1), silent = TRUE)
				if(!is(newparentresult, "try-error")) {
					tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
					if(!is(tryresult, "try-error")) {
						compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, newparentfit, method = method), silent = TRUE)
						if(!is(compresult, "try-error")) fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(newparentfit, tempfit))
					}
					waldCon[pos,] <- waldConstraint(newparentfit, newparent, waldMat, cbind(varfree[i], "~1", "", 1:ngroups))
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
			if(length(oldmarker) > 0 && oldmarker == varfree[i]) {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
			} else {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			}
			estimates[pos, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				intVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[pos, 2:ncol(estimates)] <- intVal
				totalVal <- sapply(lavaan::fitted.values(tempfit0), function(x, v) x$cov[v, v], v = pt0$lhs[runnum])
				ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
				stdIntVal <- intVal / sqrt(refTotalVal)
				stdestimates[pos,] <- stdIntVal
				stdIntVal <- stdIntVal[grouporder]
				esstd[pos,] <- stdIntVal[2:ngroups] - stdIntVal[1]

				intVal <- intVal[grouporder]
				propdiff[pos,] <- (intVal[2:ngroups] - intVal[1]) / obsdiff[pos,]
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			pos <- pos + 1
		}
		freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
		fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
		waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)
				
		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[c(indexfixvar, indexnonfixvar)]
		estimates <- cbind(estimates, stdestimates, esstd, obsdiff, propdiff)
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

		# Prop diff
		propdiff <- matrix(NA, length(varfree), ngroups - 1)
		colnames(propdiff) <- paste0("propdiff:", complab)

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
			temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
			tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
			if(!is(tryresult, "try-error")) {
				compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[i,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
			}
			listFixCon <- c(listFixCon, tryresult)
			temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			estimates[i, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[i,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				errVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[i, 2:ncol(estimates)] <- errVal
				totalVal <- sapply(lavaan::fitted.values(tempfit0), function(x, v) x$cov[v, v], v = pt0$rhs[runnum])
				ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
				stdErrVal <- errVal / sqrt(refTotalVal)
				stdestimates[i,] <- stdErrVal
				stdErrVal <- stdErrVal[grouporder]
				esstd[i,] <- stdErrVal[2:ngroups] - stdErrVal[1]
				if(any(abs(stdErrVal) > 0.9999)) warning(paste("The uniqueness of", pt0$rhs[runnum], "in some groups are over 1. The uniqueness used in arctan transformation are changed to 0.9999."))
				stdErrVal[stdErrVal > 0.9999] <- 0.9999
				zErrVal <- asin(sqrt(stdErrVal))
				esz[i,] <- zErrVal[2:ngroups] - zErrVal[1]
				
				errVal <- errVal[grouporder]
				totalVal <- totalVal[grouporder]
				errdiff <- errVal[2:ngroups] - errVal[1]
				totaldiff <- totalVal[2:ngroups] - totalVal[1]
				propdiff[i,] <- errdiff / totaldiff
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			waldCon[i,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
		}
		
		freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
		fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
		waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)
		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[index]
		estimates <- cbind(estimates, stdestimates, esstd, esz, propdiff)
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
	result <- rep(NA, 3)
	if(!any(apply(overallMat, 1, sum) != 0)) {
		try(result <- waldContrast(fit, overallMat), silent = TRUE)
	}
	return(result)
}

poolVariance <- function(var, n) {
	nm <- n - 1
	sum(var * nm) / sum(nm)
}

deltacfi <- function(parent, nested) lavaan::fitmeasures(nested)["cfi"] - lavaan::fitmeasures(parent)["cfi"]