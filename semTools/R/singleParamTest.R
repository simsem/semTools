singleParamTest <- function(model1, model2, mi = TRUE, return.fit = FALSE) {
	# Check nested models without any swaps
	if(model1@Fit@test[[1]]$df > model2@Fit@test[[1]]$df) {
        fit0 <- model1
        fit1 <- model2
    } else {
        fit0 <- model2
        fit1 <- model1
    }
	# fit0 = Nested model, fit1 = Parent model
	pt1 <- partable(fit1)
	pt0 <- partable(fit0)
	namept1 <- paramNameFromPt(pt1)
	namept0 <- paramNameFromPt(pt0)

	free1 <- (pt1$free != 0) & !(duplicated(pt1$free))
	free0 <- (pt0$free != 0) & !(duplicated(pt0$free))
	if(length(free1) != length(free0)) stop("Parameter tables in two models do not have equal lengths. This function does not work.")
	if(!all(free1[free0])) stop("Model are not nested or are not arranged in the way that this function works.")
	if(!all.equal(pt1[2:4], pt0[2:4])) stop("This function needs parameter tables of two models to have the same orders of the same parameters.")

	# Find fixed values or constraints
	differences <- !free0 & free1
	index <- which(differences)
	if(length(index) <= 0) stop("Two models are identical. No single parameter test can be done.")
	
	# Find nested model and release 1-by-1 
	freeCon <- matrix(NA, length(index), 2)
	colnames(freeCon) <- c("free.chi", "free.p")
	listFreeCon <- list()
	for(i in seq_along(index)) {
		runnum <- index[i]
		temp <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], pt0$group[runnum]) 
		previousCall <- fit0@call
		args <- as.list(previousCall[-1])
		args$model <- temp
		funcall <- as.character(previousCall[[1]])
		tryresult <- try(tempfit <- do.call(funcall[length(funcall)], args), silent = TRUE)
		if(!is(tryresult, "try-error")) {
			compresult <- try(modelcomp <- lavTestLRT(tempfit, fit0), silent = TRUE)
			if(!is(compresult, "try-error")) freeCon[i,] <- unlist(modelcomp[2, c(5, 7)])
		}
		listFreeCon <- c(listFreeCon, tryresult)
	}
	rownames(freeCon) <- names(listFreeCon) <- namept0[index]

	# Find parent model and constrain 1-by-1
	fixCon <- matrix(NA, length(index), 2)
	colnames(fixCon) <- c("fix.chi", "fix.p")
	listFixCon <- list()
	for(i in seq_along(index)) {
		runnum <- index[i]
		# Check whether is it a fixed value or a constraints
		fixed <- pt0$free[runnum] == 0
		if(fixed) {
			temp <- fixParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], pt1$group[runnum], pt0$ustart[runnum]) 
		} else {
			constrainVal <- which(pt0$free == pt0$free[runnum])
			temp <- constrainParTable(pt1, pt1$lhs[constrainVal], pt1$op[constrainVal], pt1$rhs[constrainVal], pt1$group[constrainVal])
		}
		previousCall <- fit1@call
		args <- as.list(previousCall[-1])
		args$model <- temp
		funcall <- as.character(previousCall[[1]])
		tryresult <- try(tempfit <- do.call(funcall[length(funcall)], args), silent = TRUE)
		if(!is(tryresult, "try-error")) {
			compresult <- try(modelcomp <- lavTestLRT(tempfit, fit0), silent = TRUE)
			if(!is(compresult, "try-error"))  fixCon[i,] <- unlist(modelcomp[2,c(5, 7)])
		}
		listFixCon <- c(listFixCon, tryresult)
	}
	rownames(fixCon) <- names(listFixCon) <- namept1[index]
	result <- cbind(freeCon, fixCon)
	
	if(mi) {
		# Show the modification indices of the nested model
		ngroups <- max(pt0$group)
		if(ngroups == 1) {
			target <- cbind(pt0$lhs[index], pt0$op[index], pt0$rhs[index])
		} else {
			target <- cbind(pt0$lhs[index], pt0$op[index], pt0$rhs[index], pt0$group[index])
		}
		miResult <- modindices(fit0)
		element <- apply(target, 1, matchElement, parTable=miResult)
		mi.chi <- miResult$mi[element]
		mi.p <- pchisq(mi.chi, df=1, lower.tail=FALSE)
		miCon <- cbind(mi.chi, mi.p, epc = miResult$epc[element], sepc.lv = miResult$sepc.lv[element], sepc.all = miResult$sepc.all[element], sepc.nox = miResult$sepc.nox[element])
		result <- cbind(result, miCon)
	}
	# Write out the awesome table

	if(return.fit) {
		return(invisible(list(result = result, models = list(free = listFreeCon, fix = listFixCon))))
	} else {
		return(result)
	}
}

paramNameFromPt <- function(pt) {
	ngroups <- max(pt$group)
	if(ngroups == 1) {
		return(paste0(pt$lhs, pt$op, pt$rhs))
	} else {
		grouplab <- paste0(".g", pt$group)
		grouplab[grouplab == ".g0" | grouplab == ".g1"] <- ""
		return(paste0(pt$lhs, pt$op, pt$rhs, grouplab))
	}
}
