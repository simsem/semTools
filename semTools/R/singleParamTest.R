singleParamTest <- function(model1, model2, return.fit = FALSE, method = "satorra.bentler.2001") {
	# Check nested models without any swaps
	if(lavaan::fitMeasures(model1, "df")[[1]] > lavaan::fitMeasures(model2, "df")[[1]]) {
        fit0 <- model1
        fit1 <- model2
    } else {
        fit0 <- model2
        fit1 <- model1
    }
	# fit0 = Nested model, fit1 = Parent model
	pt1 <- lavaan::partable(fit1)
	pt0 <- lavaan::partable(fit0)
	namept1 <- paramNameFromPt(pt1)
	namept0 <- paramNameFromPt(pt0)

	# Two possible constraints: fixed parameters and equality constraints
	
	free1 <- (pt1$free != 0) & !(duplicated(pt1$free))
	free0 <- (pt0$free != 0) & !(duplicated(pt0$free))
	iscon1 <- pt1$op == "=="
	iscon0 <- pt0$op == "=="
	con1 <- list(id = integer(0), lhs = character(0), op = character(0), rhs = character(0))
	con0 <- list(id = integer(0), lhs = character(0), op = character(0), rhs = character(0))
	if(any(iscon1)) con1 <- list(id = pt1$id[iscon1], lhs = pt1$lhs[iscon1], op = pt1$op[iscon1], rhs = pt1$rhs[iscon1])
	if(any(iscon0)) con0 <- list(id = pt0$id[iscon0], lhs = pt0$lhs[iscon0], op = pt0$op[iscon0], rhs = pt0$rhs[iscon0])
	
	
	if(length(free1[!iscon1]) != length(free0[!iscon0])) stop("Parameter tables in two models do not have equal lengths. This function does not work.")
	if(!all(free1[free0])) stop("Model are not nested or are not arranged in the way that this function works.")
	if(sum(iscon1) > sum(iscon0)) stop("There are equality constraints in the model with less degrees of freedom that do not exist in the model with higher degrees of freedom. Thus, two models are not nested.")

	if(!all.equal(lapply(pt1[2:4], "[", !iscon1), lapply(pt0[2:4], "[", !iscon0))) stop("This function needs parameter tables of two models to have the same orders of the same parameters.")

	# Find fixed values or constraints
	difffree <- !free0[!iscon0] & free1[!iscon1]
	textcon1 <- paste0(con1$lhs, con1$op, con1$rhs)
	textcon0 <- paste0(con0$lhs, con0$op, con0$rhs)
	indexsamecon <- match(textcon1, textcon0)
	indexdiffcon <- setdiff(seq_along(textcon0), indexsamecon)
	diffcon <- lapply(con0, "[", indexdiffcon) 
	fixval <- which(difffree)
	index <- c(fixval, diffcon$id)
	if(length(index) <= 0) stop("Two models are identical. No single parameter test can be done.")
	
	# Find nested model and release 1-by-1 
	freeCon <- matrix(NA, length(index), 2)
	colnames(freeCon) <- c("free.chi", "free.p")
	listFreeCon <- list()
	runnum <- 1
	for(i in seq_along(fixval)) {
		temp <- freeParTable(pt0, pt0$lhs[fixval[i]], pt0$op[fixval[i]], pt0$rhs[fixval[i]], pt0$group[fixval[i]]) 
		tryresult <- try(tempfit <- refit(temp, fit0), silent = TRUE)
		if(!is(tryresult, "try-error")) {
			compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit0, method = method), silent = TRUE)
			if(!is(compresult, "try-error")) freeCon[runnum,] <- unlist(modelcomp[2, c(5, 7)])
		}
		listFreeCon <- c(listFreeCon, tryresult)
		runnum <- runnum + 1
	}
	rownames(freeCon)[seq_along(fixval)] <- names(listFreeCon)[seq_along(fixval)] <- namept0[fixval]

	for(i in seq_along(diffcon$id)) {
		temp <- removeEqCon(pt0, diffcon$id[i])
		tryresult <- try(tempfit <- refit(temp, fit0), silent = TRUE)
		if(!is(tryresult, "try-error")) {
			compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit0, method = method), silent = TRUE)
			if(!is(compresult, "try-error")) freeCon[runnum,] <- unlist(modelcomp[2, c(5, 7)])
		}
		listFreeCon <- c(listFreeCon, tryresult)
		runnum <- runnum + 1
	}
	poscon <- seq_along(diffcon$id) + length(fixval)
	rownames(freeCon)[poscon] <- names(listFreeCon)[poscon] <- namept0[diffcon$id]

	
	# Find parent model and constrain 1-by-1
	fixCon <- matrix(NA, length(index), 2)
	colnames(fixCon) <- c("fix.chi", "fix.p")
	listFixCon <- list()
	runnum <- 1
	for(i in seq_along(fixval)) {
		temp <- fixParTable(pt1, pt1$lhs[fixval[i]], pt1$op[fixval[i]], pt1$rhs[fixval[i]], pt1$group[fixval[i]], pt0$ustart[fixval[i]]) 
		tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
		if(!is(tryresult, "try-error")) {
			compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
			if(!is(compresult, "try-error"))  fixCon[runnum,] <- unlist(modelcomp[2,c(5, 7)])
		}
		listFixCon <- c(listFixCon, tryresult)
		runnum <- runnum + 1
	}
	rownames(fixCon)[seq_along(fixval)] <- names(listFixCon)[seq_along(fixval)] <- namept0[fixval]
	
	for(i in seq_along(diffcon$id)) {
		temp <- patMerge(pt1, list(lhs = diffcon$lhs[i], op = diffcon$op[i], rhs = diffcon$rhs[i]))
		tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
		if(!is(tryresult, "try-error")) {
			compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
			if(!is(compresult, "try-error")) fixCon[runnum,] <- unlist(modelcomp[2, c(5, 7)])
		}
		listFixCon <- c(listFixCon, tryresult)
		runnum <- runnum + 1
	}
	poscon <- seq_along(diffcon$id) + length(fixval)
	rownames(fixCon)[poscon] <- names(listFixCon)[poscon] <- namept0[diffcon$id]
	
	result <- cbind(freeCon, fixCon)

	if(return.fit) {
		return(invisible(list(result = result, models = list(free = listFreeCon, fix = listFixCon))))
	} else {
		return(result)
	}
}

paramNameFromPt <- function(pt) {
	ngroups <- max(pt$group)
	result <- NULL
	if(ngroups == 1) {
		result <- paste0(pt$lhs, pt$op, pt$rhs)
	} else {
		grouplab <- paste0(".g", pt$group)
		grouplab[grouplab == ".g0" | grouplab == ".g1"] <- ""
		result <- paste0(pt$lhs, pt$op, pt$rhs, grouplab)
	}
	con <- pt$op == "=="
	pt$lhs[con] <- result[match(pt$lhs[con], pt$plabel)]
	pt$rhs[con] <- result[match(pt$rhs[con], pt$plabel)]
	result[con] <- paste(pt$lhs[con], pt$op[con], pt$rhs[con])
	result
}

refit <- function(pt, object, resetstart = TRUE) {
	if(resetstart && "start" %in% names(pt)) pt <- pt[-which("start" == names(pt))]
	previousCall <- lavaan::lavInspect(object, "call")
	args <- previousCall[-1]
	args$model <- pt
	funcall <- as.character(previousCall[[1]])
	tempfit <- do.call(funcall[length(funcall)], args)
}