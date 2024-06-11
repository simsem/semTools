### Sunthud Pornprasertmanit
### Last updated: 2 June 2022


##' Single Parameter Test Divided from Nested Model Comparison
##'
##' In comparing two nested models, \eqn{\Delta\chi^2} test may indicate that
##' two models are different. However, like other omnibus tests, researchers do
##' not know which fixed parameters or constraints make these two models
##' different. This function will help researchers identify the significant
##' parameter.
##'
##' This function first identifies the differences between these two models. The
##' model with more free parameters is referred to as parent model and the model
##' with fewer free parameters is referred to as nested model. Two tests are
##' implemented here:
##'
##' \enumerate{
##'  \item `free`: The nested model is used as a template. Then,
##' one parameter indicating the differences between two models is freed. The new
##' model is compared with the nested model. This process is repeated for all
##' differences between two models.
##'  \item`fix`: The parent model is used
##' as a template. Then, one parameter indicating the differences between two
##' models is fixed or constrained to be equal to other parameters. The new
##' model is then compared with the parent model. This process is repeated for
##' all differences between two models.
##'  \item`mi`: No longer available
##' because the test of modification indices is not consistent. For example, if
##' two parameters are equality constrained, the modification index from the
##' first parameter is not equal to the second parameter.
##' }
##'
##' Note that this function does not adjust for the inflated Type I error rate
##' from multiple tests.
##'
##' @param model1 Model 1.
##' @param model2 Model 2. Note that two models must be nested models. Further,
##' the order of parameters in their parameter tables are the same. That is,
##' nested models with different scale identifications may not be able to test
##' by this function.
##' @param return.fit Return the submodels fitted by this function
##' @param method The method used to calculate likelihood ratio test. See
##' [lavaan::lavTestLRT()] for available options
##' @return If `return.fit = FALSE`, the result tables are provided.
##' \eqn{\chi^2} and *p* value are provided for all methods. Note that the
##' \eqn{\chi^2} is all based on 1 *df*. Expected parameter changes
##' and their standardized forms are also provided.
##'
##' If `return.fit = TRUE`, a list with two elements are provided. The
##' first element is the tabular result. The second element is the submodels
##' used in the `free` and `fix` methods.
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##' @examples
##'
##' library(lavaan)
##'
##' # Nested model comparison by hand
##' HS.model1 <- ' visual =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6'
##' HS.model2 <- ' visual =~ a*x1 + a*x2 + a*x3
##'               textual =~ b*x4 + b*x5 + b*x6'
##'
##' m1 <- cfa(HS.model1, data = HolzingerSwineford1939, std.lv = TRUE,
##'           estimator = "MLR")
##' m2 <- cfa(HS.model2, data = HolzingerSwineford1939, std.lv = TRUE,
##'           estimator = "MLR")
##' anova(m1, m2)
##' singleParamTest(m1, m2)
##'
##' ## Nested model comparison from the measurementInvariance function
##' HW.model <- ' visual =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed =~ x7 + x8 + x9 '
##'
##' models <- measurementInvariance(model = HW.model, data = HolzingerSwineford1939,
##'                                 group = "school")
##' singleParamTest(models[[1]], models[[2]])
##'
##' ## Note that the comparison between metric (Model 2) and scalar invariance
##' ## (Model 3) cannot be done by this function because the metric invariance
##' ## model fixes factor means as 0 in Group 2 but the strong invariance model
##' ## frees the factor means in Group 2. Users may use this function to compare
##' ## scalar invariance (Model 3) to a homogeneous-means model.
##'
##' @export
singleParamTest <- function(model1, model2, return.fit = FALSE,
                            method = "satorra.bentler.2001") {
	# Check nested models without any swaps
	if(lavaan::fitMeasures(model1, "df")[[1]] > lavaan::fitMeasures(model2, "df")[[1]]) {
        fit0 <- model1
        fit1 <- model2
    } else {
        fit0 <- model2
        fit1 <- model1
    }
	# fit0 = Nested model, fit1 = Parent model
	pt1 <- parTable(fit1)
	pt0 <- parTable(fit0)
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



## ----------------
## Hidden Functions
## ----------------

paramNameFromPt <- function(pt) {
	ngroups <- max(pt$group)
	result <- NULL
	if (ngroups == 1) {
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

##' @importFrom lavaan lavInspect
refit <- function(pt, object, resetstart = TRUE) {
	if (resetstart && "start" %in% names(pt)) pt <- pt[-which("start" == names(pt))]
	previousCall <- lavInspect(object, "call")
	## Why this?
	args <- previousCall[-1]
	args$model <- pt
	funcall <- as.character(previousCall[[1]])
	do.call(funcall[length(funcall)], args)
	## instead of this?
	# previousCall$model <- pt
	# eval(previousCall)
}


