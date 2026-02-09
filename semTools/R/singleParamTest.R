### Sunthud Pornprasertmanit
### Last updated: 9 February 2026


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
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
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
##'
##' ## Nested models to test measurement invariance
##' HW.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ x7 + x8 + x9 '
##'
##' mod.config <- cfa(model = HW.model, data = HolzingerSwineford1939,
##'                   group = "school")
##' mod.metric <- cfa(model = HW.model, data = HolzingerSwineford1939,
##'                   group = "school", group.equal = "loadings")
##' singleParamTest(mod.config, mod.metric)
##'
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



## MOVED FROM lonInvariance.R when it was removed from semTools 0.5-8 (9 Feb 2026)


# rearrangeFreeElement: Rearrange the number listed in 'free' in parameter tables
rearrangeFreeElement <- function(vec) {
  vec2 <- vec
  vec <- vec[vec != 0]
  uvec <- unique(vec)
  newvec <- 1:length(unique(vec))
  vec2[vec2 != 0] <- newvec[match(vec, uvec)]
  class(vec2) <- "integer"
  vec2
}

# rearrangept: Rearrange parameter table and plabel
rearrangept <- function(pt) {

  createplabel <- function(num) {
    result <- paste0(".p", num, ".")
    result[num == 0] <- ""
    result
  }

  oldfree <- pt$free
  newfree <- rearrangeFreeElement(oldfree)
  oldplabel <- pt$plabel
  newplabel <- createplabel(seq_along(pt$op))
  eqpos <- which(pt$op == "==")
  newplabel[eqpos] <- ""
  if (length(eqpos) > 0) {
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

# freeParTable: Free elements in parameter table
# also used in partialInvariance
freeParTable <- function(parTable, lhs, op, rhs, group, ustart = NA) {
  parTable$start <- parTable$est <- parTable$se <- NULL
  target <- cbind(lhs, op, rhs, group)
  for (i in 1:nrow(target)) {
    targetElem <- matchElement(parTable = parTable, vec = target[i,])
    ptargetElem <- parTable$plabel[targetElem]
    if ((length(targetElem) == 0) || is.na(targetElem)) {
      newline <- list(lhs = as.character(target[i, 1]),
                      op = as.character(target[i, 2]),
                      rhs = as.character(target[i, 3]),
                      group = as.integer(target[i, 4]),
                      free = as.integer(max(parTable$free) + 1),
                      ustart = as.numeric(NA))
      parTable <- patMerge(pt1 = parTable, pt2 = newline)
    } else {
      if (parTable$free[targetElem] == 0) {
        parTable$ustart[targetElem] <- ustart
        parTable$user[targetElem] <- 1
        parTable$free[targetElem] <- max(parTable$free) + 1
      }
      equalelement <- which(parTable$op == "==")
      rmelem <- intersect(union(match(ptargetElem, parTable$lhs),
                                match(ptargetElem, parTable$rhs)),
                          equalelement)
      if (length(rmelem) > 0) parTable <- removeEqCon(parTable, rmelem)
    }
  }
  parTable <- rearrangept(parTable)
  parTable
}

# fixParTable: Fix elements in parameter table
# also used in partialInvariance
fixParTable <- function(parTable, lhs, op, rhs, group, ustart = NA) {
  parTable$start <- parTable$est <- parTable$se <- NULL
  target <- cbind(lhs, op, rhs, group)
  element <- apply(target, 1, matchElement, parTable=parTable)
  for (i in 1:nrow(target)) {
    ## Why was Sunthud printing warnings? (originally used warnings(), not warning()...)
    # if (parTable$free[element[i]] == 0) warning('The parameter ', lhs, op, rhs,
    #                                             ' in group ', group,
    #                                             ' is already fixed.')

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
  parTable
}


# removeEqCon: Remove equality constraints
removeEqCon <- function(pt, element) {
  pt <- lapply(pt, "[", -element)
  pt$id <- seq_along(pt$id)
  pt
}

patMerge <- function (pt1 = NULL, pt2 = NULL, remove.duplicated = FALSE,
                      fromLast = FALSE, warn = TRUE) {
  pt1 <- as.data.frame(pt1, stringsAsFactors = FALSE)
  pt2 <- as.data.frame(pt2, stringsAsFactors = FALSE)
  stopifnot(!is.null(pt1$lhs), !is.null(pt1$op), !is.null(pt1$rhs),
            !is.null(pt2$lhs), !is.null(pt2$op), !is.null(pt2$rhs))
  if (is.null(pt1$group) && is.null(pt2$group)) {
    TMP <- rbind(pt1[, c("lhs", "op", "rhs", "group")],
                 pt2[, c("lhs", "op", "rhs", "group")])
  }
  else {
    if (is.null(pt1$group) && !is.null(pt2$group)) {
      pt1$group <- rep(1L, length(pt1$lhs))
    }
    else if (is.null(pt2$group) && !is.null(pt1$group)) {
      pt2$group <- rep(1L, length(pt2$lhs))
    }
    TMP <- rbind(pt1[, c("lhs", "op", "rhs", "group")],
                 pt2[, c("lhs", "op", "rhs", "group")])
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
  if (!is.null(pt1$est)) pt1$est <- NULL
  if (!is.null(pt2$est)) pt2$est <- NULL
  if (!is.null(pt1$se)) pt1$se <- NULL
  if (!is.null(pt2$se)) pt2$se <- NULL
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



