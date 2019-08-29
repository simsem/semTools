### Terrence D. Jorgensen & Sunthud Pornprasertmanit
### Last updated: 29 August 2019
### source code for compareFit() function and FitDiff class


## -----------------
## Class and Methods
## -----------------

##' Class For Representing A Template of Model Fit Comparisons
##'
##' This class contains model fit measures and model fit comparisons among
##' multiple models
##'
##'
##' @name FitDiff-class
##' @aliases FitDiff-class show,FitDiff-method summary,FitDiff-method
##' @docType class
##'
##' @slot name \code{character}. The name of each model
##' @slot model.class \code{character}. One class to which each model belongs
##' @slot nested \code{data.frame}. Model fit comparisons between adjacently
##'   nested models that are ordered by their degrees of freedom (\emph{df})
##' @slot fit \code{data.frame}. Fit measures of all models specified in the
##'   \code{name} slot, ordered by their \emph{df}
##'
##' @section Objects from the Class: Objects can be created via the
##'  \code{\link{compareFit}} function.
##'
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##'   Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' @seealso \code{\link{compareFit}}; \code{\link{clipboard}}
##'
##' @examples
##'
##' HS.model <- ' visual =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed =~ x7 + x8 + x9 '
##'
##' out <- measurementInvariance(model = HS.model, data = HolzingerSwineford1939,
##'                              group = "school", quiet = TRUE)
##' modelDiff <- compareFit(out)
##' summary(modelDiff)
##' summary(modelDiff, fit.measures = "all")
##' summary(modelDiff, fit.measures = c("aic", "bic"))
##'
##' \dontrun{
##' ## Save results to a file
##' saveFile(modelDiff, file = "modelDiff.txt")
##'
##' ## Copy to a clipboard
##' clipboard(modelDiff)
##' }
##'
setClass("FitDiff", slots = c(name = "character",
                              model.class = "character",
                              nested = "data.frame",
                              fit = "data.frame"))



##' @rdname FitDiff-class
##' @aliases show,FitDiff-method
##' @importFrom methods getMethod
##' @export
setMethod("show", signature(object = "FitDiff"), function(object) {
  getMethod("summary", signature = "FitDiff")(object)
  invisible(object)
})



##' @rdname FitDiff-class
##' @aliases summary,FitDiff-method
##'
##' @param object object of class \code{FitDiff}
##' @param fit.measures \code{character} vector naming fit indices the user can
##'   request from \code{\link[lavaan]{fitMeasures}}. If \code{"default"}, the
##'   fit measures will be \code{c("chisq", "df", "pvalue", "cfi", "tli",
##'   "rmsea", "srmr", "aic", "bic")}. If \code{"all"}, all available fit
##'   measures will be returned.
##' @param nd number of digits printed
##'
##' @export
setMethod("summary", signature(object = "FitDiff"),
          function(object, fit.measures = "default", nd = 3) {

  out <- list()

  if (nrow(object@nested) > 0L) {
		cat("################### Nested Model Comparison #########################\n")
    out$test.statistics <- object@nested
    if (object@model.class == "lavaan") {
      print(out$test.statistics, nd = nd)
    } else {
      class(out$test.statistics) <- c("lavaan.data.frame","data.frame")
      stats::printCoefmat(out$test.statistics, P.values = TRUE, has.Pvalue = TRUE)
    }
		cat("\n")
  }


  noFit <- ncol(object@fit) == 1L && names(object@fit)[1] == "df"
  if (!noFit) {
    cat("####################### Model Fit Indices ###########################\n")
    ## this is the object to return (numeric, no printed daggers)
    out$fit.indices <- getFitSummary(object, fit.measures, return.diff = FALSE)
    class(out$fit.indices) <- c("lavaan.data.frame","data.frame")

    ## print with daggers marking each fit index's preferred model
    ## (turns "numeric" vectors into "character")
    badness <- grepl(pattern = c("chisq|rmsea|ic|rmr|ecvi|fmin"),
                     x = colnames(out$fit.indices))
    goodness <- grepl(pattern = c("cfi|tli|rfi|nfi|ifi|rni|cn|gfi|mfi"),
                      x = colnames(out$fit.indices))
    minvalue <- badness & !goodness
    minvalue[!badness & !goodness] <- NA
    fit.integer <- grepl(pattern = c("df|npar|ntotal"),
                         x = colnames(out$fit.indices))
    suppressWarnings(fitTab <- as.data.frame(mapply(tagDagger, nd = nd,
                                                    vec = out$fit.indices,
                                                    minvalue = minvalue,
                                                    print_integer = fit.integer),
                                             stringsAsFactors = FALSE))
    rownames(fitTab) <- object@name
    colnames(fitTab) <- colnames(out$fit.indices)
    print(fitTab)
    cat("\n")


    if (nrow(object@nested) > 0L) {
      cat("################## Differences in Fit Indices #######################\n")
      out$fit.diff <- getFitSummary(object, fit.measures, return.diff = TRUE)
      class(out$fit.diff) <- c("lavaan.data.frame","data.frame")
      print(out$fit.diff, nd = nd)
      cat("\n")
    }
  }


	invisible(out)
})



## "method" for saveFile() function (see "clipboard.R")
saveFileFitDiff <- function(object, file, what = "summary",
                            tableFormat = FALSE, fit.measures = "default",
                            writeArgs = list()) {
	if (tableFormat) {
	  writeArgs$file <- file
	  writeArgs$append <- TRUE
	  if (is.null(writeArgs$sep)) writeArgs$sep <- "\t"
		if (is.null(writeArgs$quote)) writeArgs$quote <- FALSE
		if (is.null(writeArgs$row.names)) writeArgs$row.names <- FALSE

		if (nrow(object@nested) > 0L) {
			cat("Nested Model Comparison\n\n", file = file, append = TRUE)
			out <- object@nested
			#out <- data.frame(model.diff = rownames(out), out)
  		writeArgs$x <- out
			do.call("write.table", writeArgs)
			cat("\n\n", file = file, append = TRUE)
		}
		out2 <- getFitSummary(object, fit.measures)
		out2 <- data.frame(model = object@name, out2)
		cat("Fit Indices Summaries\n\n", file = file, append = TRUE)
		writeArgs$x <- out2
		do.call("write.table", writeArgs)
	} else {
		write(paste(utils::capture.output(lavaan::summary(object)),
		            collapse = "\n"), file = file)
	}
}

## --------------------
## Constructor Function
## --------------------

##' Build an object summarizing fit indices across multiple models
##'
##' This function will create the template to compare fit indices across
##' multiple fitted lavaan objects. The results can be exported to a clipboard
##' or a file later.
##'
##' @importFrom lavaan lavTestLRT
##' @importMethodsFrom lavaan fitMeasures
##'
##' @param ...  fitted \code{lavaan} models or list(s) of \code{lavaan} objects.
##'   \code{\linkS4class{lavaan.mi}} objects are also accepted, but all models
##'   must belong to the same class.
##' @param nested \code{logical} indicating whether the models in \code{...} are
##'   nested. See \code{\link{net}} for an empirical test of nesting.
##' @param argsLRT \code{list} of arguments to pass to
##'   \code{\link[lavaan]{lavTestLRT}}, as well as to
##'   \code{\link{lavTestLRT.mi}} and \code{\link{fitMeasures}} when
##'   comparing \code{\linkS4class{lavaan.mi}} models.
##' @param indices \code{logical} indicating whether to return fit indices from
##'   the \code{\link[lavaan]{fitMeasures}} function.
##' @param baseline.model optional fitted \code{\linkS4class{lavaan}} model
##'   passed to \code{\link[lavaan]{fitMeasures}} to calculate incremental fit
##'   indices.
##'
##' @return A \code{\linkS4class{FitDiff}} object that saves model fit
##'   comparisons across multiple models. If the models are not nested, only
##'   fit indices for each model are returned. If the models are nested, the
##'   differences in fit indices are additionally returned, as well as test
##'   statistics comparing each sequential pair of models (ordered by their
##'   degrees of freedom).
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' @seealso \code{\linkS4class{FitDiff}}, \code{\link{clipboard}}
##'
##' @examples
##'
##' HS.model <- ' visual =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed =~ x7 + x8 + x9 '
##'
##' fit1 <- cfa(HS.model, data = HolzingerSwineford1939)
##'
##' ## non-nested model
##' m2 <- ' f1 =~ x1 + x2 + x3 + x4
##'         f2 =~ x5 + x6 + x7 + x8 + x9 '
##' fit2 <- cfa(m2, data = HolzingerSwineford1939)
##' compareFit(fit1, fit2, nested = FALSE)
##'
##'
##' ## nested model comparisons:
##' out <- measurementInvariance(model = HS.model, data = HolzingerSwineford1939,
##'                              group = "school", quiet = TRUE)
##' compareFit(out)
##'
##' \dontrun{
##' ## also applies to lavaan.mi objects (fit model to multiple imputations)
##' set.seed(12345)
##' HSMiss <- HolzingerSwineford1939[ , paste("x", 1:9, sep = "")]
##' HSMiss$x5 <- ifelse(HSMiss$x1 <= quantile(HSMiss$x1, .3), NA, HSMiss$x5)
##' HSMiss$x9 <- ifelse(is.na(HSMiss$x5), NA, HSMiss$x9)
##' HSMiss$school <- HolzingerSwineford1939$school
##' HS.amelia <- amelia(HSMiss, m = 20, noms = "school")
##' imps <- HS.amelia$imputations
##'
##' ## request robust test statistics
##' mgfit2 <- cfa.mi(HS.model, data = imps, group = "school", estimator = "mlm")
##' mgfit1 <- cfa.mi(HS.model, data = imps, group = "school", estimator = "mlm",
##'                  group.equal = "loadings")
##' mgfit0 <- cfa.mi(HS.model, data = imps, group = "school", estimator = "mlm",
##'                  group.equal = c("loadings","intercepts"))
##'
##' ## request the strictly-positive robust test statistics
##' compareFit(scalar = mgfit0, metric = mgfit1, config = mgfit2,
##'            argsLRT = list(asymptotic = TRUE,
##'                           method = "satorra.bentler.2010"))
##' }
##'
##' @export
compareFit <- function(..., nested = TRUE, argsLRT = list(),
                       indices = TRUE, baseline.model = NULL) {
  ## make sure there is something to do
  if (!(nested || indices)) {
    message('User requested neither indices nor tests of nested models.')
    return(NULL)
  }

	## separate models from lists of models
  dots <- list(...)
	idx.list <- sapply(dots, is.list)
	modLists <- dots[ idx.list]
	mods     <- dots[!idx.list]
	## capture names of any arguments passed via dots
	allnames  <- sapply(substitute(list(...))[-1], deparse)
	listnames <- allnames[ idx.list]
	modnames  <- allnames[!idx.list]

	## make sure models are named
	if (length(mods) && is.null(names(mods))) {
	  names(mods) <- modnames
	} else for (nn in seq_along(mods)) {
	  if (names(mods)[nn] == "") names(mods)[nn] <- modnames[nn]
	}
	## make sure lists are named
	if (length(modLists) && is.null(names(modLists))) {
	  names(modLists) <- listnames
	} else for (nn in seq_along(modLists)) {
	  if (names(modLists)[nn] == "") names(modLists)[nn] <- listnames[nn]
	}
	## within each list, make sure models are named
	for (i in seq_along(modLists)) {
	  if (length(modLists[[i]]) && is.null(names(modLists[[i]]))) {
	    names(modLists[[i]]) <- seq_along(modLists[[i]])
	  } else for (nn in seq_along(modLists[[i]])) {
	    if (names(modLists[[i]])[nn] == "") names(modLists[[i]])[nn] <- nn
	  }
	}

	## collapse into a single list of models
	if (length(modLists)) mods <- c(mods, unlist(modLists))

	## check for lavaan models
	not.lavaan <- !sapply(mods, inherits, what = c("lavaan","lavaanList"))
	if (any(not.lavaan)) stop("The following are not fitted lavaan models:\n",
	                          paste0(names(which(not.lavaan)), collapse = ", "))
	modClass <- unique(sapply(mods, class))
	if (length(modClass) > 1L) stop('All models must be of the same class (e.g.,',
	                                ' cannot compare lavaan objects to lavaan.mi)')


	## grab lavaan.mi options, if relevant
	if (is.null(argsLRT$pool.robust)) {
	  pool.robust <- formals(lavTestLRT.mi)$pool.robust # default value
	} else {
	  pool.robust <- argsLRT$pool.robust # user-specified value
	}
	if (is.null(argsLRT$test)) {
	  test <- eval(formals(lavTestLRT.mi)$test) # default value
	} else {
	  test <- argsLRT$test # user-specified value
	}

	## FIT INDICES
	if (indices) {
	  fitList <- lapply(mods, fitMeasures, baseline.model = baseline.model,
	                    pool.robust = pool.robust, test = test)
	  if (length(unique(sapply(fitList, length))) > 1L) {
	    warning('fitMeasures() returned vectors of different lengths for different',
	            ' models, probably because certain options are not the same. Check',
	            ' lavInspect(fit, "options")[c("estimator","test","meanstructure")]',
	            ' for each model, or run fitMeasures() on each model to investigate.')
	    indexList <- lapply(fitList, names)
	    useNames <- names(which(table(unlist(indexList)) == length(fitList)))
	    fitList <- lapply(fitList, "[", i = useNames)
	  }
	  fit <- as.data.frame(do.call(rbind, fitList))
	} else {
	  fitList <- lapply(mods, fitMeasures, fit.measures = "df",
	                    pool.robust = pool.robust, test = test)
	  ## check for scaled tests
	  nDF <- sapply(fitList, length)
	  if (any(nDF != nDF[1])) stop('Some (but not all) models have robust tests,',
	                               ' so they cannot be compared as nested models.')
	  fit <- data.frame(df = sapply(fitList, "[",
	                                i = if (any(nDF > 1L)) 2L else 1L))
	}


	## order models by increasing df (least-to-most constrained)
	ord <- order(fit$df) #FIXME: what if test == "mean.var.adjusted"?
	fit <- fit[ord, , drop = FALSE]
	mods <- mods[ord]

	## TEST STATISTICS
	if (nested) {

	  if (class(mods[[1]]) == "lavaan") {
	    argsLRT$model.names <- names(mods)
	    argsLRT$object <- mods[[1]]
	    nestedout <- do.call(lavTestLRT, c(mods[-1], argsLRT))
	  } else if (inherits(mods[[1]], "lavaan.mi")) { #FIXME: generalize to lavaan.pool

	    modsA <- mods[-1]
	    modsB <- mods[-length(mods)]
	    fitDiff <- list()

	    for (i in seq_along(modsA)) {
	      fitA <- modsA[[i]]
	      fitB <- modsB[[i]]
	      if (is.null(argsLRT$asymptotic))
	        argsLRT$asymptotic <- any(lavListInspect(fitA, "options")$test %in%
	                                    c("satorra.bentler","yuan.bentler",
	                                      "yuan.bentler.mplus","scaled.shifted",
	                                      "mean.var.adjusted","satterthwaite"))
	      tempDiff <- do.call(lavTestLRT.mi, c(list(fitA, h1 = fitB), argsLRT))

	      if (names(tempDiff)[1] == "F") {
	        statNames <- c("F", "df1", "df2", "pvalue")
	      } else statNames <- c("chisq", "df", "pvalue")
	      ## check for scaled
	      if (any(grepl(pattern = "scaled", x = names(tempDiff)))) {
	        statNames <- paste0(statNames, ".scaled")
	      }

	      diffName <- paste(names(modsA)[i], "-", names(modsB)[i])
	      fitDiff[[diffName]] <- tempDiff[statNames]
	    }

	    nestedout <- as.data.frame(do.call(rbind, fitDiff))
	  }

	  ## not nested
	} else nestedout <- data.frame()

	new("FitDiff",
	    name = names(mods), model.class = modClass,
	    nested = nestedout, fit = fit)
}



## ----------------
## Hidden Functions
## ----------------

noLeadingZero <- function(vec, fmt, nd = 3L) {
  out <- sprintf(fmt, vec)
  upper.limit <- paste0(".", paste(rep(9, nd - 1L), collapse = ""), "5")
  used <- vec < as.numeric(upper.limit) & vec >= 0
  used[is.na(used)] <- FALSE
  out[used] <- substring(out[used], 2)
  out
}

tagDagger <- function(vec, minvalue = NA, print_integer = FALSE, nd = 3L) {
  if (print_integer) {
    vec <- noLeadingZero(vec, fmt = "%.0f", nd = nd)
  } else if (is.na(minvalue)) {
    vec <- noLeadingZero(vec, fmt = paste0("%.", nd, "f"), nd = nd)
  } else {
    target <- if (minvalue) min(vec, na.rm = TRUE) else max(vec, na.rm = TRUE)
    tag <- rep(" ", length(vec))
    tag[vec == target] <- "\u2020"
    vec <- noLeadingZero(vec, fmt = paste0("%.", nd, "f"), nd = nd)
    vec <- paste0(vec, tag)
  }

  vec
}

getFitSummary <- function(object, fit.measures = "default", return.diff = FALSE) {
  if (is.null(fit.measures)) fit.measures <- colnames(object@fit)
  if ("all" %in% fit.measures) fit.measures <- colnames(object@fit)

  if (length(fit.measures) == 1 && fit.measures == "default") {
    ## robust or scaled test statistics?
    if (is.null(object@fit$cfi.scaled)) {
      fit.measures <- c("chisq","df","pvalue","rmsea","cfi","tli","srmr")
    } else if (all(!is.na(object@fit$cfi.robust)) && !is.null(object@fit$cfi.robust)) {
      fit.measures <- c("chisq.scaled","df.scaled","pvalue.scaled",
                        "rmsea.robust","cfi.robust","tli.robust","srmr")
    } else {
      fit.measures <- c("chisq.scaled","df.scaled","pvalue.scaled",
                        "rmsea.scaled","cfi.scaled","tli.scaled","srmr")
    }

    if ("aic" %in% colnames(object@fit)) {
      fit.measures <- c(fit.measures, "aic", "bic")
    }
  }

  ## chi-squared difference test already reported, so remove (diff in p-value)
  if (return.diff) {
    fit.measures <- fit.measures[!grepl(pattern = "chisq|pvalue|ntotal",
                                        x = fit.measures)]
  }

  ## return numeric values
  fitTab <- object@fit[ , colnames(object@fit) %in% fit.measures]
  if (!return.diff) return(fitTab)

  ## or return differences in fit indices
  diffTab <- as.data.frame(do.call(cbind, lapply(fitTab, diff)))
  rownames(diffTab) <- paste(object@name[-1], "-", object@name[-length(object@name)])
  diffTab
}


