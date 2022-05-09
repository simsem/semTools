### Terrence D. Jorgensen & Sunthud Pornprasertmanit
### Last updated: 9 May 2022
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
##' @slot fit.diff \code{data.frame}. Sequential differences in fit measures in
##'   the \code{fit} slot
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
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ x7 + x8 + x9 '
##' fit.config <- cfa(HS.model, data = HolzingerSwineford1939, group = "school")
##' ## invariance constraints
##' fit.metric <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
##'                   group.equal = "loadings")
##' fit.scalar <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
##'                   group.equal = c("loadings","intercepts"))
##' fit.strict <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
##'                   group.equal = c("loadings","intercepts","residuals"))
##' measEqOut <- compareFit(fit.config, fit.metric, fit.scalar, fit.strict)
##' summary(measEqOut)
##' summary(measEqOut, fit.measures = "all")
##' summary(measEqOut, fit.measures = c("aic", "bic"))
##'
##' \dontrun{
##' ## Save results to a file
##' saveFile(measEqOut, file = "measEq.txt")
##'
##' ## Copy to a clipboard
##' clipboard(measEqOut)
##' }
##'
setClass("FitDiff", slots = c(name = "character", # list of model names
                              model.class = "character", # lavaan or lavaan.mi
                              nested = "data.frame",     # anova() table
                              fit = "data.frame",        # fitMeasures() output
                              fit.diff = "data.frame"))  # index differences


##' @rdname FitDiff-class
##' @aliases show,FitDiff-method
##' @importFrom methods getMethod
##' @export
setMethod("show", signature(object = "FitDiff"), function(object) {
  cat("The following", object@model.class, "models were compared:\n    ")
  cat(object@name, sep = "\n    ")
  cat("To view results, assign the compareFit() output to an object and ",
      "use the summary() method; see the class?FitDiff help page.\n")
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
##' @param tag single \code{character} used to flag the model preferred by each
##'   fit index. To omit tags, set to \code{NULL} or \code{NA}.
##'
##' @export
setMethod("summary", signature(object = "FitDiff"),
          function(object, fit.measures = "default", nd = 3, tag = "\u2020") {

  if (nrow(object@nested) > 0L) {
		cat("################### Nested Model Comparison #########################\n")
    test.statistics <- object@nested
    if (object@model.class == "lavaan") {
      print(test.statistics, nd = nd)
    } else {
      class(test.statistics) <- c("lavaan.data.frame","data.frame")
      stats::printCoefmat(test.statistics, P.values = TRUE, has.Pvalue = TRUE)
    }
		cat("\n")
  }


  noFit <- ncol(object@fit) == 1L && names(object@fit)[1] == "df"
  if (!noFit) {
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


    cat("####################### Model Fit Indices ###########################\n")
    ## this is the object to return (numeric, no printed daggers)
    fit.indices <- object@fit[ , fit.measures , drop = FALSE]

    ## print with daggers marking each fit index's preferred model
    ## (turns "numeric" vectors into "character")
    badness <- grepl(pattern = c("chisq|rmsea|ic|rmr|ecvi|fmin|hqc"),
                     x = colnames(fit.indices))
    goodness <- grepl(pattern = c("cfi|tli|rfi|nfi|ifi|rni|cn|gfi|mfi|Hat"),
                      x = colnames(fit.indices))
    minvalue <- badness & !goodness
    minvalue[!badness & !goodness] <- NA
    fit.integer <- grepl(pattern = c("df|npar|ntotal"),
                         x = colnames(fit.indices))
    suppressWarnings(fitTab <- as.data.frame(mapply(tagCharacter, nd = nd,
                                                    char = tag,
                                                    vec = fit.indices,
                                                    minvalue = minvalue,
                                                    print_integer = fit.integer),
                                             stringsAsFactors = FALSE))
    rownames(fitTab) <- object@name
    colnames(fitTab) <- colnames(fit.indices)
    class(fitTab) <- c("lavaan.data.frame","data.frame")
    print(fitTab, nd = nd)
    cat("\n")


    if (nrow(object@nested) > 0L) {
      fit.diff.measures <- fit.measures[!grepl(pattern = "chisq|pvalue|ntotal",
                                               x = fit.measures)]
      cat("################## Differences in Fit Indices #######################\n")
      fit.diff <- object@fit.diff[ , fit.diff.measures, drop = FALSE]
      class(fit.diff) <- c("lavaan.data.frame","data.frame")
      print(fit.diff, nd = nd)
      cat("\n")
    }
  }


	invisible(object)
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
		write(paste(utils::capture.output(getMethod("summary",
		                                            signature = "FitDiff")(object)),
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
##'   the \code{\link[lavaan]{fitMeasures}} function. Selecting particular
##'   indices is controlled in the \code{summary} method; see
##'   \code{\linkS4class{FitDiff}}.
##' @param moreIndices \code{logical} indicating whether to return fit indices
##'   from the \code{\link{moreFitIndices}} function. Selecting particular
##'   indices is controlled in the \code{summary} method; see
##'   \code{\linkS4class{FitDiff}}.
##' @param baseline.model optional fitted \code{\linkS4class{lavaan}} model
##'   passed to \code{\link[lavaan]{fitMeasures}} to calculate incremental fit
##'   indices.
##' @param nPrior passed to \code{\link{moreFitIndices}}, if relevant
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
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ x7 + x8 + x9 '
##'
##' ## non-nested models
##' fit1 <- cfa(HS.model, data = HolzingerSwineford1939)
##'
##' m2 <- ' f1 =~ x1 + x2 + x3 + x4
##'         f2 =~ x5 + x6 + x7 + x8 + x9 '
##' fit2 <- cfa(m2, data = HolzingerSwineford1939)
##'
##' (out1 <- compareFit(fit1, fit2, nested = FALSE))
##' summary(out1)
##'
##'
##' ## nested model comparisons: measurement equivalence/invariance
##' fit.config <- cfa(HS.model, data = HolzingerSwineford1939, group = "school")
##' fit.metric <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
##'                   group.equal = "loadings")
##' fit.scalar <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
##'                   group.equal = c("loadings","intercepts"))
##' fit.strict <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
##'                   group.equal = c("loadings","intercepts","residuals"))
##'
##' measEqOut <- compareFit(fit.config, fit.metric, fit.scalar, fit.strict,
##'                         moreIndices = TRUE) # include moreFitIndices()
##' summary(measEqOut)
##' summary(measEqOut, fit.measures = "all")
##' summary(measEqOut, fit.measures = c("aic", "bic", "sic", "ibic"))
##'
##'
##' \dontrun{
##' ## also applies to lavaan.mi objects (fit model to multiple imputations)
##' set.seed(12345)
##' HSMiss <- HolzingerSwineford1939[ , paste("x", 1:9, sep = "")]
##' HSMiss$x5 <- ifelse(HSMiss$x1 <= quantile(HSMiss$x1, .3), NA, HSMiss$x5)
##' HSMiss$x9 <- ifelse(is.na(HSMiss$x5), NA, HSMiss$x9)
##' HSMiss$school <- HolzingerSwineford1939$school
##'
##' library(Amelia)
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
##' out2 <- compareFit(scalar = mgfit0, metric = mgfit1, config = mgfit2,
##'                    argsLRT = list(asymptotic = TRUE,
##'                                   method = "satorra.bentler.2010"))
##' ## note that moreFitIndices() does not work for lavaan.mi objects, but the
##' ## fitMeasures() method for lavaan.mi objects already returns gammaHat(s)
##' summary(out2, fit.measures = c("ariv","fmi","df","crmr","srmr",
##'                                "cfi.robust","tli.robust",
##'                                "adjGammaHat.scaled","rmsea.ci.lower.robust",
##'                                "rmsea.robust","rmsea.ci.upper.robust"))
##' }
##'
##' @export
compareFit <- function(..., nested = TRUE, argsLRT = list(), indices = TRUE,
                       moreIndices = FALSE, baseline.model = NULL, nPrior = 1) {
  ## make sure there is something to do
  if (!(nested || indices || moreIndices)) {
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
	not.lavaan <- !sapply(mods, inherits, what = c("lavaan","lavaan.mi"))
	if (any(not.lavaan)) stop("The following are not fitted lavaan(.mi) models:\n",
	                          paste0(names(which(not.lavaan)), collapse = ", "))
	modClass <- unique(sapply(mods, class))
	if (length(modClass) > 1L) stop('All models must be of the same class (e.g.,',
	                                ' cannot compare lavaan objects to lavaan.mi)')
	if (inherits(mods[[1]], "lavaan")) {
	  nonConv <- !sapply(mods, lavInspect, what = "converged")
	} else if (inherits(mods[[1]], "lavaan.mi")) {
	  nonConv <- !sapply(mods, function(fit) {
	    any(sapply(fit@convergence, "[", i = "converged"))
	  })
	}
	if (all(nonConv)) {
	  stop('No models converged')
	} else if (any(nonConv)) {
	  message('The following models did not converge, so they are ignored:\n',
	          paste(names(nonConv)[nonConv], collapse = ",\t"))
	  mods <- mods[which(!nonConv)]
	}


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
	if (indices || moreIndices) {
	  fitList <- lapply(mods, fitMeasures, baseline.model = baseline.model,
	                    pool.robust = pool.robust, test = test)
	  if (moreIndices && modClass == "lavaan") {
	    moreFitList <- lapply(mods, moreFitIndices, nPrior = nPrior)
	    fitList <- mapply(c, fitList, moreFitList, SIMPLIFY = FALSE)
	  }

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

	  if (inherits(mods[[1]], "lavaan")) {
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

	## DIFFERENCES IN FIT INDICES
	if (indices && length(names(mods)) > 1L) {
	  fitSubset <-  colnames(fit)[!grepl(pattern = "chisq|pvalue|ntotal",
	                                     x =  colnames(fit))]
	  fitTab <- fit[ , fitSubset, drop = FALSE]
	  diffTab <- as.data.frame(do.call(cbind, lapply(fitTab, diff)))
	  rownames(diffTab) <- paste(names(mods)[-1], "-", names(mods)[-length(names(mods))])

	} else diffTab <- data.frame(df = diff(fit))

	new("FitDiff", name = names(mods), model.class = modClass,
	    nested = nestedout, fit = fit, fit.diff = diffTab)
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

tagCharacter <- function(vec, char = "\u2020", minvalue = NA,
                         print_integer = FALSE, nd = 3L) {
  char <- if (is.null(char)) as.character(NA) else as.character(char)[1]
  if (nchar(char) != 1L) {
    message('Only a single character can be used to tag= preferred models, so ',
            'no tags were added.  To omit tags, specify tag=NULL or tag=NA.')
    char <- as.character(NA)
  }
  if (print_integer) {
    vec <- noLeadingZero(vec, fmt = "%.0f", nd = nd)
  } else if (is.na(minvalue)) {
    vec <- noLeadingZero(vec, fmt = paste0("%.", nd, "f"), nd = nd)
  } else {
    target <- if (minvalue) min(vec, na.rm = TRUE) else max(vec, na.rm = TRUE)
    tag <- rep(" ", length(vec))
    if (!is.na(char)) tag[vec == target] <- char
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
  fitTab <- object@fit[ , colnames(object@fit) %in% fit.measures, drop = FALSE]
  if (!return.diff) return(fitTab)

  ## or return differences in fit indices
  diffTab <- as.data.frame(do.call(cbind, lapply(fitTab, diff)))
  rownames(diffTab) <- paste(object@name[-1], "-", object@name[-length(object@name)])
  diffTab
}




