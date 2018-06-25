### Sunthud Pornprasertmanit
### Last updated: 25 June 2018
### source code for compareFit() function and FitDiff class


## -----------------
## Class and Methods
## -----------------

#' Class For Representing A Template of Model Fit Comparisons
#'
#' This class contains model fit measures and model fit comparisons among
#' multiple models
#'
#'
#' @name FitDiff-class
#' @aliases FitDiff-class show,FitDiff-method summary,FitDiff-method
#' @docType class
#' @slot name The name of each model
#' @slot nested Model fit comparisons between adjacent nested models that are
#'  ordered based on their degrees of freedom (\emph{df})
#' @slot ordernested The order of nested models regarding to their \emph{df}
#' @slot fit Fit measures of all models specified in the \code{name} slot
#' @section Objects from the Class: Objects can be created via the
#'  \code{\link{compareFit}} function.
#' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#' @seealso \code{\link{compareFit}}; \code{\link{clipboard}}
#' @examples
#'
#' HW.model <- ' visual =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed =~ x7 + x8 + x9 '
#'
#' out <- measurementInvariance(model = HW.model, data = HolzingerSwineford1939,
#'                              group = "school", quiet = TRUE)
#' modelDiff <- compareFit(out)
#' summary(modelDiff)
#' summary(modelDiff, fit.measures = "all")
#' summary(modelDiff, fit.measures = c("aic", "bic"))
#'
#' \dontrun{
#' ## Save results to a file
#' saveFile(modelDiff, file = "modelDiff.txt")
#'
#' ## Copy to a clipboard
#' clipboard(modelDiff)
#' }
#'
setClass("FitDiff", representation(name = "vector",
                                   nested = "data.frame",
                                   ordernested = "vector",
                                   fit = "data.frame"))

#' @rdname FitDiff-class
#' @export
setMethod("show", signature(object = "FitDiff"), function(object) {
    summary(object)
})

#' @rdname FitDiff-class
#' @param object object of class \code{FitDiff}
#' @param fit.measures \code{character} vector naming fit indices the user can
#'  request from \code{\link[lavaan]{fitMeasures}}. If \code{"default"}, the
#'  fit measures will be \code{c("chisq", "df", "pvalue", "cfi", "tli",
#'  "rmsea", "srmr", "aic", "bic")}. If \code{"all"}, all available fit measures
#'  will be returned.
#' @export
setMethod("summary", signature(object = "FitDiff"),
          function(object, fit.measures = "default") {
  if (isNested(object)) {
		cat("################### Nested Model Comparison #########################\n")
		print(getNestedTable(object))
		cat("\n")
	}
	cat("#################### Fit Indices Summaries ##########################\n")
	print(getFitSummary(object, fit.measures))
})

getNestedTable <- function(object) {
	ord <- object@ordernested
	nameDiff <- paste(object@name[ord[-1]], "-", object@name[ord[-length(ord)]])
	pprint <- noLeadingZero(object@nested$p, "%4.3f")
	pprint[object@nested$p < 0.001] <- " <.001"
	nestedTab <- object@nested
	nestedTab$chi <- sprintf("%.2f", object@nested$chi)
	nestedTab$p <- pprint
	nestedTab$delta.cfi <- sprintf("%.4f", object@nested$delta.cfi)
	nestedTab <- data.frame(nestedTab)
	rownames(nestedTab) <- nameDiff
	nestedTab[nrow(nestedTab):1,]
}

getFitSummary <- function(object, fit.measures = "default") {
	if (is.null(fit.measures)) fit.measures <- "all"
	if (length(fit.measures) == 1) {
		if (fit.measures == "default") {
			fit.measures <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr", "aic", "bic")
		} else if (fit.measures == "all") {
			fit.measures <- colnames(object@fit)
		}
	}
	fitTab <- object@fit
	orderThing <- rep(NA, ncol(fitTab))
	orderThing[colnames(fitTab) %in% c("rmsea","aic","bic","bic2","srmr","rmr",
	                                   "srmr_nomean","rmr_nomean","ecvi")] <- TRUE
	orderThing[colnames(fitTab) %in% c("pvalue","cfi","tli","nnfi","rfi","nfi",
	                                   "pnfi","ifi","rni","cn_05","cn_01","gfi",
	                                   "agfi","pgfi","mfi")] <- FALSE
	isDF <- rep(FALSE, ncol(fitTab))
	isDF[grep("df", colnames(fitTab))] <- TRUE
	suppressWarnings(fitTab <- as.data.frame(mapply(tagDagger, fitTab, orderThing, is.df = isDF)))
	rownames(fitTab) <- object@name
	fitTab[ , colnames(fitTab) %in% fit.measures]
}

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

		if (isNested(object)) {
			cat("Nested Model Comparison\n\n", file = file, append = TRUE)
			out <- getNestedTable(object)
			out <- data.frame(model.diff = rownames(out), out)
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

#' Build an object summarizing fit indices across multiple models
#'
#' This function will create the template to compare fit indices across
#' multiple fitted lavaan objects. The results can be exported to a clipboard
#' or a file later.
#'
#'
#' @param ...  fitted \code{lavaan} models or list(s) of \code{lavaan} objects
#' @param nested \code{logical} indicating whether the models in \code{...} are
#'  nested. See the \code{\link{net}} function for an empirical test of nesting.
#' @return A \code{\linkS4class{FitDiff}} object that saves model fit
#' comparisons across multiple models. If the output is not assigned as an
#' object, the output is printed in two parts: (1) nested model comparison (if
#' models are nested) and (2) summary of fit indices. In the fit indices
#' summaries, daggers are tagged to the model with the best fit according to
#' each fit index.
#' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#' @seealso \code{\linkS4class{FitDiff}}, \code{\link{clipboard}}
#' @examples
#'
#' m1 <- ' visual  =~ x1 + x2 + x3
#'         textual =~ x4 + x5 + x6
#'         speed   =~ x7 + x8 + x9 '
#'
#' fit1 <- cfa(m1, data = HolzingerSwineford1939)
#'
#' m2 <- ' f1  =~ x1 + x2 + x3 + x4
#'         f2 =~ x5 + x6 + x7 + x8 + x9 '
#' fit2 <- cfa(m2, data = HolzingerSwineford1939)
#' compareFit(fit1, fit2, nested = FALSE)
#'
#' HW.model <- ' visual =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed =~ x7 + x8 + x9 '
#'
#' out <- measurementInvariance(model = HW.model, data = HolzingerSwineford1939,
#'                              group = "school", quiet = TRUE)
#' compareFit(out)
#'
#' @export
compareFit <- function(..., nested = TRUE) {
	arg <- match.call()
	mods <- input <- list(...)
	if(any(sapply(mods, is, "list"))) {
		temp <- list()
		for(i in seq_along(mods)) {
			if(!is(mods[[i]], "list")) { temp <- c(temp, list(mods[[i]])) } else { temp <- c(temp, mods[[i]]) }
		}
		mods <- temp
	}
	if(any(!sapply(mods, is, "lavaan"))) stop("Some models specified here are not lavaan outputs or list of lavaan outputs")
	nameMods <- NULL
	tempname <- as.list(arg)[-1]
	if(!is.null(names(tempname))) tempname <- tempname[!(names(tempname) %in% "nested")]
	tempname <- lapply(tempname, as.character)
	for(i in seq_along(input)) {
		if(is(input[[i]], "list")) {
			if(length(tempname[[i]]) == 1) {
				temp2 <- paste0(tempname[[i]], "[[", seq_along(input[[i]]), "]]")
				if(!is.null(names(input[[i]]))) temp2 <- names(input[[i]])
				nameMods <- c(nameMods, temp2)
			} else {
				temp2 <- tempname[[i]][tempname[[i]] != "list"]
				nameMods <- c(nameMods, temp2)
			}
		} else {
			nameMods <- c(nameMods, tempname[[i]])
		}
	}
	nestedout <- data.frame()
	ord <- NA
	if(nested) {
		dfs <- sapply(mods, function(x) lavaan::fitMeasures(x)["df"])
		ord <- order(dfs, decreasing = TRUE)
		modsTemp <- mods[ord]
		modsA <- modsTemp[-1]
		modsB <- modsTemp[-length(mods)]
		chisqdiff <- NULL
		dfdiff <- NULL
		pdiff <- NULL
		cfidiff <- NULL
		# Need the for loop because the mapply function does not work.
		for(i in seq_along(modsA)) {
			fitA <- modsA[[i]]
			fitB <- modsB[[i]]
			fitDiff <- lavaan::anova(fitA, fitB)
			cfidiff <- c(cfidiff, lavaan::fitMeasures(fitA)["cfi"] - lavaan::fitMeasures(fitB)["cfi"])
			chisqdiff <- c(chisqdiff, fitDiff["Chisq diff"][2, 1])
			dfdiff <- c(dfdiff, fitDiff["Df diff"][2, 1])
			pdiff <- c(pdiff, fitDiff["Pr(>Chisq)"][2, 1])
		}
		nestedout <- data.frame(chi = chisqdiff, df = dfdiff, p = pdiff, delta.cfi = cfidiff)
	}
	fit <- as.data.frame(t(sapply(mods, lavaan::fitMeasures)))
	new("FitDiff", name = nameMods, nested = nestedout, ordernested = ord, fit = fit)
}



## ----------------
## Hidden Functions
## ----------------

isNested <- function(object) length(object@ordernested) > 1 || !is.na(object@ordernested)

noLeadingZero <- function(vec, fmt) {
  out <- sprintf(fmt, vec)
  used <- vec < 1 & vec >= 0
  used[is.na(used)] <- FALSE
  out[used] <- substring(out[used], 2)
  out
}

tagDagger <- function(vec, minvalue = NA, is.df = FALSE) {
  if(is.na(minvalue)) {
    if(is.df) {
      vec <- noLeadingZero(vec, fmt="%.0f")
    } else {
      vec <- noLeadingZero(vec, fmt="%.3f")
    }
  } else  {
    target <- max(vec, na.rm=TRUE)
    if (minvalue) {
      target <- min(vec, na.rm=TRUE)
    }
    tag <- rep(" ", length(vec))
    tag[vec == target] <- "\u2020"
    vec <- noLeadingZero(vec, fmt="%.3f")
    vec <- paste0(vec, tag)
  }
  vec
}


