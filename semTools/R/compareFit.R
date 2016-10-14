### Sunthud Pornprasertmanit
### Last updated: 14 October 2016


setClass("FitDiff", representation(name = "vector", nested = "data.frame", ordernested = "vector", fit="data.frame"))

isNested <- function(object) length(object@ordernested) > 1 || !is.na(object@ordernested)

noLeadingZero <- function(vec, fmt) {
	out <- sprintf(fmt, vec)
	used <- vec < 1 & vec >= 0
	used[is.na(used)] <- FALSE
	out[used] <- substring(out[used], 2)
	out
}

setMethod("show", signature(object = "FitDiff"), function(object) {
    summary(object)
}) 

setMethod("summary", signature(object = "FitDiff"), function(object, fit.measures = "default") {
    if(isNested(object)) {
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
	if(is.null(fit.measures)) fit.measures <- "all"
	if(length(fit.measures) == 1) {
		if(fit.measures == "default") {
			fit.measures <- c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr", "aic", "bic")
		} else if (fit.measures == "all") {
			fit.measures <- colnames(object@fit)
		}
	}
	fitTab <- object@fit
	orderThing <- rep(NA, ncol(fitTab))
	orderThing[colnames(fitTab) %in% c("rmsea", "aic", "bic", "bic2", "srmr", "srmr_nomean", "rmr", "rmr_nomean", "ecvi")] <- TRUE 
	orderThing[colnames(fitTab) %in% c("pvalue", "cfi", "tli", "nnfi", "rfi", "nfi", "pnfi", "ifi", "rni", "cn_05", "cn_01", "gfi", "agfi", "pgfi", "mfi")] <- FALSE 
	isDF <- rep(FALSE, ncol(fitTab))
	isDF[grep("df", colnames(fitTab))] <- TRUE
	suppressWarnings(fitTab <- as.data.frame(mapply(tagDagger, fitTab, orderThing, is.df=isDF)))
	rownames(fitTab) <- object@name
	fitTab[,colnames(fitTab) %in% fit.measures]
}

saveFileFitDiff <- function(object, filewrite, what="summary", tableFormat=FALSE, fit.measures = "default") {
	if(tableFormat) {
		filetemplate <- file(filewrite, 'w')
		if(isNested(object)) {
			cat("Nested Model Comparison\n\n", file=filetemplate)
			out <- getNestedTable(object)
			out <- data.frame(model.diff = rownames(out), out)
			write.table(out, file=filetemplate, sep="\t", quote=FALSE, row.names=FALSE)
			cat("\n\n", file=filetemplate)
		}
		out2 <- getFitSummary(object, fit.measures)
		out2 <- data.frame(model = object@name, out2)
		cat("Fit Indices Summaries\n\n", file=filetemplate)
		write.table(out2, file=filetemplate, sep="\t", quote=FALSE, row.names=FALSE)
		close(filetemplate)
	} else {
		write(paste(capture.output(lavaan::summary(object)), collapse="\n"), file=filewrite)
	}
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

