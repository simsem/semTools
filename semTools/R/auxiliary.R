## Title: Automatically accounts for auxiliary variable in full information maximum likelihood
## Author: Sunthud Pornprasertmanit
# Description: Automatically accounts for auxiliary variable in full information maximum likelihood

setClass("lavaanStar", contains = "lavaan", representation(nullfit = "vector", imputed="list"), prototype(nullfit=c(chi=0,df=0), imputed=list()))

setMethod("inspect", "lavaanStar",
function(object, what="free") {
	what <- tolower(what)
	if(what == "fit" ||
              what == "fitmeasures" ||
              what == "fit.measures" ||
              what == "fit.indices") {
		fitMeasuresLavaanStar(object)
	} else if(what == "imputed" ||
              what == "impute") {
		result <- object@imputed
		if(length(result) > 0) {
			return(result)
		} else {
			stop("This method did not made by multiple imputation.")
		}
	} else {
		getMethod("inspect", "lavaan")(object, what=what)
	}
})

setMethod("summary", "lavaanStar",
function(object, fit.measures=FALSE, ...) {
	getMethod("summary", "lavaan")(object, fit.measures=FALSE, ...)
	if(fit.measures) {
		cat("Because the original method to find the baseline model does not work, \nplease do not use any fit measures relying on baseline model, including CFI and TLI. \nTo find the correct one, please use the inspect function: inspect(object, what='fit').\n")
	}
})
	
# auxiliary: Automatically accounts for auxiliary variable in full information maximum likelihood

cfa.auxiliary <- function(model, aux, ...) {
	auxiliary(model = model, aux = aux, fun = "cfa", ...)
}

sem.auxiliary <- function(model, aux, ...) {
	auxiliary(model = model, aux = aux, fun = "sem", ...)
}

growth.auxiliary <- function(model, aux, ...) {
	auxiliary(model = model, aux = aux, fun = "growth", ...)
}

lavaan.auxiliary <- function(model, aux, ...) {
	auxiliary(model = model, aux = aux, fun = "lavaan", ...)
}

auxiliary <- function(model, aux, fun, ...) {
	args <- list(...)
	args$fixed.x <- FALSE
	
	if(is(model, "lavaan")) {
		model <- model@ParTable
	}
	
	# Check whether the model is a parameter table
	if(!(is.list(model) && ("lhs" %in% names(model)))) {
		fit <- do.call(fun, c(list(model=model, do.fit=FALSE), args))
		model <- fit@ParTable
	}
	
	if(any(model$exo == 1)) {
		stop("All exogenous variables (covariates) must be treated as endogenous variables by the 'auxiliary' function (fixed.x = FALSE).")
	}
	
	facName <- NULL
	indName <- NULL
	singleIndicator <- NULL
	facName <- unique(model$lhs[model$op == "=~"])
	indName <- setdiff(unique(model$rhs[model$op == "=~"]), facName)
	singleIndicator <- setdiff(unique(c(model$lhs, model$rhs)), c(facName, indName, ""))
	facSingleIndicator <- paste0("f", singleIndicator)
	for(i in seq_along(singleIndicator)) {
		model$lhs <- gsub(singleIndicator[i], facSingleIndicator[i], model$lhs)
		model$rhs <- gsub(singleIndicator[i], facSingleIndicator[i], model$rhs)
	}
		
	ngroups <- max(model$group)
		if(!is.null(singleIndicator) && length(singleIndicator) != 0) model <- attachPT(model, facSingleIndicator, "=~", singleIndicator, ngroups, fixed = TRUE, ustart = 1, expand = FALSE)
		if(!is.null(singleIndicator) && length(singleIndicator) != 0) model <- attachPT(model, singleIndicator, "~~", singleIndicator, ngroups, fixed = TRUE, ustart = 0, expand = FALSE)
	if(is.null(indName) || length(indName) == 0) {
		faux <- paste0("f", aux)
		model <- attachPT(model, faux, "=~", aux, ngroups, fixed = TRUE, ustart = 1, expand = FALSE)
		model <- attachPT(model, aux, "~~", aux, ngroups, fixed = TRUE, ustart = 0, expand = FALSE)
		model <- attachPT(model, facSingleIndicator, "~~", faux, ngroups)
		model <- attachPT(model, faux, "~~", faux, ngroups, symmetric=TRUE)
		if(any(model$op == "~1")) model <- attachPT(model, faux, "~1", "", ngroups)
	} else {
		if(!is.null(indName) && length(indName) != 0) model <- attachPT(model, indName, "~~", aux, ngroups)
		model <- attachPT(model, aux, "~~", aux, ngroups, symmetric=TRUE, useUpper=TRUE)
		if(!is.null(singleIndicator) && length(singleIndicator) != 0) model <- attachPT(model, facSingleIndicator, "=~", aux, ngroups)
		if(any(model$op == "~1")) model <- attachPT(model, aux, "~1", "", ngroups)
	}
	args$model <- model
	result <- do.call(fun, args)
	
	codeNull <- nullAuxiliary(aux, union(indName, singleIndicator), NULL, any(model$op == "~1"), max(model$group))
	resultNull <- lavaan(codeNull, ...)
	result <- as(result, "lavaanStar")
	fit <- fitMeasures(resultNull)
	name <- names(fit)
	fit <- as.vector(fit)
	names(fit) <- name
	result@nullfit <- fit
	return(result)
}

attachPT <- function(pt, lhs, op, rhs, ngroups, symmetric=FALSE, exo=FALSE, fixed=FALSE, useUpper=FALSE, ustart = NA, expand = TRUE, diag = TRUE) {
	if(expand) {
		element <- expand.grid(lhs, rhs, stringsAsFactors = FALSE)
	} else {
		element <- cbind(lhs, rhs)
	}
	if(symmetric) { 
		if(useUpper) {
			element <- element[as.vector(upper.tri(diag(length(lhs)), diag=diag)),]
		} else {
			element <- element[as.vector(lower.tri(diag(length(lhs)), diag=diag)),]
		}
	}
	num <- nrow(element) * ngroups
	pt$id <- c(pt$id, (max(pt$id)+1):(max(pt$id)+num))
	pt$lhs <- c(pt$lhs, rep(element[,1], ngroups))
	pt$op <- c(pt$op, rep(op, num))
	pt$rhs <- c(pt$rhs, rep(element[,2], ngroups))
	pt$user <- c(pt$user, rep(1, num))
	pt$group <- c(pt$group, rep(1:ngroups, each=nrow(element)))
	free <- (max(pt$free)+1):(max(pt$free)+num)
	if(fixed) free <- rep(0L, num)
	pt$free <- c(pt$free, free)
	pt$ustart <- c(pt$ustart, rep(ustart, num))
	pt$exo <- c(pt$exo, rep(as.numeric(exo), num))
	pt$label <- c(pt$label, rep("", num))
	pt$eq.id <- c(pt$eq.id, rep(0L, num))
	unco <- (max(pt$unco)+1):(max(pt$unco)+num)
	if(fixed) unco <- rep(0L, num)
	pt$unco <- c(pt$unco, unco)
	return(pt)
}

nullAuxiliary <- function(aux, indName, covName=NULL, meanstructure, ngroups) {
	covName <- rev(covName)
	pt <- list()
	num <- length(indName) * ngroups
	if(meanstructure) num <- num*2
	pt$id <- 1:num
	pt$lhs <- rep(indName, ngroups)
	pt$op <- rep("~~", num)
	pt$rhs <- rep(indName, ngroups)
	pt$user <- rep(1, num)
	pt$group <- rep(1:ngroups, each=length(indName))
	pt$free <- 1:num
	pt$ustart <- rep(NA, num)
	pt$exo <- rep(0, num)
	pt$label <- rep("", num)
	pt$eq.id <- rep(0, num)
	pt$unco <- 1:num
	if(meanstructure) {
		pt$lhs <- rep(rep(indName, ngroups), 2)
		pt$op <- rep(c("~~", "~1"), each=num/2)
		pt$rhs <- c(rep(indName, ngroups), rep("", num/2))
		pt$group <- rep(rep(1:ngroups, each=length(indName)), 2)
	}
	pt <- attachPT(pt, aux, "~~", aux, ngroups, symmetric=TRUE)
	pt <- attachPT(pt, indName, "~~", aux, ngroups)
	if(meanstructure) pt <- attachPT(pt, aux, "~1", "", ngroups)
	if(!is.null(covName) && length(covName) != 0) {
		pt <- attachPT(pt, aux, "~~", covName, ngroups)	
		pt <- attachPT(pt, covName, "~~", covName, ngroups, symmetric=TRUE, useUpper=TRUE)
		if(meanstructure) pt <- attachPT(pt, covName, "~1", "", ngroups)
	}
	return(pt)
}


fitMeasuresLavaanStar <- function(object) {
	result <- getMethod("inspect", "lavaan")(object, what="fit")
	result[c("baseline.chisq", "baseline.df", "baseline.pvalue")] <- object@nullfit[c("chisq", "df", "pvalue")]
		
    if(object@Options$test %in% c("satorra.bentler", "yuan.bentler")) {
        scaled <- TRUE
    } else {
        scaled <- FALSE
    }
    
	if(scaled) {
		result[c("baseline.chisq.scaled", "baseline.df.scaled", "baseline.pvalue.scaled", "baseline.chisq.scaling.factor")] <- object@nullfit[c("chisq.scaled", "df.scaled", "pvalue.scaled", "chisq.scaling.factor")]
	}
	
	
	X2.null <- object@nullfit["chisq"]
	df.null <- object@nullfit["df"]
	X2 <- result["chisq"]
	df <- result["df"]
	result["cfi"] <- ( 1 - max(c(X2 - df,0)) / 
                                    max( c(X2-df, X2.null-df.null, 0) ) 
                              )
	if(df > 0) {
		result["tli"] <- (X2.null/df.null - X2/df)/(X2.null/df.null - 1)
	} else {
		result["tli"] <- 1
	}
	if(scaled) {
		X2.scaled <- result["chisq.scaled"]
		df.scaled <- result["df.scaled"]
		X2.null.scaled <- object@nullfit["chisq.scaled"]
		df.null.scaled <- object@nullfit["df.scaled"]
		result["cfi.scaled"] <- 
			( 1 - max( c(X2.scaled - df.scaled,0)) / 
				  max( c(X2.scaled - df.scaled, 
						 X2.null.scaled - df.null.scaled, 0) ) 
			)
		if(df > 0) {
			result["tli.scaled"] <- (X2.null.scaled/df.null.scaled - 
					X2.scaled/df.scaled) / 
				   (X2.null.scaled/df.null.scaled - 1)
		} else {
			result["tli.scaled"] <- 1
		}
	} 
	return(result)
		
}
