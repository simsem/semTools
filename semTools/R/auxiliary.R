## Title: Automatically accounts for auxiliary variable in full information maximum likelihood
## Author: Sunthud Pornprasertmanit
# Description: Automatically accounts for auxiliary variable in full information maximum likelihood

setClass("lavaanStar", contains = "lavaan", representation(nullfit = "vector"), prototype(nullfit=c(chi=0,df=0)))

setMethod("inspect", "lavaanStar",
function(object, what="free") {
	if(what == "fit" ||
              what == "fitmeasures" ||
              what == "fit.measures" ||
              what == "fit.indices") {
		fitMeasuresLavaanStar(object)
	} else {
		getMethod("inspect", "lavaan")(object, what=what)
	}
})

setMethod("summary", "lavaanStar",
function(object, fit.measures=FALSE, ...) {
	getMethod("summary", "lavaan")(object, fit.measures=FALSE, ...)
	if(fit.measures) {
		cat("Because the original method to find the baseline model does not work, \n
		    please do not use any fit measures relying on baseline model, including CFI and TLI. \n
			To find the correct one, please use the inspect function: inspect(object, what='fit').\n")
	}
})
	
# auxiliary: Automatically accounts for auxiliary variable in full information maximum likelihood
auxiliary <- function(object, aux, ...) {
	if(is(object, "lavaan")) {
		object <- object@ParTable
	}
	facName <- NULL
	indName <- NULL
	covName <- NULL
	if("=~" %in% object$op) {
		facName <- unique(object$lhs[object$op == "=~"])
		indName <- unique(object$rhs[object$op == "=~"])
		covName <- setdiff(unique(object$rhs), c(facName, indName, ""))
	} else {
		indName <- unique(object$lhs[object$op == "~"])
		covName <- setdiff(unique(object$rhs), c(indName, ""))
	}
	
	code <- auxiliaryCode(object, aux, correlate=indName, extradv=covName)
	result <- lavaan(code, ...)
	codeNull <- nullAuxiliary(aux, indName, covName, any(object$op == "~1"), max(object$group))
	resultNull <- lavaan(codeNull, ...)
	result <- as(result, "lavaanStar")
	fit <- fitMeasures(resultNull)
	name <- names(fit)
	fit <- as.vector(fit)
	names(fit) <- name
	result@nullfit <- fit
	return(result)
}

auxiliaryCode <- function(object, aux, correlate=NULL, extradv=NULL) {
	ngroups <- max(object$group)
	object <- attachPT(object, aux, "~~", aux, ngroups, symmetric=TRUE)
	if(!is.null(correlate) && length(correlate) != 0) object <- attachPT(object, correlate, "~~", aux, ngroups)
	if(!is.null(extradv) && length(extradv) != 0) object <- attachPT(object, aux, "~", extradv, ngroups)
	if(any(object$op == "~1")) object <- attachPT(object, aux, "~1", "", ngroups)
	return(object)
}

attachPT <- function(pt, lhs, op, rhs, ngroups, symmetric=FALSE, exo=FALSE, fixed=FALSE, useUpper=FALSE) {
	element <- expand.grid(lhs, rhs, stringsAsFactors = FALSE)
	if(symmetric) { 
		if(useUpper) {
			element <- element[as.vector(upper.tri(diag(length(lhs)), diag=TRUE)),]
		} else {
			element <- element[as.vector(lower.tri(diag(length(lhs)), diag=TRUE)),]
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
	if(fixed) free <- rep(0, num)
	pt$free <- c(pt$free, free)
	pt$ustart <- c(pt$ustart, rep(NA, num))
	pt$exo <- c(pt$exo, rep(as.numeric(exo), num))
	pt$label <- c(pt$label, rep("", num))
	pt$eq.id <- c(pt$eq.id, rep(0, num))
	unco <- (max(pt$unco)+1):(max(pt$unco)+num)
	if(fixed) unco <- rep(0, num)
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

# In case we will need it later.
# reorderPT <- function(pt, aux, indName, facName) {
	# pt1 <- lapply(pt, function(x, select) x[select], select=(pt$op == "=~"))
	# pt2 <- lapply(pt, function(x, select) x[select], select=(pt$op == "~"))
	# pt31 <- lapply(pt, function(x, select) x[select], select=((pt$op == "~~") & (pt$lhs %in% aux)))
	# pt32 <- lapply(pt, function(x, select) x[select], select=((pt$op == "~~") & (pt$lhs %in% indName)))
	# pt33 <- lapply(pt, function(x, select) x[select], select=((pt$op == "~~") & (pt$lhs %in% facName)))
	# pt4 <- lapply(pt, function(x, select) x[select], select=(pt$op == "~1"))
	# pt <- mapply(pt1, pt2, pt31, pt32, pt33, pt4, FUN=c, SIMPLIFY=FALSE)
	# pt$id <- 1:max(pt$id)
	# pt$free <- sort(unique(pt$free))[sapply(as.list(pt$free), function(x, new) which(x == new), new=unique(pt$free))]
	# pt$eq.id <- sort(unique(pt$eq.id))[sapply(as.list(pt$eq.id), function(x, new) which(x == new), new=unique(pt$eq.id))]
	# pt$unco <- sort(unique(pt$unco))[sapply(as.list(pt$unco), function(x, new) which(x == new), new=unique(pt$unco))]
	
	
	# return(pt)
# }
