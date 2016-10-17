### Title: Automatically accounts for auxiliary variable in full information maximum likelihood
### Author: Sunthud Pornprasertmanit
### Last updated: 17 October 2016
### Description: Automatically accounts for auxiliary variable in full information maximum likelihood

setClass("lavaanStar", contains = "lavaan", representation(nullfit = "vector", imputed="list", imputedResults="list", auxNames="vector"), prototype(nullfit=c(chi=0,df=0), imputed=list(), imputedResults=list(), auxNames = ""))

setMethod("inspect", "lavaanStar",  ## FIXME: get rid of lavaanStar object
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
	} else if(what == "aux" ||
              what == "auxiliary") {
		print(object@auxNames)
	} else {
		getMethod("inspect", "lavaan")(object, what=what)  ## FIXME: don't set a new inspect method
	}
})

setMethod("summary", "lavaanStar",
function(object, fit.measures=FALSE, ...) {
	getMethod("summary", "lavaan")(object, fit.measures=FALSE, ...)
	if(fit.measures) {
		cat("Because the original method to find the baseline model does not work, \nplease do not use any fit measures relying on baseline model, including CFI and TLI. \nTo find the correct one, please use the function: lavInspect(object, what='fit').\n")
	}
})

setMethod("anova", signature(object = "lavaanStar"),
function(object, ...) {
	imputed <- object@imputed
	if(length(imputed) > 0) {
		dots <- list(...)
		if(length(dots) > 1) stop("Multiple Imputed Results: Cannot compare more than two objects")
		object2 <- dots[[1]]
		imputed2 <- object2@imputed
		if(length(imputed) == 0) stop("The second object must come from multiple imputation.")
		listlogl1 <- imputed[["indivlogl"]]
		listlogl2 <- imputed2[["indivlogl"]]
		df1 <- lavaan::lavInspect(object, "fit")["df"]
		df2 <- lavaan::lavInspect(object2, "fit")["df"]
		if(df2 < df1) {
			templogl <- listlogl1
			listlogl1 <- listlogl2
			listlogl2 <- templogl
		}
		dfdiff <- df2 - df1
		anovaout <- mapply(lavaan::anova, object@imputedResults, object2@imputedResults, SIMPLIFY = FALSE)
		chidiff <- sapply(anovaout, function(u) u[2, "Chisq diff"])
		dfdiff2 <- mean(sapply(anovaout, function(u) u[2, "Df diff"]))
		fit.altcc <- mean(chidiff)
		naive <- c(fit.altcc, dfdiff2, 1 - pchisq(fit.altcc, dfdiff2))
		names(naive) <- c("chisq", "df", "pvalue")	
		lmrr <- lmrrPooledChi(chidiff, dfdiff2)
		mr <- NULL
		mplus <- NULL
		if(!is.null(listlogl1[["loglmod"]]) | !is.null(listlogl2[["loglmod"]])) {
			logl1 <- listlogl1[["loglmod"]]
			logl2 <- listlogl2[["loglmod"]]
			chimean <- mean((logl1 - logl2)*2)
			m <- length(logl1)
			ariv <- ((m+1)/((m-1)*dfdiff))*(fit.altcc-chimean)

			mplus <- mplusPooledChi(chimean, dfdiff, ariv)
			mr <- mrPooledChi(chimean, m, dfdiff, ariv)
		}
		result <- list(naive = naive, lmrr = lmrr, mr = mr, mplus = mplus)
		return(result)
	} else {
		return(getMethod("anova", "lavaan")(object, ...))
	}
})

setMethod("vcov", "lavaanStar",
function(object, ...) {
	result <- object@imputed
	if(length(result) == 0) {
		return(getMethod("vcov", "lavaan")(object, ...))
	} else {
		out <- object@vcov$vcov
		rownames(out) <- colnames(out) <- lavaan::lav_partable_labels(lavaan::partable(object), type="free")
		return(out)
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
	args$missing <- "fiml"
	
	if(is(model, "lavaan")) {
		if(!lavaan::lavInspect(model, "options")$meanstructure) stop("The lavaan fitted model must evaluate the meanstructure. Please re-fit the lavaan object again with 'meanstructure=TRUE'")
		model <- lavaan::parTable(model)
	} else if(!(is.list(model) && ("lhs" %in% names(model)))) {
		fit <- do.call(fun, c(list(model=model, do.fit=FALSE), args))
		model <- lavaan::parTable(fit)
	}
	model <- model[setdiff(1:length(model), which(names(model) == "start"))]

	if(any(model$exo == 1)) {
		stop("All exogenous variables (covariates) must be treated as endogenous variables by the 'auxiliary' function (fixed.x = FALSE).")
	}

	auxResult <- craftAuxParTable(model = model, aux = aux, ...)
	if(checkOrdered(args$data, auxResult$indName, ...)) {
		stop("The analysis model or the analysis data have ordered categorical variables. The auxiliary variable feature is not available for the models for categorical variables with the weighted least square approach.")
	}
	
	args$model <- auxResult$model
	result <- do.call(fun, args)
	
	codeNull <- nullAuxiliary(aux, auxResult$indName, NULL, any(model$op == "~1"), max(model$group))
	resultNull <- lavaan::lavaan(codeNull, ...)
	result <- as(result, "lavaanStar")
	fit <- lavaan::fitMeasures(resultNull)
	name <- names(fit)
	fit <- as.vector(fit)
	names(fit) <- name
	result@nullfit <- fit
	result@auxNames <- aux
	return(result)
}

checkOrdered <- function(dat, varnames, ...) {
	ord <- list(...)$ordered
	if(is.null(ord)) { 
		ord <- FALSE
	} else {
		ord <- TRUE
	}
	if(is.null(dat)) {
		orderedVar <- FALSE
	} else {
		orderedVar <- sapply(dat[,varnames], function(x) "ordered" %in% is(x))
	}
	any(c(ord, orderedVar))
}

craftAuxParTable <- function(model, aux, ...) {
	constraintLine <- model$op %in% c("==", ":=", ">", "<")
	modelConstraint <- lapply(model, "[", constraintLine)
	model <- lapply(model, "[", !constraintLine)
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
	if(!is.null(singleIndicator) && length(singleIndicator) != 0) model <- attachPT(model, singleIndicator, "~1", "", ngroups, fixed = TRUE, ustart = 0, expand = FALSE)
	if(is.null(indName) || length(indName) == 0) {
		faux <- paste0("f", aux)
		model <- attachPT(model, faux, "=~", aux, ngroups, fixed = TRUE, ustart = 1, expand = FALSE)
		model <- attachPT(model, aux, "~~", aux, ngroups, fixed = TRUE, ustart = 0, expand = FALSE)
		model <- attachPT(model, facSingleIndicator, "~~", faux, ngroups)
		model <- attachPT(model, faux, "~~", faux, ngroups, symmetric=TRUE)
		if(any(model$op == "~1")) {
			model <- attachPT(model, faux, "~1", "", ngroups)
			model <- attachPT(model, aux, "~1", "", ngroups, fixed = TRUE, ustart = 0, expand = FALSE)
		}
	} else {
		if(!is.null(indName) && length(indName) != 0) model <- attachPT(model, indName, "~~", aux, ngroups)
		model <- attachPT(model, aux, "~~", aux, ngroups, symmetric=TRUE, useUpper=TRUE)
		if(!is.null(singleIndicator) && length(singleIndicator) != 0) model <- attachPT(model, facSingleIndicator, "=~", aux, ngroups)
		if(any(model$op == "~1")) model <- attachPT(model, aux, "~1", "", ngroups)
	}
	model <- attachConstraint(model, modelConstraint)

	list(model = model, indName = union(indName, singleIndicator))
}

attachConstraint <- function(pt, con) {
	len <- length(con$id)
	if(len > 0) {
		pt$id <- c(pt$id, (max(pt$id)+1):(max(pt$id)+len))
		pt$lhs <- c(pt$lhs, con$lhs)
		pt$op <- c(pt$op, con$op)
		pt$rhs <- c(pt$rhs, con$rhs)
		pt$user <- c(pt$user, con$user)
		pt$group <- c(pt$group, con$group)
		pt$free <- c(pt$free, con$free)
		pt$ustart <- c(pt$ustart, con$ustart)
		pt$exo <- c(pt$exo, con$exo)
		pt$label <- c(pt$label, con$label)
		pt$plabel <- c(pt$plabel, con$plabel)
		pt$start <- c(pt$start, con$start)
		pt$est <- c(pt$est, con$est)
		pt$se <- c(pt$se, con$se)
	}
	pt
}

attachPT <- function(pt, lhs, op, rhs, ngroups, symmetric=FALSE, exo=FALSE, fixed=FALSE, useUpper=FALSE, ustart = NA, expand = TRUE, diag = TRUE) {
	pt$start <- pt$est <- pt$se <- NULL
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
	pt$plabel <- c(pt$plabel, rep("", num))
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
	pt$plabel <- rep("", num)
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
	notused <- capture.output(result <- suppressWarnings(getMethod("inspect", "lavaan")(object, what="fit"))) ## FIXME: don't set a new inspect method
	result[c("baseline.chisq", "baseline.df", "baseline.pvalue")] <- object@nullfit[c("chisq", "df", "pvalue")]
		
    if(lavaan::lavInspect(object, "options")$test %in% c("satorra.bentler", "yuan.bentler", 
                   "mean.var.adjusted", "scaled.shifted")) {
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
	
	if(df.null == 0) {
		result["cfi"] <- NA
		result["tli"] <- NA
		result["nnfi"] <- NA
		result["rfi"] <- NA
		result["nfi"] <- NA
		result["pnfi"] <- NA
		result["ifi"] <- NA
		result["rni"] <- NA
	} else {
		# CFI
		if("cfi" %in% names(result)) {
			t1 <- max( c(X2 - df, 0) )
			t2 <- max( c(X2 - df, X2.null - df.null, 0) )
			if(t1 == 0 && t2 == 0) {
				result["cfi"] <- 1
			} else {
				result["cfi"] <- 1 - t1/t2
			}
		}
		
		# TLI
		if("tli" %in% names(result)) {
			if(df > 0) {
				t1 <- X2.null/df.null - X2/df
				t2 <- X2.null/df.null - 1
				# note: TLI original formula was in terms of fx/df, not X2/df
				# then, t1 <- fx_0/df.null - fx/df
				# t2 <- fx_0/df.null - 1/N (or N-1 for wishart)
				if(t1 < 0 && t2 < 0) {
					TLI <- 1
				} else {
					TLI <- t1/t2
				}
			} else {
			   TLI <- 1
			}
			result["tli"] <- result["nnfi"] <- TLI
		}

		# RFI
		if("rfi" %in% names(result)) {
			if(df > 0) {
				t1 <- X2.null/df.null - X2/df
				t2 <- X2.null/df.null
				if(t1 < 0 || t2 < 0) {
					RLI <- 1
				} else {
					RLI <- t1/t2
				}
			} else {
			   RLI <- 1
			}
			result["rfi"] <- RLI
		}
		
		# NFI
		if("nfi" %in% names(result)) {
			t1 <- X2.null - X2
			t2 <- X2.null
			NFI <- t1/t2
			result["nfi"] <- NFI
		}
		
		# PNFI
		if("pnfi" %in% names(result)) {
			t1 <- X2.null - X2
			t2 <- X2.null
			PNFI <- (df/df.null) * t1/t2
			result["pnfi"] <- PNFI
		}
		
		# IFI
		if("ifi" %in% names(result)) {
			t1 <- X2.null - X2
			t2 <- X2.null - df
			if(t2 < 0) {
				IFI <- 1
			} else {
				IFI <- t1/t2
			}
			result["ifi"] <- IFI
		}
		
		# RNI
		if("rni" %in% names(result)) {
			t1 <- X2 - df
			t2 <- X2.null - df.null
			if(df.null == 0) {
				RNI <- NA
			} else if(t1 < 0 || t2 < 0) {
				RNI <- 1
			} else {
				RNI <- 1 - t1/t2
			}
			result["rni"] <- RNI
		}
	}
	
	if(scaled) {
		X2.scaled <- result["chisq.scaled"]
		df.scaled <- result["df.scaled"]
		X2.null.scaled <- object@nullfit["chisq.scaled"]
		df.null.scaled <- object@nullfit["df.scaled"]
		
		if(df.null.scaled == 0) {
			result["cfi.scaled"] <- NA
			result["tli.scaled"] <- result["nnfi.scaled"] <- NA
			result["rfi.scaled"] <- NA
			result["nfi.scaled"] <- NA
			result["pnfi.scaled"] <- NA
			result["ifi.scaled"] <- NA
			result["rni.scaled"] <- NA
		} else {
			if("cfi.scaled" %in% names(result)) {
				t1 <- max( c(X2.scaled - df.scaled, 0) )
				t2 <- max( c(X2.scaled - df.scaled,
							 X2.null.scaled - df.null.scaled, 0) )
				if(t1 == 0 && t2 == 0) {
					result["cfi.scaled"] <- 1
				} else {
					result["cfi.scaled"] <- 1 - t1/t2
				}
			}
			
			if("tli.scaled" %in% names(result)) {
				if(df > 0) {
					t1 <- X2.null.scaled/df.null.scaled - X2.scaled/df.scaled
					t2 <- X2.null.scaled/df.null.scaled - 1
					if(t1 < 0 && t2 < 0) {
						TLI <- 1
					} else {
						TLI <- t1/t2
					}
				} else {
					TLI <- 1
				}
				result["tli.scaled"] <- result["nnfi.scaled"] <- TLI
			}
			
			if("rfi.scaled" %in% names(result)) {
				if(df > 0) {
					t1 <- X2.null.scaled/df.null.scaled - X2.scaled/df.scaled
					t2 <- X2.null.scaled/df.null.scaled
					if(t1 < 0 || t2 < 0) {
						RLI <- 1
					} else {
						RLI <- t1/t2
					}
				} else {
				   RLI <- 1
				}
				result["rfi.scaled"] <- RLI
			}
			
			if("nfi.scaled" %in% names(result)) {
				t1 <- X2.null.scaled - X2.scaled
				t2 <- X2.null.scaled
				NFI <- t1/t2
				result["nfi.scaled"] <- NFI
			}
			
			if("pnfi.scaled" %in% names(result)) {
				t1 <- X2.null.scaled - X2.scaled
				t2 <- X2.null.scaled
				PNFI <- (df/df.null) * t1/t2
				result["pnfi.scaled"] <- PNFI
			}
			
			if("ifi.scaled" %in% names(result)) {
				t1 <- X2.null.scaled - X2.scaled
				t2 <- X2.null.scaled
				if(t2 < 0) {
					IFI <- 1
				} else {
					IFI <- t1/t2
				}
				result["ifi.scaled"] <- IFI
			}
			
			if("rni.scaled" %in% names(result)) {
				t1 <- X2.scaled - df.scaled
				t2 <- X2.null.scaled - df.null.scaled
				t2 <- X2.null - df.null
				if(t1 < 0 || t2 < 0) {
					RNI <- 1
				} else {
					RNI <- 1 - t1/t2
				}
				result["rni.scaled"] <- RNI
			}
		}
	} 
	
	#logl
	imputed <- object@imputed
	if(length(imputed) > 0) {
		loglikval <- unlist(imputed[["logl"]])
		npar <- result["npar"]
		result["unrestricted.logl"] <- loglikval["unrestricted.logl"]
		result["logl"] <- loglikval["logl"]
		result["aic"] <-  -2*loglikval["logl"] + 2*npar
        result["bic"] <- -2*loglikval["logl"] + npar*log(result["ntotal"])
		N.star <- (result["ntotal"] + 2) / 24
		result["bic2"] <- -2*loglikval["logl"] + npar*log(N.star)
		result <- result[-which("fmin" == names(result))]
	}
	result	
}


