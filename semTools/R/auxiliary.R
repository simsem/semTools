### Terrence D. Jorgensen
### Last updated: 8 March 2018
### new auxiliary function does NOT create a lavaanStar-class instance

#' Implement Saturated Correlates with FIML
#'
#' Automatically add auxiliary variables to a lavaan model when using full
#' information maximum likelihood (FIML) to handle missing data
#'
#' These functions are wrappers around the corresponding lavaan functions.
#' You can use them the same way you use \code{\link[lavaan]{lavaan}}, but you
#' \emph{must} pass your full \code{data.frame} to the \code{data} argument.
#' Because the saturated-correlates approaches (Enders, 2008) treates exogenous
#' variables as random, \code{fixed.x} must be set to \code{FALSE}. Because FIML
#' requires continuous data (although nonnormality corrections can still be
#' requested), no variables in the model nor auxiliary variables specified in
#' \code{aux} can be declared as \code{ordered}.
#'
#' @aliases auxiliary lavaan.auxiliary cfa.auxiliary sem.auxiliary growth.auxiliary
#' @importFrom lavaan lavInspect parTable
#' @importFrom stats cov quantile
#'
#' @param model The analysis model can be specified with 1 of 2 objects:
#'   \enumerate{
#'     \item  lavaan \code{\link[lavaan]{model.syntax}} specifying a hypothesized
#'            model \emph{without} mention of auxiliary variables in \code{aux}
#'     \item  a parameter table, as returned by \code{\link[lavaan]{parTable}},
#'            specifying the target model \emph{without} auxiliary variables.
#'            This option requires these columns (and silently ignores all others):
#'            \code{c("lhs","op","rhs","user","group","free","label","plabel","start")}
#'   }
#' @param data \code{data.frame} that includes auxiliary variables as well as
#'   any observed variables in the \code{model}
#' @param aux \code{character}. Names of auxiliary variables to add to \code{model}
#' @param fun \code{character}. Name of a specific lavaan function used to fit
#'   \code{model} to \code{data} (i.e., \code{"lavaan"}, \code{"cfa"},
#'   \code{"sem"}, or \code{"growth"}). Only required for \code{auxiliary}.
#' @param ... additional arguments to pass to \code{\link[lavaan]{lavaan}}.
#'
#' @author
#' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
#'
#' @references Enders, C. K. (2008). A note on the use of missing auxiliary
#'   variables in full information maximum likelihood-based structural equation
#'   models. \emph{Structural Equation Modeling, 15}(3), 434-448.
#'   doi:10.1080/10705510802154307
#'
#' @return a fitted \code{\linkS4class{lavaan}} object.  Additional
#'   information is stored as a \code{list} in the \code{\@external} slot:
#'   \itemize{
#'     \item \code{baseline.model}. a fitted \code{\linkS4class{lavaan}}
#'           object. Results of fitting an appropriate independence model for
#'           the calculation of incremental fit indices (e.g., CFI, TLI) in
#'           which the auxiliary variables remain saturated, so only the target
#'           variables are constrained to be orthogonal. See Examples for how
#'           to send this baseline model to \code{\link[lavaan]{fitMeasures}}.
#'     \item \code{aux}. The character vector of auxiliary variable names.
#'   }
#'
#' @examples
#' dat1 <- lavaan::HolzingerSwineford1939
#' set.seed(12345)
#' dat1$z <- rnorm(nrow(dat1))
#' dat1$x5 <- ifelse(dat1$z < quantile(dat1$z, .3), NA, dat1$x5)
#' dat1$x9 <- ifelse(dat1$z > quantile(dat1$z, .8), NA, dat1$x9)
#'
#' targetModel <- "
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#' "
#'
#' ## works just like cfa(), but with an extra "aux" argument
#' fitaux1 <- cfa.auxiliary(targetModel, data = dat1, aux = "z",
#'                          missing = "fiml", estimator = "mlr")
#'
#' ## with multiple auxiliary variables and multiple groups
#' fitaux2 <- cfa.auxiliary(targetModel, data = dat1, aux = c("z","ageyr","grade"),
#'                          group = "school", group.equal = "loadings")
#'
#' ## calculate correct incremental fit indices (e.g., CFI, TLI) using
#' ## the internally stored baseline model
#' fitMeasures(fitaux2, fit.measures = c("cfi","tli"),
#'             baseline.model = fitaux2@external$baseline.model)
#' @export
auxiliary <- function(model, data, aux, fun, ...) {
  lavArgs <- list(...)
  ## check for constraints on observed variables that affect auxiliaries
  if (!is.null(lavArgs$group.equal)) {
    if (length(intersect(lavArgs$group.equal,
                         c("intercepts","residuals","residual.covariances")))) {
      warning('Using the "group.equal" argument to constrain intercepts will',
              ' constrain all observed-variable intercepts, including',
              ' auxiliary-variables means. Likewise, using the "group.equal"',
              ' argument to constrain residual (co)variances will apply ',
              ' those constraints to all observed variables, including',
              ' exogenous and auxiliary variables. In order not to fit an',
              ' overly constrained model, it is safer to apply equality',
              ' constraints across groups manually using labels in the syntax.')
    }
  }

  if (missing(aux))
    stop("Please provide a character vector with names of auxiliary variables")
  if (missing(data))
    stop("Please provide a data.frame that includes modeled and auxiliary variables")
  if (!all(sapply(data[aux], function(x) class(x[1]) %in% c("numeric","integer"))))
    stop("missing = 'FIML' is unavailable for categorical data")

  ## Easiest scenario: model is a character string
  if (is.character(model)) {
    varnames <- lavaan::lavNames(lavaan::lavaanify(model), type = "ov")
    if (length(intersect(varnames, aux))) stop('modeled variable declared as auxiliary')
    ## concatenate saturated auxiliaries
    covstruc <- outer(aux, aux, function(x, y) paste(x, "~~", y))
    lavArgs$model <- c(model, paste(aux, "~ 1"),
                       covstruc[lower.tri(covstruc, diag = TRUE)],
                       outer(aux, varnames, function(x, y) paste(x, "~~", y)))
    lavArgs$data <- data
    lavArgs$fixed.x <- FALSE
    lavArgs$missing <- "fiml"
    lavArgs$meanstructure <- TRUE
    lavArgs$ordered <- NULL
#FIXME   if (lavArgs$group.equal != "")
    result <- do.call(fun, lavArgs)
    ## specify, fit, and attach an appropriate independence model
    baseArgs <- list()
    baseArgs$data                <- data
    baseArgs$missing             <- "fiml"
    baseArgs$group               <- lavArgs$group
    baseArgs$group.label         <- lavArgs$group.label
    baseArgs$cluster             <- lavArgs$cluster
    baseArgs$sample.cov.rescale  <- lavArgs$sample.cov.rescale
    baseArgs$information         <- lavArgs$information
    baseArgs$se                  <- lavArgs$se
    baseArgs$test                <- lavArgs$test
    baseArgs$bootstrap           <- lavArgs$bootstrap
    baseArgs$control             <- lavArgs$control
    baseArgs$optim.method        <- lavArgs$optim.method
    baseArgs$model <- c(paste(varnames, "~~", varnames),
                        paste(varnames, "~ 1"), paste(aux, "~ 1"),
                        covstruc[lower.tri(covstruc, diag = TRUE)],
                        outer(aux, varnames, function(x, y) paste(x, "~~", y)))
    result@external$baseline.model <- do.call(lavaan::lavaan, baseArgs)
    result@external$aux <- aux
    return(result)
  }

  ## otherwise, only accept a parameter table
  PTcols <- c("lhs","op","rhs","user","block","group","free","label","plabel","start")

	if (is.list(model)) {
	  if (any(model$exo == 1))
	    stop("All exogenous variables (covariates) must be treated as endogenous",
	         " by the 'auxiliary' function. Please set 'fixed.x = FALSE'")

	  if (is.null(model$start)) {
	    startArgs <- lavArgs
	    startArgs$model <- model
	    startArgs$data <- data
	    startArgs$fixed.x <- FALSE
	    startArgs$missing <- "fiml"
	    startArgs$meanstructure <- TRUE
	    startArgs$do.fit <- FALSE
	    model$start <- parTable(do.call(fun, startArgs))$start
	  }

	  missingCols <- setdiff(PTcols, names(model))
	  if (length(missingCols)) stop("If the 'model' argument is a parameter table",
	                                " it must also include these columns: \n",
	                                paste(missingCols, collapse = ", "))
	  PT <- as.data.frame(model, stringsAsFactors = FALSE)[PTcols]
	} else stop("The 'model' argument must be a character vector of",
	            " lavaan syntax or a parameter table")

  ## separately store rows with constraints or user-defined parameters
  conRows <- PT$op %in% c("==","<",">",":=")
  if (any(conRows)) {
    CON <- PT[  conRows, ]
    PT <-  PT[ !conRows, ]
  } else CON <- data.frame(NULL)

  ## variable names
  varnames <- lavaan::lavNames(PT, type = "ov")
	if (length(intersect(varnames, aux))) stop('modeled variable declared as auxiliary')

	## specify a saturated model for auxiliaries
	covstruc <- outer(aux, aux, function(x, y) paste(x, "~~", y))
	satMod <- c(covstruc[lower.tri(covstruc, diag = TRUE)], paste(aux, "~ 1"), # among auxiliaries
	            outer(aux, varnames, function(x, y) paste(x, "~~", y)))        # between aux and targets
	satPT <- lavaan::lavaanify(satMod,
	                           ngroups = max(PT$group))[c("lhs","op","rhs","user","block","group")]

	## after omitting duplicates, check number of added parameters, add columns
	mergedPT <- lavaan::lav_partable_merge(PT, satPT, remove.duplicated = TRUE, warn = FALSE)
	nAuxPar <- nrow(mergedPT) - nrow(PT)
	newRows <- 1L:nAuxPar + nrow(PT)
	##FIXME:  mergedPT$user[newRows] <- 2L (list as constraints to omit printing?) or new code (9L)?
	mergedPT$free[newRows] <- 1L:nAuxPar + max(PT$free)
	mergedPT$plabel[newRows] <- paste0(".p", 1L:nAuxPar + nrow(PT), ".")
	## calculate sample moments as starting values (recycle over groups)
	# if (is.null(lavArgs$group)) {
	#   auxCov <- cov(data[aux], use = "pairwise.complete.obs")
	#   auxM <- colMeans(data[aux], na.rm = TRUE)
	#   auxTarget <- cov(data[c(aux, varnames)],
	#                    use = "pairwise.complete.obs")[aux, varnames]
	#   ## match order of parameters in syntax above
	#   mergedPT$start[newRows] <- c(auxCov[lower.tri(auxCov, diag = TRUE)],
	#                                auxM, as.numeric(auxTarget))
	# } else {
	#   auxCovs <- list()
	#   auxMs <- list()
	#   auxTargets <- list()
	#   startVals <- numeric(0)
	#   for (g in unique(data[ , lavArgs$group])) {
	#     auxCovs[[g]] <- cov(data[data[ , lavArgs$group] == g, aux],
	#                         use = "pairwise.complete.obs")
	#     auxMs[[g]] <- colMeans(data[data[ , lavArgs$group] == g, aux], na.rm = TRUE)
	#     auxTargets[[g]] <- cov(data[data[ , lavArgs$group] == g, c(aux, varnames)],
	#                            use = "pairwise.complete.obs")[aux, varnames]
	#     startVals <- c(startVals, auxCovs[[g]][lower.tri(auxCovs[[g]], diag = TRUE)],
	#                    auxMs[[g]], as.numeric(auxTargets[[g]]))
	#   }
	#   ## match order of parameters in syntax above
	#   mergedPT$start[newRows] <- startVals
	# }
	finalPT <- lavaan::lav_partable_complete(rbind(mergedPT, CON))

  lavArgs$model <- finalPT
  lavArgs$data <- data
  lavArgs$fixed.x <- FALSE
  lavArgs$missing <- "fiml"
  lavArgs$meanstructure <- TRUE
  lavArgs$ordered <- NULL
  result <- do.call(fun, lavArgs)

  ## specify, fit, and attach an appropriate independence model
  baseArgs <- list()
  baseArgs$model                 <- lavaan::lav_partable_complete(satPT)
  baseArgs$data                  <- data
  baseArgs$group                 <- lavArgs$group
  baseArgs$group.label           <- lavArgs$group.label
  baseArgs$missing               <- "fiml"
  baseArgs$cluster               <- lavArgs$cluster
  baseArgs$sample.cov.rescale    <- lavArgs$sample.cov.rescale
  baseArgs$information           <- lavArgs$information
  baseArgs$se                    <- lavArgs$se
  baseArgs$test                  <- lavArgs$test
  baseArgs$bootstrap             <- lavArgs$bootstrap
  baseArgs$control               <- lavArgs$control
  baseArgs$optim.method          <- lavArgs$optim.method
  result@external$baseline.model <- do.call(lavaan::lavaan, baseArgs)
  result@external$aux <- aux

	result
}

#' @rdname auxiliary
#' @aliases lavaan.auxiliary
#' @export
lavaan.auxiliary <- function(model, data, aux, ...) {
	auxiliary(model = model, data = data, aux = aux, fun = "lavaan", ...)
}

#' @rdname auxiliary
#' @aliases cfa.auxiliary
#' @export
cfa.auxiliary <- function(model, data, aux, ...) {
	auxiliary(model = model, data = data, aux = aux, fun = "cfa", ...)
}

#' @rdname auxiliary
#' @aliases sem.auxiliary
#' @export
sem.auxiliary <- function(model, data, aux, ...) {
	auxiliary(model = model, data = data, aux = aux, fun = "sem", ...)
}

#' @rdname auxiliary
#' @aliases growth.auxiliary
#' @export
growth.auxiliary <- function(model, data, aux, ...) {
	auxiliary(model = model, data = data, aux = aux, fun = "growth", ...)
}



#################################################################
## After enough time passes, delete everything below this line ##
#################################################################

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
	notused <- utils::capture.output(result <- suppressWarnings(getMethod("inspect", "lavaan")(object, what="fit"))) ## FIXME: don't set a new inspect method
	result[c("baseline.chisq", "baseline.df", "baseline.pvalue")] <- object@nullfit[c("chisq", "df", "pvalue")]

    if(lavInspect(object, "options")$test %in% c("satorra.bentler", "yuan.bentler",
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


