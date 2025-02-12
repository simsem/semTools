### Terrence D. Jorgensen
### Last updated: 12 February 2025

##' Implement Saturated Correlates with FIML
##'
##' Automatically add auxiliary variables to a lavaan model when using full
##' information maximum likelihood (FIML) to handle missing data
##'
##' These functions are wrappers around the corresponding lavaan functions.
##' You can use them the same way you use [lavaan::lavaan()], but you
##' *must* pass your full `data.frame` to the `data` argument.
##' Because the saturated-correlates approaches (Enders, 2008) treats exogenous
##' variables as random, `fixed.x` must be set to `FALSE`. Because FIML
##' requires continuous data (although nonnormality corrections can still be
##' requested), no variables in the model nor auxiliary variables specified in
##' `aux` can be declared as `ordered`.
##'
##' @aliases auxiliary lavaan.auxiliary cfa.auxiliary sem.auxiliary growth.auxiliary
##' @importFrom lavaan lavInspect parTable lavNames
##' @importFrom stats cov quantile
##'
##' @param model The analysis model can be specified with 1 of 2 objects:
##'   \enumerate{
##'     \item  lavaan [lavaan::model.syntax()] specifying a hypothesized
##'            model *without* mention of auxiliary variables in `aux`
##'     \item  a parameter table, as returned by [lavaan::parTable()],
##'            specifying the target model *without* auxiliary variables.
##'            This option requires these columns (and silently ignores all others):
##'            `c("lhs","op","rhs","user","group","free","label","plabel","start")`
##'   }
##' @param data `data.frame` that includes auxiliary variables as well as
##'   any observed variables in the `model`
##' @param aux `character`. Names of auxiliary variables to add to `model`
##' @param fun `character`. Name of a specific lavaan function used to fit
##'   `model` to `data` (i.e., `"lavaan"`, `"cfa"`,
##'   `"sem"`, or `"growth"`). Only required for `auxiliary`.
##' @param ... Additional arguments to pass to `fun=`.
##' @param envir Passed to [do.call()].
##' @param return.syntax `logical` indicating whether to return a
##'   `character` string of [lavaan::model.syntax()] that can be
##'   added to a target `model=` that is also a `character` string.
##'   This can be advantageous, for example, to use add saturated correlates to
##'   a \pkg{blavaan} model.
##'
##' @author
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Enders, C. K. (2008). A note on the use of missing auxiliary variables in
##'   full information maximum likelihood-based structural equation models.
##'   *Structural Equation Modeling, 15*(3), 434--448.
##'   \doi{10.1080/10705510802154307}
##'
##' @return a fitted [lavaan::lavaan-class] object.  Additional
##'   information is stored as a `list` in the `@@external` slot:
##'   \itemize{
##'     \item `baseline.model`. a fitted [lavaan::lavaan-class]
##'           object. Results of fitting an appropriate independence model for
##'           the calculation of incremental fit indices (e.g., CFI, TLI) in
##'           which the auxiliary variables remain saturated, so only the target
##'           variables are constrained to be orthogonal. See Examples for how
##'           to send this baseline model to [lavaan::fitMeasures()].
##'     \item `aux`. The character vector of auxiliary variable names.
##'     \item `baseline.syntax`. A character vector generated within the
##'           `auxiliary` function, specifying the `baseline.model`
##'           syntax.
##'   }
##'
##' @examples
##' dat1 <- lavaan::HolzingerSwineford1939
##' set.seed(12345)
##' dat1$z <- rnorm(nrow(dat1))
##' dat1$x5 <- ifelse(dat1$z < quantile(dat1$z, .3), NA, dat1$x5)
##' dat1$x9 <- ifelse(dat1$z > quantile(dat1$z, .8), NA, dat1$x9)
##'
##' targetModel <- "
##'   visual  =~ x1 + x2 + x3
##'   textual =~ x4 + x5 + x6
##'   speed   =~ x7 + x8 + x9
##' "
##'
##' ## works just like cfa(), but with an extra "aux" argument
##' fitaux1 <- cfa.auxiliary(targetModel, data = dat1, aux = "z",
##'                          missing = "fiml", estimator = "mlr")
##'
##' ## with multiple auxiliary variables and multiple groups
##' fitaux2 <- cfa.auxiliary(targetModel, data = dat1, aux = c("z","ageyr","grade"),
##'                          group = "school", group.equal = "loadings")
##'
##' ## calculate correct incremental fit indices (e.g., CFI, TLI)
##' fitMeasures(fitaux2, fit.measures = c("cfi","tli"))
##' ## NOTE: lavaan will use the internally stored baseline model, which
##' ##       is the independence model plus saturated auxiliary parameters
##' lavInspect(fitaux2@external$baseline.model, "free")
##'
##'
##' @export
auxiliary <- function(model, data, aux, fun, ...,
                      envir = getNamespace("lavaan"), return.syntax = FALSE) {
  lavArgs <- list(...)
  lavArgs$data <- substitute(data)
  lavArgs$fixed.x <- FALSE
  if (fun %in% c("lavaan","cfa","sem","growth")) {
    #TODO: allow fun = "efa"?
    lavArgs$missing <- "fiml"
    lavArgs$meanstructure <- TRUE
  }
  lavArgs$ordered <- NULL

  if (missing(aux))
    stop("Please provide a character vector with names of auxiliary variables")
  if (missing(data))
    stop("Please provide a data.frame that includes modeled and auxiliary variables")
  if (fun %in% c("lavaan","cfa","sem","growth") && !all(sapply(data[aux], is.numeric)))
    stop("missing = 'FIML' is unavailable for categorical data")


  PTcols <- c("lhs","op","rhs","user","block","group","free","label","plabel","start")
  #TODO: add any? if fun %in% c("blavaan","bcfa","bsem","bgrowth")
  ## check parameter table, or create one from syntax
  if (is.list(model)) {
	  if (any(model$exo == 1))
	    stop("All exogenous variables (covariates) must be treated as endogenous",
	         " by the 'auxiliary' function. Please set 'fixed.x = FALSE'")

    if (!is.null(lavArgs$group.equal))
      warning("The 'group.equal' argument is ignored when 'model' is a parameter table.")

	  if (is.null(model$start)) {
	    startArgs <- lavArgs
	    startArgs$model <- model
	    startArgs$do.fit <- FALSE
	    model$start <- parTable(do.call(fun, startArgs, envir = envir))$start
	  }

	  missingCols <- setdiff(PTcols, names(model))
	  if (length(missingCols)) stop("If the 'model' argument is a parameter table",
	                                " it must also include these columns: \n",
	                                paste(missingCols, collapse = ", "))
	  PT <- as.data.frame(model, stringsAsFactors = FALSE)[PTcols]
	} else if (is.character(model)) {
	  ptArgs <- lavArgs
	  ptArgs$model <- model
	  ptArgs$do.fit <- FALSE #FIXME: does this work with blavaan?
	  PT <- parTable(do.call(fun, ptArgs, envir = envir))[PTcols]
	} else stop("The 'model' argument must be a character vector of",
	            " lavaan syntax or a parameter table")


  ## separately store rows with constraints or user-defined parameters
  conRows <- PT$op %in% c("==","<",">",":=")
  if (any(conRows)) {
    CON <- PT[  conRows, ]
    PT <-  PT[ !conRows, ]
  } else CON <- data.frame(NULL)

  ## variable names
  varnames <- lavNames(PT, type = "ov")
	if (length(intersect(varnames, aux))) stop('modeled variable declared as auxiliary')

	## specify a saturated model for auxiliaries
	covstruc <- outer(aux, aux, function(x, y) paste(x, "~~", y))
	satMod <- c(covstruc[lower.tri(covstruc, diag = TRUE)], paste(aux, "~ 1"), # among auxiliaries
	            outer(aux, varnames, function(x, y) paste(x, "~~", y)))        # between aux and targets
	if (return.syntax) return(satMod)
	satPT <- lavaan::lavaanify(satMod, ngroups = max(PT$group))[c("lhs","op","rhs",
	                                                              "user","block","group")]

	## after omitting duplicates, check number of added parameters, add columns
	mergedPT <- lavaan::lav_partable_merge(PT, satPT, remove.duplicated = TRUE, warn = FALSE)
	nAuxPar <- nrow(mergedPT) - nrow(PT)
	newRows <- 1L:nAuxPar + nrow(PT)
	#TODO:  mergedPT$user[newRows] <- 2L (list as constraints to omit printing?) or new code (9L)?
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
	lavArgs$model <- lavaan::lav_partable_complete(rbind(mergedPT, CON)) #FIXME for blavaan
  result <- do.call(fun, lavArgs, envir = envir)

  if (fun %in% c("lavaan","cfa","sem","growth")) {
    ## specify, fit, and attach an appropriate independence model
    baseArgs <- list()
    baseArgs$model                  <- lavaan::lav_partable_complete(satPT)
    baseArgs$data                   <- data
    baseArgs$group                  <- lavArgs$group
    baseArgs$group.label            <- lavArgs$group.label
    baseArgs$missing                <- "fiml"
    baseArgs$cluster                <- lavArgs$cluster
    baseArgs$sample.cov.rescale     <- lavArgs$sample.cov.rescale
    baseArgs$estimator              <- lavArgs$estimator
    baseArgs$information            <- lavArgs$information
    baseArgs$se                     <- lavArgs$se
    baseArgs$test                   <- lavArgs$test
    baseArgs$bootstrap              <- lavArgs$bootstrap
    baseArgs$control                <- lavArgs$control
    baseArgs$optim.method           <- lavArgs$optim.method

    result@external$baseline.model  <- do.call("lavaan", baseArgs, envir = envir)
    result@external$aux             <- aux
    result@external$baseline.syntax <- satMod
  } # else something similar for blavaan?  Better to fit manually for blavFitIndices

	result
}

##' @rdname auxiliary
##' @aliases lavaan.auxiliary
##' @export
lavaan.auxiliary <- function(model, data, aux, ...,
                             envir = getNamespace("lavaan")) {
  mc <- match.call(expand.dots = TRUE)
  mc$fun <- "lavaan"
  mc[[1L]] <- quote(semTools::auxiliary)
  eval(mc, parent.frame())
}

##' @rdname auxiliary
##' @aliases cfa.auxiliary
##' @export
cfa.auxiliary <- function(model, data, aux, ...,
                          envir = getNamespace("lavaan")) {
  mc <- match.call(expand.dots = TRUE)
  mc$fun <- "cfa"
  mc[[1L]] <- quote(semTools::auxiliary)
  eval(mc, parent.frame())
}

##' @rdname auxiliary
##' @aliases sem.auxiliary
##' @export
sem.auxiliary <- function(model, data, aux, ...,
                          envir = getNamespace("lavaan")) {
  mc <- match.call(expand.dots = TRUE)
  mc$fun <- "sem"
  mc[[1L]] <- quote(semTools::auxiliary)
  eval(mc, parent.frame())
}

##' @rdname auxiliary
##' @aliases growth.auxiliary
##' @export
growth.auxiliary <- function(model, data, aux, ...,
                             envir = getNamespace("lavaan")) {
  mc <- match.call(expand.dots = TRUE)
  mc$fun <- "growth"
  mc[[1L]] <- quote(semTools::auxiliary)
  eval(mc, parent.frame())
}


