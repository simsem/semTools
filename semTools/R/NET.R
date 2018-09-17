### Terrence D. Jorgensen
### Last updated: 17 September 2018
### semTools functions for Nesting and Equivalence Testing


## -----------------
## Class and Methods
## -----------------

##' Class For the Result of Nesting and Equivalence Testing
##'
##' This class contains the results of nesting and equivalence testing among
##' multiple models
##'
##'
##' @name Net-class
##' @aliases Net-class show,Net-method summary,Net-method
##' @docType class
##'
##' @slot test Logical \code{matrix} indicating nesting/equivalence among models
##' @slot df The degrees of freedom of tested models
##'
##' @section Objects from the Class: Objects can be created via the
##' \code{\link{net}} function.
##'
##' @param object An object of class \code{Net}.
##'
##' @return
##' \item{show}{\code{signature(object = "Net")}: prints the logical matrix of
##'   test results.}
##' \item{summary}{\code{signature(object = "Net")}: prints a narrative
##'   description of results. The original \code{object} is invisibly returned.}
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @seealso \code{\link{net}}
##'
##' @examples
##'
##' # See the example in the net function.
##'
setClass("Net", representation(test = "matrix", df = "vector"))


##' @rdname Net-class
##' @aliases show,Net-method
##' @export
setMethod("show", "Net",
function(object) {
  if (length(object@test)) {
    m <- as.matrix(unclass(object@test))
    m[upper.tri(m, diag = TRUE)] <- ""
    cat("
        If cell [R, C] is TRUE, the model in row R is nested within column C.

        If the models also have the same degrees of freedom, they are equivalent.

        NA indicates the model in column C did not converge when fit to the
        implied means and covariance matrix from the model in row R.

        The hidden diagonal is TRUE because any model is equivalent to itself.
        The upper triangle is hidden because for models with the same degrees
        of freedom, cell [C, R] == cell [R, C].  For all models with different
        degrees of freedom, the upper diagonal is all FALSE because models with
        fewer degrees of freedom (i.e., more parameters) cannot be nested
        within models with more degrees of freedom (i.e., fewer parameters).
        \n")
    print(m, quote = FALSE)
  } else {
    cat(data.class(object@test), "(0)\n", sep = "")
  }
  invisible(object)
})


##' @rdname Net-class
##' @aliases summary,Net-method
##' @export
setMethod("summary", "Net",
function(object) {
  DFs <- object@df
  x <- object@test
  mods <- colnames(x)
  for (R in 2:nrow(x)) {
    for (C in (R - 1):1) {
      ## if model didn't converge (logical value is missing), go to next iteration
      if (is.na(x[R, C])) next
      ## if the models are not nested, go to next iteration
      if (!x[R, C]) next
      ## choose message based on whether models are equivalent or nested
      if (identical(DFs[R], DFs[C])) {
        rel <- "equivalent to"
      } else {
        rel <- "nested within"
      }
      cat("Model \"", mods[R], "\" is ", rel, " model \"", mods[C], "\"\n", sep = "")
    }
  }
  invisible(object)
})



## --------------------
## Constructor Function
## --------------------

##' Nesting and Equivalence Testing
##'
##' This test examines whether models are nested or equivalent based on Bentler
##' and Satorra's (2010) procedure.
##'
##' The concept of nesting/equivalence should be the same regardless of
##' estimation method. However, the particular method of testing
##' nesting/equivalence (as described in Bentler & Satorra, 2010) employed by
##' the net function analyzes summary statistics (model-implied means and
##' covariance matrices, not raw data). In the case of robust methods like MLR,
##' the raw data is only utilized for the robust adjustment to SE and chi-sq,
##' and the net function only checks the unadjusted chi-sq for the purposes of
##' testing nesting/equivalence.  This method does not apply to models that
##' estimate thresholds for categorical data, so an error message will be issued
##' if such a model is provided.
##'
##'
##' @importFrom lavaan lavInspect
##'
##' @param \dots The \code{lavaan} objects used for test of nesting and
##'   equivalence
##' @param crit The upper-bound criterion for testing the equivalence of models.
##'   Models are considered nested (or equivalent) if the difference between
##'   their \eqn{\chi^2} fit statistics is less than this criterion.
##'
##' @return The \linkS4class{Net} object representing the outputs for nesting
##'   and equivalent testing, including a logical matrix of test results and a
##'   vector of degrees of freedom for each model.
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##' Bentler, P. M., & Satorra, A. (2010). Testing model nesting and equivalence.
##' \emph{Psychological Methods, 15}(2), 111--123. doi:10.1037/a0019625
##'
##' @examples
##'
##' \dontrun{
##' m1 <- ' visual  =~ x1 + x2 + x3
##' 	       textual =~ x4 + x5 + x6
##' 	       speed   =~ x7 + x8 + x9 '
##'
##'
##' m2 <- ' f1  =~ x1 + x2 + x3 + x4
##' 	       f2 =~ x5 + x6 + x7 + x8 + x9 '
##'
##' m3 <- ' visual  =~ x1 + x2 + x3
##' 	       textual =~ eq*x4 + eq*x5 + eq*x6
##' 	       speed   =~ x7 + x8 + x9 '
##'
##' fit1 <- cfa(m1, data = HolzingerSwineford1939)
##' fit1a <- cfa(m1, data = HolzingerSwineford1939, std.lv = TRUE) # Equivalent to fit1
##' fit2 <- cfa(m2, data = HolzingerSwineford1939) # Not equivalent to or nested in fit1
##' fit3 <- cfa(m3, data = HolzingerSwineford1939) # Nested in fit1 and fit1a
##'
##' tests <- net(fit1, fit1a, fit2, fit3)
##' tests
##' summary(tests)
##' }
##'
##' @export
net <- function(..., crit = .0001) {
  ## put fitted objects in a list
  fitList <- list(...)
  nFits <- length(fitList)

  ## check that they are all lavaan objects
  notLavaan <- sapply(fitList, class) != "lavaan"
  if (any(notLavaan)) {
    fitNames <- sapply(as.list(substitute(list(...)))[-1], deparse)
    stop(paste("The following arguments are not fitted lavaan objects:\n",
               paste(fitNames[notLavaan], collapse = "\t")))
  }

  ## check for meanstructure
  meanstructure <- sapply(fitList, function(x) lavInspect(x, "options")$meanstructure)
  if (!(all(meanstructure) || !any(meanstructure)))
    stop('Some (but not all) fitted lavaan objects include a mean structure. ',
         'Please re-fit all models with the argument meanstructure=TRUE.')

  ## check whether any models include categorical outcomes
  catMod <- sapply(fitList, function(x) lavInspect(x, "options")$categorical)
  if (any(catMod)) stop("This method only applies to continuous outcomes.")

  ## get degrees of freedom for each model
  DFs <- sapply(fitList, function(x) lavInspect(x, "fit")["df"])

  ## name according to named objects, with DF in parentheses
  fitNames <- names(fitList)
  dotNames <- sapply(as.list(substitute(list(...)))[-1], deparse)
  if (is.null(names(fitList))) {
    fitNames <- dotNames
  } else {
    noName <- which(fitNames == "")
    fitNames[noName] <- dotNames[noName]
  }
  names(fitList) <- paste(fitNames, " (df = ", DFs, ")", sep = "")

  ## sort list according to DFs
  fitList <- fitList[order(DFs)]
  fitNames <- fitNames[order(DFs)]
  orderedDFs <- DFs[order(DFs)]

  ## create structure for sequence of tests (logical matrix), FALSE by default
  nestMat <- matrix(FALSE, nFits, nFits, dimnames = list(names(fitList), fitNames))
  diag(nestMat) <- TRUE # every model is equivalent with itself

  ## Loop through sorted models in sequence of most to least restricted model
  for (R in 2:nrow(nestMat)) {
    for (C in (R - 1):1) {
      ## test for nesting/equivalence
      nestMat[R, C] <- x.within.y(x = fitList[[R]], y = fitList[[C]], crit = crit)
      ## if models are equivalent, set above-diagonal value to TRUE
      if (identical(orderedDFs[R], orderedDFs[C])) nestMat[C, R] <- nestMat[R, C]
      if (C == 1) next # to prevent the next 2 tests from returning an error
      ## if model didn't converge (logical value is missing), go to next iteration
      if (is.na(nestMat[R, C]) | is.na(nestMat[R - 1, C - 1])) next
      ## check whether nesting is implied, to skip unnecessary tests
      if (nestMat[R, C] & nestMat[R - 1, C - 1]) {
        nestMat[R, C - 1] <- TRUE
        next
      }
    }
  }
  out <- new("Net", test = nestMat, df = orderedDFs)
  out
}


## --------------------------------------------------------------------
## Hidden Function to test whether model "x" is nested within model "y"
## --------------------------------------------------------------------

#' @importFrom lavaan lavInspect
x.within.y <- function(x, y, crit = .0001) {
  if (length(c(lavaan::lavNames(x, "ov.ord"), lavaan::lavNames(y, "ov.ord"))))
    stop("The net() function is not available for categorical-data estimators.")

  exoX <- lavInspect(x, "options")$fixed.x & length(lavaan::lavNames(x, "ov.x"))
  exoY <- lavInspect(y, "options")$fixed.x & length(lavaan::lavNames(y, "ov.x"))
  if (exoX | exoY) {
    stop(c("The net() function does not work with exogenous variables.\n",
           "Fit the model again with 'fixed.x = FALSE'"))
  }
  ## variable names
  Xnames <- lavaan::lavNames(x)
  Ynames <- lavaan::lavNames(y)
  if (!identical(sort(Xnames), sort(Ynames)))
    stop("Models do not contain the same variables")

  ## check that the analyzed data matches
  xData <- lavInspect(x, "data")
  if (is.list(xData)) xData <- do.call(rbind, xData)
  xData <- xData[ , order(Xnames)]
  yData <- lavInspect(y, "data")
  if (is.list(yData)) yData <- do.call(rbind, yData)
  yData <- yData[ , order(Ynames)]
  if (!identical(xData, yData)) stop("Models must apply to the same data")

  ## check degrees of freedom support nesting structure
  if (lavInspect(x, "fit")["df"] < lavInspect(y, "fit")["df"])
    stop("x cannot be nested within y because y is more restricted than x")

  ## model-implied moments
  Sigma <- lavInspect(x, "cov.ov")
  if (lavInspect(x, "options")$meanstructure) {
    Mu <- lavInspect(x, "mean.ov")
  } else {
    Mu <- NULL
  }
  N <- lavInspect(x, "nobs")

  ## fit model and check that chi-squared < crit

  suppressWarnings(try(newFit <- lavaan::update(y, data = NULL,
                                                sample.cov = Sigma,
                                                sample.mean = Mu,
                                                sample.nobs = N,
                                                estimator = "ML",
                                                se = "none", # to save time
                                                test = "standard")))
  if(!lavInspect(newFit, "converged")) return(NA) else {
    result <- lavInspect(newFit, "fit")["chisq"] < crit
    names(result) <- NULL
    if (lavInspect(x, "fit")["df"] ==
        lavInspect(y, "fit")["df"]) return(c(Equivalent = result))
  }
  c(Nested = result)
}


