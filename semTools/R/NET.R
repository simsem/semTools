### Terrence D. Jorgensen
### Last updated: 14 October 2016
### semTools function for Nesting and Equivalence Testing

setClass("Net", representation(test = "matrix", df = "vector"))


## function to test whether model "x" is nested within model "y"
x.within.y <- function(x, y, crit = crit) {
  if (length(c(lavaan::lavNames(x, "ov.ord"), lavaan::lavNames(y, "ov.ord"))))
    stop("The net() function is not available for categorical-data estimators.")

  exoX <- lavaan::lavInspect(x, "options")$fixed.x & length(lavaan::lavNames(x, "ov.x"))
  exoY <- lavaan::lavInspect(y, "options")$fixed.x & length(lavaan::lavNames(y, "ov.x"))
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
  xData <- lavaan::lavInspect(x, "data")
  if (is.list(xData)) xData <- do.call(rbind, xData)
  xData <- xData[ , rank(Xnames)]
  yData <- lavaan::lavInspect(y, "data")
  if (is.list(yData)) yData <- do.call(rbind, yData)
  yData <- yData[ , rank(Ynames)]
  if (!identical(xData, yData)) stop("Models must apply to the same data")
  ##############################################################################

  ## check degrees of freedom support nesting structure
  if (lavaan::lavInspect(x, "fit")["df"] < lavaan::lavInspect(y, "fit")["df"])
    stop("x cannot be nested within y because y is more restricted than x")

  ## model-implied moments
  Sigma <- lavaan::lavInspect(x, "cov.ov")
  Mu <- lavaan::lavInspect(x, "mean.ov")
  N <- lavaan::lavInspect(x, "nobs")

  ## fit model and check that chi-squared < crit

  suppressWarnings(try(newFit <- lavaan::update(y, data = NULL,
                                                sample.cov = Sigma,
                                                sample.mean = Mu,
                                                sample.nobs = N)))
  if(!lavaan::lavInspect(newFit, "converged")) return(NA) else {
    result <- lavaan::lavInspect(newFit, "fit")["chisq"] < crit
    names(result) <- NULL
    if (lavaan::lavInspect(x, "fit")["df"] ==
        lavaan::lavInspect(y, "fit")["df"]) return(c(Equivalent = result))
  }
  c(Nested = result)
}

## generic function that utilizes "x.within.y" to test a set of models
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

  ## check whether any models include categorical outcomes
  catMod <- sapply(fitList, function(x) lavaan::lavInspect(x, "options")$categorical)
  if (any(catMod)) stop("This method only applies to continuous outcomes.")

  ## get degrees of freedom for each model
  DFs <- sapply(fitList, function(x) lavaan::lavInspect(x, "fit")["df"])

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

    # class(nestMat) <- c("Net", class(nestMat))
  # attr(nestMat, "df") <- orderedDFs
    out <- new("Net",
               test         = nestMat,
               df       = orderedDFs
              )
  out
}

setMethod("show", "Net",
function(object) {
  if (length(object@test)) {
    m <- as.matrix(unclass(object@test))
    m[upper.tri(m, diag = TRUE)] <- ""
    cat("
     If cell [R, C] is TRUE, the model in row R is nested within column C.

     If cell [R, C] is TRUE and the models have the same degrees of freedom,
     they are equivalent models.  See Bentler & Satorra (2010) for details.

     If cell [R, C] is NA, then the model in column C did not converge when
     fit to the implied means and covariance matrix from the model in row R.

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

