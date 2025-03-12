### Terrence D. Jorgensen
### Last updated: 12 March 2025
### runMI creates OLDlavaan.mi object, inherits from lavaanList class

### DEPRECATED: 16 June 2024
### supplanted by lavaan.mi package

## -------------
## Main function
## -------------


##' Fit a lavaan Model to Multiple Imputed Data Sets
##'
##' This function fits a lavaan model to a list of imputed data sets, and can
##' also implement multiple imputation for a single `data.frame` with
##' missing observations, using either the Amelia package or the mice package.
##'
##'
##' @aliases runMI-deprecated lavaan.mi-deprecated cfa.mi-deprecated sem.mi-deprecated growth.mi-deprecated
##' @importFrom lavaan lavInspect parTable
##' @importFrom methods getMethod
##'
##' @param model The analysis model can be specified using
##'   [lavaan::model.syntax()] or a [lavaan::parTable()]
##' @param data A `data.frame` with missing observations, or a `list`
##'   of imputed data sets (if data are imputed already). If `runMI()` has
##'   already been called, then imputed data sets are stored in the
##'   `@@DataList` slot, so `data=` can also be an `OLDlavaan.mi` object
##'   from which the same imputed data will be used for additional analyses.
##' @param fun `character`. Name of a specific lavaan function used to fit
##'   `model=` to `data=` (i.e., `"lavaan"`, `"cfa"`, `"sem"`, or `"growth"`).
##'   Only required for `runMI()`.
##' @param \dots additional arguments to pass to [lavaan::lavaan()] or
##'   [lavaan::lavaanList()]. See also [lavaan::lavOptions()].
##'   Note that `lavaanList` provides parallel computing options, as well as
##'   a `FUN` argument so the user can extract custom output after the model
##'   is fitted to each imputed data set (see **Examples**).  TIP: If a
##'   custom `FUN` is used *and* `parallel = "snow"` is requested,
##'   the user-supplied function should explicitly call `library` or use
##'   \code{\link[base]{::}} for any functions not part of the base distribution.
##' @param m `integer`. Request the number of imputations. Ignored if `data=` is
##'   already a `list` of imputed data sets or an `OLDlavaan.mi` object.
##' @param miArgs Addition arguments for the multiple-imputation function
##'   (`miPackage`). The arguments should be put in a list (see example
##'   below). Ignored if `data=` is already a `list` of imputed data
##'   sets or an `OLDlavaan.mi` object.
##' @param miPackage Package to be used for imputation. Currently these
##'   functions only support `"Amelia"` or `"mice"` for imputation.
##'   Ignored if `data` is already a `list` of imputed data sets or an
##'   `OLDlavaan.mi` object.
##' @param seed `integer`. Random number seed to be set before imputing the
##'   data. Ignored if `data` is already a `list` of imputed data sets
##'   or an `OLDlavaan.mi` object.
##'
##' @return A [semTools::OLDlavaan.mi-class] object
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Enders, C. K. (2010). *Applied missing data analysis*. New
##'   York, NY: Guilford.
##'
##'   Rubin, D. B. (1987). *Multiple imputation for nonresponse in surveys*.
##'   New York, NY: Wiley.
##'
##' @examples
##'
##' ## See the new lavaan.mi package
##'
##' @name runMI-deprecated
##' @usage
##' runMI(model, data, fun = "lavaan", ...,
##'       m, miArgs = list(), miPackage = "Amelia", seed = 12345)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL


##' @rdname semTools-deprecated
##' @section `runMI()` and the `lavaan.mi` class functionality:
##' The `runMI()` function and support for `lavaan.mi-class` objects became
##' such a large part of semTools that it made sense to move that functionality
##' to its own package.  The \pkg{lavaan.mi} package is now available for users
##' to fit `lavaan` models to their multiply imputed data.  The new package
##' already fixes many bugs and provides many new features that make the
##' semTools `OLDlavaan.mi-class` obsolete.  Please immediately discontinue
##' your dependence on `semTools::runMI()` amd transition to the new
##' \pkg{lavaan.mi} package, which also provides a more similar user interface
##' as the \pkg{lavaan} package provides for a single `data=` set. The README
##' on <https://github.com/TDJorgensen/lavaan.mi> provides a list of analogous
##' functionality in \pkg{lavaan} and \pkg{lavaan.mi}, and the NEWS file
##' documents new features and other differences from the deprecated
##' `semTools::runMI()` functionality.
##'
##' @export
runMI <- function(model, data, fun = "lavaan", ...,
                  m, miArgs = list(), miPackage = "Amelia", seed = 12345) {

  .Deprecated(msg = c("\nThe runMI() function and lavaan.mi-class have been ",
                      "deprecated and will cease to be included in future ",
                      "versions of semTools.\n\nSupport is still provided for ",
                      "analyzing lavaan.mi-class objects (e.g., compRelSEM() ",
                      "can estimate reliability using pooled results), which ",
                      "can now be created using the lavaan.mi package.\n\nThe ",
                      "deprecated runMI() function now creates an object of ",
                      "class OLDlavaan.mi, which can be analyzed using the ",
                      "deprecated functions in semTools, like lavTestLRT.mi(),",
                      " that have been updated and improved in the lavaan.mi ",
                      "package.\n\nFind more details help('semTools-deprecated)"))

  CALL <- match.call()
  dots <- list(...)

  ## check for (Bollen-Stine) bootstrap request
  if (all(!is.null(dots$test),
          tolower(dots$test) %in% c("boot","bootstrap","bollen.stine")) ||
      all(!is.null(dots$se), tolower(dots$se) %in% c("boot","bootstrap"))) {
    stop('Bootstraping unavailable (and not recommended) in combination with ',
         'multiple imputations. For robust confidence intervals of indirect',
         ' effects, see the ?semTools::monteCarloCI help page. To bootstrap ',
         'within each imputation, users can pass a custom function to the ',
         'FUN= argument (see ?lavaanList) to save bootstrap distributions in ',
         'the @funList slot, then manually combine afterward.')
  }

  seed <- as.integer(seed[1])
  ## Create (or acknowledge) list of imputed data sets
  imputedData <- NULL
  if (missing(data)) {
    #TODO: check for summary statistics
    #TODO: make lavaanList() accept lists of summary stats
    #TODO: Add argument to implement Li Cai's pool-polychorics first, pass
    #      to lavaan for DWLS with pooled WLS.V= and NACOV=, return(lavaan).

  } else if (is.data.frame(data)) {
    if (miPackage[1] == "Amelia") {
      requireNamespace("Amelia")
      if (!"package:Amelia" %in% search()) attachNamespace("Amelia")
      imputeCall <- c(list(Amelia::amelia, x = data, m = m, p2s = 0), miArgs)
      set.seed(seed)
      imputedData <- unclass(eval(as.call(imputeCall))$imputations)
    } else if (miPackage[1] == "mice") {
      requireNamespace("mice")
      if (!"package:mice" %in% search()) attachNamespace("mice")
      imputeCall <- c(list(mice::mice, data = data, m = m, diagnostics = FALSE,
                           printFlag = FALSE), miArgs)
      set.seed(seed)
      miceOut <- eval(as.call(imputeCall))
      imputedData <- list()
      for (i in 1:m) {
        imputedData[[i]] <- mice::complete(data = miceOut, action = i, include = FALSE)
      }
    } else stop("Currently runMI only supports imputation by Amelia or mice")
  } else if (is.list(data)) {
    ## check possibility that it is a mids object (inherits from list)
    if (requireNamespace("mice", quietly = TRUE)) {
      if (mice::is.mids(data)) {
        m <- data$m
        imputedData <- list()
        for (i in 1:m) {
          imputedData[[i]] <- mice::complete(data, action = i, include = FALSE)
        }
        imputeCall <- list()
      } else {
        seed <- integer(length = 0)
        imputeCall <- list()
        imputedData <- data
        m <- length(data)
        class(imputedData) <- "list" # override inheritance (e.g., "mi" if Amelia)
      }
    } else {
      ## can't check for mice, so probably isn't mids
      seed <- integer(length = 0)
      imputeCall <- list()
      imputedData <- data
      m <- length(data)
      class(imputedData) <- "list" # override inheritance (e.g., "mi" if Amelia)
    }
  } else if (is(data, "OLDlavaan.mi")) {
    seed <- data@seed
    imputeCall <- data@imputeCall
    imputedData <- data@DataList
    m <- length(imputedData)
  } else stop("data is not a valid input type: a partially observed data.frame,",
              " a list of imputed data.frames, or previous OLDlavaan.mi object")

  ## Function to get custom output for OLDlavaan.mi object
  ## NOTE: Need "lavaan::" to allow for parallel computations
  .getOutput. <- function(obj) {
    converged <- lavaan::lavInspect(obj, "converged")
    if (converged) {
      se <- lavaan::parTable(obj)$se
      se.test <- all(!is.na(se)) & all(se >= 0) & any(se != 0)
      if (lavaan::lavInspect(obj, "ngroups") == 1L && lavaan::lavInspect(obj, "nlevels") == 1L) {
        Heywood.lv <- det(lavaan::lavInspect(obj, "cov.lv")) <= 0
        Heywood.ov <- det(lavaan::lavInspect(obj, "theta")) <= 0
      } else {
        Heywood.lv <- !all(sapply(lavaan::lavInspect(obj, "cov.lv"), det) > 0)
        Heywood.ov <- !all(sapply(lavaan::lavInspect(obj, "theta"), det) > 0)
      }
    } else {
      se.test <- Heywood.lv <- Heywood.ov <- NA
    }
    list(sampstat = lavaan::lavInspect(obj, "sampstat"),
         coefMats = lavaan::lavInspect(obj, "est"),
         satPT = data.frame(lavaan::lav_partable_unrestricted(obj),
                            #FIXME: do starting values ALWAYS == estimates?
                            stringsAsFactors = FALSE),
         modindices = try(lavaan::modindices(obj), silent = TRUE),
         cov.lv = lavaan::lavInspect(obj, "cov.lv"), #TODO: calculate from pooled estimates for reliability()
         converged = converged, SE = se.test,
         Heywood.lv = Heywood.lv, Heywood.ov = Heywood.ov)
  }

  ## fit model using lavaanList
  lavListCall <- list(lavaan::lavaanList, model = model, dataList = imputedData,
                      cmd = fun)
  lavListCall <- c(lavListCall, dots)
  lavListCall$store.slots <- c("partable","vcov","test","h1","baseline")
  lavListCall$FUN <- if (is.null(dots$FUN)) .getOutput. else function(obj) {
    temp1 <- .getOutput.(obj)
    temp2 <- dots$FUN(obj)
    if (!is.list(temp2)) temp2 <- list(userFUN1 = temp2)
    if (is.null(names(temp2))) names(temp2) <- paste0("userFUN", 1:length(temp2))
    duplicatedNames <- which(sapply(names(temp2), function(x) {
      x %in% c("sampstat","coefMats","satPT","modindices","converged",
               "SE","Heywood.lv","Heywood.ov","cov.lv")
    }))
    for (i in duplicatedNames) names(temp2)[i] <- paste0("userFUN", i)
    c(temp1, temp2)
  }
  fit <- eval(as.call(lavListCall))
  ## Store custom @DataList and @SampleStatsList
  fit@SampleStatsList <- lapply(fit@funList, "[[", i = "sampstat")
  fit@DataList <- imputedData
  ## add parameter table to @h1List
  for (i in 1:m) fit@h1List[[i]] <- c(fit@h1List[[i]],
                                      list(PT = fit@funList[[i]]$satPT))
  ## assign class and add new slots
  fit <- as(fit, "OLDlavaan.mi")
  fit@coefList <- lapply(fit@funList, "[[", i = "coefMats")
  fit@miList <- lapply(fit@funList, "[[", i = "modindices")
  fit@phiList <- lapply(fit@funList, "[[", i = "cov.lv")
  fit@seed <- seed
  fit@call <- CALL
  fit@lavListCall <- lavListCall
  fit@imputeCall <- imputeCall
  convList <- lapply(fit@funList, "[", i = c("converged","SE",
                                             "Heywood.lv","Heywood.ov"))
  nonConv <- which(sapply(convList, is.null))
  if (length(nonConv)) for (i in nonConv) {
    convList[[i]] <- list(converged = FALSE, SE = NA, Heywood.lv = NA, Heywood.ov = NA)
  }
  fit@convergence <- lapply(convList, function(x) do.call(c, x))
  conv <- which(sapply(fit@convergence, "[", i = "converged"))
  if (!length(conv)) warning('The model did not converge for any imputed data sets.')

  ## keep any remaining funList slots (if allowing users to supply custom FUN)
  funNames <- names(fit@funList[[1]])
  keepIndex <- which(!sapply(funNames, function(x) {
    x %in% c("sampstat","coefMats","satPT","modindices","converged",
             "SE","Heywood.lv","Heywood.ov","cov.lv")
  }))
  if (length(keepIndex)) {
    fit@funList <- lapply(fit@funList, "[", i = keepIndex)
    if (length(keepIndex) > 1L) {
      keepNames <- funNames[keepIndex]
      noNames <- which(keepNames == "")
      for (i in seq_along(noNames)) keepNames[ noNames[i] ] <- paste0("userFUN", i)
      fit@funList <- lapply(fit@funList, "names<-", value = keepNames)
    }
  } else fit@funList <- list()

  NewStartVals <- try(getMethod("coef", "OLDlavaan.mi")(fit, type = "user",
                                                        labels = FALSE),
                      silent = TRUE)
  if (!inherits(NewStartVals, "try-error")) fit@ParTable$start <- NewStartVals
  fit
}




##' @name runMI-deprecated
##' @usage
##' lavaan.mi(model, data, ...,
##'           m, miArgs = list(), miPackage = "Amelia", seed = 12345)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL

##' @rdname semTools-deprecated
##' @export
lavaan.mi <- function(model, data, ...,
                      m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  # runMI(model = model, data = data, fun = "lavaan", ...,
  #       m = m, miArgs = miArgs, miPackage = miPackage, seed = seed)
  mc <- match.call(expand.dots = TRUE)
  mc$fun <- "lavaan"
  mc[[1L]] <- quote(semTools::runMI)
  eval(mc, parent.frame())
}

##' @name runMI-deprecated
##' @usage
##' cfa.mi(model, data, ...,
##'        m, miArgs = list(), miPackage = "Amelia", seed = 12345)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL

##' @rdname semTools-deprecated
##' @export
cfa.mi <- function(model, data, ...,
                   m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  # runMI(model = model, data = data, fun = "cfa", ...,
  #       m = m, miArgs = miArgs, miPackage = miPackage, seed = seed)
  mc <- match.call(expand.dots = TRUE)
  mc$fun <- "cfa"
  mc[[1L]] <- quote(semTools::runMI)
  eval(mc, parent.frame())
}

##' @name runMI-deprecated
##' @usage
##' sem.mi(model, data, ...,
##'        m, miArgs = list(), miPackage = "Amelia", seed = 12345)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL

##' @rdname semTools-deprecated
##' @export
sem.mi <- function(model, data, ...,
                   m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  # runMI(model = model, data = data, fun = "sem", ...,
  #       m = m, miArgs = miArgs, miPackage = miPackage, seed = seed)
  mc <- match.call(expand.dots = TRUE)
  mc$fun <- "sem"
  mc[[1L]] <- quote(semTools::runMI)
  eval(mc, parent.frame())
}

##' @name runMI-deprecated
##' @usage
##' growth.mi(model, data, ...,
##'           m, miArgs = list(), miPackage = "Amelia", seed = 12345)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL

##' @rdname semTools-deprecated
##' @export
growth.mi <- function(model, data, ...,
                      m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  # runMI(model = model, data = data, fun = "growth", ...,
  #       m = m, miArgs = miArgs, miPackage = miPackage, seed = seed)
  mc <- match.call(expand.dots = TRUE)
  mc$fun <- "growth"
  mc[[1L]] <- quote(semTools::runMI)
  eval(mc, parent.frame())
}



## -----------------
## Utility functions
## -----------------

##' Calculate the "D2" statistic
##'
##' This is a utility function used to calculate the "D2" statistic for pooling
##' test statistics across multiple imputations. This function is called by
##' several functions used for [OLDlavaan.mi-class] objects, such as
##' [lavTestLRT.mi()], [lavTestWald.mi()], and
##' [lavTestScore.mi()]. But this function can be used for any general
##' scenario because it only requires a vector of \eqn{\chi^2} statistics (one
##' from each imputation) and the degrees of freedom for the test statistic.
##' See Li, Meng, Raghunathan, & Rubin (1991) and Enders (2010, chapter 8) for
##' details about how it is calculated.
##'
##' @importFrom stats var pf pchisq
##'
##' @param w `numeric` vector of Wald \eqn{\chi^2} statistics. Can also
##'   be Wald *z* statistics, which will be internally squared to make
##'   \eqn{\chi^2} statistics with one *df* (must set `DF = 0L`).
##' @param DF degrees of freedom (*df*) of the \eqn{\chi^2} statistics.
##'   If `DF = 0L` (default), `w` is assumed to contain *z*
##'   statistics, which will be internally squared.
##' @param asymptotic `logical`. If `FALSE` (default), the pooled test
##'   will be returned as an *F*-distributed statistic with numerator
##'   (`df1`) and denominator (`df2`) degrees of freedom.
##'   If `TRUE`, the pooled *F* statistic will be multiplied by its
##'   `df1` on the assumption that its `df2` is sufficiently large
##'   enough that the statistic will be asymptotically \eqn{\chi^2} distributed
##'   with `df1`.
##'
##' @return A `numeric` vector containing the test statistic, *df*,
##'   its *p* value, and 2 missing-data diagnostics: the relative invrease
##'   in variance (RIV, or average for multiparameter tests: ARIV) and the
##'   fraction missing information (FMI = ARIV / (1 + ARIV)).
##'
##' @seealso [lavTestLRT.mi()], [lavTestWald.mi()], [lavTestScore.mi()]
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Enders, C. K. (2010). *Applied missing data analysis*. New
##'   York, NY: Guilford.
##'
##'   Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
##'   Significance levels from repeated *p*-values with multiply-imputed
##'   data. *Statistica Sinica, 1*(1), 65--92. Retrieved from
##'   <https://www.jstor.org/stable/24303994>
##'
##' @examples
##' ## generate a vector of chi-squared values, just for example
##' DF <- 3 # degrees of freedom
##' M <- 20 # number of imputations
##' CHI <- rchisq(M, DF)
##'
##' ## pool the "results"
##' calculate.D2(CHI, DF) # by default, an F statistic is returned
##' calculate.D2(CHI, DF, asymptotic = TRUE) # asymptotically chi-squared
##'
##' ## generate standard-normal values, for an example of Wald z tests
##' Z <- rnorm(M)
##' calculate.D2(Z) # default DF = 0 will square Z to make chisq(DF = 1)
##' ## F test is equivalent to a t test with the denominator DF
##'
##' @name calculate.D2-deprecated
##' @usage calculate.D2(w, DF = 0L, asymptotic = FALSE)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL


##' @rdname semTools-deprecated
##'
##' @export
calculate.D2 <- function(w, DF = 0L, asymptotic = FALSE) {

  .Deprecated(msg = c("\nThe calculate.D2() function has been deprecated from ",
                      "semTools and moved to the new lavaan.mi package, along ",
                      "with other multiple-imputation functionality related to",
                      " runMI().\nSee help('semTools-deprecated)"))

  if (length(w) == 0L) return(NA)
  w <- as.numeric(w)
  DF <- as.numeric(DF)

  nImps <- sum(!is.na(w))
  if (nImps == 0) return(NA)

  if (DF <= 0L) {
    ## assume they are Z scores
    w <- w^2
    DF <- 1L
  }

  ## pool test statistics
  if (length(w) > 1L) {
    w_bar <- mean(w, na.rm = TRUE)
    ariv <- (1 + 1/nImps) * var(sqrt(w), na.rm = TRUE)
    test.stat <- (w_bar/DF - ((nImps + 1) * ariv / (nImps - 1))) / (1 + ariv)
  } else {
    warning('There was only 1 non-missing value to pool, leading to zero ',
            'variance, so D2 cannot be calculated.')
    test.stat <- ariv <- NA
  }
  if (test.stat < 0) test.stat <- 0
  if (asymptotic) {
    out <- c("chisq" = test.stat * DF, df = DF,
             pvalue = pchisq(test.stat * DF, df = DF, lower.tail = FALSE),
             ariv = ariv, fmi = ariv / (1 + ariv))
  } else {
    v3 <- DF^(-3 / nImps) * (nImps - 1) * (1 + (1 / ariv))^2
    out <- c("F" = test.stat, df1 = DF, df2 = v3,
             pvalue = pf(test.stat, df1 = DF, df2 = v3, lower.tail = FALSE),
             ariv = ariv, fmi = ariv / (1 + ariv))
  }
  class(out) <- c("lavaan.vector","numeric")
  out
}

