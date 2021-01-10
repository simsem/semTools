### Terrence D. Jorgensen
### Last updated: 10 January 2021
### runMI creates lavaan.mi object, inherits from lavaanList class


## -------------
## Main function
## -------------


##' Fit a lavaan Model to Multiple Imputed Data Sets
##'
##' This function fits a lavaan model to a list of imputed data sets, and can
##' also implement multiple imputation for a single \code{data.frame} with
##' missing observations, using either the Amelia package or the mice package.
##'
##'
##' @aliases runMI lavaan.mi cfa.mi sem.mi growth.mi
##' @importFrom lavaan lavInspect parTable
##' @importFrom methods getMethod
##'
##' @param model The analysis model can be specified using lavaan
##'   \code{\link[lavaan]{model.syntax}} or a parameter table (as returned by
##'   \code{\link[lavaan]{parTable}}).
##' @param data A \code{data.frame} with missing observations, or a \code{list}
##'   of imputed data sets (if data are imputed already). If \code{runMI} has
##'   already been called, then imputed data sets are stored in the
##'   \code{@@DataList} slot, so \code{data} can also be a \code{lavaan.mi} object
##'   from which the same imputed data will be used for additional analyses.
##' @param fun \code{character}. Name of a specific lavaan function used to fit
##'   \code{model} to \code{data} (i.e., \code{"lavaan"}, \code{"cfa"},
##'   \code{"sem"}, or \code{"growth"}). Only required for \code{runMI}.
##' @param \dots additional arguments to pass to \code{\link[lavaan]{lavaan}} or
##'   \code{\link[lavaan]{lavaanList}}. See also \code{\link[lavaan]{lavOptions}}.
##'   Note that \code{lavaanList} provides parallel computing options, as well as
##'   a \code{FUN} argument so the user can extract custom output after the model
##'   is fitted to each imputed data set (see \strong{Examples}).  TIP: If a
##'   custom \code{FUN} is used \emph{and} \code{parallel = "snow"} is requested,
##'   the user-supplied function should explicitly call \code{library} or use
##'   \code{\link[base]{::}} for any functions not part of the base distribution.
##' @param m \code{integer}. Request the number of imputations. Ignored if
##'   \code{data} is already a \code{list} of imputed data sets or a
##'   \code{lavaan.mi} object.
##' @param miArgs Addition arguments for the multiple-imputation function
##'   (\code{miPackage}). The arguments should be put in a list (see example
##'   below). Ignored if \code{data} is already a \code{list} of imputed data
##'   sets or a \code{lavaan.mi} object.
##' @param miPackage Package to be used for imputation. Currently these
##'   functions only support \code{"Amelia"} or \code{"mice"} for imputation.
##'   Ignored if \code{data} is already a \code{list} of imputed data sets or a
##'   \code{lavaan.mi} object.
##' @param seed \code{integer}. Random number seed to be set before imputing the
##'   data. Ignored if \code{data} is already a \code{list} of imputed data sets
##'   or a \code{lavaan.mi} object.
##'
##' @return A \code{\linkS4class{lavaan.mi}} object
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Enders, C. K. (2010). \emph{Applied missing data analysis}. New
##'   York, NY: Guilford.
##'
##'   Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}.
##'   New York, NY: Wiley.
##'
##' @examples
##'  \dontrun{
##' ## impose missing data for example
##' HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
##'                                       "ageyr","agemo","school")]
##' set.seed(12345)
##' HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
##' age <- HSMiss$ageyr + HSMiss$agemo/12
##' HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
##'
##' ## specify CFA model from lavaan's ?cfa help page
##' HS.model <- '
##'   visual  =~ x1 + x2 + x3
##'   textual =~ x4 + x5 + x6
##'   speed   =~ x7 + x8 + x9
##' '
##'
##' ## impute data within runMI...
##' out1 <- cfa.mi(HS.model, data = HSMiss, m = 20, seed = 12345,
##'                miArgs = list(noms = "school"))
##'
##' ## ... or impute missing data first
##' library(Amelia)
##' set.seed(12345)
##' HS.amelia <- amelia(HSMiss, m = 20, noms = "school", p2s = FALSE)
##' imps <- HS.amelia$imputations
##' out2 <- cfa.mi(HS.model, data = imps)
##'
##' ## same results (using the same seed results in the same imputations)
##' cbind(impute.within = coef(out1), impute.first = coef(out2))
##'
##' summary(out1, fit.measures = TRUE)
##' summary(out1, ci = FALSE, fmi = TRUE, output = "data.frame")
##' summary(out1, ci = FALSE, stand = TRUE, rsq = TRUE)
##'
##' ## model fit. D3 includes information criteria
##' anova(out1)
##' ## equivalently:
##' lavTestLRT.mi(out1)
##' ## request D2
##' anova(out1, test = "D2")
##' ## request fit indices
##' fitMeasures(out1)
##'
##'
##' ## fit multigroup model without invariance constraints
##' mgfit.config <- cfa.mi(HS.model, data = imps, estimator = "mlm",
##'                        group = "school")
##' ## add invariance constraints, and use previous fit as "data"
##' mgfit.metric <- cfa.mi(HS.model, data = mgfit.config, estimator = "mlm",
##'                        group = "school", group.equal = "loadings")
##' mgfit.scalar <- cfa.mi(HS.model, data = mgfit.config, estimator = "mlm",
##'                        group = "school",
##'                        group.equal = c("loadings","intercepts"))
##'
##' ## compare fit of 2 models to test metric invariance
##' ## (scaled likelihood ratio test)
##' lavTestLRT.mi(mgfit.metric, h1 = mgfit.config)
##' ## To compare multiple models, you must use anova()
##' anova(mgfit.config, mgfit.metric, mgfit.scalar)
##' ## or compareFit(), which also includes fit indices for comparison
##' ## (optional: name the models)
##' compareFit(config = mgfit.config, metric = mgfit.metric,
##'            scalar = mgfit.scalar,
##'            argsLRT = list(test = "D2", method = "satorra.bentler.2010"))
##'
##' ## correlation residuals to investigate local misfit
##' resid(mgfit.scalar, type = "cor.bentler")
##' ## modification indices for fixed parameters, to investigate local misfit
##' modindices.mi(mgfit.scalar)
##' ## or lavTestScore.mi for modification indices about equality constraints
##' lavTestScore.mi(mgfit.scalar)
##'
##' ## Wald test of whether latent means are == (fix 3 means to zero in group 2)
##' eq.means <- ' .p70. == 0
##'               .p71. == 0
##'               .p72. == 0 '
##' lavTestWald.mi(mgfit.scalar, constraints = eq.means)
##'
##'
##'
##' ## ordered-categorical data
##' data(datCat)
##' lapply(datCat, class) # indicators already stored as ordinal
##' ## impose missing values
##' set.seed(123)
##' for (i in 1:8) datCat[sample(1:nrow(datCat), size = .1*nrow(datCat)), i] <- NA
##'
##' ## impute ordinal missing data using mice package
##' library(mice)
##' set.seed(456)
##' miceImps <- mice(datCat)
##' ## save imputations in a list of data.frames
##' impList <- list()
##' for (i in 1:miceImps$m) impList[[i]] <- complete(miceImps, action = i)
##'
##' ## fit model, save zero-cell tables and obsolete "WRMR" fit indices
##' catout <- cfa.mi(' f =~ 1*u1 + 1*u2 + 1*u3 + 1*u4 ', data = impList,
##'                  FUN = function(fit) {
##'                    list(wrmr = lavaan::fitMeasures(fit, "wrmr"),
##'                         zeroCells = lavaan::lavInspect(fit, "zero.cell.tables"))
##'                  })
##' summary(catout)
##' lavTestLRT.mi(catout, test = "D2", pool.robust = TRUE)
##' fitMeasures(catout, fit.measures = c("rmsea","srmr","cfi"),
##'             test = "D2", pool.robust = TRUE)
##'
##' ## extract custom output
##' sapply(catout@funList, function(x) x$wrmr) # WRMR for each imputation
##' catout@funList[[1]]$zeroCells # zero-cell tables for first imputation
##' catout@funList[[2]]$zeroCells # zero-cell tables for second imputation ...
##'
##' }
##'
##' @export
runMI <- function(model, data, fun = "lavaan", ...,
                  m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  CALL <- match.call()
  dots <- list(...)

  ## check for (Bollen-Stine) bootstrap request
  if (all(!is.null(dots$test),
          tolower(dots$test) %in% c("boot","bootstrap","bollen.stine")) ||
      all(!is.null(dots$se), tolower(dots$se) %in% c("boot","bootstrap"))) {
    stop('Bootstraping unavailable (and not recommended) in combination with ',
         'multiple imputations. For robust confidence intervals of indirect',
         ' effects, see the ?semTools::monteCarloMed help page. To bootstrap ',
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
  } else if (is(data, "lavaan.mi")) {
    seed <- data@seed
    imputeCall <- data@imputeCall
    imputedData <- data@DataList
    m <- length(imputedData)
  } else stop("data is not a valid input type: a partially observed data.frame,",
              " a list of imputed data.frames, or previous lavaan.mi object")

  ## Function to get custom output for lavaan.mi object
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
  fit <- as(fit, "lavaan.mi")
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

  NewStartVals <- try(getMethod("coef", "lavaan.mi")(fit, type = "user",
                                                          labels = FALSE),
                           silent = TRUE)
  if (!inherits(NewStartVals, "try-error")) fit@ParTable$start <- NewStartVals
  fit
}

##' @rdname runMI
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

##' @rdname runMI
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

##' @rdname runMI
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

##' @rdname runMI
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
##' several functions used for \code{\linkS4class{lavaan.mi}} objects, such as
##' \code{\link{lavTestLRT.mi}}, \code{\link{lavTestWald.mi}}, and
##' \code{\link{lavTestScore.mi}}. But this function can be used for any general
##' scenario because it only requires a vector of \eqn{\chi^2} statistics (one
##' from each imputation) and the degrees of freedom for the test statistic.
##' See Li, Meng, Raghunathan, & Rubin (1991) and Enders (2010, chapter 8) for
##' details about how it is calculated.
##'
##' @importFrom stats var pf pchisq
##'
##' @param w \code{numeric} vector of Wald \eqn{\chi^2} statistics. Can also
##'   be Wald \emph{z} statistics, which will be internally squared to make
##'   \eqn{\chi^2} statistics with one \emph{df} (must set \code{DF = 0L}).
##' @param DF degrees of freedom (\emph{df}) of the \eqn{\chi^2} statistics.
##'   If \code{DF = 0L} (default), \code{w} is assumed to contain \emph{z}
##'   statistics, which will be internally squared.
##' @param asymptotic \code{logical}. If \code{FALSE} (default), the pooled test
##'   will be returned as an \emph{F}-distributed statistic with numerator
##'   (\code{df1}) and denominator (\code{df2}) degrees of freedom.
##'   If \code{TRUE}, the pooled \emph{F} statistic will be multiplied by its
##'   \code{df1} on the assumption that its \code{df2} is sufficiently large
##'   enough that the statistic will be asymptotically \eqn{\chi^2} distributed
##'   with \code{df1}.
##'
##' @return A \code{numeric} vector containing the test statistic, \emph{df},
##'   its \emph{p} value, and 2 missing-data diagnostics: the relative invrease
##'   in variance (RIV, or average for multiparameter tests: ARIV) and the
##'   fraction missing information (FMI = ARIV / (1 + ARIV)).
##'
##' @seealso \code{\link{lavTestLRT.mi}}, \code{\link{lavTestWald.mi}},
##'   \code{\link{lavTestScore.mi}}
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Enders, C. K. (2010). \emph{Applied missing data analysis}. New
##'   York, NY: Guilford.
##'
##'   Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
##'   Significance levels from repeated \emph{p}-values with multiply-imputed
##'   data. \emph{Statistica Sinica, 1}(1), 65--92. Retrieved from
##'   \url{https://www.jstor.org/stable/24303994}
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
##'
##' @export
calculate.D2 <- function(w, DF = 0L, asymptotic = FALSE) {
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

