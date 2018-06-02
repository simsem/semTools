### Terrence D. Jorgensen
### Last updated: 3 June 2018
### runMI creates lavaan.mi object, inherits from lavaanList class


#' Fit a lavaan Model to Multiple Imputed Data Sets
#'
#' This function fits a lavaan model to a list of imputed data sets, and can
#' also implement multiple imputation for a single \code{data.frame} with
#' missing observations, using either the Amelia package or the mice package.
#'
#'
#' @aliases runMI lavaan.mi cfa.mi sem.mi growth.mi
#' @importFrom lavaan lavInspect parTable
#'
#' @param model The analysis model can be specified using lavaan
#' \code{\link[lavaan]{model.syntax}} or a parameter table (as returned by
#' \code{\link[lavaan]{parTable}}).
#' @param data A \code{data.frame} with missing observations, or a \code{list}
#' of imputed data sets (if data are imputed already). If \code{runMI} has
#' already been called, then imputed data sets are stored in the
#' \code{@DataList} slot, so \code{data} can also be a \code{lavaan.mi} object
#' from which the same imputed data will be used for additional analyses.
#' @param fun \code{character}. Name of a specific lavaan function used to fit
#' \code{model} to \code{data} (i.e., \code{"lavaan"}, \code{"cfa"},
#' \code{"sem"}, or \code{"growth"}). Only required for \code{runMI}.
#' @param \dots additional arguments to pass to \code{\link[lavaan]{lavaan}} or
#' \code{\link[lavaan]{lavaanList}}. See also \code{\link[lavaan]{lavOptions}}.
#' Note that \code{lavaanList} provides parallel computing options, as well as
#' a \code{FUN} argument so the user can extract custom output after the model
#' is fitted to each imputed data set (see \strong{Examples}).  TIP: If a
#' custom \code{FUN} is used \emph{and} \code{parallel = "snow"} is requested,
#' the user-supplied function should explicitly call \code{library} or use
#' \code{\link[base]{::}} for any functions not part of the base distribution.
#' @param m \code{integer}. Request the number of imputations. Ignored if
#' \code{data} is already a \code{list} of imputed data sets or a
#' \code{lavaan.mi} object.
#' @param miArgs Addition arguments for the multiple-imputation function
#' (\code{miPackage}). The arguments should be put in a list (see example
#' below). Ignored if \code{data} is already a \code{list} of imputed data sets
#' or a \code{lavaan.mi} object.
#' @param miPackage Package to be used for imputation. Currently these
#' functions only support \code{"Amelia"} or \code{"mice"} for imputation.
#' Ignored if \code{data} is already a \code{list} of imputed data sets or a
#' \code{lavaan.mi} object.
#' @param seed \code{integer}. Random number seed to be set before imputing the
#'  data. Ignored if \code{data} is already a \code{list} of imputed data sets
#'  or a \code{lavaan.mi} object.
#' @return A \code{\linkS4class{lavaan.mi}} object
#' @author Terrence D. Jorgensen (University of Amsterdam;
#' \email{TJorgensen314@@gmail.com})
#' @references Enders, C. K. (2010). \emph{Applied missing data analysis}. New
#' York, NY: Guilford.
#'
#' Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}.
#' New York, NY: Wiley.
#' @examples
#'  \dontrun{
#' ## impose missing data for example
#' HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
#'                                       "ageyr","agemo","school")]
#' set.seed(12345)
#' HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
#' age <- HSMiss$ageyr + HSMiss$agemo/12
#' HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
#'
#' ## specify CFA model from lavaan's ?cfa help page
#' HS.model <- '
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#' '
#'
#' ## impute data within runMI...
#' out1 <- cfa.mi(HS.model, data = HSMiss, m = 20, seed = 12345,
#'                miArgs = list(noms = "school"))
#'
#' ## ... or impute missing data first
#' library(Amelia)
#' set.seed(12345)
#' HS.amelia <- amelia(HSMiss, m = 20, noms = "school", p2s = FALSE)
#' imps <- HS.amelia$imputations
#' out2 <- cfa.mi(HS.model, data = imps)
#'
#' ## same results (using the same seed results in the same imputations)
#' cbind(impute.within = coef(out1), impute.first = coef(out2))
#'
#' summary(out1)
#' summary(out1, ci = FALSE, fmi = TRUE, add.attributes = FALSE)
#' summary(out1, ci = FALSE, stand = TRUE, rsq = TRUE)
#'
#' ## model fit. D3 includes information criteria
#' anova(out1)
#' anova(out1, test = "D2", indices = TRUE) # request D2 and fit indices
#'
#'
#'
#' ## fit multigroup model without invariance constraints
#' mgfit1 <- cfa.mi(HS.model, data = imps, estimator = "mlm", group = "school")
#' ## add invariance constraints, and use previous fit as "data"
#' mgfit0 <- cfa.mi(HS.model, data = mgfit1, estimator = "mlm", group = "school",
#'                  group.equal = c("loadings","intercepts"))
#'
#' ## compare fit (scaled likelihood ratio test)
#' anova(mgfit0, h1 = mgfit1)
#'
#' ## correlation residuals
#' resid(mgfit0, type = "cor.bentler")
#'
#'
#' ## use D1 to test a parametrically nested model (whether latent means are ==)
#' anova(mgfit0, test = "D1", constraints = '
#'       .p70. == 0
#'       .p71. == 0
#'       .p72. == 0')
#'
#'
#'
#' ## ordered-categorical data
#' data(datCat)
#' lapply(datCat, class)
#' ## impose missing values
#' set.seed(123)
#' for (i in 1:8) datCat[sample(1:nrow(datCat), size = .1*nrow(datCat)), i] <- NA
#'
#' catout <- cfa.mi(' f =~ u1 + u2 + u3 + u4 ', data = datCat,
#'                  m = 3, seed = 456,
#'                  miArgs = list(ords = paste0("u", 1:8), noms = "g"),
#'                  FUN = function(fit) {
#'                    list(wrmr = lavaan::fitMeasures(fit, "wrmr"),
#'                         zeroCells = lavaan::lavInspect(fit, "zero.cell.tables"))
#'                  })
#' summary(catout)
#' anova(catout, indices = "all") # note the scaled versions of indices, too
#'
#' ## extract custom output
#' sapply(catout@funList, function(x) x$wrmr) # WRMR for each imputation
#' catout@funList[[1]]$zeroCells # zero-cell tables for first imputation
#' catout@funList[[2]]$zeroCells # zero-cell tables for second imputation ...
#'
#' }
#'
#' @export
runMI <- function(model, data, fun = "lavaan", ...,
                  m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  CALL <- match.call()
  dots <- list(...)
  if (!is.null(dots$fixed.x)) {
    if (dots$fixed.x) warning('fixed.x set to FALSE')
  }
  if (!is.null(dots$conditional.x)) {
    if (dots$conditional.x) warning('conditional.x set to FALSE')
  }
  dots$fixed.x <- dots$conditional.x <- FALSE

  seed <- as.integer(seed[1])
  ## Create (or acknowledge) list of imputed data sets
  imputedData <- NULL
  if (is.data.frame(data)) {
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
        imputedData[[i]] <- mice::complete(x = miceOut, action = i, include = FALSE)
      }
    } else stop("Currently runMI only supports imputation by Amelia or mice")
  } else if (is.list(data)) {
    seed <- integer(length = 0)
    imputeCall <- list()
    imputedData <- data
    m <- length(data)
    class(imputedData) <- "list" # override inheritance (e.g., "mi" if Amelia)
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
      if (lavaan::lavInspect(obj, "ngroups") == 1L) {
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
         modindices = try(lavaan::modindices(obj), silent = TRUE),
         GLIST = obj@Model@GLIST, # FIXME: @Model slot may disappear; need GLIST for std.all and in simsem
         converged = converged, SE = se.test,
         Heywood.lv = Heywood.lv, Heywood.ov = Heywood.ov)
  }

  ## fit model using lavaanList
  lavListCall <- list(lavaan::lavaanList, model = model, dataList = imputedData,
                      cmd = fun)
  lavListCall <- c(lavListCall, dots)
  lavListCall$store.slots <- c("partable","vcov","test")
  lavListCall$FUN <- if (is.null(dots$FUN)) .getOutput. else function(obj) {
    temp1 <- .getOutput.(obj)
    temp2 <- dots$FUN(obj)
    if (!is.list(temp2)) temp2 <- list(userFUN1 = temp2)
    if (is.null(names(temp2))) names(temp2) <- paste0("userFUN", 1:length(temp2))
    duplicatedNames <- which(sapply(names(temp2), function(x) {
      x %in% c("sampstat","coefMats","modindices","converged",
               "SE","Heywood.lv","Heywood.ov","GLIST")
    }))
    for (i in duplicatedNames) names(temp2)[i] <- paste0("userFUN", i)
    c(temp1, temp2)
  }
  fit <- eval(as.call(lavListCall))
  ## Store custom @DataList and @SampleStatsList
  fit@SampleStatsList <- lapply(fit@funList, "[[", i = "sampstat")
  fit@DataList <- imputedData
  ## assign class and add new slots
  fit <- as(fit, "lavaan.mi")
  fit@coefList <- lapply(fit@funList, "[[", i = "coefMats")
  fit@miList <- lapply(fit@funList, "[[", i = "modindices")
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
  if (length(conv)) {
    firstConv <- conv[1]
    fit@GLIST <- list()
    ## loop over GLIST elements
    for (mat in seq_along(fit@funList[[firstConv]][["GLIST"]])) {
      matList <- lapply(fit@funList[conv], function(x) x$GLIST[[mat]])
      fit@GLIST[[mat]] <- Reduce("+", matList) / length(matList)
    }
    names(fit@GLIST) <- names(fit@funList[[firstConv]][["GLIST"]])
  } else {
    fit@GLIST <- list()
    warning('The model did not converge for any imputed data sets.')
  }

  ## keep any remaining funList slots (if allowing users to supply custom FUN)
  funNames <- names(fit@funList[[1]])
  keepIndex <- which(!sapply(funNames, function(x) {
    x %in% c("sampstat","coefMats","modindices","converged",
             "SE","Heywood.lv","Heywood.ov","GLIST")
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

  fit@ParTable$start <- getMethod("coef", "lavaan.mi")(fit, type = "user", labels = FALSE)
  fit
}

#' @rdname runMI
#' @export
lavaan.mi <- function(model, data, ...,
                      m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  runMI(model = model, data = data, fun = "lavaan", ...,
        m = m, miArgs = miArgs, miPackage = miPackage, seed = seed)
}

#' @rdname runMI
#' @export
cfa.mi <- function(model, data, ...,
                   m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  runMI(model = model, data = data, fun = "cfa", ...,
        m = m, miArgs = miArgs, miPackage = miPackage, seed = seed)
}

#' @rdname runMI
#' @export
sem.mi <- function(model, data, ...,
                   m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  runMI(model = model, data = data, fun = "sem", ...,
        m = m, miArgs = miArgs, miPackage = miPackage, seed = seed)
}

#' @rdname runMI
#' @export
growth.mi <- function(model, data, ...,
                      m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  runMI(model = model, data = data, fun = "growth", ...,
        m = m, miArgs = miArgs, miPackage = miPackage, seed = seed)
}

