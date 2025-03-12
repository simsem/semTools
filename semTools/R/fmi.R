### Mauricio Garnier Villarreal & Terrence D. Jorgensen
### Last updated: 12 March 2025
### This function estimates the Fraction of Missing Information for means and
### (co)variances of each variable in a partially observed data set or from
### a list of multiple imputed data sets

##' Fraction of Missing Information.
##'
##' This function estimates the Fraction of Missing Information (FMI) for
##' summary statistics of each variable, using either an incomplete data set or
##' a list of imputed data sets.
##'
##' The function estimates a saturated model with [lavaan::lavaan()] for a
##' single incomplete data set using FIML, or with [lavaan.mi::lavaan.mi()]
##' for a list of imputed data sets. If method = `"saturated"`, FMI will be
##' estiamted for all summary statistics, which could take a lot of time with
##' big data sets. If method = `"null"`, FMI will only be estimated for
##' univariate statistics (e.g., means, variances, thresholds). The saturated
##' model gives more reliable estimates, so it could also help to request a
##' subset of variables from a large data set.
##'
##'
##' @importFrom lavaan lavListInspect lavInspect
##'
##' @param data Either a single `data.frame` with incomplete observations,
##'   or a `list` of imputed data sets.
##' @param method character. If `"saturated"` or `"sat"` (default),
##'   the model used to estimate FMI is a freely estimated covariance matrix and
##'   mean vector for numeric variables, and/or polychoric correlations and
##'   thresholds for ordered categorical variables, for each group (if
##'   applicable). If `"null"`, only means and variances are estimated for
##'   numeric variables, and/or thresholds for ordered categorical variables
##'   (i.e., covariances and/or polychoric/polyserial correlations are
##'   constrained to zero). See **Details** for more information.
##' @param group `character`. The optional name of a grouping variable, to
##'   request FMI in each group.
##' @param ords Optional `character` vector naming ordered-categorical
##'   variables, if they are not already stored as class `ordered` in `data`.
##' @param varnames Optional `character` vector of variable names, to calculate
##'   FMI for a subset of variables in `data`. By default, all numeric and
##'   `ordered=` variables will be included, unless `data=` is a single
##'   incomplete `data.frame`, in which case only numeric variables can be
##'   used with FIML estimation. Other variable types will be removed.
##' @param exclude Optional `character` vector naming variables to exclude from
##'   the analysis.
##' @param return.fit logical. If `TRUE`, the fitted [lavaan::lavaan-class] or
##'   [lavaan.mi::lavaan.mi-class] model is returned, so FMI can be found from
##'   `summary(..., fmi=TRUE)`.
##'
##' @return `fmi()` returns a list with at least 2 of the following:
##'
##' \item{Covariances}{A list of symmetric matrices: (1) the estimated/pooled
##'   covariance matrix, or a list of group-specific matrices (if applicable)
##'   and (2) a matrix of FMI, or a list of group-specific matrices (if
##'   applicable). Only available if `method = "saturated"`.  When
##'   `method="cor"`, this element is replaced by `Correlations`.}
##' \item{Variances}{The estimated/pooled variance for each numeric variable.
##'   Only available if `method = "null"` (otherwise, it is on the diagonal
##'   of Covariances).}
##' \item{Means}{The estimated/pooled mean for each numeric variable.}
##' \item{Thresholds}{The estimated/pooled threshold(s) for each
##'   ordered-categorical variable.}
##'
##' @author
##' Mauricio Garnier Villarreal (Vrije Universiteit Amsterdam; \email{m.garniervillarreal@@vu.nl})
##'
##' Terrence Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##' Rubin, D. B. (1987). *Multiple imputation for nonresponse in surveys*.
##' New York, NY: Wiley.
##'
##' Savalei, V. & Rhemtulla, M. (2012). On obtaining estimates of the fraction
##' of missing information from full information maximum likelihood.
##' *Structural Equation Modeling, 19*(3), 477--494.
##' \doi{10.1080/10705511.2012.687669}
##'
##' Wagner, J. (2010). The fraction of missing information as a tool for
##' monitoring the quality of survey data. *Public Opinion Quarterly,
##' 74*(2), 223--243. \doi{10.1093/poq/nfq007}
##'
##' @examples
##'
##' HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
##'                                       "ageyr","agemo","school")]
##' set.seed(12345)
##' HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
##' age <- HSMiss$ageyr + HSMiss$agemo/12
##' HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
##'
##' ## calculate FMI (using FIML, provide partially observed data set)
##' (out1 <- fmi(HSMiss, exclude = "school"))
##' (out2 <- fmi(HSMiss, exclude = "school", method = "null"))
##' (out3 <- fmi(HSMiss, varnames = c("x5","x6","x7","x8","x9")))
##' (out4 <- fmi(HSMiss, method = "cor", group = "school")) # correlations by group
##'
##' ## significance tests in lavaan(.mi) object
##' out5 <- fmi(HSMiss, method = "cor", return.fit = TRUE)
##' summary(out5) # factor loading == SD, covariance = correlation
##'
##' if(requireNamespace("lavaan.mi")){
##'   ## ordered-categorical data
##'   data(binHS5imps, package = "lavaan.mi")
##'
##'   ## calculate FMI, using list of imputed data sets
##'   fmi(binHS5imps, group = "school")
##' }
##'
##' @export
fmi <- function(data, method = "saturated", group = NULL, ords = NULL,
                varnames = NULL, exclude = NULL, return.fit = FALSE) {
  fiml <- is.data.frame(data)
  ## check for single data set or list of imputed data sets
  data1 <- if (fiml) data else data[[1]]
  ## select user-specified variables
  vars <- if (is.character(varnames)) varnames else colnames(data1)
  ## remove grouping variable and user-specified exclusions, if applicable
  vars <- setdiff(vars, c(group, exclude))
  ## check classes
  ordvars <- vars[sapply(data1[vars], is.ordered)]
  if (!is.null(ords)) ordvars <- c(ordvars, ords)
  numvars <- vars[sapply(data1[vars], is.numeric)]
  vars <- union(numvars, ordvars)
  numvars <- setdiff(vars, ordvars)
  if (fiml) {
    #TODO: enable estimator = "PML"
    #      pass missing= option as another fmi() argument?
    if (length(ordvars)) message(c("By providing a single data set, only the ",
                                   "FIML option is available to calculate FMI,",
                                   " which requires continuous variables. The ",
                                   "following variables were removed: ",
                                   paste(ordvars, collapse = ", ")))
    if (!length(numvars)) stop("No numeric variables were provided.")
    vars <- numvars
  }

  ## construct model
  if (method == "saturated" | method == "sat") {
    covstruc <- outer(vars, vars, function(x, y) paste(x, "~~", y))
    diag(covstruc)[which(vars %in% ordvars)] <- ""
    model <- covstruc[lower.tri(covstruc, diag = TRUE)]

  } else if (method == "null") {
    covstruc <- outer(vars, vars, function(x, y) paste(x, "~~", y))
    diag(covstruc)[which(vars %in% ordvars)] <- ""
    model <- diag(covstruc)

  } else if (method == "cor") {
    # phantoms <- paste0(".", vars, ".")
    model <- c(paste0(".", vars, ". =~ ", ifelse(vars %in% ordvars, "1*", "NA*"),
                      vars), # loadings = SDs (fixed to 1 when ordinal)
               paste0(vars, " ~~ 0*", vars))     # "residual" variances = 0

  } else stop('Invalid method= argument: "', method, '"')
  if (length(numvars)) model <- c(model, paste(numvars, "~ 1"))

  ## fit model
  if (fiml) {
    if (method == "cor") {
      fit <- lavaan::cfa(model, data = data, missing = "fiml", group = group,
                         std.lv = TRUE)
    } else {
      fit <- lavaan::lavaan(model, data = data, missing = "fiml", group = group)
    }
    if (return.fit) return(fit)

    comb.results <- lavaan::parameterEstimates(fit, fmi = TRUE, zstat = FALSE,
                                               pvalue = FALSE, ci = FALSE)
    nG <- lavInspect(fit, "ngroups")
    if (nG == 1L) comb.results$group <- 1L
    group.label <- lavInspect(fit, "group.label")


  } else {
    ## list of imputations
    if (!"package:lavaan.mi" %in% search()) attachNamespace("lavaan.mi")

    if (method == "cor") {
      fit <- lavaan.mi::cfa.mi(model, data = data, group = group,
                               ordered = ordvars, std.lv = TRUE)
    } else {
      fit <- lavaan.mi::lavaan.mi(model, data = data, group = group,
                                  ordered = ordvars, auto.th = TRUE)
    }
    if (return.fit) return(fit)

    comb.results <- lavaan.mi::parameterEstimates.mi(fit, fmi = TRUE,
                                                     zstat = FALSE,
                                                     pvalue = FALSE, ci = FALSE)
    nG <- lavListInspect(fit, "ngroups")
    if (nG == 1L) comb.results$group <- 1L
    group.label <- lavListInspect(fit, "group.label")
    #FIXME: also return RIV?  Or make it an argument (FMI or RIV)
    comb.results$riv <- NULL
  }

  ## Variances from null model, if applicable
  if (method == "null") {
    if (length(numvars)) {
      Variances <- comb.results[comb.results$lhs == comb.results$rhs,
                                c("lhs","group","est","fmi")]
      colnames(Variances)[c(1, 3)] <- c("variable","coef")
      if (nG > 1L) Variances$group <- group.label[Variances$group]
      class(Variances) <- c("lavaan.data.frame","data.frame")
      attr(Variances, "header") <- paste("Null-model estimates may not be as",
                                         "accurate as saturated-model estimates.")
      ## start list of results
      results <- list(Variances = Variances)
    } else results <- list()
  } else {
    ## covariances from saturated model, including polychorics (if applicable)
    if (fiml) {
      if (nG == 1L) {
        covmat <- lavInspect(fit, "est")[[ifelse(method == "cor", "psi", "theta")]]
        covmat <- list(covmat)
      } else {
        covmat <- sapply(lavInspect(fit, "est"),    "[[",
                         i = ifelse(method == "cor", "psi", "theta"),
                         simplify = FALSE)
      }

    } else {
      useImps <- sapply(fit@convergence, "[[", "converged")
      m <- sum(useImps)
      if (nG == 1L) {
        CovList <- lapply(fit@coefList[useImps],
                            function(x) x[[ifelse(method == "cor", "psi", "theta")]])
        covmat <- list(Reduce("+", CovList) / m)
      } else {
        covmat <- list()
        for (i in group.label) {
          groupList <- lapply(fit@coefList[useImps],"[[", i)
          CovList <- lapply(groupList, function(x) x[[ifelse(method == "cor", "psi", "theta")]])
          covmat[[i]] <- Reduce("+", CovList) / m
        }
      }
    }

    fmimat <- covmat
    covars <- comb.results[comb.results$op == "~~", c("lhs","rhs","group","est","fmi")]
    for (i in 1:nG) {
      theseCovars <- covars[covars$group == i,]
      if (method == "cor") {
        phantomRows <- !(theseCovars$lhs %in% vars)
        theseCovars <- theseCovars[phantomRows,]
      }
      fmimat[[i]][as.matrix(theseCovars[, 1:2])] <- theseCovars$fmi
      fmimat[[i]][as.matrix(theseCovars[, 2:1])] <- theseCovars$fmi
      if (method == "cor") {
        ## reset variable names
        dimnames(fmimat[[i]]) <- dimnames(covmat[[i]]) <- list(vars, vars)
      }
    }
    if (nG == 1L) {
      Covariances <- list(coef = covmat[[1]], fmi = fmimat[[1]])
    } else Covariances <- list(coef = covmat, fmi = fmimat)
    ## start list of results
    results <- setNames(list(Covariances),
                        nm = ifelse(method == "cor", "Correlations", "Covariances"))
  }

  ## Means, if applicable
  if (length(numvars)) {
    results$Means <- comb.results[comb.results$op == "~1" & comb.results$lhs %in% numvars,
                                  c("lhs","group","est","fmi")]
    colnames(results$Means)[c(1, 3)] <- c("variable","coef")
    if (nG > 1L) results$Means$group <- group.label[results$Means$group]
    class(results$Means) <- c("lavaan.data.frame","data.frame")
  }
  ## Thresholds, if applicable
  if (length(ordvars)) {
    results$Thresholds <- comb.results[comb.results$op == "|",
                                  c("lhs","rhs","group","est","fmi")]
    colnames(results$Thresholds)[c(1, 2, 4)] <- c("variable","threshold","coef")
    if (nG > 1L) results$Thresholds$group <- group.label[results$Thresholds$group]
    class(results$Thresholds) <- c("lavaan.data.frame","data.frame")
  }

  results
}


