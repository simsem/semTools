### Mauricio Garnier Villarreal & Terrence D. Jorgensen
### Last updated: 5 July 2019
### This function estimates the Fraction of Missing Information for means and
### (co)variances of each variable in a partially observed data set or from
### a list of multiple imputed data sets

#' Fraction of Missing Information.
#'
#' This function estimates the Fraction of Missing Information (FMI) for
#' summary statistics of each variable, using either an incomplete data set or
#' a list of imputed data sets.
#'
#' The function estimates a saturated model with \code{\link[lavaan]{lavaan}}
#' for a single incomplete data set using FIML, or with \code{\link{lavaan.mi}}
#' for a list of imputed data sets. If method = \code{"saturated"}, FMI will be
#' estiamted for all summary statistics, which could take a lot of time with
#' big data sets. If method = \code{"null"}, FMI will only be estimated for
#' univariate statistics (e.g., means, variances, thresholds). The saturated
#' model gives more reliable estimates, so it could also help to request a
#' subset of variables from a large data set.
#'
#'
#' @importFrom lavaan lavListInspect lavInspect
#'
#' @param data Either a single \code{data.frame} with incomplete observations,
#' or a \code{list} of imputed data sets.
#' @param method character. If \code{"saturated"} or \code{"sat"} (default),
#' the model used to estimate FMI is a freely estimated covariance matrix and
#' mean vector for numeric variables, and/or polychoric correlations and
#' thresholds for ordered categorical variables, for each group (if
#' applicable). If \code{"null"}, only means and variances are estimated for
#' numeric variables, and/or thresholds for ordered categorical variables
#' (i.e., covariances and/or polychoric correlations are constrained to zero).
#' See Details for more information.
#' @param group character. The optional name of a grouping variable, to request
#' FMI in each group.
#' @param ords character. Optional vector of names of ordered-categorical
#' variables, which are not already stored as class \code{ordered} in
#' \code{data}.
#' @param varnames character. Optional vector of variable names, to calculate
#' FMI for a subset of variables in \code{data}. By default, all numeric and
#' ordered variables will be included, unless \code{data} is a single
#' incomplete \code{data.frame}, in which case only numeric variables can be
#' used with FIML estimation. Other variable types will be removed.
#' @param exclude character. Optional vector of variable names to exclude from
#' the analysis.
#' @param fewImps logical. If \code{TRUE}, use the estimate of FMI that applies
#' a correction to the estimated between-imputation variance. Recommended when
#' there are few imputations; makes little difference when there are many
#' imputations. Ignored when \code{data} is not a list of imputed data sets.
#' @return \code{fmi} returns a list with at least 2 of the following:
#' \item{Covariances}{A list of symmetric matrices: (1) the estimated/pooled
#' covariance matrix, or a list of group-specific matrices (if applicable) and
#' (2) a matrix of FMI, or a list of group-specific matrices (if applicable).
#' Only available if \code{method = "saturated"}.} \item{Variances}{The
#' estimated/pooled variance for each numeric variable. Only available if
#' \code{method = "null"} (otherwise, it is on the diagonal of Covariances).}
#' \item{Means}{The estimated/pooled mean for each numeric variable.}
#' \item{Thresholds}{The estimated/pooled threshold(s) for each
#' ordered-categorical variable.} \item{message}{A message indicating caution
#' when the null model is used.}
#' @author Mauricio Garnier Villarreal (University of Kansas;
#' \email{mauricio.garniervillarreal@@marquette.edu}) Terrence Jorgensen
#' (University of Amsterdam; \email{TJorgensen314@@gmail.com})
#' @references Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse
#' in surveys}. New York, NY: Wiley.
#'
#' Savalei, V. & Rhemtulla, M. (2012). On obtaining estimates of the fraction
#' of missing information from full information maximum likelihood.
#' \emph{Structural Equation Modeling, 19}(3), 477--494.
#' doi:10.1080/10705511.2012.687669
#'
#' Wagner, J. (2010). The fraction of missing information as a tool for
#' monitoring the quality of survey data. \emph{Public Opinion Quarterly,
#' 74}(2), 223--243. doi:10.1093/poq/nfq007
#' @examples
#'
#' HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
#'                                       "ageyr","agemo","school")]
#' set.seed(12345)
#' HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
#' age <- HSMiss$ageyr + HSMiss$agemo/12
#' HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
#'
#' ## calculate FMI (using FIML, provide partially observed data set)
#' (out1 <- fmi(HSMiss, exclude = "school"))
#' (out2 <- fmi(HSMiss, exclude = "school", method = "null"))
#' (out3 <- fmi(HSMiss, varnames = c("x5","x6","x7","x8","x9")))
#' (out4 <- fmi(HSMiss, group = "school"))
#'
#' \dontrun{
#' ## ordered-categorical data
#' data(datCat)
#' lapply(datCat, class)
#' ## impose missing values
#' set.seed(123)
#' for (i in 1:8) datCat[sample(1:nrow(datCat), size = .1*nrow(datCat)), i] <- NA
#' ## impute data m = 3 times
#' library(Amelia)
#' set.seed(456)
#' impout <- amelia(datCat, m = 3, noms = "g", ords = paste0("u", 1:8), p2s = FALSE)
#' imps <- impout$imputations
#' ## calculate FMI, using list of imputed data sets
#' fmi(imps, group = "g")
#' }
#'
#' @export
fmi <- function(data, method = "saturated", group = NULL, ords = NULL,
                varnames = NULL, exclude = NULL, fewImps = FALSE) {
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
    if (length(ordvars)) message(c("By providing a single data set, only the ",
                                   "FIML option is available to calculate FMI,",
                                   " which requires continuous variables. The ",
                                   "following variables were removed: ",
                                   paste(ordvars, collapse = ", ")))
    if (!length(numvars)) stop("No numeric variables were provided.")
    vars <- numvars
  }

  ## construct model
  covstruc <- outer(vars, vars, function(x, y) paste(x, "~~", y))
  if (method == "saturated" | method == "sat") {
    diag(covstruc)[which(ordvars %in% vars)] <- ""
    model <- covstruc[lower.tri(covstruc, diag = TRUE)]
  } else if (method == "null") model <- diag(covstruc)
  if (length(numvars)) model <- c(model, paste(numvars, "~1"))

  ## fit model
  if (fiml) {
    fit <- lavaan::lavaan(model, data = data, missing = "fiml", group = group)
    comb.results <- lavaan::parameterEstimates(fit, fmi = TRUE, zstat = FALSE,
                                               pvalue = FALSE, ci = FALSE)
    nG <- lavInspect(fit, "ngroups")
    if (nG == 1L) comb.results$group <- 1L
    group.label <- lavInspect(fit, "group.label")
  } else {
    fit <- lavaan.mi(model, data, group = group, ordered = ordvars, auto.th = TRUE)
    comb.results <- getMethod("summary","lavaan.mi")(fit, fmi = TRUE, ci = FALSE,
                                                     output = "data.frame")
    nG <- lavListInspect(fit, "ngroups")
    if (nG == 1L) comb.results$group <- 1L
    group.label <- lavListInspect(fit, "group.label")
    if (fewImps) {
      comb.results["fmi1"] <- NULL
      names(comb.results)[names(comb.results) == "fmi2"] <- "fmi"
    } else {
      comb.results["fmi2"] <- NULL
      names(comb.results)[names(comb.results) == "fmi1"] <- "fmi"
    }
    for (i in c("t","df","pvalue","riv")) comb.results[i] <- NULL
  }

  ## Variances from null model, if applicable
  if (method == "null") {
    if (length(numvars)) {
      Variances <- comb.results[comb.results$lhs == comb.results$rhs,
                                c("lhs","group","est","fmi")]
      colnames(Variances)[c(1, 3)] <- c("variable","coef")
      if (nG > 1L) Variances$group <- group.label[Variances$group]
      class(Variances) <- c("lavaan.data.frame","data.frame")
      ## start list of results
      results <- list(Variances = Variances)
    } else results <- list()
  } else {
    ## covariances from saturated model, including polychorics (if applicable)
    if (fiml) {
      covmat <- lavInspect(fit, "theta")
      if (nG == 1L) covmat <- list(covmat)
    } else {
      useImps <- sapply(fit@convergence, "[[", "converged")
      m <- sum(useImps)
      if (nG == 1L) {
        ThetaList <- lapply(fit@coefList[useImps], function(x) x$theta)
        covmat <- list(Reduce("+", ThetaList) / m)
      } else {
        covmat <- list()
        for (i in group.label) {
          groupList <- lapply(fit@coefList[useImps],"[[", i)
          ThetaList <- lapply(groupList, function(x) x$theta)
          covmat[[i]] <- Reduce("+", ThetaList) / m
        }
      }
    }

    fmimat <- covmat
    covars <- comb.results[comb.results$op == "~~", c("lhs","rhs","group","est","fmi")]
    for (i in 1:nG) {
      fmimat[[i]][as.matrix(covars[covars$group == i, 1:2])] <- covars$fmi[covars$group == i]
      fmimat[[i]][as.matrix(covars[covars$group == i, 2:1])] <- covars$fmi[covars$group == i]
    }
    if (nG == 1L) {
      Covariances <- list(coef = covmat[[1]], fmi = fmimat[[1]])
    } else Covariances <- list(coef = covmat, fmi = fmimat)
    ## start list of results
    results <- list(Covariances = Covariances)
  }

  ## Means, if applicable
  if (length(numvars)) {
    results$Means <- comb.results[comb.results$op == "~1" & comb.results$lhs %in% numvars,
                                  c("lhs","group","est","fmi")]
    colnames(results$Means)[c(1, 3)] <- c("variable","coef")
    if (nG > 1L) results$Means$group <- group.label[results$Means$group]
    class(results$Means) <- c("lavaan.data.frame","data.frame")
  }
  ## Thresholds, if  applicable
  if (length(ordvars)) {
    results$Thresholds <- comb.results[comb.results$op == "|",
                                  c("lhs","rhs","group","est","fmi")]
    colnames(results$Thresholds)[c(1, 2, 4)] <- c("variable","threshold","coef")
    if (nG > 1L) results$Thresholds$group <- group.label[results$Thresholds$group]
    class(results$Thresholds) <- c("lavaan.data.frame","data.frame")
  }

  ## return results, with message if using null model
  if (method == "null")
    results$message <- paste("Null-model estimates may not be as",
                             "precise as saturated-model estimates.")
  results
}


