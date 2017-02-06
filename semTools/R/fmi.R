### Mauricio Garnier Villarreal & Terrence D. Jorgensen
### Last updated: 14 January 2017
### This function estimates the Fraction of Missing Information for means and
### (co)variances of each variable in a partially observed data set or from
### a list of multiple imputed data sets

## data = either a single incomplete data set of a list of imputed data sets
## method = the model used for the estimation
## group = the name of a grouping variable
## ords = ordinal variable names
## varnames = variable names to include in the analysis
## exclude = variable names to exclude from the analysis
## fewImps = If TRUE, apply small-M correction to between-imps variance
fmi <- function(data, method = "saturated", group = NULL, ords = NULL,
                varnames = NULL, exclude = NULL, fewImps = TRUE) {
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
                                   "FIML option is available to calculate FMI, ",
                                   "which requires continuous variables. The ",
                                   "following variables were removed: ",
                                   paste(ordvars, collapse = ", ")))
    if (!length(numvars)) stop("No numeric variables were provided.")
    vars <- numvars
  }

  ## construct model
  covstruc <- outer(vars, vars, function(x, y) paste(x, "~~", y))
  if (method == "saturated" | method == "sat") {
    model <- covstruc[lower.tri(covstruc, diag = TRUE)]
  } else if(method == "null") model <- diag(covstruc)
  if (length(numvars)) model <- c(model, paste(numvars, "~1"))

  ## fit model
  if (fiml) {
    fit <- lavaan::lavaan(model, data = data, missing = "fiml", group = group)
    comb.results <- lavaan::parameterEstimates(fit, fmi = TRUE, zstat = FALSE,
                                               pvalue = FALSE, ci = FALSE)
    nG <- lavaan::lavInspect(fit, "ngroups")
    if(nG == 1L) comb.results$group <- 1L
    group.label <- lavaan::lavInspect(fit, "group.label")
  } else {
    fit <- lavaan.mi(model, data, group = group, ordered = ordvars, auto.th = TRUE)
    comb.results <- getMethod("summary","lavaan.mi")(fit, fmi = TRUE, ci = FALSE,
                                                     add.attributes = FALSE)
    nG <- lavaan::lavListInspect(fit, "ngroups")
    if(nG == 1L) comb.results$group <- 1L
    group.label <- lavaan::lavListInspect(fit, "group.label")
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
      covmat <- lavaan::lavInspect(fit, "theta")
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


