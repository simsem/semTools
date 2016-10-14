### Terrence D. Jorgensen
### Last updated: 24 April 2016
### permutation randomization test for measurement equivalence and DIF


## create s4 class for result object
setClass("permuteMeasEq", slots = c(PT = "data.frame",
                                    modelType = "character",
                                    ANOVA = "vector",
                                    AFI.obs = "vector",
                                    AFI.dist = "data.frame",
                                    AFI.pval = "vector",
                                    MI.obs = "data.frame",
                                    MI.dist = "vector",
                                    extra.obs = "vector",
                                    extra.dist = "data.frame",
                                    n.Permutations = "integer",
                                    n.Converged = "integer",
                                    n.nonConverged = "vector",
                                    n.Sparse = "vector",
                                    oldSeed = "integer"))

## function to check validity of arguments to permuteMeasEq()
checkPermArgs <- function(nPermute, modelType, con, uncon, null,
                          param, freeParam, covariates, AFIs, moreAFIs,
                          maxSparse, maxNonconv, showProgress, warn,
                          datafun, extra, parallelType, ncpus, cl, iseed) {
  fixedCall <- as.list(match.call())[-1]

  fixedCall$nPermute <- as.integer(nPermute[1])
  fixedCall$modelType <- modelType[1]
  if (!fixedCall$modelType %in% c("mgcfa","mimic","long"))
    stop('modelType must be one of c("mgcfa","mimic","long")')
  if (fixedCall$modelType == "long") stop('modelType "long" is not yet available.')
  if (fixedCall$modelType == "mgcfa" && lavaan::lavInspect(con, "ngroups") == 1L)
    stop('modelType = "mgcfa" applies only to multigroup models.')
  if (fixedCall$modelType == "mimic") {
    uncon <- NULL
    fixedCall$uncon <- NULL
    fixedCall <- c(fixedCall, list(uncon = NULL))
  }
  ## strip white space
  if (is.list(param)) {
    fixedCall$param <- lapply(param, function(cc) gsub("[[:space:]]+", "", cc))
  } else if (!is.null(param)) fixedCall$param <- gsub("[[:space:]]+", "", param)
  if (!is.null(freeParam)) fixedCall$freeParam <- gsub("[[:space:]]+", "", freeParam)
  if (fixedCall$modelType == "mimic") {
    # PT <- lavaan::lavaanify(fixedCall$param)
    # checkCovs <- unique(PT$rhs[PT$op == "~"])
    # if (is.null(covariates)) covariates <- checkCovs
    # if (length(setdiff(covariates, checkCovs)))
    #   warning('Argument "covariates" includes predictors not in argument "param"')
    ##### ordVars <- lavaan::lavNames(con, type = "ov.ord")
    fixedCall$covariates <- as.character(covariates)
  }
  fixedCall$maxSparse <- as.integer(maxSparse[1])
  fixedCall$maxNonconv <- as.integer(maxNonconv[1])
  fixedCall$showProgress <- as.logical(showProgress[1])
  fixedCall$warn <- as.integer(warn[1])
  fixedCall$oldSeed <- as.integer(NULL)
  parallelType <- as.character(parallelType[1])
  if (!parallelType %in% c("none","multicore","snow")) parallelType <- "none"
  if (!is.null(cl)) {
    if (!is(cl, "cluster")) stop("Invalid cluster object.  Check class(cl)")
    parallelType <- "snow"
    ncpus <- length(cl)
  }
  if (parallelType == "multicore" && .Platform$OS.type == "windows") {
    parallelType <- "snow"
    message("'multicore' option unavailable on Windows. Using 'snow' instead.")
  }
  ## parallel settings, adapted from boot::boot()
  if (parallelType != "none") {
    if (is.null(ncpus) || ncpus > parallel::detectCores()) {
      ncpus <- parallel::detectCores() - 1
    }
    if (ncpus <= 1L) {
      parallelType <- "none"
    } else {
      fixedCall$showProgress <- FALSE
      fixedCall$old_RNG <- RNGkind()
      fixedCall$oldSeed <- .Random.seed
      if (fixedCall$old_RNG[1] != "L'Ecuyer-CMRG") {
        RNGkind("L'Ecuyer-CMRG")
        message("Your RNGkind() was changed from ", fixedCall$old_RNG[1],
                " to L'Ecuyer-CMRG, which is required for reproducibility ",
                " in parallel jobs.  Your RNGkind() has been returned to ",
                fixedCall$old_RNG[1], " but the seed has not been set. ",
                " The state of your previous RNG is saved in the slot ",
                " named 'oldSeed', if you want to restore it using ",
                " the syntax:\n",
                ".Random.seed[-1] <- permuteMeasEqObject@oldSeed[-1]")
      }
      fixedCall$iseed <- as.integer(iseed[1])
      if (is.na(fixedCall$iseed)) fixedCall$iseed <- 12345
    }
  }
  fixedCall$parallelType <- parallelType
  if (is.null(ncpus)) {
    fixedCall$ncpus <- NULL
    fixedCall <- c(fixedCall, list(ncpus = NULL))
  } else fixedCall$ncpus <- ncpus

  ## check that "param" is NULL if uncon is NULL, and check for lavaan class
  notLavaan <- "Non-NULL 'con', 'uncon', or 'null' must be fitted lavaan object."
  if (is.null(uncon)) {
    if (!is.null(fixedCall$param) && fixedCall$modelType == "mgcfa") {
      message(c(" When 'uncon = NULL', only configural invariance is tested.",
                "\n So the 'param' argument was changed to NULL."))
      fixedCall$param <- NULL
      fixedCall <- c(fixedCall, list(param = NULL))
    }
    if (class(con) != "lavaan") stop(notLavaan)
  } else {
    if (class(con) != "lavaan") stop(notLavaan)
    if (class(uncon) != "lavaan") stop(notLavaan)
  }
  if (!is.null(null)) {
    if (class(null) != "lavaan") stop(notLavaan)
  }

  ############ FIXME: check that lavaan::lavInspect(con, "options")$conditional.x = FALSE (find defaults for continuous/ordered indicators)
  if (!is.null(fixedCall$param)) {
    ## Temporarily warn about testing thresholds without necessary constraints.   FIXME: check for binary indicators
    if ("thresholds" %in% fixedCall$param | any(grepl("\\|", fixedCall$param))) {
      warning(c("This function is not yet optimized for testing thresholds.\n",
                "Necessary identification contraints might not be specified."))
    }
    ## collect parameter types for "mgcfa"
    if (fixedCall$modelType != "mimic") {
      ## save all estimates from constrained model
      PT <- lavaan::parTable(con)[ , c("lhs","op","rhs","group","plabel")]
      ## extract parameters of interest
      paramTypes <- c("loadings","intercepts","thresholds","residuals","means",
                      "residual.covariances","lv.variances","lv.covariances")
      params <- PT[paste0(PT$lhs, PT$op, PT$rhs) %in% setdiff(fixedCall$param,
                                                              paramTypes), ]
      ## add parameters by type, if any are specified
      types <- intersect(fixedCall$param, paramTypes)
      ov.names <- lavaan::lavNames(con, "ov")
      isOV <- PT$lhs %in% ov.names
      lv.names <- con@pta$vnames$lv[[1]]
      isLV <- PT$lhs %in% lv.names & PT$rhs %in% lv.names
      if ("loadings" %in% types) params <- rbind(params, PT[PT$op == "=~", ])
      if ("intercepts" %in% types) {
        params <- rbind(params, PT[isOV & PT$op == "~1", ])
      }
      if ("thresholds" %in% types) params <- rbind(params, PT[PT$op == "|", ])
      if ("residuals" %in% types) {
        params <- rbind(params, PT[isOV & PT$lhs == PT$rhs & PT$op == "~~", ])
      }
      if ("residual.covariances" %in% types) {
        params <- rbind(params, PT[isOV & PT$lhs != PT$rhs & PT$op == "~~", ])
      }
      if ("means" %in% types) {
        params <- rbind(params, PT[PT$lhs %in% lv.names & PT$op == "~1", ])
      }
      if ("lv.variances" %in% types) {
        params <- rbind(params, PT[isLV & PT$lhs == PT$rhs & PT$op == "~~", ])
      }
      if ("lv.covariances" %in% types) {
        params <- rbind(params, PT[isLV & PT$lhs != PT$rhs & PT$op == "~~", ])
      }
      ## remove parameters specified by "freeParam" argument
      params <- params[!paste0(params$lhs, params$op, params$rhs) %in% fixedCall$freeParam, ]
      fixedCall$param <- paste0(params$lhs, params$op, params$rhs)
    }
  }


  if (is.null(AFIs) & is.null(moreAFIs)) {
    message("No AFIs were selected, so only chi-squared will be permuted.\n")
    fixedCall$AFIs <- "chisq"
    AFIs <- "chisq"
  }
  if ("ecvi" %in% AFIs & lavaan::lavInspect(con, "ngroups") > 1L)
    stop("ECVI is not available for multigroup models.")

  ## check estimators
  leastSq <- grepl("LS", lavaan::lavInspect(con, "options")$estimator)
  if (!is.null(uncon)) {
    if (uncon@Options$estimator != lavaan::lavInspect(con, "options")$estimator)
      stop("Models must be fit using same estimator.")
  }
  if (!is.null(null)) {
    if (lavaan::lavInspect(null, "options")$estimator != lavaan::lavInspect(con, "options")$estimator)
      stop("Models must be fit using same estimator.")
  }

  ## check extra functions, if any
  restrictedArgs <- c("con","uncon","null","param","freeParam","covariates",
                      "AFIs","moreAFIs","maxSparse","maxNonconv","iseed")
  if (!missing(datafun)) {
    if (!is.function(datafun)) stop('Argument "datafun" must be a function.')
    extraArgs <- formals(datafun)
    if (!all(names(extraArgs) %in% c(restrictedArgs, "data")))
      stop('The user-supplied function "datafun" can only have any among the ',
           'following arguments:\n', paste(restrictedArgs, collapse = ", "))
  }
  if (!missing(extra)) {
    if (!is.function(extra)) stop('Argument "extra" must be a function.')
    extraArgs <- formals(extra)
    if (!all(names(extraArgs) %in% restrictedArgs))
      stop('The user-supplied function "extra" can only have any among the ',
           'following arguments:\n', paste(restrictedArgs, collapse = ", "))
  }

  ## return evaluated list of other arguments
  lapply(fixedCall, eval)
}


## function to extract fit measures
getAFIs <- function(...) {
  dots <- list(...)

  AFI1 <- list()
  AFI0 <- list()
  leastSq <- grepl("LS", lavaan::lavInspect(dots$con, "options")$estimator)
  ## check validity of user-specified AFIs, save output
  if (!is.null(dots$AFIs)) {
    IC <- grep("ic|logl", dots$AFIs, value = TRUE)
    if (leastSq & length(IC)) {
      stop(paste("Argument 'AFIs' includes invalid options:",
                 paste(IC, collapse = ", "),
                 "Information criteria unavailable for least-squares estimators.",
                 sep = "\n"))
    }
    if (!is.null(dots$uncon))
      AFI1[[1]] <- lavaan::fitMeasures(dots$uncon, fit.measures = dots$AFIs,
                                       baseline.model = dots$null)
    AFI0[[1]] <- lavaan::fitMeasures(dots$con, fit.measures = dots$AFIs,
                                     baseline.model = dots$null)
  }
  ## check validity of user-specified moreAFIs
  if (!is.null(dots$moreAFIs)) {
    IC <- grep("ic|hqc", dots$moreAFIs, value = TRUE)
    if (leastSq & length(IC)) {
      stop(paste("Argument 'moreAFIs' includes invalid options:",
                 paste(IC, collapse = ", "),
                 "Information criteria unavailable for least-squares estimators.",
                 sep = "\n"))
    }
    if (!is.null(dots$uncon))
      AFI1[[2]] <- moreFitIndices(dots$uncon, fit.measures = dots$moreAFIs)
    AFI0[[2]] <- moreFitIndices(dots$con, fit.measures = dots$moreAFIs)
  }

  ## save observed AFIs or delta-AFIs
  if (is.null(dots$uncon)) {
    AFI.obs <- unlist(AFI0)
  } else {
    AFI.obs <- unlist(AFI0) - unlist(AFI1)
  }
  AFI.obs
}

## Function to extract modification indices for equality constraints
getMIs <- function(...) {
  dots <- list(...)

  if (dots$modelType == "mgcfa") {
    ## save all estimates from constrained model
    PT <- lavaan::parTable(dots$con)[ , c("lhs","op","rhs","group","plabel")]
    ## extract parameters of interest
    params <- PT[paste0(PT$lhs, PT$op, PT$rhs) %in% dots$param, ]
    ## return modification indices for specified constraints (param)
    MIs <- lavaan::lavTestScore(dots$con)$uni
    MI.obs <- MIs[MIs$lhs %in% params$plabel, ]
  } else if (dots$modelType == "mimic") {
    if (is.list(dots$param)) {
      MI <- lapply(dots$param, function(x) lavaan::lavTestScore(dots$con, add = x)$test)
      MI.obs <- do.call(rbind, MI)
    } else MI.obs <- lavaan::lavTestScore(dots$con, add = dots$param)$uni
  } else if (dots$modelType == "long") {
    ## coming soon
  }

  MI.obs
}

## Functions to find delta-AFIs & maximum modification index in one permutation
permuteOnce.mgcfa <- function(i, d, G, con, uncon, null, param, freeParam,
                              covariates, AFIs, moreAFIs, maxSparse, maxNonconv,
                              iseed, warn, extra = NULL, datafun = NULL) {
  old_warn <- options()$warn
  options(warn = warn)
  ## save arguments from call
  argNames <- names(formals(permuteOnce.mgcfa))
  availableArgs <- lapply(argNames, function(x) eval(as.name(x)))
  names(availableArgs) <- argNames

  nSparse <- 0L
  nTries <- 1L
  while ( (nSparse <= maxSparse) & (nTries <= maxNonconv) ) {
    ## permute grouping variable
    d[ , G] <- sample(d[ , G])
    ## transform data?
    if (!is.null(datafun)) {
      extraArgs <- formals(datafun)
      neededArgs <- intersect(names(extraArgs), names(availableArgs))
      extraArgs <- do.call(c, lapply(neededArgs, function(nn) availableArgs[nn]))
      extraArgs$data <- d
      originalNames <- colnames(d)
      d <- do.call(datafun, extraArgs)
      ## coerce extraOut to data.frame
      if (!is.data.frame(d)) stop('Argument "datafun" did not return a data.frame')
      if (!all(originalNames %in% colnames(d)))
        stop('The data.frame returned by argument "datafun" did not contain ',
             'column names required by the model:\n',
             paste(setdiff(originalNames, colnames(d)), collapse = ", "))
    }

    ## for ordered indicators, check that groups have same observed categories
    ordVars <- lavaan::lavNames(con, type = "ov.ord")
    if (length(ordVars) > 0) {
      try(onewayTables <- lavaan::lavTables(d, dimension = 1L,
                                            categorical = ordVars, group = G),
          silent = TRUE)
      if (exists("onewayTables")) {
        if (any(onewayTables$obs.prop == 1)) {
          nSparse <- nSparse + 1L
          next
        }
      } else {
        ## no "onewayTables" probably indicates empty categories in 1+ groups
        nSparse <- nSparse + 1L
        next
      }
    }
    ## fit null model, if it exists
    if (!is.null(null)) {
      out.null <- lavaan::update(null, data = d, group.label = lavaan::lavInspect(con, "group.label"))
    }

    ## fit constrained model, check for convergence
    try(out0 <- lavaan::update(con, data = d, group.label = lavaan::lavInspect(con, "group.label")))
    if (!exists("out0")) {
      nTries <- nTries + 1L
      next
    }
    if (!lavaan::lavInspect(out0, "converged")) {
      nTries <- nTries + 1L
      next
    }

    ## fit unconstrained model (unless NULL), check for convergence
    if (!is.null(uncon)) {
      try(out1 <- lavaan::update(uncon, data = d, group.label = lavaan::lavInspect(con, "group.label")))
      if (!exists("out1")) {
        nTries <- nTries + 1L
        next
      }
      if (!lavaan::lavInspect(out1, "converged")) {
        nTries <- nTries + 1L
        next
      }

    }
    ## If you get this far, everything converged, so break WHILE loop
    break
  }
  ## if WHILE loop ended before getting results, return NA
  if ( (nSparse == maxSparse) | (nTries == maxNonconv) ) {
    allAFIs <- c(AFIs, moreAFIs)
    AFI <- rep(NA, sum(!is.na(allAFIs)))
    names(AFI) <- allAFIs[!is.na(allAFIs)]
    MI <- if (is.null(param)) NULL else NA
    extra.obs <- NA
    nTries <- nTries + 1L
  } else {
    availableArgs$con <- out0
    if (exists("out1")) availableArgs$uncon <- out1
    if (exists("out.null")) availableArgs$null <- out.null
    AFI <- do.call(getAFIs, availableArgs)
    ## save max(MI) if !is.null(param)
    if (is.null(param)) {
      MI <- NULL
    } else {
      MI <- max(do.call(getMIs, c(availableArgs, modelType = "mgcfa"))$X2)
    }
    ## anything extra?
    if (!is.null(extra)) {
      extraArgs <- formals(extra)
      neededArgs <- intersect(names(extraArgs), names(availableArgs))
      extraArgs <- do.call(c, lapply(neededArgs, function(nn) availableArgs[nn]))
      extraOut <- do.call(extra, extraArgs)
      ## coerce extraOut to data.frame
      if (!is.list(extraOut)) extraOut <- as.list(extraOut)
      extra.obs <- data.frame(extraOut)
    } else extra.obs <- data.frame(NULL)
  }
  options(warn = old_warn)
  list(AFI = AFI, MI = MI, extra = extra.obs,
       n.nonConverged = nTries - 1L, n.Sparse = nSparse)
}

permuteOnce.mimic <- function(i, d, G, con, uncon, null, param, freeParam,
                              covariates, AFIs, moreAFIs, maxSparse, maxNonconv,
                              iseed, warn, extra = NULL, datafun = NULL) {
  old_warn <- options()$warn
  options(warn = warn)
  ## save arguments from call
  argNames <- names(formals(permuteOnce.mimic))
  availableArgs <- lapply(argNames, function(x) eval(as.name(x)))
  names(availableArgs) <- argNames

  nTries <- 1L
  while (nTries <= maxNonconv) {
    ## permute covariate(s) within each group
    if (length(G)) {
      for (gg in lavaan::lavInspect(con, "group.label")) {
        dG <- d[ d[[G]] == gg, ]
        N <- nrow(dG)
        newd <- dG[sample(1:N, N), covariates, drop = FALSE]
        for (COV in covariates) d[d[[G]] == gg, COV] <- newd[ , COV]
      }
    } else {
      N <- nrow(d)
      newd <- d[sample(1:N, N), covariates, drop = FALSE]
      for (COV in covariates) d[ , COV] <- newd[ , COV]
    }
    ## transform data?
    if (!is.null(datafun)) {
      extraArgs <- formals(datafun)
      neededArgs <- intersect(names(extraArgs), names(availableArgs))
      extraArgs <- do.call(c, lapply(neededArgs, function(nn) availableArgs[nn]))
      extraArgs$data <- d
      originalNames <- colnames(d)
      d <- do.call(datafun, extraArgs)
      ## coerce extraOut to data.frame
      if (!is.data.frame(d)) stop('Argument "datafun" did not return a data.frame')
      if (!all(originalNames %in% colnames(d)))
        stop('The data.frame returned by argument "datafun" did not contain ',
             'column names required by the model:\n',
             paste(setdiff(originalNames, colnames(d)), collapse = ", "))
    }


    ## fit null model, if it exists
    if (!is.null(null)) {
      out.null <- lavaan::update(null, data = d, group.label = lavaan::lavInspect(con, "group.label"))
    }

    ## fit constrained model
    try(out0 <- lavaan::update(con, data = d, group.label = lavaan::lavInspect(con, "group.label")))
    ## check for convergence
    if (!exists("out0")) {
      nTries <- nTries + 1L
      next
    }
    if (!lavaan::lavInspect(out0, "converged")) {
      nTries <- nTries + 1L
      next
    }
    ## If you get this far, everything converged, so break WHILE loop
    break
  }
  ## if WHILE loop ended before getting results, return NA
  if (nTries == maxNonconv) {
    allAFIs <- c(AFIs, moreAFIs)
    AFI <- rep(NA, sum(!is.na(allAFIs)))
    names(AFI) <- allAFIs[!is.na(allAFIs)]
    MI <- if (is.null(param)) NULL else NA
    extra.obs <- NA
    nTries <- nTries + 1L
  } else {
    availableArgs$con <- out0
    if (exists("out.null")) availableArgs$null <- out.null
    AFI <- do.call(getAFIs, availableArgs)
    if (is.null(param)) {
      MI <- NULL
    } else {
      MI <- max(do.call(getMIs, c(availableArgs, modelType = "mimic"))$X2)
    }
    ## anything extra?
    if (!is.null(extra)) {
      extraArgs <- formals(extra)
      neededArgs <- intersect(names(extraArgs), names(availableArgs))
      extraArgs <- do.call(c, lapply(neededArgs, function(nn) availableArgs[nn]))
      extraOut <- do.call(extra, extraArgs)
      ## coerce extraOut to data.frame
      if (!is.list(extraOut)) extraOut <- as.list(extraOut)
      extra.obs <- data.frame(extraOut)
    } else extra.obs <- data.frame(NULL)
  }
  options(warn = old_warn)
  list(AFI = AFI, MI = MI, extra = extra.obs,
       n.nonConverged = nTries - 1L, n.Sparse = integer(length = 0))
}


## Function to permute difference in fits
permuteMeasEq <- function(nPermute, modelType = c("mgcfa","mimic"),
                          con, uncon = NULL, null = NULL,
                          param = NULL, freeParam = NULL, covariates = NULL,
                          AFIs = NULL, moreAFIs = NULL,
                          maxSparse = 10, maxNonconv = 10, showProgress = TRUE,
                          warn = -1, datafun, extra,
                          parallelType = c("none","multicore","snow"),
                          ncpus = NULL, cl = NULL, iseed = 12345) {

  ## save arguments from call
  availableArgs <- as.list(formals(permuteMeasEq))
  argNames <- names(availableArgs)
  if (missing(datafun)) argNames <- setdiff(argNames, "datafun")
  if (missing(extra)) argNames <- setdiff(argNames, "extra")
  for (aa in argNames) {
    if (!is.null(eval(as.name(aa))))
      suppressWarnings(availableArgs[[aa]] <- eval(as.name(aa)))
  }
  ## check and return them
  fullCall <- do.call(checkPermArgs, availableArgs)
  ## assign them to workspace (also adds old_RNG & oldSeed to workspace)
  for (aa in names(fullCall)) assign(aa, fullCall[[aa]])

  ###################### SAVE OBSERVED RESULTS ##########################
  AFI.obs <- do.call(getAFIs, fullCall)
  ## save modification indices if !is.null(param)
  if (is.null(param)) {
    MI.obs <- data.frame(NULL)
  } else MI.obs <- do.call(getMIs, fullCall)

  ## anything extra?
  if (!missing(extra)) {
    extraArgs <- formals(extra)
    neededArgs <- intersect(names(extraArgs), names(fullCall))
    extraArgs <- do.call(c, lapply(neededArgs, function(nn) fullCall[nn]))
    extraOut <- do.call(extra, extraArgs)
    ## check that extra() returns a named list of scalars
    if (!is.list(extraOut)) extraOut <- as.list(extraOut)
    wrongFormat <- paste('Function "extra" must return a numeric vector or a',
                         'list of scalars, with each element named.')
    if (!all(sapply(extraOut, is.numeric))) stop(wrongFormat)
    if (!all(sapply(extraOut, length) == 1L)) stop(wrongFormat)
    if (is.null(names(extraOut)) | any(names(extraOut) == "")) stop(wrongFormat)
    extra.obs <- do.call(c, extraOut)
  } else extra.obs <- numeric(length = 0L)

  ######################### PREP DATA ##############################
  argList <- fullCall[c("con","uncon","null","param","freeParam","covariates",
                        "AFIs","moreAFIs","maxSparse","maxNonconv","warn","iseed")]
  argList$G <- lavaan::lavInspect(con, "group")
    ## check for categorical variables
    # catVars <- lavaan::lavNames(con, type = "ov.ord")
    # numVars <- lavaan::lavNames(con, type = "ov.num")
    # latentVars <- lavaan::lavNames(con, type = "lv.regular")
  ## assemble data to which the models were fit
  if (length(argList$G)) {
    dataList <- mapply(FUN = function(x, g, n) {
      y <- data.frame(as.data.frame(x), g, stringsAsFactors = FALSE)
      names(y) <- c(n, argList$G)
      y
    }, SIMPLIFY = FALSE,
    x = lavaan::lavInspect(con, "data"), g = lavaan::lavInspect(con, "group.label"),
    n = lavaan::lavNames(con, type = "ov",
                         group = seq_along(lavaan::lavInspect(con, "group.label"))))
    argList$d <- do.call(rbind, dataList)
  } else {
    argList$d <- as.data.frame(lavaan::lavInspect(con, "data"))
    names(argList$d) <- lavaan::lavNames(con, type = "ov")
  }
  ## check that covariates are actual variables
  if (modelType == "mimic") {
    if (length(covariates) && !all(covariates %in% names(argList$d)))
      stop('These specified covariates are not columns in the data.frame:\n',
           paste(setdiff(covariates, names(argList$d)), collapse = ", "))
  }
  ## anything extra?
  if (!missing(extra)) argList$extra <- extra
  if (!missing(datafun)) argList$datafun <- datafun

  ###################### PERMUTED RESULTS ###########################
  ## permute and return distributions of (delta)AFIs, largest MI, and extras
  if (showProgress) {
    mypb <- txtProgressBar(min = 1, max = nPermute, initial = 1, char = "=",
                           width = 50, style = 3, file = "")
    permuDist <- list()
    for (j in 1:nPermute) {
      permuDist[[j]] <- do.call(paste("permuteOnce", modelType, sep = "."),
                                args = c(argList, i = j))
      setTxtProgressBar(mypb, j)
    }
    close(mypb)
  } else if (parallelType == "multicore") {
    if (length(iseed)) set.seed(iseed)
    argList$FUN <- paste("permuteOnce", modelType, sep = ".")
    argList$X <- 1:nPermute
    argList$mc.cores <- ncpus
    argList$mc.set.seed <- TRUE
	pmcl <- function(...) { parallel::mclapply(...) }
    permuDist <- do.call(pmcl, args = argList)
    ## restore old RNG type
    if (fullCall$old_RNG[1] != "L'Ecuyer-CMRG") RNGkind(fullCall$old_RNG[1])
  } else if (parallelType == "snow") {
    stopTheCluster <- FALSE
    if (is.null(cl)) {
      stopTheCluster <- TRUE
      cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
    }
    parallel::clusterSetRNGStream(cl, iseed = iseed)
    # clusterExport(cl, c("getAFIs","getMIs","permuteOnce.mgcfa","permuteOnce.mimic"))
    argList$cl <- cl
    argList$X <- 1:nPermute
    argList$fun <- paste("permuteOnce", modelType, sep = ".")
	tempppl <- function(...) { parallel::parLapply(...) }
    permuDist <- do.call(tempppl, args = argList)
    if (stopTheCluster) parallel::stopCluster(cl)
    ## restore old RNG type
    if (fullCall$old_RNG[1] != "L'Ecuyer-CMRG") RNGkind(fullCall$old_RNG[1])
  } else {
    argList$X <- 1:nPermute
    argList$FUN <- paste("permuteOnce", modelType, sep = ".")
    permuDist <- do.call(lapply, args = argList)
  }

  ## extract AFI distribution
  if (length(AFI.obs) > 1) {
    AFI.dist <- as.data.frame(t(sapply(permuDist, function(x) x$AFI)))
  }
  if (length(AFI.obs) == 1L) {
    AFI.dist <- data.frame(sapply(permuDist, function(x) x$AFI))
    colnames(AFI.dist) <- names(AFI.obs)
  }
  ## identify badness-of-fit measures
  badness <- grepl(pattern = "fmin|chi|aic|bic|rmr|rmsea|cn|sic|hqc",
                   x = names(AFI.obs), ignore.case = TRUE)
  ## calculate all one-directional p-values
  AFI.pval <- mapply(FUN = function(x, y, b) {
      if (b) return(mean(x >= y, na.rm = TRUE))
      mean(x <= y, na.rm = TRUE)
    }, x = unclass(AFI.dist), y = AFI.obs, b = badness)

  ## extract distribution of maximum modification indices
  MI.dist <- as.numeric(unlist(lapply(permuDist, function(x) x$MI)))
  ## calculate Tukey-adjusted p values for modification indices
  if (!is.null(param)) {
    MI.obs$tukey.p.value <- sapply(MI.obs$X2,
                                   function(i) mean(i <= MI.dist, na.rm = TRUE))
    MI.obs <- as.data.frame(unclass(MI.obs))
    rownames(MI.obs) <- names(param)
  }

  ## anything extra?
  if (!missing(extra)) {
    extra.dist <- do.call(rbind, lapply(permuDist, function(x) x$extra))
  } else extra.dist <- data.frame(NULL)

  ## save parameter table for show/summary methods
  PT <- as.data.frame(lavaan::parTable(con))
  PT$par <- paste0(PT$lhs, PT$op, PT$rhs)
  if (length(lavaan::lavInspect(con, "group")))
    PT$group.label[PT$group > 0] <- lavaan::lavInspect(con, "group.label")[PT$group[PT$group > 0] ]

  ## return observed results, permutation p values, and ANOVA results
  if (is.null(uncon)) {
    delta <- lavaan::anova(con)
  } else {
    delta <- lavaan::anova(uncon, con)
  }
  ANOVA <- sapply(delta[,c("Chisq diff","Df diff","Pr(>Chisq)")], function(x) x[2])
  out <- new("permuteMeasEq", PT = PT, modelType = modelType, ANOVA = ANOVA,
             AFI.obs = AFI.obs, AFI.dist = AFI.dist, AFI.pval = AFI.pval,
             MI.obs = MI.obs, MI.dist = MI.dist,
             extra.obs = extra.obs, extra.dist = extra.dist,
             n.Permutations = nPermute, n.Converged = sum(!is.na(AFI.dist[,1])),
             n.nonConverged = sapply(permuDist, function(x) x$n.nonConverged),
             n.Sparse = sapply(permuDist, function(x) x$n.Sparse),
             oldSeed = fullCall$oldSeed)
  out
}

## methods
setMethod("show", "permuteMeasEq", function(object) {
  ## print warning if there are nonConverged permutations
  if (object@n.Permutations != object@n.Converged) {
    warning(paste("Only", object@n.Converged, "out of",
                  object@n.Permutations, "models converged within",
                  max(object@n.nonConverged), "attempts per permutation.\n\n"))
  }
  ## print ANOVA
  cat("Omnibus p value based on parametric chi-squared difference test:\n\n")
  print(round(object@ANOVA, digits = 3))
  ## print permutation results
  cat("\n\nOmnibus p values based on nonparametric permutation method: \n\n")
  AFI <- data.frame(AFI.Difference = object@AFI.obs, p.value = object@AFI.pval)
  class(AFI) <- c("lavaan.data.frame","data.frame")
  print(AFI, nd = 3)
  invisible(object)
})

setMethod("summary", "permuteMeasEq", function(object, alpha = .05, nd = 3,
                                               extra = FALSE) {
  ## print warning if there are nonConverged permutations
  if (object@n.Permutations != object@n.Converged) {
    warning(paste("Only", object@n.Converged, "out of",
                  object@n.Permutations, "models converged within",
                  max(object@n.nonConverged), "attempts per permutation.\n\n"))
  }
  ## print ANOVA
  cat("Omnibus p value based on parametric chi-squared difference test:\n\n")
  print(round(object@ANOVA, digits = nd))
  ## print permutation results
  cat("\n\nOmnibus p values based on nonparametric permutation method: \n\n")
  AFI <- data.frame(AFI.Difference = object@AFI.obs, p.value = object@AFI.pval)
  class(AFI) <- c("lavaan.data.frame","data.frame")
  print(AFI, nd = nd)

  ## print extras or DIF test results, if any were requested
  if (extra && length(object@extra.obs)) {
    cat("\n\nUnadjusted p values of extra statistics,\n",
        "based on permutation distribution of each statistic: \n\n")
    MI <- data.frame(Statistic = object@extra.obs)
    class(MI) <- c("lavaan.data.frame","data.frame")
    MI$p.value <- sapply(names(object@extra.dist), function(nn) {
      mean(abs(object@extra.dist[,nn]) >= abs(object@extra.obs[nn]), na.rm = TRUE)
    })
    MI$flag <- ifelse(MI$p.value < alpha, "*   ", "")
    print(MI, nd = nd)
  } else if (length(object@MI.dist)) {
    cat("\n\n Modification indices for equality constrained parameter estimates,\n",
        "with unadjusted 'p.value' based on chi-squared distribution and\n",
        "adjusted 'tukey.p.value' based on permutation distribution of the\n",
        "maximum modification index per iteration: \n\n")
    MI <- do.call(paste("summ", object@modelType, sep = "."),
                  args = list(object = object, alpha = alpha))
    print(MI, nd = nd)

    ## print messages about potential DIF
    if (all(MI$tukey.p.value > alpha)) {
      cat("\n\n No equality constraints were flagged as significant.\n\n")
      return(invisible(MI))
    }
    if (object@modelType == "mgcfa") {
      cat("\n\nThe following equality constraints were flagged as significant:\n\n")
      for (i in which(MI$tukey.p.value < alpha)) {
        cat("Parameter '", MI$parameter[i], "' may differ between Groups '",
            MI$group.lhs[i], "' and '", MI$group.rhs[i], "'.\n", sep = "")
      }
      cat("\nUse lavTestScore(..., epc = TRUE) on your constrained model to",
          "display expected parameter changes for these equality constraints\n\n")
    }

  } else return(invisible(object))

  invisible(MI)
})

summ.mgcfa <- function(object, alpha) {
  MI <- object@MI.obs
  class(MI) <- c("lavaan.data.frame","data.frame")
  PT <- object@PT
  eqPar <- rbind(PT[PT$plabel %in% MI$lhs, ], PT[PT$plabel %in% MI$rhs, ])
  MI$flag <- ""
  MI$parameter <- ""
  MI$group.lhs <- ""
  MI$group.rhs <- ""
  for (i in 1:nrow(MI)) {
    par1 <- eqPar$par[ eqPar$plabel == MI$lhs[i] ]
    par2 <- eqPar$par[ eqPar$plabel == MI$rhs[i] ]
    MI$parameter[i] <- par1
    MI$group.lhs[i] <- eqPar$group.label[ eqPar$plabel == MI$lhs[i] ]
    MI$group.rhs[i] <- eqPar$group.label[ eqPar$plabel == MI$rhs[i] ]
    if (par1 != par2) {
      myMessage <- paste0("Constraint '", MI$lhs[i], "==", MI$rhs[i],
                          "' refers to different parameters: \n'",
                          MI$lhs[i], "' is '", par1, "' in group '",
                          MI$group.lhs[i], "'\n'",
                          MI$rhs[i], "' is '", par2, "' in group '",
                          MI$group.rhs[i], "'\n")
      warning(myMessage)
    }
    if (MI$tukey.p.value[i] < alpha) MI$flag[i] <- "*  -->"
  }
  MI
}

summ.mimic <- function(object, alpha) {
  MI <- object@MI.obs
  class(MI) <- c("lavaan.data.frame","data.frame")
  MI$flag <- ifelse(MI$tukey.p.value < alpha, "*   ", "")
  MI
}

setMethod("hist", "permuteMeasEq", function(x, ..., AFI, alpha = .05, nd = 3,
                                            printLegend = TRUE,
                                            legendArgs = list(x = "topleft")) {
  histArgs <- list(...)
  histArgs$x <- x@AFI.dist[[AFI]]
  if (is.null(histArgs$col)) histArgs$col <- "grey69"
  histArgs$freq <- !grepl("chi", AFI)
  histArgs$ylab <- if (histArgs$freq) "Frequency" else "Probability Density"

  if (printLegend) {
    if (is.null(legendArgs$box.lty)) legendArgs$box.lty <- 0
    if (nd < length(strsplit(as.character(1 / alpha), "")[[1]]) - 1) {
      warning(paste0("The number of digits argument (nd = ", nd ,
                     ") is too low to display your p value at the",
                     " same precision as your requested alpha level (alpha = ",
                     alpha, ")"))
    }
    if (x@AFI.pval[[AFI]] < (1 / 10^nd)) {
      pVal <- paste(c("< .", rep(0, nd - 1),"1"), collapse = "")
    } else {
      pVal <- paste("=", round(x@AFI.pval[[AFI]], nd))
    }
  }

  delta <- length(x@MI.dist) > 0L && x@modelType == "mgcfa"
  if (grepl("chi", AFI)) {   ####################################### Chi-squared
    ChiSq <- x@AFI.obs[AFI]
    DF <- x@ANOVA[2]
    histArgs$xlim <- range(c(ChiSq, x@AFI.dist[[AFI]], qchisq(c(.01, .99), DF)))
    xVals <- seq(histArgs$xlim[1], histArgs$xlim[2], by = .1)
    theoDist <- dchisq(xVals, df = DF)
    TheoCrit <- round(qchisq(p = alpha, df = DF, lower.tail = FALSE), 2)
    Crit <- quantile(histArgs$x, probs = 1 - alpha)
    if (ChiSq > histArgs$xlim[2]) histArgs$xlim[2] <- ChiSq
    if (delta) {
      histArgs$main <- expression(Permutation~Distribution~of~Delta*chi^2)
      histArgs$xlab <- expression(Delta*chi^2)
      if (printLegend) {
        legendArgs$legend <- c(bquote(Theoretical~Delta*chi[Delta*.(paste("df =", DF))]^2 ~ Distribution),
                               bquote(Critical~chi[alpha~.(paste(" =", alpha))]^2 == .(round(TheoCrit, nd))),
                               bquote(.(paste("Permuted Critical Value =", round(Crit, nd)))),
                               bquote(Observed~Delta*chi^2 == .(round(ChiSq, nd))),
                               expression(paste("")),
                               bquote(Permuted~italic(p)~.(pVal)))
      }
    } else {
      histArgs$main <- expression(Permutation~Distribution~of~chi^2)
      histArgs$xlab <- expression(chi^2)
      if (printLegend) {
        legendArgs$legend <- c(bquote(Theoretical~chi[.(paste("df =", DF))]^2 ~ Distribution),
                               bquote(Critical~chi[alpha~.(paste(" =", alpha))]^2 == .(round(TheoCrit, nd))),
                               bquote(.(paste("Permuted Critical Value =", round(Crit, nd)))),
                               bquote(Observed~chi^2 == .(round(ChiSq, nd))),
                               expression(paste("")),
                               bquote(Permuted~italic(p)~.(pVal)))
      }
    }
    H <- do.call(hist, c(histArgs["x"], plot = FALSE))
    histArgs$ylim <- c(0, max(H$density, theoDist))
    if (printLegend) {
      legendArgs <- c(legendArgs, list(lty = c(2, 2, 1, 1, 0, 0),
                                       lwd = c(2, 2, 2, 3, 0, 0),
                                       col = c("black","black","black","red","","")))
    }
  } else {        ################################################### other AFIs
    badness <- grepl(pattern = "fmin|aic|bic|rmr|rmsea|cn|sic|hqc",
                     x = AFI, ignore.case = TRUE)
    if (badness) {
      Crit <- quantile(histArgs$x, probs = 1 - alpha)
    } else {
      Crit <- quantile(histArgs$x, probs = alpha)
    }
    histArgs$xlim <- range(histArgs$x, x@AFI.obs[AFI])
    if (delta) {
      histArgs$main <- bquote(~Permutation~Distribution~of~Delta*.(toupper(AFI)))
      histArgs$xlab <- bquote(~Delta*.(toupper(AFI)))
      if (printLegend) {
        legendArgs$legend <- c(bquote(Critical~Delta*.(toupper(AFI))[alpha~.(paste(" =", alpha))] == .(round(Crit, nd))),
                               bquote(Observed~Delta*.(toupper(AFI)) == .(round(x@AFI.obs[AFI], nd))),
                               expression(paste("")),
                               bquote(Permuted~italic(p)~.(pVal)))

      }
    } else {
      histArgs$main <- paste("Permutation Distribution of", toupper(AFI))
      histArgs$xlab <- toupper(AFI)
      if (printLegend) {
        legendArgs$legend <- c(bquote(Critical~.(toupper(AFI))[alpha~.(paste(" =", alpha))] == .(round(Crit, nd))),
                               bquote(Observed~.(toupper(AFI)) == .(round(x@AFI.obs[AFI], nd))),
                               expression(paste("")),
                               bquote(Permuted~italic(p)~.(pVal)))

      }
    }
    if (printLegend) {
      legendArgs <- c(legendArgs, list(lty = c(1, 1, 0, 0),
                                       lwd = c(2, 3, 0, 0),
                                       col = c("black","red","","")))
    }
  }
  ## print histogram (and optionally, print legend)
  suppressWarnings({
    do.call(hist, histArgs)
    if (grepl("chi", AFI)) {
      lines(x = xVals, y = theoDist, lwd = 2, lty = 2)
      abline(v = TheoCrit, col = "black", lwd = 2, lty = 2)
    }
    abline(v = Crit, col = "black", lwd = 2)
    abline(v = x@AFI.obs[AFI], col = "red", lwd = 3)
    if (printLegend) do.call(legend, legendArgs)
  })
  ## return arguments to create histogram (and optionally, legend)
  invisible(list(hist = histArgs, legend = legendArgs))
})

