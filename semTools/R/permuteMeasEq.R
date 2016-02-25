### Terrence D. Jorgensen
### Last updated: 24 February 2016
### permutation randomization test for measurement equivalence and DIF


## create s4 class for result object
setClass("permuteMeasEq", slots = c(PT = "data.frame",
                                    ANOVA = "vector",
                                    AFI.obs = "vector",
                                    AFI.dist = "data.frame",
                                    AFI.pval = "vector",
                                    MI.obs = "data.frame",
                                    MI.dist = "vector",
                                    n.Permutations = "integer",
                                    n.Converged = "integer",
                                    n.nonConverged = "vector",
                                    n.Sparse = "vector"))

## Function to extract modification indices for equality constraints
getMIs <- function(con, param) {
  if (class(con) != "lavaan") stop("This function only applies to fitted lavaan models.")
  if (con@Data@ngroups == 1L) stop("This function only applies to multiple-group models.")
  ## save all estimates from less constrained model
  PT <- parTable(con)[ , c("lhs","op","rhs","group","plabel")]
  ## extract parameters of interest
  if (param[1] == "loadings") params <- PT[PT$op == "=~", ]
  if (param[1] == "intercepts") params <- PT[PT$lhs %in% con@Data@ov$name & PT$op == "~1", ]
  if (param[1] == "thresholds") params <- PT[PT$op == "|", ]
  if (param[1] == "residuals") params <- PT[(PT$lhs %in% con@Data@ov$name) &
                                            (PT$lhs == PT$rhs) & PT$op == "~~", ]
  if (!param[1] %in% c("loadings","intercepts","thresholds","residuals")) {
    allNames <- paste0(PT$lhs, PT$op, PT$rhs)
    params <- PT[allNames %in% param, ]
  }
  ## return modification indices for specified constraints (param)
  MIs <- lavTestScore(con)$uni
  MIs[MIs$lhs %in% params$plabel, ]
}

## Function to find delta-AFIs AND maximum modification index in one permutation
permuteOnce <- function(i, d, con, uncon, null, param, G, AFIs, moreAFIs,
                        maxSparse = 10, maxNonconv = 10) {
  nSparse <- 0L
  nTries <- 1L
  while ( (nSparse <= maxSparse) & (nTries <= maxNonconv) ) {
    ## permute grouping variable
    d[ , G] <- sample(d[ , G])

    ## for ordered indicators, check that groups have same observed categories
    ordVars <- con@Data@ov$name[con@Data@ov$type == "ordered"]
    if (length(ordVars) > 0) {
      onewayTables <- lavTables(d, dimension = 1L,
                                categorical = ordVars, group = G)
      if (any(onewayTables$obs.prop == 1)) {
        nSparse <- nSparse + 1L
        next
      }
    }
    ## fit null model, if it exists
    if (!is.null(null)) {
      out.null <- lavaan(data = d, group = G, slotParTable = null@ParTable,
                         slotOptions = null@Options)
    } else out.null <- NULL

    ## fit constrained model, check for convergence
    suppressWarnings(try(out0 <- lavaan(data = d, group = G,
                                        slotParTable = con@ParTable,
                                        slotOptions = con@Options)))
    if (!exists("out0")) {
      nTries <- nTries + 1L
      next
    }
    if (!inspect(out0, "converged")) {
      nTries <- nTries + 1L
      next
    }

    ## fit unconstrained model (unless NULL), check for convergence
    if (!is.null(uncon)) {
    suppressWarnings(try(out1 <- lavaan(data = d, group = G,
                                        slotParTable = uncon@ParTable,
                                        slotOptions = uncon@Options)))
    if (!exists("out1")) {
      nTries <- nTries + 1L
      next
    }
    if (!inspect(out1, "converged")) {
      nTries <- nTries + 1L
      next
    }

    }

    fit1 <- list()
    fit0 <- list()
    if (!is.null(AFIs[1])) {
      if (!is.null(uncon)) fit1[[1]] <- fitMeasures(out1, fit.measures = AFIs,
                                                    baseline.model = out.null)
      fit0[[1]] <- fitMeasures(out0, fit.measures = AFIs,
                               baseline.model = out.null)
    }
    if (!is.null(moreAFIs[1])) {
      if (!is.null(uncon)) fit1[[2]] <- moreFitIndices(out1, fit.measures = moreAFIs)
      fit0[[2]] <- moreFitIndices(out0, fit.measures = moreAFIs)
    }

    break
  } ## end WHILE loop
  ## if loop ended before getting results, return NA
  if ( (nSparse == maxSparse) | (nTries == maxNonconv) ) {
    allAFIs <- c(AFIs, moreAFIs)
    AFI <- rep(NA, sum(!is.na(allAFIs)))
    names(AFI) <- allAFIs[!is.na(allAFIs)]
    MI <- if (is.null(uncon)) NULL else NA
    nTries <- nTries + 1L
  } else {
    ## calculate AFI for configural, otherwise delta-AFI & max(MI)
    if (is.null(uncon)) {
      AFI <- unlist(fit0)
      MI <- NULL
    } else {
      AFI <- unlist(fit1) - unlist(fit0)
      MI <- max(getMIs(out0, param)$X2)
    }
  }
  list(AFI = AFI, MI = MI, n.nonConverged = nTries - 1L, n.Sparse = nSparse)
}

## Function to permute difference in fits
permuteMeasEq <- function(nPermute, con, uncon = NULL, null = NULL,
                          param = NULL, AFIs = NULL, moreAFIs = NULL,
                          maxSparse = 10, maxNonconv = 10, showProgress = TRUE) {
  library(lavaan)
  nPermute <- as.integer(nPermute[1])

  ## check that "param" is NULL if uncon is NULL, and check for lavaan class
  notLavaan <- "Non-NULL 'con', 'uncon', or 'null' must be fitted lavaan object."
  if (is.null(uncon)) {
    if (!is.null(param)) {
      message(c(" When 'uncon = NULL', only configural invariance is tested.",
                "\n So the 'param' argument was changed to NULL."))
      param <- NULL
    }
    if (class(con) != "lavaan") stop(notLavaan)
  } else {
    if (class(con) != "lavaan") stop(notLavaan)
    if (class(uncon) != "lavaan") stop(notLavaan)
  }
  if (!is.null(null)) {
    if (class(null) != "lavaan") stop(notLavaan)
  }

  ## FIXME: Temporarily warn about testing thresholds without necessary constraints
  if (!is.null(param)) {
    if (param == "thresholds" | any(grepl("\\|", param))) {
      warning(c("This function is not yet optimized for testing thresholds.\n",
                "Necessary identification contraints might not be specified."))
    }
  }
  if (is.null(AFIs) & is.null(moreAFIs)) {
    message("No AFIs were selected, so only chi-squared will be permuted.")
    AFIs <- "chisq"
  }
  if ("ecvi" %in% AFIs) stop("ECVI is not available for multigroup models.")

  ## check estimators
  leastSq <- grepl("LS", con@Options$estimator)
  if (!is.null(uncon)) {
    if (uncon@Options$estimator != con@Options$estimator)
      stop("Models must be fit using same estimator.")
  }

  ###################### OBSERVED RESULTS ##########################
  AFI1 <- list()
  AFI0 <- list()
  ## check validity of user-specified AFIs, save output
  if (!is.null(AFIs)) {
    IC <- grep("ic|logl", AFIs, value = TRUE)
    if (leastSq & length(IC)) {
      stop(paste("Argument 'AFIs' includes invalid options:",
                 paste(IC, collapse = ", "),
                 "Information criteria unavailable for least-squares estimators.",
                 sep = "\n"))
    }
    if (!is.null(uncon)) AFI1[[1]] <- fitMeasures(uncon, fit.measures = AFIs,
                                                  baseline.model = null)
    AFI0[[1]] <- fitMeasures(con, fit.measures = AFIs, baseline.model = null)
  }
  ## check validity of user-specified moreAFIs
  if (!is.null(moreAFIs)) {
    IC <- grep("ic|hqc", moreAFIs, value = TRUE)
    if (leastSq & length(IC)) {
      stop(paste("Argument 'moreAFIs' includes invalid options:",
                 paste(IC, collapse = ", "),
                 "Information criteria unavailable for least-squares estimators.",
                 sep = "\n"))
    }
    if (!is.null(uncon)) AFI1[[2]] <- moreFitIndices(uncon, fit.measures = moreAFIs)
    AFI0[[2]] <- moreFitIndices(con, fit.measures = moreAFIs)
  }

  ## save observed delta-AFIs, and modification indices if !is.null(uncon)
  if (is.null(uncon)) {
    AFI.obs <- unlist(AFI0)
    MI.obs <- data.frame(NULL)
  } else {
    AFI.obs <- unlist(AFI1) - unlist(AFI0)
    MI.obs <- getMIs(con, param = param)
  }

  ######################### PREP DATA ##############################
  ## save name of grouping variable
  G <- con@Data@group
  ## check for categorical variables
  # catVars <- which(con@Data@ov$type[!con@Data@ov$exo] == "ordered")
  # numVars <- which(con@Data@ov$type[!con@Data@ov$exo] != "ordered")
  # latentVars <- con@pta$vnames$lv[[1]]

  ## assemble data to which the models were fit
  dataList <- mapply(FUN = function(x, g, n) {
    y <- data.frame(as.data.frame(x), g, stringsAsFactors = FALSE)
    names(y) <- c(n, con@Data@group)
    y
  }, x = con@Data@X, g = con@Data@group.label,
     n = con@Data@ov.names, SIMPLIFY = FALSE)
  allData <- do.call(rbind, dataList)

  ###################### PERMUTED RESULTS ###########################
  ## permute groups and return distributions of delta-AFIs and largest MI
  if (showProgress) {
    mypb <- txtProgressBar(min = 1, max = nPermute, initial = 1, char = "=",
                           width = 50, style = 3, file = "")
    permuDist <- list()
    for (j in 1:nPermute) {
      permuDist[[j]] <- permuteOnce(j, d = allData, con = con, uncon = uncon,
                                    null = null, param = param, G = G,
                                    AFIs = AFIs, moreAFIs = moreAFIs,
                                    maxSparse = maxSparse, maxNonconv = maxNonconv)
      setTxtProgressBar(mypb, j)
    }
    close(mypb)
  } else {
    permuDist <- lapply(1:nPermute, permuteOnce, d = allData, con = con,
                        uncon = uncon, null = null, param = param, G = G,
                        AFIs = AFIs, moreAFIs = moreAFIs,
                        maxSparse = maxSparse, maxNonconv = maxNonconv)
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
      if (b) return(mean(x <= y, na.rm = TRUE))
      mean(x >= y, na.rm = TRUE)
    }, x = unclass(AFI.dist), y = AFI.obs, b = badness)

  ## extract distribution of maximum modification indices
  MI.dist <- as.numeric(unlist(lapply(permuDist, function(x) x$MI)))
  ## calculate Tukey-adjusted p values for modification indices
  if (!is.null(uncon)) {
    MI.obs$tukey.p.value <- sapply(MI.obs$X2, function(i) mean(i <= MI.dist))
  }

  ## save parameter table for show/summary methods
  PT <- as.data.frame(parTable(con))
  PT$par <- paste0(PT$lhs, PT$op, PT$rhs)
  PT$group.label[PT$group > 0] <- con@Data@group.label[PT$group]

  ## return observed results, permutation p values, and ANOVA results
  if (is.null(uncon)) {
    delta <- anova(con)
  } else {
    delta <- anova(uncon, con)
  }
  ANOVA <- sapply(delta[,c("Chisq diff","Df diff","Pr(>Chisq)")], function(x) x[2])
  out <- new("permuteMeasEq", PT = PT, ANOVA = ANOVA,
             AFI.obs = AFI.obs, AFI.dist = AFI.dist, AFI.pval = AFI.pval,
             MI.obs = as.data.frame(unclass(MI.obs)), MI.dist = MI.dist,
             n.Permutations = nPermute, n.Converged = sum(!is.na(AFI.dist[,1])),
             n.nonConverged = sapply(permuDist, function(x) x$n.nonConverged),
             n.Sparse = sapply(permuDist, function(x) x$n.Sparse))
  out
}

## methods
setMethod("show", "permuteMeasEq", function(object) {
  cat("Omnibus p value based on parametric chi-squared difference test:\n\n")
  print(round(object@ANOVA, digits = 3))

  cat("\n\nOmnibus p values based on nonparametric permutation method: \n\n")
  AFI <- data.frame(AFI.Difference = object@AFI.obs, p.value = object@AFI.pval)
  class(AFI) <- c("lavaan.data.frame","data.frame")
  print(AFI, nd = 3)
  ## print warning if there are nonConverged permutations
  if (object@n.Permutations != object@n.Converged) {
    warning(paste("Only", object@n.Converged, "out of",
                  object@n.Permutations, "models converged within",
                  max(object@n.nonConverged), "attempts per permutation."))
  }
  invisible(object)
})

setMethod("summary", "permuteMeasEq", function(object, alpha = .05, nd = 3) {
  cat("Omnibus p value based on parametric chi-squared difference test:\n\n")
  print(round(object@ANOVA, digits = nd))

  cat("\n\nOmnibus p values based on nonparametric permutation method: \n\n")
  AFI <- data.frame(AFI.Difference = object@AFI.obs, p.value = object@AFI.pval)
  class(AFI) <- c("lavaan.data.frame","data.frame")
  print(AFI, nd = nd)

  ## unless testing configural, print requested DIF test results
  if (length(object@MI.dist)) {
    cat("\n\n Modification indices for equality constrained parameter estimates,\n",
        "with unadjusted 'p.value' based on chi-squared distribution and\n",
        "adjusted 'tukey.p.value' based on permutation distribution of the\n",
        "maximum modification index per iteration: \n\n")
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
    print(MI, nd = nd)

    ## print messages about potential DIF
    if (all(MI$tukey.p.value > alpha)) {
      cat("\n\n No equality constraints were flagged as significant.\n\n")
      return(invisible(MI))
    }
    cat("\n\nThe following equality constraints were flagged as significant:\n\n")
    for (i in which(MI$tukey.p.value < alpha)) {
      cat("Parameter '", MI$parameter[i], "' may differ between Groups '",
          MI$group.lhs[i], "' and '", MI$group.rhs[i], "'.\n", sep = "")
    }
    cat("\nUse lavTestScore(..., epc = TRUE) on your constrained model to",
        "display expected parameter changes for these equality constraints\n\n")
    return(invisible(MI))
  }

  ## print warning if there are nonConverged permutations
  if (object@n.Permutations != object@n.Converged) {
    warning(paste("Only", object@n.Converged, "out of",
                  object@n.Permutations, "models converged within",
                  max(object@n.nonConverged), "attempts per permutation."))
  }
  invisible(object)
})

