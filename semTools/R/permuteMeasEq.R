### Terrence D. Jorgensen
### Last updated: 12 April 2016
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
getMIs <- function(con, param, freeParam) {
  if (class(con) != "lavaan") stop("This function applies to fitted lavaan models.")
  if (con@Data@ngroups == 1L) stop("This function applies to multigroup models.")
  ## strip white space
  param <- gsub("[[:space:]]+", "", param)
  freeParam <- gsub("[[:space:]]+", "", freeParam)
  ## save all estimates from constrained model
  PT <- lavaan::parTable(con)[ , c("lhs","op","rhs","group","plabel")]
  ## extract parameters of interest
  paramTypes <- c("loadings","intercepts","thresholds","residuals","means",
                  "residual.covariances","lv.variances","lv.covariances")
  params <- PT[paste0(PT$lhs, PT$op, PT$rhs) %in% setdiff(param, paramTypes), ]
  ## add parameters by type, if any are specified
  types <- intersect(param, paramTypes)
  ov.names <- con@Data@ov$name
  isOV <- PT$lhs %in% ov.names
  lv.names <- unique(PT$lhs[PT$op == "=~"])
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
  params <- params[!paste0(params$lhs, params$op, params$rhs) %in% freeParam, ]
  ## return modification indices for specified constraints (param)
  MIs <- lavaan::lavTestScore(con)$uni
  MIs[MIs$lhs %in% params$plabel, ]
}

## Function to find delta-AFIs AND maximum modification index in one permutation
permuteOnce <- function(i, d, con, uncon, null, param, freeParam, G,
                        AFIs, moreAFIs, maxSparse = 10, maxNonconv = 10) {
  nSparse <- 0L
  nTries <- 1L
  while ( (nSparse <= maxSparse) & (nTries <= maxNonconv) ) {
    ## permute grouping variable
    d[ , G] <- sample(d[ , G])

    ## for ordered indicators, check that groups have same observed categories
    ordVars <- con@Data@ov$name[con@Data@ov$type == "ordered"]
    if (length(ordVars) > 0) {
      onewayTables <- lavaan::lavTables(d, dimension = 1L,
                                categorical = ordVars, group = G)
      if (any(onewayTables$obs.prop == 1)) {
        nSparse <- nSparse + 1L
        next
      }
    }
    ## fit null model, if it exists
    if (!is.null(null)) {
      out.null <- lavaan::lavaan(data = d, group = G, slotParTable = null@ParTable,
                         slotOptions = null@Options)
    } else out.null <- NULL

    ## fit constrained model, check for convergence
    suppressWarnings(try(out0 <- lavaan::lavaan(data = d, group = G,
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
    suppressWarnings(try(out1 <- lavaan::lavaan(data = d, group = G,
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
      if (!is.null(uncon)) fit1[[1]] <- lavaan::fitMeasures(out1, fit.measures = AFIs,
                                                    baseline.model = out.null)
      fit0[[1]] <- lavaan::fitMeasures(out0, fit.measures = AFIs,
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
      AFI <- unlist(fit0) - unlist(fit1)
      MI <- max(getMIs(out0, , param = param, freeParam = freeParam)$X2)
    }
  }
  list(AFI = AFI, MI = MI, n.nonConverged = nTries - 1L, n.Sparse = nSparse)
}

## Function to permute difference in fits
permuteMeasEq <- function(nPermute, con, uncon = NULL, null = NULL, param = NULL,
                          freeParam = NULL, AFIs = NULL, moreAFIs = NULL,
                          maxSparse = 10, maxNonconv = 10, showProgress = TRUE) {
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
    if ("thresholds" %in% param | any(grepl("\\|", param))) {
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
    if (!is.null(uncon)) AFI1[[1]] <- lavaan::fitMeasures(uncon, fit.measures = AFIs,
                                                  baseline.model = null)
    AFI0[[1]] <- lavaan::fitMeasures(con, fit.measures = AFIs, baseline.model = null)
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
    AFI.obs <- unlist(AFI0) - unlist(AFI1)
    MI.obs <- getMIs(con, param = param, freeParam = freeParam)
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
                                    null = null, G = G,
                                    param = param, freeParam = freeParam,
                                    AFIs = AFIs, moreAFIs = moreAFIs,
                                    maxSparse = maxSparse, maxNonconv = maxNonconv)
      setTxtProgressBar(mypb, j)
    }
    close(mypb)
  } else {
    permuDist <- lapply(1:nPermute, permuteOnce, d = allData, con = con,
                        uncon = uncon, null = null, G = G, param = param,
                        freeParam = freeParam, AFIs = AFIs, moreAFIs = moreAFIs,
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
      if (b) return(mean(x >= y, na.rm = TRUE))
      mean(x <= y, na.rm = TRUE)
    }, x = unclass(AFI.dist), y = AFI.obs, b = badness)

  ## extract distribution of maximum modification indices
  MI.dist <- as.numeric(unlist(lapply(permuDist, function(x) x$MI)))
  ## calculate Tukey-adjusted p values for modification indices
  if (!is.null(uncon)) {
    MI.obs$tukey.p.value <- sapply(MI.obs$X2, function(i) mean(i <= MI.dist))
  }

  ## save parameter table for show/summary methods
  PT <- as.data.frame(lavaan::parTable(con))
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

  delta <- length(x@MI.dist) > 0L
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

