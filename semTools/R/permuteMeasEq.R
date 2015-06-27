### Terrence D. Jorgensen
### 6 April 2015
### permutation randomization test for measurement equivalence and DIF


## create s4 classes for result object
setClass("MeasEq.observed", representation(AFI = "vector",
                                           DIF = "matrix"))
setClass("MeasEq.p.values", representation(AFI = "vector",
                                           DIF.ss = "vector",
                                           DIF.all = "matrix",
                                           DIF.each = "matrix",
                                           DIF.pairs = "matrix"))
setClass("permuteMeasEq", representation(observed = "MeasEq.observed",
                                         p.values = "MeasEq.p.values",
                                         ANOVA = "vector",
                                         n.Permutations = "integer",
                                         n.Converged = "integer",
                                         n.nonConverged = "vector",
                                         n.Sparse = "vector",
                                         DIF.dist = "matrix",
                                         AFI.dist = "data.frame"))

## function to calculate DIF in specified parameters from an unconstrained model
calculateDIF <- function(uncon, param) {
  if (class(uncon) != "lavaan") stop("This function only applies to fitted lavaan models.")
  if (uncon@Data@ngroups == 1L) stop("This function only applies to multiple-group models.")
  ## save all estimates from less constrained model
  allCoefs <- lavaan::parameterEstimates(uncon)
  ## extract parameters of interest
  if (param[1] == "loadings") params <- allCoefs[allCoefs$op == "=~", c("lhs","op","rhs","group","est")]
  if (param[1] == "intercepts") params <- allCoefs[allCoefs$lhs %in% uncon@Data@ov$name & allCoefs$op == "~1",
                                                   c("lhs","op","rhs","group","est")]
  if (param[1] == "thresholds") params <- allCoefs[allCoefs$op == "|", c("lhs","op","rhs","group","est")]
  if (param[1] == "residuals") params <- allCoefs[allCoefs$lhs %in% uncon@Data@ov$name &
                                                    allCoefs$lhs == allCoefs$rhs & allCoefs$op == "~~",
                                                  c("lhs","op","rhs","group","est")]
  if (!param[1] %in% c("loadings","intercepts","thresholds","residuals")) {
    allNames <- mapply(function(x, y, z) paste0(x, y, z), x = allCoefs$lhs,
                       y = allCoefs$op, z = allCoefs$rhs)
    params <- allCoefs[allNames %in% param, ]
  }
  ## store estimates for each group in a list
  g1coefs <- params[params$group == 1L, ]
  g1.names <- mapply(function(x, y, z) paste0(x, y, z), x = g1coefs$lhs,
                     y = g1coefs$op, z = g1coefs$rhs)
  parList <- list(g1coefs$est)
  for (i in 2:length(table(allCoefs$group))) {
    parList <- c(parList, list(params$est[params$group == i]))
  }
  names(parList) <- uncon@Data@group.label

    ## all possible group comparisons, sorted so permutations are always in the same order
  comp <- combn(sort(uncon@Data@group.label), 2)
  ## function to compare one pair of groups
  parDiff <- function(groups = "", parList) {
    if (length(groups) != 2) stop("Must compare exactly 2 vectors of parameters")
    (parList[[ groups[1] ]] - parList[[ groups[2] ]])
  }
  ## Calculate observed DIF estimates for all comparisons from less constrained model
  diffs <- apply(comp, 2, parDiff, parList = parList)
  colnames(diffs) <- apply(comp, 2, paste, collapse = "-")
  rownames(diffs) <- g1.names
  diffs
}

## function to get the maximum sum-of-squared DIF across items
getMaxSS <- function(difmat) {
  if (is.null(nrow(difmat))) return(sum(difmat^2))
  max(apply(difmat, 1, function(x) sum(x^2)))
}

## function to find delta-AFIs AND maximum DIF in one permutation
permuteOnce <- function(i, d, uncon, con, null = NULL, param, G, diffs,
                        AFIs, moreAFIs, maxSparse = 10, maxNonconv = 10) {
  nSparse <- 0L
  nTries <- 1L
  while ( (nSparse <= maxSparse) & (nTries <= maxNonconv) ) {
    ## permute grouping variable
    d[ , G] <- sample(d[ , G])

    ## for ordered indicators, check that groups have same observed categories
    ordVars <- uncon@Data@ov$name[uncon@Data@ov$type == "ordered"]
    if (length(ordVars) > 0) {
      try(passTest <- lavaan::lavTables(d, dim = 1L, categorical = ordVars, group = G))
      if (!exists("passTest")) {
        nSparse <- nSparse + 1
        next
      }
    }
    ## fit null model, if it exists
    if (!is.null(null)) {
      out.null <- lavaan::lavaan(data = d, group = G, slotParTable = null@ParTable, slotOptions = null@Options)
    } else out.null <- NULL

    ## fit other models, check for convergence
    try(out1 <- lavaan::lavaan(data = d, group = G, slotParTable = uncon@ParTable, slotOptions = uncon@Options))
    if (! exists("out1")) {
      nTries <- nTries + 1
      next
    }
    if (!inspect(out1, "converged")) {
      nTries <- nTries + 1
      next
    }

    try(out0 <- lavaan::lavaan(data = d, group = G, slotParTable = con@ParTable, slotOptions = con@Options))
    if (! exists("out0")) {
      nTries <- nTries + 1
      next
    }
    if (!inspect(out0, "converged")) {
      nTries <- nTries + 1
      next
    }

    fit1 <- list()
    fit0 <- list()
    if (!is.na(AFIs[1])) {
      fit1[[1]] <- lavaan::fitMeasures(out1, fit.measures = AFIs, baseline.model = out.null)
      fit0[[1]] <- lavaan::fitMeasures(out0, fit.measures = AFIs, baseline.model = out.null)
    }
    if (!is.na(moreAFIs[1])) {
      fit1[[2]] <- moreFitIndices(out1, fit.measures = moreAFIs)
      fit0[[2]] <- moreFitIndices(out0, fit.measures = moreAFIs)
    }

    break
  } ## end WHILE loop
  ## if loop ended before getting results, return NA
  if ( (nSparse == maxSparse) | (nTries == maxNonconv) ) {
    allAFIs <- c(AFIs, moreAFIs)
    AFI <- rep(NA, sum(!is.na(allAFIs)))
    names(AFI) <- allAFIs[!is.na(allAFIs)]
    temp <- calculateDIF(uncon, param)[,1]
    temp[1:length(temp)] <- NA
    return(list(AFI = AFI, DIF = c(temp, all = NA, ss = NA), n.nonConverged = nTries, n.Sparse = nSparse))
  }

  permuted <- calculateDIF(out1, param)
  list(AFI = unlist(fit1) - unlist(fit0),
       DIF = c(apply(abs(permuted), 1, max),
               all = max(abs(permuted)),
               ss = getMaxSS(permuted)),
       DIF0 = abs(permuted) >= abs(diffs),
       n.nonConverged = nTries - 1L, n.Sparse = nSparse)
}

## function to permute difference in fits
permuteMeasEq <- function(nPermute, uncon, con, null = NULL, AFIs = NULL, moreAFIs = NULL,
                          param = "loadings", maxSparse = 10, maxNonconv = 10) {
  ## FIXME: Temporarily warn about testing thresholds; learn more about necessary constraints
  if (param == "thresholds" | any(grepl("\\|", param))) {
    stop("This function is not yet optimized for testing thresholds.
       Necessary identification contraints might not be specified.")
  }
  if (all(is.na(AFIs[1]), is.na(moreAFIs[1]))) warning("No AFIs were selected, so only the chi-squared statistic will be permuted.")

  nPermute <- as.integer(nPermute[1])
  ## logical check that models were fit to the same data
  # if (!all.equal(con@Data, uncon@Data)) stop("Models must be fit to the same groups")
  ## check for least-squares estimators
  leastSq <- any(c(uncon@Options$estimator, con@Options$estimator) %in% c("ULS","GLS","WLS", "DWLS"))

  ## specify default AFIs if NULL, or check validity of user-specified AFIs
  if (is.null(AFIs)) {
    AFIs <- c("cfi","rni","tli","rmsea","srmr","mfi")
    if (!leastSq) AFIs <- c(AFIs,"aic","bic","bic2")
  } else {
    IC <- grep("ic|ecvi|logl", AFIs, value = TRUE)
    if (leastSq & length(IC)) stop(paste("Argument 'AFIs' includes invalid options:",
                                         paste(IC, collapse = ", "),
                                         "Information criteria unavailable for least-squares estimators.",
                                         sep = "\n"))
  }
  ## specify default moreAFIs if NULL, or check validity of user-specified moreAFIs
  if (is.null(moreAFIs)) {
    moreAFIs <- c("gammaHat","adjGammaHat")
    if (!leastSq) moreAFIs <- c(moreAFIs,"aic.smallN","bic.priorN")
  } else {
    IC <- grep("ic|hqc", moreAFIs, value = TRUE)
    if (leastSq & length(IC)) stop(paste("Argument 'moreAFIs' includes invalid options:",
                                         paste(IC, collapse = ", "),
                                         "Information criteria unavailable for least-squares estimators.",
                                         sep = "\n"))
  }

  ###################### OBSERVED RESULTS ##########################
  ## save observed delta-AFIs
  AFI1 <- list()
  AFI0 <- list()
  if (!is.na(AFIs[1])) {
    AFI1[[1]] <- lavaan::fitMeasures(uncon, fit.measures = AFIs, baseline.model = null)
    AFI0[[1]] <- lavaan::fitMeasures(con, fit.measures = AFIs, baseline.model = null)
  }
  if (!is.na(moreAFIs[1])) {
    AFI1[[2]] <- moreFitIndices(uncon, fit.measures = moreAFIs)
    AFI0[[2]] <- moreFitIndices(con, fit.measures = moreAFIs)
  }
  AFI.diff <- unlist(AFI1) - unlist(AFI0)
  ## calculate observed DIF estimates
  diffs <- calculateDIF(uncon, param = param)
  obs <- new("MeasEq.observed", AFI = AFI.diff, DIF = diffs)

  ######################### PREP DATA ##############################
  ## save name of grouping variable
  G <- con@Data@group
  ## check for categorical variables
  # catVars <- which(con@Data@ov$type[!con@Data@ov$exo] == "ordered")
  numVars <- which(con@Data@ov$type[!con@Data@ov$exo] != "ordered")
  ## check whether mean-str or cov-str measurement parameters are being tested
  centerVars <- param %in% c("loadings","residuals") | any(grepl("~~", param)) | any(grepl("=~", param))
  scaleVars <- param == "intercepts" | any(grepl("~1", param))
  ## assemble data to which the models were fit
  dataList <- mapply(FUN = function(x, g, n, numVars, centerVars, scaleVars) {
    y <- data.frame(as.data.frame(x), g, stringsAsFactors = FALSE)
    names(y) <- c(n, con@Data@group)
    for (i in numVars) {
      if (centerVars) y[ , i] <- scale(y[ , i], scale = FALSE)[,1]
      if (scaleVars) y[ , i] <- scale(y[ , i], center = FALSE)[,1]
    }
    y
  }, x = con@Data@X, g = con@Data@group.label, n = con@Data@ov.names,
     numVars = numVars, centerVars = centerVars, scaleVars = scaleVars, SIMPLIFY = FALSE)
  allData <- do.call(rbind, dataList)

  ############ FIXME: remove the next 5 lines if standardizing does not work
  ## fit standardized data again
  # uncon <- lavaan(data = allData, group = uncon@Data@group, slotParTable = uncon@ParTable, slotOptions = uncon@Options)
  ## calculate observed DIF estimates
  # diffs <- calculateDIF(uncon, param = param)
  # obs <- new("MeasEq.observed", AFI = AFI.diff, DIF = diffs)

  ###################### PERMUTED RESULTS ###########################
  ## permute groups and return distributions of delta-AFIs and largest DIF
  permuDist <- lapply(1:nPermute, permuteOnce, d = allData, G = G, diffs = diffs,
                      AFIs = AFIs, moreAFIs = moreAFIs, uncon = uncon, con = con, null = null,
                      param = param, maxSparse = maxSparse, maxNonconv = maxNonconv)

  ## extract AFI distribution
  if (length(AFI.diff) > 1) {
    AFI.dist <- as.data.frame(t(sapply(permuDist, function(x) x$AFI)))
  }
  if (length(AFI.diff) == 1L) {
    AFI.dist <- data.frame(sapply(permuDist, function(x) x$AFI))
    colnames(AFI.dist) <- names(AFI.diff)
  }
  ## reverse-score badness-of-fit indices
  badness <- grep(pattern = "fmin|chi|aic|bic|rmr|rmsea|cn|sic|hqc",
                  x = names(AFI.diff), value = TRUE)
  for (i in badness) {
    AFI.diff[i] <- -AFI.diff[i]
    AFI.dist[ , i] <- -AFI.dist[ , i]
  }
  ## calculate all one-directional p-values
  p.AFI <- mapply(FUN = function(x, y) mean(x >= y, na.rm = TRUE),
                  x = unclass(AFI.dist), y = AFI.diff)

  ## extract distribution of maximum absolute DIF
  DIF.dist <- sapply(permuDist, function(x) x$DIF)

  ## calculate p values for DIF estimates
  p.all <- matrix(NA, nrow(diffs), ncol(diffs), dimnames = dimnames(diffs))
  p.each <- matrix(NA, nrow(diffs), ncol(diffs), dimnames = dimnames(diffs))
  p.ss <- rep(NA, times = nrow(diffs))
  names(p.ss) <- rownames(diffs)
  obs.ss <- apply(diffs, 1, function(x) sum(x^2))
  for (i in 1:nrow(diffs)) {
    p.ss[i] <- mean(DIF.dist["ss", ] >= obs.ss[i], na.rm = TRUE)
    for (j in 1:ncol(diffs)) {
      p.all[i, j] <- mean(DIF.dist["all", ] >= abs(diffs[i, j]), na.rm = TRUE)
      p.each[i, j] <- mean(DIF.dist[i, ] >= abs(diffs[i, j]), na.rm = TRUE)
    }
  }
  DIF.pairs <- Reduce("+", lapply(permuDist, function(x) x$DIF0)) / length(permuDist)
  #apply(simplify2array(lapply(tempDist, function(x) x$DIF0)), 1:2, mean)

  ## save all p-values
  pVals <- new("MeasEq.p.values", AFI = p.AFI, DIF.all = p.all,
               DIF.each = p.each, DIF.pairs = DIF.pairs, DIF.ss = p.ss)
  ## return observed results, permutation p values, and ANOVA results
  delta <- anova(uncon, con)
  ANOVA <- sapply(delta[,c("Chisq diff","Pr(>Chisq)")], function(x) x[2])
  out <- new("permuteMeasEq", observed = obs, ANOVA = ANOVA, p.values = pVals,
             n.Permutations = nPermute, n.Converged = sum(!is.na(AFI.dist[ , 1])),
             n.nonConverged = sapply(permuDist, function(x) x$n.nonConverged),
             n.Sparse = sapply(permuDist, function(x) x$n.Sparse),
             DIF.dist = DIF.dist, AFI.dist = AFI.dist)
  out
}


## methods
setMethod("show", "permuteMeasEq", function(object) {
  cat("Omnibus p values based on nonparametric permutation method: \n\n")
  print(data.frame(AFI_diff = object@observed@AFI, p.value = object@p.values@AFI))
  cat("\n\nOmnibus p value based on parametric method (chi-sq difference test): \n\n")
  print(object@ANOVA)
  cat("\n\nItemwise omnibus p values based on permutation distribution of maximum-SS-DIF per item: \n\n")
  print(object@p.values@DIF.ss)
  ## print warning if there are nonConverged permutations
  if (object@n.Permutations != object@n.Converged) {
    warning(paste("Only", object@n.Converged, "out of",
                  object@n.Permutations, "models converged within",
                  max(object@n.nonConverged), "attempts per permutation."))
  }
  invisible(object)
})

setMethod("summary", "permuteMeasEq", function(object, alpha = .05, digits = 3,
                                               type = c("all","each","pairs","step-up")) {
  printMeasEq <- function(mat, crit, digits = 3) {
    printMat <- matrix("", nrow(mat), ncol(mat), dimnames = dimnames(mat))
    printMat[crit] <- round(mat[crit], digits = digits)
    print(printMat, quote = FALSE)
    invisible(NULL)
  }

  type <- type[1]

  ## print matching omnibus test (for type == "pairs", print both)
  if (type %in% c("all","pairs","step-up")) {
    cat("\nOmnibus p values based on nonparametric permutation method: \n\n")
    print(data.frame(AFI_diff = object@observed@AFI, p.value = object@p.values@AFI))
  }
  if (type %in% c("each","pairs")) {
    cat("\nItemwise omnibus p values based on permutation distribution of maximum-SS-DIF per item: \n\n")
    print(object@p.values@DIF.ss)
  }

  ## extract requested DIF test results
  if (type == "step-up") {
    DIF <- object@p.values@DIF.pairs
    rankedAlpha <- matrix((alpha * rank(DIF)) / length(DIF),
                          nrow(DIF), ncol(DIF), dimnames = dimnames(DIF))
    sig <- which(DIF < rankedAlpha, arr.ind = TRUE)
  } else {
    DIF <- slot(object@p.values, paste0("DIF.", type))
    sig <- which(DIF < alpha, arr.ind = TRUE)
  }

  ## print requested DIF test results
  if (nrow(sig) == 0) {
    cat("\n\nThere are no significant differences between groups. \n\n")
    return(invisible(DIF < alpha))
  }
  cat("\n\nThere are significant differences between the following groups' parameters: \n\n")
  for (i in 1:nrow(sig)) {
    print(c(Groups = colnames(DIF)[[ sig[i, 2] ]],
            Parameter = rownames(DIF)[[ sig[i, 1] ]]))
    cat("\n")
  }
  cat("Significant DIF estimates:\n\n")
  printMeasEq(mat = object@observed@DIF, crit = sig, digits = digits)
  ## print warning if there are nonConverged permutations
  if (object@n.Permutations != object@n.Converged) {
    warning(paste("Only", object@n.Converged, "out of",
                  object@n.Permutations, "models converged within",
                  max(object@n.nonConverged), "attempts per permutation."))
  }
  invisible(DIF < alpha)
})

