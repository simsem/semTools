### Terrence D. Jorgensen
### Last updated: 24 February 2017
### source code for new runMI, extending lavaanList class instead of lavaanStar

cfa.mi <- function(model, data, ...,
                   m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  runMI(model = model, data = data, fun = "cfa", ...,
        m = m, miArgs = miArgs, miPackage = miPackage, seed = seed)
}

sem.mi <- function(model, data, ...,
                   m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  runMI(model = model, data = data, fun = "sem", ...,
        m = m, miArgs = miArgs, miPackage = miPackage, seed = seed)
}

growth.mi <- function(model, data, ...,
                      m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  runMI(model = model, data = data, fun = "growth", ...,
        m = m, miArgs = miArgs, miPackage = miPackage, seed = seed)
}

lavaan.mi <- function(model, data, ...,
                      m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
  runMI(model = model, data = data, fun = "lavaan", ...,
        m = m, miArgs = miArgs, miPackage = miPackage, seed = seed)
}

runMI <- function(model, data, fun = "lavaan", ...,
                  m, miArgs = list(), miPackage = "Amelia", seed = 12345) {
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
  getOutput <- function(obj) {
    converged <- lavaan::lavInspect(obj, "converged")
    if (converged) {
      se <- lavaan::parTable(obj)$se
      se.test <- all(!is.na(se)) & all(se >= 0) & any(se != 0)
      if (lavaan::lavInspect(obj, "ngroups") == 1L) {
        Heywood.lv <- det(lavaan::lavInspect(obj, "cov.lv")) > 0
        Heywood.ov <- det(lavaan::lavInspect(obj, "theta")) > 0
      } else {
        Heywood.lv <- all(sapply(lavaan::lavInspect(obj, "cov.lv"), det) > 0)
        Heywood.ov <- all(sapply(lavaan::lavInspect(obj, "theta"), det) > 0)
      }
    } else {
      se.test <- Heywood.lv <- Heywood.ov <- NA
    }
    list(sampstat = lavaan::lavInspect(obj, "sampstat"),
         coefMats = lavaan::lavInspect(obj, "coef"),
         GLIST = obj@Model@GLIST, # FIXME: @Model slot may disappear; need GLIST for std.all
         converged = converged, SE = se.test,
         Heywood.lv = Heywood.lv, Heywood.ov = Heywood.ov)
  }
  ## FIXME: in case of user-supplied FUN for lavaanList, combine with getOutput

  ## fit model using lavaanList
  lavListCall <- list(lavaan::lavaanList, model = model, dataList = imputedData,
                      cmd = fun)
  lavListCall <- c(lavListCall, dots)
  lavListCall$store.slots <- c("partable","vcov","test")
  lavListCall$FUN <- getOutput
  fit <- eval(as.call(lavListCall))
  ## Store custom @DataList and @SampleStatsList
  fit@SampleStatsList <- lapply(fit@funList, "[[", i = "sampstat")
  fit@DataList <- imputedData
  ## assign class and add new slots
  fit <- as(fit, "lavaan.mi")
  fit@coefList <- lapply(fit@funList, "[[", i = "coefMats")
  fit@seed <- seed
  fit@imputeCall <- imputeCall
  convList <- lapply(fit@funList, "[", i = c("converged","SE",
                                             "Heywood.lv","Heywood.ov"))
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
  keepNames <- funNames[ ! funNames %in% c("sampstat","coefMats","converged",
                                           "SE","Heywood.lv","Heywood.ov","GLIST")]
  fit@funList <- if (length(keepNames)) {
    lapply(fit@funList, "[[", i = keepNames) # FIXME: add user FUN to list as "other"?
    #lapply(fit@funList, "[[", i = "other")
  } else list()

  fit@ParTable$start <- getMethod("coef", "lavaan.mi")(fit, type = "user", labels = FALSE)
  fit
}


#############
## Methods ##
#############

## create s4 class for result object
setClass("lavaan.mi", contains = "lavaanList",
         slots = c(coefList = "list",   # coefficients in matrix format
                   GLIST = "list",      # list of pooled coefs in GLIST format
                   seed = "integer",    # seed set before running imputations
                   imputeCall = "list", # store call from imputation, if used
                   convergence = "list")) # also check SEs and Heywood cases
## lavaan.mi replaces lavaanStar

setMethod("show", "lavaan.mi", function(object) {
  nData <- object@meta$ndat
  nConverged <- sum(do.call(c, lapply(object@convergence, "[", "converged")))
  SE <- do.call(c, lapply(object@convergence, "[", "SE"))
  Heywood.ov <- do.call(c, lapply(object@convergence, "[", "Heywood.ov"))
  Heywood.lv <- do.call(c, lapply(object@convergence, "[", "Heywood.lv"))

  cat('lavaan.mi object based on ', nData, ' imputed data sets. \n',
      'See class?lavaan.mi help page for available methods. \n\n',
      'Convergence information:\n', 'The model converged on ',
      nConverged, ' imputed data sets \n\n', sep = "")

  if (!all(SE)) cat('Standard errors could not be computed for data set(s)',
                    paste(which(!SE), collapse = ", "), '\nTry fitting the',
                    'model to the individual data set(s) to diagnose',
                    'problems. If they cannot be fixed, try inspecting the',
                    'imputations. It may be necessary to reimpute the data',
                    'with some restrictions imposed. \n\n')

  if (!all(Heywood.ov & Heywood.lv))
    cat('Heywood cases detected for data set(s)',
        paste(which(!Heywood.ov | !Heywood.lv), collapse = ", "),
        '\nThese are not necessarily a cause for concern, unless a pooled',
        'estimate is also a Heywood case. \n\n')

  object
})
setMethod("summary", "lavaan.mi", function(object, se = TRUE, ci = TRUE, level = .95,
                                           standardized = FALSE, rsquare = FALSE,
                                           fmi = FALSE, add.attributes = TRUE) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  m <- sum(useImps)
  ## extract parameter table with attributes for printing
  PT <- lavaan::parTable(object)
  PE <- PT[ , c("lhs","op","rhs","block","group")]
  free <- PT$free > 0L | PT$op == ":="
  STDs <- !(PT$op %in% c("==","<",">"))
  PE$est <- rowMeans(sapply(object@ParTableList, "[[", i = "est"))

  if (lavaan::lavListInspect(object, "options")$se == "none") {
    warning('pooled variances and tests unavailable when se="none" is requested')
    se <- FALSE
  }
  if (!se) fmi <- FALSE
  messPool <- paste0("Rubin's (1987) rules were used to pool point",
                     if (se) " and SE",
                     " estimates across ", m, " imputed data sets",
                     if (se) ", and to calculate degrees of freedom for each",
                     if (se) " parameter's t test and CI.",
                     "\n")
  if (se) {
    PE$se <- lavaan::lav_model_vcov_se(object@Model, lavpartable = object@ParTable,
                                       VCOV = getMethod("vcov","lavaan.mi")(object))
    PE$t[free] <- PE$est[free] / PE$se[free]
    ## calculate df for t test
    W <- rowMeans(sapply(object@ParTableList, "[[", i = "se")^2)
    B <- apply(sapply(object@ParTableList, "[[", i = "est"), 1, var)
    Bm <- B + B/m
    Tot <- W + Bm
    ## can't do finite-sample correction because Wald z tests have no df (see Enders, 2010, p. 231, eq. 8.13 & 8.14)
    PE$df[free] <- (m - 1) * (1 + W[free] / Bm[free])^2
    ## if DF are obscenely large, set them to infinity for pretty printing
    PE$df <- ifelse(PE$df > 9999, Inf, PE$df)
    PE$pvalue <- pt(-abs(PE$t), df = PE$df)*2
    if (ci) {
      crit <- qt(1 - (1 - level) / 2, df = PE$df)
      PE$ci.lower <- PE$est - crit * PE$se
      PE$ci.upper <- PE$est + crit * PE$se
      PE$ci.lower[!free] <- PE$ci.upper[!free] <- PE$est[!free]
    }
  }

  if (is.logical(standardized)) {
    if (standardized) {
      PE$std.lv[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                      type = "std.lv",
                                                      GLIST = object@GLIST,
                                                      est = PE$est)$est.std
      PE$std.all[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                       type = "std.all",
                                                       GLIST = object@GLIST,
                                                       est = PE$est)$est.std
    }
  } else if (as.character(standardized)[1] == "std.lv") {
    PE$std.lv[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                    type = "std.lv",
                                                    GLIST = object@GLIST,
                                                    est = PE$est)$est.std
  } else if (as.character(standardized)[1] == "std.all") {
    PE$std.all[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                     type = "std.all",
                                                     GLIST = object@GLIST,
                                                     est = PE$est)$est.std
  }
  if (fmi) {
    PE$fmi1[free] <- Bm[free] / Tot[free]
    PE$fmi2[free] <- (Bm[free] + 2 / (PE$df[free] + 3)) / Tot[free]
    PE$riv[free] <- Bm[free] / W[free] # (Enders, 2010, p. 226, eq. 8.10)
    # == PE$riv[free] <- PE$fmi1[free] / (1 - PE$fmi1[free])
    messFMI <- paste("FMI(2) adjusts for using few imputations, but can",
                     "exceed 1 when the df are very large, making it",
                     "uninterpretable. The RIV will exceed 1 whenever",
                     "between-imputation variance exceeds",
                     "within-imputation variance (when FMI(1) > 50%).\n\n")
  }
  ## fancy or not?
  if (add.attributes) {
    PE$label <- PT$label
    PE$exo <- 0L # because PT$exo must be when !fixed.x
    class(PE) <- c("lavaan.parameterEstimates","lavaan.data.frame","data.frame")
    attr(PE, "information") <- lavaan::lavListInspect(object, "options")$information
    attr(PE, "se") <- lavaan::lavListInspect(object, "options")$se
    attr(PE, "group.label") <- lavaan::lavListInspect(object, "group.label")
    attr(PE, "missing") <- lavaan::lavListInspect(object, "options")$missing
    cat(messPool)
    if (fmi) cat("\n", messFMI, sep = "")
  } else {
    ## if not, attach header
    class(PE) <- c("lavaan.data.frame","data.frame")
    attr(PE, "header") <- if (fmi) c(messPool, "\n", messFMI) else messPool
  }
  ## requested R-squared?
  endoNames <- c(lavaan::lavNames(object, "ov.nox"),
                 lavaan::lavNames(object, "lv.nox"))
  if (rsquare & length(endoNames)) {
    isEndo <- sapply(PE$lhs, function(x) x %in% endoNames)
    rsqPE <- PE[PE$lhs == PE$rhs & PE$op == "~~" & isEndo, ]
    rsqPE$op <- "r2"
    for (i in which(!sapply(colnames(PE),
                            function(x) x %in% c("lhs","op","rhs","block","group","est","exo")))) {
      rsqPE[ , i] <- NA
    }
    STD <- lavaan::standardizedSolution(object, se = FALSE, type = "std.all",
                                        GLIST = object@GLIST, est = PE$est)
    isEndoSTD <- sapply(STD$lhs, function(x) x %in% endoNames)
    std.all <- STD$est.std[STD$lhs == STD$rhs & STD$op == "~~" & isEndoSTD]
    rsqPE$est <- ifelse(std.all < 0, NA, 1 - std.all) # negative variances
    if (add.attributes) rsqPE$label <- ""
    PE <- rbind(PE, rsqPE)
  }

  getMethod("show", "lavaan.mi")(object)
  if (!add.attributes) PE <- PE[!(PE$op %in% c("==","<",">")), ]
  rownames(PE) <- NULL
  PE
})

setMethod("nobs", "lavaan.mi", function(object, total = TRUE) {
  if (total) return(lavaan::lavListInspect(object, "ntotal"))
  N <- lavaan::lavListInspect(object, "norig")
  if (length(N) > 1L) names(N) <- lavaan::lavListInspect(object, "group.label")
  N
})

setMethod("coef", "lavaan.mi", function(object, type = "free", labels = TRUE) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  PT <- lavaan::parTable(object)
  if(type == "user" || type == "all") {
    type <- "user"
    idx <- 1:length(PT$lhs)
  } else if(type == "free") {
    ## FIXME: duplicated leftover from old way of handling EQ constraints?
    idx <- which(PT$free > 0L & !duplicated(PT$free))
  }
  ## extract coefficients for converged models
  coefList <- lapply(object@ParTableList[useImps], "[[", i = "est")
  out <- colMeans(do.call(rbind, coefList))[idx]
  ## attach names, set class
  if(labels) names(out) <- lavaan::lav_partable_labels(PT, type = type)
  class(out) <- c("lavaan.vector","numeric")
  out
})

setMethod("vcov", "lavaan.mi", function(object, type = c("pooled","between","within")) {
  if (lavaan::lavListInspect(object, "options")$se == "none") {
    warning('requested se="none", so only between-imputation (co)variance can',
            ' be computed')
    type <- "between"
  }
  PT <- lavaan::parTable(object)
  npar <- max(PT$free) - sum(PT$op == "==")
  useImps <- sapply(object@convergence, "[[", i = "converged")
  m <- sum(useImps)
  type <- tolower(type[1])

  useSE <- sapply(object@convergence, "[[", i = "SE")
  coefList <- lapply(object@ParTableList[useImps], "[[", i = "est")
  B <- cov(do.call(rbind, coefList)[ , PT$free > 0L & !duplicated(PT$free)])
  class(B) <- c("lavaan.matrix.symmetric","matrix")
  rownames(B) <- colnames(B) <- lavaan::lav_partable_labels(PT, type = "free")
  if (type == "between") return(B)

  W <- Reduce("+", lapply(object@vcovList[useSE], function(x) x$vcov)) / sum(useSE)
  class(W) <- c("lavaan.matrix.symmetric","matrix")
  dimnames(W) <- dimnames(B)
  if (type == "within") {
    return(W)
  } else if (type != "pooled") stop("'", type, "' is not a valid option for 'type'")

  if (!all(useImps == useSE))
    warning('Between-imputation covariance matrix based on estimated parameters',
            ' from ', m, ' converged solutions, but the mean within-imputation',
            ' covariance matrix based on ', sum(useSE), ' solutions for which',
            ' standard errors could be calculated.  Pooled total covariance',
            ' matrix is therefore based on different imputed data sets.')

  ## check whether equality constraints prevent inversion of W
  if (max(PT$free) == npar) {
    ## relative increase in variance due to missing data
    r <- (1 + 1/m)/npar * sum(diag(B %*% solve(W))) # Enders (2010, p. 235) eqs. 8.20-21
    Total <- (1 + r) * W # FIXME: asked Yves for a hack, says it can't be inverted back to infoMat
  } else {
    ## less reliable, but constraints prevent inversion of W
    Total <- W + B + (1/m)*B ## Enders (2010, p. 235) eq. 8.19
  }
  ## return pooled variance
  Total
})

## "borrowed" lavTestWald()
D1 <- function(object, constraints = NULL, asymptotic = FALSE, verbose = FALSE) {
  nImps <- sum(sapply(object@convergence, "[[", i = "converged"))
  if (nImps == 1L) stop("model did not converge on any imputations")
  if (is.null(constraints) || nchar(constraints) == 0L) stop("constraints are empty")

  # remove == constraints from parTable, save as list
  PT <- lavaan::parTable(object)
  partable <- as.list(PT[PT$op != "==", ])
  if (sum(PT$op == "==") > 0L) {
    message("When the unrestricted model already has equality constraints,",
            " D1 requires 'asymptotic = TRUE', so an approximate chi-squared",
            " will be returned instead of an F test statistic. \n")
    asymptotic <- TRUE
  }

  # parse constraints
  FLAT <- lavaan::lavParseModelString( constraints )
  CON <- attr(FLAT, "constraints")
  LIST <- list()
  if (length(CON) > 0L) {
    lhs <- unlist(lapply(CON, "[[", i = "lhs"))
    op <- unlist(lapply(CON, "[[", i = "op"))
    rhs <- unlist(lapply(CON, "[[", i = "rhs"))
    LIST$lhs <- c(LIST$lhs, lhs) # FIXME: why concatenate with NULL?
    LIST$op  <- c(LIST$op,  op)
    LIST$rhs <- c(LIST$rhs, rhs)
  } else stop("no equality constraints found in constraints argument")

  # theta = free parameters only (equality-constrained allowed)
  theta <- getMethod("coef", "lavaan.mi")(object) #object@optim$x

  # build constraint function
  ceq.function <- lavaan::lav_partable_constraints_ceq(partable = partable,
                                                       con = LIST, debug = FALSE)
  # compute jacobian restrictions
  JAC <- try(lavaan::lav_func_jacobian_complex(func = ceq.function, x = theta),
             silent = TRUE)
  if (inherits(JAC, "try-error")) { # eg. pnorm()
    JAC <- lavaan::lav_func_jacobian_simple(func = ceq.function, x = theta)
  }
  if (verbose) {cat("Restriction matrix (jacobian):\n"); print(JAC); cat("\n")}

  # linear restriction
  theta.r <- ceq.function( theta )
  if (verbose) {cat("Restricted theta values:\n"); print(theta.r); cat("\n")}

  # get VCOV
  VCOV <- getMethod("vcov", "lavaan.mi")(object)

  # restricted vcov
  info.r  <- JAC %*% VCOV %*% t(JAC)

  # Wald test statistic
  test.stat <- as.numeric(t(theta.r) %*% solve( info.r ) %*% theta.r)

  # number of constraints (k in Enders (2010, p. 235) eqs. 8.23-25)
  DF <- nrow(JAC)

  if (asymptotic) {
    out <- c("chisq" = test.stat * DF, df = DF,
             pvalue = pchisq(test.stat * DF, df = DF, lower.tail = FALSE))
  } else {
    npar <- max(PT$free) - sum(PT$op == "==")
    W <- getMethod("vcov", "lavaan.mi")(object, type = "within")
    B <- getMethod("vcov", "lavaan.mi")(object, type = "between")
    ## relative increase in variance due to missing data
    ariv <- (1 + 1/nImps) * sum(diag(B %*% solve(W))) / npar
    ########### FIXME: can't invert with equality constraints.
    ##                 Asked Yves for a hack, says it can't be done

    ## calculate denominator DF for F statistic
    a <- DF*(nImps - 1)
    if (a > 4) {
      v2 <- 4 + (a - 4) * (1 + (1 - 2/a)*(1 / ariv))^2 # Enders (eq. 8.24)
    } else {
      v2 <- a*(1 + 1/DF) * (1 + 1/ariv)^2 / 2 # Enders (eq. 8.25)
    }
    out <- c("F" = test.stat, df1 = DF, df2 = v2,
             pvalue = pf(test.stat, df1 = DF, df2 = v2, lower.tail = FALSE))
  }

  class(out) <- c("lavaan.vector","numeric")
  out
}
D2 <- function(object, h1 = NULL, asymptotic = FALSE,
               robust = FALSE, scaleshift = FALSE) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  nImps <- sum(useImps)
  ## check for robust test
  test <- if (robust) 2L else 1L ## only included for simulation studies
  #test <- ifelse(lavListInspect(object, "options")$test == "standard", 1L, 2L)

  ## pool Wald tests
  if (is.null(h1)) {
    DF <- mean(sapply(object@testList[useImps], function(x) x[[test]][["df"]]))
    w <- sapply(object@testList[useImps], function(x) x[[test]][["stat"]])
  } else {
    DF0 <- mean(sapply(object@testList[useImps], function(x) x[[test]][["df"]]))
    DF1 <- mean(sapply(h1@testList[useImps], function(x) x[[test]][["df"]]))
    DF <- DF0 - DF1
    w0 <- sapply(object@testList[useImps], function(x) x[[test]][["stat"]])
    w1 <- sapply(h1@testList[useImps], function(x) x[[test]][["stat"]])
    w <- w0 - w1
  }
  w_bar <- mean(w)
  ariv <- (1 + 1/nImps) * var(sqrt(w))
  test.stat <- (w_bar/DF - ((nImps + 1) * ariv / (nImps - 1))) / (1 + ariv)
  if (test.stat < 0) test.stat <- 0
  if (asymptotic) {
    out <- c("chisq" = test.stat * DF, df = DF,
             pvalue = pchisq(test.stat * DF, df = DF, lower.tail = FALSE))
  } else {
    v3 <- DF^(-3 / nImps) * (nImps - 1) * (1 + (1 / ariv))^2
    out <- c("F" = test.stat, df1 = DF, df2 = v3,
             pvalue = pf(test.stat, df1 = DF, df2 = v3, lower.tail = FALSE))
  }
  if (is.null(h1)) {
    PT <- lavaan::parTable(object)
    out <- c(out, npar = max(PT$free) - sum(PT$op == "=="),
             ntotal = lavaan::lavListInspect(object, "ntotal"))
  }

  class(out) <- c("lavaan.vector","numeric")
  out
}
getLLs <- function(object) {
  meanstructure <- lavaan::lavListInspect(object, "meanstructure")
  nG <- lavaan::lavListInspect(object, "ngroups")
  group <- lavaan::lavListInspect(object, "group")
  useImps <- sapply(object@convergence, "[[", i = "converged")
  nImps <- sum(useImps)
  m <- length(useImps)
  implied <- getMethod("fitted", "lavaan.mi")(object)
  ## Multiple groups?
  if (nG > 1L) {
    group.label <- lavaan::lavListInspect(object, "group.label")
    varnames <- lavaan::lavNames(object, group = 1:nG) # in case order differs?
    names(varnames) <- group.label
    S <- lapply(implied, "[[", i = "cov")
    if (meanstructure) {
      M <- lapply(implied, "[[", i = "mean")
    } else {
      M <- list()
      for (g in group.label) {
        M[[g]] <- Reduce("+", lapply(lapply(object@SampleStatsList[useImps],
                                            "[[", i = g),   # Group g's list of means
                                     "[[", i = "mean")) / nImps # average them
      }
    }
    LL <- numeric(length = m)
    for (i in 1:m) {
      if (!useImps[i]) next
      LLg <- numeric(length = nG)
      names(LLg) <- group.label
      dd <- object@DataList[[i]]
      for (g in group.label) {
        LLg[g] <- sum(apply(as.matrix(dd[ dd[,group] == g, varnames[[g]]]),
                            MARGIN = 1, FUN = mnormt::dmnorm, log = TRUE,
                            mean = M[[g]], varcov = unclass(S[[g]])))
      }
      LL[i] <- sum(LLg)
    }
  } else {
    varnames <- lavaan::lavNames(object)
    S <- implied$cov
    M <- if (meanstructure) implied$mean else {
      Reduce("+", lapply(object@SampleStatsList[useImps], "[[", i = "mean")) / nImps
    }
    LL <- numeric(length = m)
    for (i in 1:m) {
      if (!useImps[i]) next
      LL[i] <- sum(apply(as.matrix(object@DataList[[i]][ , varnames]),
                         MARGIN = 1, FUN = mnormt::dmnorm,
                         mean = M, varcov = unclass(S), log = TRUE))
    }
  }
  LL
}
D3 <- function(object, h1 = NULL, asymptotic = FALSE) {
  N <- lavaan::lavListInspect(object, "ntotal")
  nG <- lavaan::lavListInspect(object, "ngroups")
  group <- lavaan::lavListInspect(object, "group")
  useImps <- sapply(object@convergence, "[[", i = "converged")
  nImps <- sum(useImps)
  m <- length(object@testList)
  if (is.null(h1)) {
    DF <- object@testList[[ which(useImps)[1] ]][[1]][["df"]]
  } else {
    DF1 <- h1@testList[[ which(useImps)[1] ]][[1]][["df"]]
    DF0 <- object@testList[[ which(useImps)[1] ]][[1]][["df"]]
    DF <- DF0 - DF1
  }

  ## calculate m log-likelihoods under pooled H0 estimates
  LL0 <- getLLs(object)

  ## calculate m log-likelihoods under pooled H1 estimates
  if (is.null(h1)) {
    ## calculate log-likelihood under saturated model as alternative (H1)
    LL1 <- numeric(length = m)
    if (nG > 1L) {
      group.label <- lavaan::lavListInspect(object, "group.label")
      varnames <- lavaan::lavNames(object, group = 1:nG) # in case order changes?
      names(varnames) <- group.label
      Ns <- lavaan::lavListInspect(object, "nobs") # group sample sizes
      names(Ns) <- group.label
      S1 <- M1 <- list()
      for (g in group.label) {
        S1[[g]] <- Reduce("+", lapply(lapply(object@SampleStatsList[useImps],
                                             "[[", i = g),  # Group g's cov list
                                      "[[", i = "cov")) / nImps  # average them
        M1[[g]] <- Reduce("+", lapply(lapply(object@SampleStatsList[useImps],
                                             "[[", i = g),  # Group g's list of means
                                      "[[", i = "mean")) / nImps # average them
        if (lavaan::lavListInspect(object, "options")$sample.cov.rescale) {
          S1[[g]] <- S1[[g]] * (Ns[g] - 1) / Ns[g]
        }
      }
      ## within 1:m, iterate over 1:nG
      for (i in 1:m) {
        if (!useImps[i]) next
        LL1g <- numeric(length = nG)
        names(LL1g) <- group.label
        dd <- object@DataList[[i]]
        for (g in group.label) {
          LL1g[g] <- sum(apply(as.matrix(dd[ dd[,group] == g, varnames[[g]]]),
                               MARGIN = 1, FUN = mnormt::dmnorm, log = TRUE,
                               mean = M1[[g]], varcov = unclass(S1[[g]])))
        }
        LL1[i] <- sum(LL1g)
      }
    } else {
      varnames <- lavaan::lavNames(object)
      S1 <- Reduce("+", lapply(object@SampleStatsList[useImps], "[[", i = "cov")) / nImps
      if (lavaan::lavListInspect(object, "options")$sample.cov.rescale) S1 <- S1 * (N - 1) / N
      M1 <- Reduce("+", lapply(object@SampleStatsList[useImps], "[[", i = "mean")) / nImps
      for (i in 1:m) {
        if (!useImps[i]) next
        LL1[i] <- sum(apply(as.matrix(object@DataList[[i]][ , varnames]),
                            MARGIN = 1, FUN = mnormt::dmnorm,
                            mean = M1, varcov = unclass(S1), log = TRUE))
      }
    }
  } else LL1 <- getLLs(h1)

  ## calculate average of m LRTs
  LRT_con <- mean(-2*(LL0[useImps] - LL1[useImps]))
  ## average chisq across imputations
  LRT_bar <- mean(sapply(object@testList[useImps], function(x) x[[1]]$stat))
  ## calculate average relative increase in variance
  a <- DF*(nImps - 1)
  ariv <- ((nImps + 1) / a) * (LRT_bar - LRT_con)
  test.stat <- LRT_con / (DF*(1 + ariv))
  if (test.stat < 0) {
    message('Negative test statistic set to zero \n')
    test.stat <- 0
  }
  if (asymptotic) {
    out <- c("chisq" = test.stat * DF, df = DF,
             pvalue = pchisq(test.stat * DF, df = DF, lower.tail = FALSE))
  } else {
    ## F statistic
    if (a > 4) {
      v4 <- 4 + (a - 4) * (1 + (1 - (2 / a))*(1 / ariv))^2 # Enders (eq. 8.34)
    } else {
      v4 <- a*(1 + 1/DF)*(1 + 1/ariv)^2 / 2 # Enders (eq. 8.35)
      # v4 <- (DF + 1)*(m - 1)*(1 + (1 / ariv))^2 / 2 # Grund et al. (eq. 9)
    }
    out <- c("F" = test.stat, df1 = DF, df2 = v4,
             pvalue = pf(test.stat, df1 = DF, df2 = v4, lower.tail = FALSE))
  }
  ## add log-likelihood and AIC/BIC for target model
  if (is.null(h1)) {
    PT <- lavaan::parTable(object)
    npar <- max(PT$free) - sum(PT$op == "==")
    if (lavaan::lavListInspect(object, "options")$sample.cov.rescale) N <- N - nG
    out <- c(out, npar = npar, ntotal = lavaan::lavListInspect(object, "ntotal"),
             logl = mean(LL0), unrestricted.logl = mean(LL1),
             aic = -2*mean(LL0) + 2*npar, bic = -2*mean(LL0) + npar*log(N),
             bic2 = -2*mean(LL0) + npar*log((N + 2) / 24))
    ## NOTE: Mplus reports the average of m likelihoods evaluated at the
    ##       m point estimates, not evaluated at the pooled point estimates.
    ##       Mplus also uses those to calcluate AIC and BIC.
  }

  class(out) <- c("lavaan.vector","numeric")
  out
}
robustify <- function(ChiSq, object, h1 = NULL) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  scaleshift <- lavaan::lavListInspect(object, "options")$test == "scaled.shifted"

  d0 <- mean(sapply(object@testList[useImps], function(x) x[[2]][["df"]]))
  c0 <- mean(sapply(object@testList[useImps],
                    function(x) x[[2]][["scaling.factor"]]))
  if (!is.null(h1)) {
    d1 <- mean(sapply(h1@testList[useImps], function(x) x[[2]][["df"]]))
    c1 <- mean(sapply(h1@testList[useImps],
                      function(x) x[[2]][["scaling.factor"]]))
    delta_c <- (d0*c0 - d1*c1) / (d0 - d1)
    ChiSq["chisq.scaled"] <- ChiSq[["chisq"]] / delta_c
    ChiSq["df.scaled"] <- d0 - d1
    ChiSq["pvalue.scaled"] <- pchisq(ChiSq[["chisq.scaled"]],
                                     df = ChiSq[["df.scaled"]],
                                     lower.tail = FALSE)
    ChiSq["chisq.scaling.factor"] <- delta_c
  } else {
    ChiSq["chisq.scaled"] <- ChiSq[["chisq"]] / c0
    ChiSq["df.scaled"] <- d0
    if (scaleshift) {
      ## add shift parameter here (copy from below) or below
      shift <- mean(sapply(object@testList[useImps],
                           function(x) x[[2]][["shift.parameter"]]))
      ChiSq["chisq.scaled"] <- ChiSq[["chisq.scaled"]] + shift
      ChiSq["pvalue.scaled"] <- pchisq(ChiSq[["chisq.scaled"]],
                                       df = ChiSq[["df.scaled"]],
                                       lower.tail = FALSE)
      ChiSq["chisq.scaling.factor"] <- c0
      ChiSq["chisq.shift.parameter"] <- shift
    } else {
      ChiSq["pvalue.scaled"] <- pchisq(ChiSq[["chisq.scaled"]],
                                       df = ChiSq[["df.scaled"]],
                                       lower.tail = FALSE)
      ChiSq["chisq.scaling.factor"] <- c0
    }
  }
  ChiSq
}
setMethod("anova", "lavaan.mi", function(object, h1 = NULL,
                                         test = c("D3","D2","D1"),
                                         asymptotic = FALSE, constraints = NULL,
                                         indices = FALSE, baseline = NULL) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  nImps <- sum(useImps)
  ## check class
  if (!is(object, "lavaan.mi")) stop("object is not class 'lavaan.mi'")
  if (!is.null(h1) & !is(object, "lavaan.mi")) stop("h1 is not class 'lavaan.mi'")

  ## Everything else obsolete if test = "D1"
  if (toupper(test[1]) == "D1") {
    if (!asymptotic) asymptotic <- TRUE ## FIXME: until W can be inverted with eq. constraints
    out <- D1(object = object, constraints = constraints, asymptotic = asymptotic)
    message('D1 (Wald test) calculated using pooled "',
            lavaan::lavListInspect(object, "options")$se,
            '" asymptotic covariance matrix of model parameters')
    return(out)
  }

  ## check for robust
  robust <- lavaan::lavListInspect(object, "options")$test != "standard"
  if (robust) asymptotic <- TRUE
  scaleshift <- lavaan::lavListInspect(object, "options")$test == "scaled.shifted"
  if (scaleshift & !is.null(h1)) stop("Robust correction unavailable for model",
                                      " comparison when test = 'scaled.shifted'")
  ################### FIXME: unless possible to mimic DIFFTEST behavior?


  ## check request for fit indices
  incremental <- c("cfi","tli","nnfi","rfi","nfi","pnfi","ifi","rni")
  if (is.logical(indices)) {
    moreFit <- is.null(h1) & indices
    if (moreFit) indices <- c("cfi","tli","rmsea","srmr")
  } else if (is.character(indices)) {
    indices <- tolower(indices)
    moreFit <- is.null(h1) & indices %in% c(incremental, "all","mfi","rmsea",
                                            "gammaHat","rmr","srmr")
    if (moreFit & any(indices == "all")) {
      indices <- c(incremental, "mfi","rmsea","gammaHat","rmr","srmr")
    }
  } else indices <- moreFit <- FALSE
  ## fit baseline model if necessary
  if (moreFit & any(indices %in% incremental)) {
    if (is.null(baseline)) {
      PTb <- lavaan::lav_partable_independence(lavdata = object@Data,
                         lavoptions = lavaan::lavListInspect(object, "options"))
      baseFit <- runMI(model = PTb, data = object@DataList[useImps],
                       group = lavaan::lavListInspect(object, "group"),
                       se = "none", # to save time
                       test = lavaan::lavListInspect(object, "options")$test,
                       estimator = lavaan::lavListInspect(object, "options")$estimator,
                       ordered = lavaan::lavListInspect(object, "ordered"),
                       parameterization = lavaan::lavListInspect(object,
                                                                 "parameterization"))
    } else if (!is(baseline, "lavaan.mi")) {
      stop('User-supplied baseline model must be "lavaan.mi" class fit',
           ' to the same imputed data')
    } else baseFit <- baseline
    baseImps <- sapply(baseFit@convergence, "[[", i = "converged")
    if (!all(baseImps)) warning('baseline model did not converge for data set(s): ',
                                which(useImps)[!baseImps])
  }

  ## check DF
  DF0 <- object@testList[[ which(useImps)[1] ]][[1]][["df"]]
  if (!is.null(h1)) {
    if (!is(h1, "lavaan.mi")) stop("h1 is not class 'lavaan.mi'")
    DF1 <- h1@testList[[ which(useImps)[1] ]][[1]][["df"]]
    if (DF0 == DF1) stop("models have the equal degrees of freedom")
    if (DF0 < DF1) {
      H0 <- h1
      h1 <- object
      object <- H0
      H0 <- DF1
      DF1 <- DF0
      DF0 <- H0
    }
    DF <- DF0 - DF1
  } else DF <- DF0
  if (DF == 0) indices <- moreFit <- FALSE # arbitrary perfect fit, no indices
  if (moreFit) asymptotic <- TRUE

  ## check test options, backward compatibility?
  if (tolower(test[1]) == "mplus") {
    test <- "D3"
    asymptotic <- TRUE
  }
  if (tolower(test[1]) %in% c("mr","meng.rubin","likelihood","lrt")) test <- "D3"
  if (tolower(test[1]) %in% c("lmrr","li.et.al","pooled.wald")) test <- "D2"
  if (toupper(test[1]) == "D3" & !lavaan::lavListInspect(object, "options")$estimator %in% c("ML","PML","FML")) {
    message('"D3" only available using maximum likelihood estimation. ',
            'Changed test to "D2".')
    test <- "D2"
  }
  ## calculate pooled test
  if (toupper(test[1]) == "D3") {
    ## check estimator
    if (lavaan::lavListInspect(object, "options")$estimator != "ML")
      stop("D3 is only available using ML estimation")
    out <- D3(object = object, h1 = h1, asymptotic = asymptotic)
    if (any(indices %in% incremental)) baseOut <- D3(baseFit, asymptotic = TRUE)
  } else if (toupper(test[1]) == "D2") {
    out <- D2(object = object, h1 = h1, asymptotic = asymptotic)
    if (any(indices %in% incremental)) baseOut <- D2(baseFit, asymptotic = TRUE)
  }
  ## If test statistic is negative, return without any indices or robustness
  if (asymptotic & (moreFit | robust)) {
    if (out[["chisq"]] == 0) {
      message('Negative test statistic set to zero, so fit will appear to be ',
              'arbitrarily perfect.  Robust corrections and additional fit ',
              'indices are not returned because they are uninformative.\n')
      class(out) <- c("lavaan.vector","numeric")
      return(out)
    }
  }

  ## add robust statistics
  if (robust) {
    out <- robustify(ChiSq = out, object, h1)
    if (scaleshift) {
      extraWarn <- ' and shift parameter'
    } else if (lavaan::lavListInspect(object, "options")$test == "mean.var.adjusted") {
      extraWarn <- ' and degrees of freedom'
    } else extraWarn <- ''
    message('Robust corrections are made to the naive (pooled) chi-squared',
            ' test statistic using the mean scaling factor', extraWarn,
            ' across ', nImps, ' imputations for which the model converged. \n')
  }

  ## add fit indices for single model
  if (moreFit) {
    X2 <- out[["chisq"]]
    if (robust) {
      X2.sc <- out[["chisq.scaled"]]
      DF.sc <- out[["df.scaled"]] ## for mean.var.adjusted, mean DF across imputations
      ch <- out[["chisq.scaling.factor"]] ## mean c_hat across imputations
      if (X2 < .Machine$double.eps && DF == 0) ch <- 0
      ## for RMSEA
      if ("rmsea" %in% indices) {
        d <- mean(sapply(object@testList[useImps],
                         function(x) sum(x[[2]][["trace.UGamma"]])))
        if (is.na(d) || d == 0) d <- NA # FIXME: only relevant when scaleshift?
      }
    }
    ## for CFI, TLI, etc.
    if (any(indices %in% incremental)) {
      bX2 <- baseOut[["chisq"]]
      bDF <- baseOut[["df"]]
      out <- c(out, baseline.chisq = bX2, baseline.df = bDF,
               baseline.pvalue = baseOut[["pvalue"]])
      if (robust) {
        baseOut <- robustify(ChiSq = baseOut, object = baseFit)
        out["baseline.chisq.scaled"] <- bX2.sc <- baseOut[["chisq.scaled"]]
        out["baseline.df.scaled"]    <- bDF.sc <- baseOut[["df.scaled"]]
        out["baseline.pvalue.scaled"] <- baseOut[["pvalue.scaled"]]
        cb <- baseOut[["chisq.scaling.factor"]]
        out["baseline.chisq.scaling.factor"] <- cb
      }
    }
  }
  if ("cfi" %in% indices) {
    t1 <- max(X2 - DF, 0)
    t2 <- max(X2 - DF, bX2 - bDF, 0)
    out["cfi"] <- if(t1 == 0 && t2 == 0) 1 else 1 - t1/t2
    if (robust) {
      ## scaled
      t1 <- max(X2.sc - DF.sc, 0)
      t2 <- max(X2.sc - DF.sc, bX2.sc - bDF.sc, 0)
      if (is.na(t1) || is.na(t2)) {
        out["cfi.scaled"] <- NA
      } else if (t1 == 0 && t2 == 0) {
        out["cfi.scaled"] <- 1
      } else out["cfi.scaled"] <- 1 - t1/t2
      ## Brosseau-Liard & Savalei MBR 2014, equation 15
      if (lavaan::lavListInspect(object, "options")$test %in%
          c("satorra.bentler","yuan.bentler")) {
        t1 <- max(X2 - ch*DF, 0)
        t2 <- max(X2 - ch*DF, bX2 - cb*bDF, 0)
        if (is.na(t1) || is.na(t2)) {
          out["cfi.robust"] <- NA
        } else if (t1 == 0 && t2 == 0) {
          out["cfi.robust"] <- 1
        } else out["cfi.robust"] <- 1 - t1/t2
      }
    }
  }
  if ("rni" %in% indices) {
    t1 <- X2 - DF
    t2 <- bX2 - bDF
    out["rni"] <- if (t2 == 0) NA else 1 - t1/t2
    if (robust) {
      ## scaled
      t1 <- X2.sc - DF.sc
      t2 <- bX2.sc - bDF.sc
      if (is.na(t1) || is.na(t2)) {
        out["rni.scaled"] <- NA
      } else if (t2 == 0) {
        out["rni.scaled"] <- NA
      } else out["rni.scaled"] <- 1 - t1/t2
      ## Brosseau-Liard & Savalei MBR 2014, equation 15
      if (lavaan::lavListInspect(object, "options")$test %in%
          c("satorra.bentler","yuan.bentler")) {
        t1 <- X2 - ch*DF
        t2 <- bX2 - cb*bDF
        if (is.na(t1) || is.na(t2)) {
          out["rni.robust"] <- NA
        } else if (t1 == 0 && t2 == 0) {
          out["rni.robust"] <- NA
        } else out["rni.robust"] <- 1 - t1/t2
      }
    }
  }
  if (any(indices %in% c("tli","nnfi"))) {
    t1 <- (X2 - DF)*bDF
    t2 <- (bX2 - bDF)*DF
    out["tli"] <- out["nnfi"] <- if (DF > 0) 1 - t1/t2 else 1
    if (robust) {
      ## scaled
      t1 <- (X2.sc - DF.sc)*bDF.sc
      t2 <- (bX2.sc - bDF.sc)*DF.sc
      if (is.na(t1) || is.na(t2)) {
        out["tli.scaled"] <- out["nnfi.scaled"] <- NA
      } else if (DF > 0 && t2 != 0) {
        out["tli.scaled"] <- out["nnfi.scaled"] <- 1 - t1/t2
      } else {
        out["tli.scaled"] <- out["nnfi.scaled"] <- 1
      }
      ## Brosseau-Liard & Savalei MBR 2014, equation 15
      if (lavaan::lavListInspect(object, "options")$test %in%
          c("satorra.bentler","yuan.bentler")) {
        t1 <- (X2 - ch*DF)*bDF
        t2 <- (bX2 - cb*bDF)*DF
        if (is.na(t1) || is.na(t2)) {
          out["tli.robust"] <- out["nnfi.robust"] <- NA
        } else if (t1 == 0 && t2 == 0) {
          out["tli.robust"] <- out["nnfi.robust"] <- 1 - t1/t2
        } else out["tli.robust"] <- out["nnfi.robust"] <- 1
      }
    }
  }
  if ("rfi" %in% indices) {
    if (DF > 0) {
      t2 <- bX2 / bDF
      t1 <- t2 - X2/DF
      out["rfi"] <- if (t1 < 0 || t2 < 0) 1 else t1/t2
    } else out["rfi"] <- 1
    if (robust) {
      if (DF > 0) {
        t2 <- bX2.sc / bDF.sc
        t1 <- t2 - X2.sc/DF.sc
        out["rfi.scaled"] <- if (t1 < 0 || t2 < 0) 1 else t1/t2
      } else out["rfi.scaled"] <- 1
    }
  }
  if ("nfi" %in% indices) {
    if (DF > 0) {
      t1 <- bX2 - X2
      t2 <- bX2
      out["nfi"] <- t1 / t2
    } else out["nfi"] <- 1
    if (robust) out["nfi.scaled"] <- (bX2.sc - X2.sc) / bX2.sc
  }
  if ("pnfi" %in% indices) {
    t1 <- bX2 - X2
    t2 <- bX2
    out["pnfi"] <- (DF / bDF) * t1/t2
    if (robust) {
      t1 <- bX2.sc - X2.sc
      t2 <- bX2.sc
      out["pnfi.scaled"] <- (DF / bDF) * t1/t2
    }
  }
  if ("ifi" %in% indices) {
    t1 <- bX2 - X2
    t2 <- bX2 - DF
    out["ifi"] <- if (t2 < 0) 1 else t1/t2
    if (robust) {
      t1 <- bX2.sc - X2.sc
      t2 <- bX2.sc - DF.sc
      if (is.na(t2)) {
        out["ifi.scaled"] <- NA
      } else if (t2 < 0) {
        out["ifi.scaled"] <- 1
      } else out["ifi.scaled"] <- t1/t2
    }
  }

  N <- lavaan::lavListInspect(object, "ntotal")
  Ns <- lavaan::lavListInspect(object, "nobs")
  nG <- lavaan::lavListInspect(object, "ngroups")
  nVars <- length(lavaan::lavNames(object))
  if (!(lavaan::lavListInspect(object, "options")$likelihood == "normal" |
        lavaan::lavListInspect(object, "options")$estimator %in% c("ML","PML","FML"))) {
    N <- N - nG
    Ns <- Ns - 1
  }

  if ("mfi" %in% indices) {
    out["mfi"] <- exp(-0.5 * (X2 - DF) / N)
  }

  if ("rmsea" %in% indices) {
    N.RMSEA <- max(N, X2*4) # FIXME: good strategy??

    if (is.na(X2) || is.na(DF)) {
      out["rmsea"] <- as.numeric(NA)
    } else if (DF > 0) {
      getLambda <- function(lambda, chi, df, p) pchisq(chi, df, ncp=lambda) - p

      out["rmsea"] <- sqrt( max(0, (X2/N)/DF - 1/N) ) * sqrt(nG)
      ## lower confidence limit
      if (getLambda(0, X2, DF, .95) < 0.0) out["rmsea.ci.lower"] <- 0 else {
        lambda.l <- try(uniroot(f = getLambda, chi = X2, df = DF, p = .95,
                                lower = 0, upper = X2)$root, silent = TRUE)
        if (inherits(lambda.l, "try-error")) lambda.l <- NA
        out["rmsea.ci.lower"] <- sqrt( lambda.l/(N*DF) ) * sqrt(nG)
      }
      ## upper confidence limit
      if (getLambda(N.RMSEA, X2, DF, .05) > 0 || getLambda(0, X2, DF, .05) < 0) {
        out["rmsea.ci.upper"] <- 0
      } else {
        lambda.u <- try(uniroot(f = getLambda, chi = X2, df = DF, p = .05,
                                lower = 0, upper = N.RMSEA)$root, silent = TRUE)
        if (inherits(lambda.u, "try-error")) lambda.u <- NA
        out["rmsea.ci.upper"] <- sqrt( lambda.u/(N*DF) ) * sqrt(nG)
      }
      ## p value
      out["rmsea.pvalue"] <- pchisq(X2, DF, ncp = N*DF*0.05^2/nG,
                                    lower.tail = FALSE)

      ## Scaled versions (naive and robust)
      if (robust & !scaleshift) {
        ## naive
        out["rmsea.scaled"] <- sqrt( max(0, (X2/N)/d - 1/N) ) * sqrt(nG)
        ## lower confidence limit
        if (DF.sc < 1 | getLambda(0, X2, DF.sc, .95) < 0.0) {
          out["rmsea.ci.lower.scaled"] <- 0
         } else {
          lambda.l <- try(uniroot(f = getLambda, chi = X2, df = DF.sc, p = .95,
                                  lower = 0, upper = X2)$root, silent = TRUE)
          if (inherits(lambda.l, "try-error")) lambda.l <- NA
          out["rmsea.ci.lower.scaled"] <- sqrt( lambda.l/(N*DF) ) * sqrt(nG)
        }
        ## upper confidence limit
        if (DF.sc < 1 | getLambda(N.RMSEA, X2, DF.sc, .05) > 0.0) {
          out["rmsea.ci.upper.scaled"] <- 0
        } else {
          lambda.u <- try(uniroot(f = getLambda, chi = X2, df = DF.sc, p = .05,
                                  lower = 0, upper = N.RMSEA)$root, silent = TRUE)
          if (inherits(lambda.u, "try-error")) lambda.u <- NA
          out["rmsea.ci.upper.scaled"] <- sqrt( lambda.u/(N*DF) ) * sqrt(nG)
        }
        ## p value
        out["rmsea.pvalue.scaled"] <- pchisq(X2, DF.sc, ncp = N*DF.sc*0.05^2/nG,
                                             lower.tail = FALSE)

        if (object@Options$test %in% c("satorra.bentler","yuan.bentler")) {
          ## robust
          out["rmsea.robust"] <- sqrt( max(0, (X2/N)/DF - ch/N ) ) * sqrt(nG)
          ## lower confidence limit
          if (DF.sc < 1 | getLambda(0, X2.sc, DF.sc, .95) < 0.0) {
            out["rmsea.ci.lower.robust"] <- 0
          } else {
            lambda.l <- try(uniroot(f = getLambda, chi = X2.sc, df = DF.sc, p = .95,
                                    lower = 0, upper = X2)$root, silent = TRUE)
            if (inherits(lambda.l, "try-error")) lambda.l <- NA
            out["rmsea.ci.lower.robust"] <- sqrt( (ch*lambda.l)/(N*DF.sc) ) * sqrt(nG)
          }
          ## upper confidence limit
          if (DF.sc < 1 | getLambda(N.RMSEA, X2.sc, DF.sc, .05) > 0.0) {
            out["rmsea.ci.upper.robust"] <- 0
          } else {
            lambda.u <- try(uniroot(f = getLambda, chi = X2.sc, df = DF.sc, p = .05,
                                    lower = 0, upper = N.RMSEA)$root, silent = TRUE)
            if (inherits(lambda.u, "try-error")) lambda.u <- NA
            out["rmsea.ci.upper.robust"] <- sqrt( (ch*lambda.u)/(N*DF.sc) ) * sqrt(nG)
          }
          ## p value
          ########## To be discovered?
        }
      } else if (scaleshift) {
        ## naive only
        out["rmsea.scaled"] <- sqrt( max(0, (X2.sc/N)/DF - 1/N) ) * sqrt(nG)
        ## lower confidence limit
        if (DF.sc < 1 | getLambda(0, X2.sc, DF.sc, .95) < 0.0) {
          out["rmsea.ci.lower.scaled"] <- 0
        } else {
          lambda.l <- try(uniroot(f = getLambda, chi = X2.sc, df = DF.sc, p = .95,
                                  lower = 0, upper = X2.sc)$root, silent = TRUE)
          if (inherits(lambda.l, "try-error")) lambda.l <- NA
          out["rmsea.ci.lower.scaled"] <- sqrt( lambda.l/(N*DF.sc) ) * sqrt(nG)
        }
        ## upper confidence limit
        if (DF.sc < 1 | getLambda(N.RMSEA, X2.sc, DF.sc, .05) > 0.0) {
          out["rmsea.ci.upper.scaled"] <- 0
        } else {
          lambda.u <- try(uniroot(f = getLambda, chi = X2.sc, df = DF.sc, p = .05,
                                  lower = 0, upper = N.RMSEA)$root, silent = TRUE)
          if (inherits(lambda.u, "try-error")) lambda.u <- NA
          out["rmsea.ci.upper.scaled"] <- sqrt( lambda.u/(N*DF.sc) ) * sqrt(nG)
        }
        ## p value
        out["rmsea.pvalue.scaled"] <- pchisq(X2.sc, DF.sc, ncp = N*DF.sc*0.05^2/nG,
                                             lower.tail = FALSE)
      }
    }
  }

  if ("gammaHat" %in% indices) {
    out["gammaHat"] <- nVars / (nVars + 2*((X2 - DF) / N))
    out["adjGammaHat"] <- 1 - (((nG * nVars * (nVars + 1)) / 2) / DF) * (1 - out["gammaHat"])
    if (robust) {
      out["gammaHat.scaled"] <- nVars / (nVars + 2*((X2.sc - DF.sc) / N))
      out["adjGammaHat.scaled"] <- 1 - (((nG * nVars * (nVars + 1)) / 2) / DF.sc) * (1 - out["gammaHat.scaled"])
    }
  }

  getSRMR <- function(object, type) {
    vv <- lavaan::lavNames(object, type = "ov.num")
    R <- getMethod("resid", "lavaan.mi")(object, type = type)
    index <- if (type == "raw") "cov" else "cor"
    if (nG > 1L) {
      RR <- list()
      for (g in 1:nG) {
        RR[[g]] <- c(R[[g]][[index]][lower.tri(R[[g]][[index]], diag = FALSE)]^2,
                     diag(R[[g]][[index]])[vv]^2)
      }
    } else RR <- c(R[[index]][lower.tri(R[[index]], diag = FALSE)]^2,
                   diag(R[[index]])[vv]^2)

    if (lavaan::lavListInspect(object, "meanstructure")) {
      if (nG > 1L) {
        for (g in 1:nG) RR[[g]] <- c(RR[[g]], R[[g]]$mean[vv]^2)
      } else RR <- c(RR, R$mean[vv]^2)
    }

    SS <- if (nG > 1L) sqrt(sapply(RR, mean)) else sqrt(mean(RR))
    as.numeric( (lavaan::lavListInspect(object, "nobs") %*% SS) / lavaan::lavListInspect(object, "ntotal") )
  }
  if("rmr" %in% indices) out["rmr"] <- getSRMR(object, type = "raw")
  if("srmr" %in% indices) {
    out["srmr_bollen"] <- getSRMR(object, type = "cor.bollen")
    out["srmr_bentler"] <- getSRMR(object, type = "cor.bentler")
  }

  class(out) <- c("lavaan.vector","numeric")
  out # FIXME: in future, accept more than 2 models, arrange sequentially by DF
})

## function to pool each group's list of sample stats
sampstat.lavaan.mi <- function(lst, means = FALSE, categ = FALSE, m = m) {
  ## average sample stats across imputations
  out <- list(cov = Reduce("+", lapply(lst, "[[", i = "cov")) / m)
  if (means) out$mean <- Reduce("+", lapply(lst, "[[", i = "mean")) / m
  if (categ) out$th <- Reduce("+", lapply(lst, "[[", i = "th")) / m
  out
}
fitted.lavaan.mi <- function(object) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  m <- sum(useImps)
  meanstructure <- lavaan::lavListInspect(object, "meanstructure")
  categ <- lavaan::lavListInspect(object, "categorical")
  nG <- lavaan::lavListInspect(object, "ngroups")
  ov.names <- lavaan::lavNames(object)

  est <- getMethod("coef", "lavaan.mi")(object)
  imp <- lavaan::lav_model_implied(lavaan::lav_model_set_parameters(object@Model,
                                                                    x = est))
  out <- list()
  if (nG > 1L) {
    group.label <- lavaan::lavListInspect(object, "group.label")
    for (i in seq_along(imp)) names(imp[[i]]) <- group.label
    for (g in group.label) {
      out[[g]]$cov <- imp$cov[[g]]
      dimnames(out[[g]]$cov) <- list(ov.names, ov.names)
      class(out[[g]]$cov) <- c("lavaan.matrix.symmetric","matrix")
      if (meanstructure) {
        out[[g]]$mean <- as.numeric(imp$mean[[g]])
        names(out[[g]]$mean) <- ov.names
        class(out[[g]]$mean) <- c("lavaan.vector","numeric")
      } else {
        out[[g]]$mean <- sampstat.lavaan.mi(lapply(object@SampleStatsList[useImps], "[[", g),
                                            means = TRUE, categ = categ, m = m)$mean
      }
      if (categ) {
        out[[g]]$th <- imp$th[[g]]
        names(out[[g]]$th) <- lavaan::lavNames(object, "th")
        class(out[[g]]$th) <- c("lavaan.vector","numeric")
      }
    }
  } else {
    out$cov <- imp$cov[[1]]
    dimnames(out$cov) <- list(ov.names, ov.names)
    class(out$cov) <- c("lavaan.matrix.symmetric","matrix")
    if (meanstructure) {
      out$mean <- as.numeric(imp$mean[[1]])
      names(out$mean) <- ov.names
      class(out$mean) <- c("lavaan.vector","numeric")
    } else {
      out$mean <- sampstat.lavaan.mi(object@SampleStatsList[useImps],
                                     means = TRUE, categ = categ, m = m)$mean
    }
    if (categ) {
      out$th <- imp$th[[1]]
      names(out$th) <- lavaan::lavNames(object, "th")
      class(out$th) <- c("lavaan.vector","numeric")
    }
  }
  out
}
setMethod("fitted", "lavaan.mi", fitted.lavaan.mi)
setMethod("fitted.values", "lavaan.mi", fitted.lavaan.mi)

## function to calculate residuals for one group
gp.resid.lavaan.mi <- function(Observed, N, Implied, type,
                               means = FALSE, categ = FALSE, m) {
  obsMats <- sampstat.lavaan.mi(Observed, means = means, categ = categ, m = m)
  ## average sample stats across imputations
  S_mean <- if (is.null(N)) obsMats$cov else (obsMats$cov * ((N - 1L) / N))
  if (means) M_mean <- obsMats$mean
  if (categ) Th_mean <- obsMats$th

  if (type == "raw") {
    out <- list(cov = S_mean - Implied$cov)
    if (means) out$mean <- M_mean - Implied$mean else {
      out$mean <- rep(0, nrow(out$cov))
      names(out$mean) <- rownames(out$cov)
    }
    if (categ) out$th <- Th_mean - Implied$th
    return(out)
  } else if (type == "cor.bollen") {
    out <- list(cor = cov2cor(S_mean) - cov2cor(Implied$cov))
    if (!means) {
      out$mean <- rep(0, nrow(out$cor))
      names(out$mean) <- rownames(out$cor)
    } else {
      std.obs.M <- M_mean / sqrt(diag(S_mean))
      std.mod.M <- Implied$mean / sqrt(diag(Implied$cov))
      out$mean <- std.obs.M - std.mod.M
    }
  } else if (type == "cor.bentler") {
    SDs <- diag(sqrt(diag(S_mean)))
    dimnames(SDs) <- dimnames(S_mean)
    out <- list(cor = solve(SDs) %*% (S_mean - Implied$cov) %*% solve(SDs))
    class(out$cor) <- c("lavaan.matrix.symmetric","matrix")
    if (!means) {
      out$mean <- rep(0, nrow(out$cor))
      names(out$mean) <- rownames(out$cor)
    } else out$mean <- (M_mean - Implied$mean) / diag(SDs)
  } else stop("argument 'type' must be 'raw', 'cor', 'cor.bollen', ",
              "or 'cor.bentler'.")
  if (categ) out$th <- Th_mean - Implied$th
  out
}
resid.lavaan.mi <- function(object, type = c("raw","cor")) {
  ## @SampleStatsList is (for each imputation) output from:
  ##    getSampStats <- function(obj) lavaan::lavInspect(obj, "sampstat")
  useImps <- sapply(object@convergence, "[[", i = "converged")
  m <- sum(useImps)
  rescale <- lavaan::lavListInspect(object, "options")$sample.cov.rescale
  meanstructure <- lavaan::lavListInspect(object, "meanstructure")
  categ <- lavaan::lavListInspect(object, "categorical")
  type <- tolower(type[1])
  ## check for type = "cor" ("cor.bollen") or "cor.bentler"
  if (type == "cor") type <- "cor.bollen"
  ## model-implied moments, already pooled
  Implied <- getMethod("fitted", "lavaan.mi")(object)
  ## Calculate residuals
  nG <- lavaan::lavListInspect(object, "ngroups")
  N <- lavaan::lavListInspect(object, "nobs")
  if (nG > 1L) {
    group.label <- names(Implied)
    if (is.null(group.label)) group.label <- 1:length(Implied) else names(N) <- group.label
    out <- list()
    for (g in group.label) {
      out[[g]] <- gp.resid.lavaan.mi(Observed = lapply(object@SampleStatsList[useImps], "[[", g),
                                     N = if (rescale) N[g] else NULL,
                                     Implied = Implied[[g]], type = type,
                                     means = meanstructure, m = m, categ = categ)
    }
  } else {
    out <- gp.resid.lavaan.mi(Observed = object@SampleStatsList[useImps],
                              N = if (rescale) N else NULL,
                              Implied = Implied, type = type,
                              means = meanstructure, m = m, categ = categ)
  }
  out
}
setMethod("residuals", "lavaan.mi", resid.lavaan.mi)
setMethod("resid", "lavaan.mi", resid.lavaan.mi)


