## Terrence D. Jorgensen
### Last updated: 14 October 2016
### semTools function to implement 2-stage ML

setClass("twostage",
         slots = c(saturated = "lavaan", target = "lavaan", baseline = "lavaan",
                   auxNames = "character"))

cfa.2stage <- function(..., aux = NULL, baseline.model = NULL) {
  twostage(..., aux = aux, fun = "cfa", baseline.model = baseline.model)
}

sem.2stage <- function(..., aux = NULL, baseline.model = NULL) {
  twostage(..., aux = aux, fun = "sem", baseline.model = baseline.model)
}

growth.2stage <- function(..., aux = NULL, baseline.model = NULL) {
  twostage(..., aux = aux, fun = "growth", baseline.model = baseline.model)
}

lavaan.2stage <- function(..., aux = NULL, baseline.model = NULL) {
  twostage(..., aux = aux, fun = "lavaan", baseline.model = baseline.model)
}

twostage <- function(..., aux, fun, baseline.model = NULL) {
  if (all(aux == "")) aux <- NULL
  dots <- list(...)
  if (is.null(dots$model)) stop("lavaan model syntax argument must be named 'model'.")
  lavaanifyArgs <- dots[intersect(names(dots), names(formals(lavaan::lavaanify)))]
  funArgs <- dots[intersect(names(dots), names(formals(lavaan::lavaan)))]
  ## set some non-optional lavaan arguments
  funArgs$meanstructure <- TRUE
  funArgs$conditional.x <- FALSE
  funArgs$fixed.x <- FALSE
  funArgs$missing <- "fiml"
  funArgs$estimator <- "ML"
  funArgs$test <- "standard"
  if (is.null(funArgs$information)) funArgs$information <- "observed"

  if (funArgs$information == "expected") {
    message("If data are MAR, only the observed information matrix is consistent.")
    if (!is.null(aux)) {
      funArgs$information <- "observed"
      message(c("Using auxiliary variables implies assuming that data are MAR. ",
                "The lavaan argument 'information' was set to 'observed'."))
    }
    if (!is.null(funArgs$se)) if(funArgs$se != "standard") {
      funArgs$information <- "observed"
      message(c("The lavaan argument 'information' was set to 'observed' ",
                "because adjusting SEs for non-normality requires it."))
    }
  }
  funArgs$NACOV <- NULL
  funArgs$do.fit <- NULL

  ## STAGE 1:
  ## fit saturated model
  if (!is.null(funArgs$group))
    lavaanifyArgs$ngroups <- length(table(funArgs$data[ , funArgs$group]))
  targetNames <- lavaan::lavNames(do.call(lavaan::lavaanify, lavaanifyArgs))
  varnames <- c(targetNames, aux)
  covstruc <- outer(varnames, varnames, function(x, y) paste(x, "~~", y))
  satArgs <- funArgs
  satArgs$constraints <- NULL
  satArgs$group.equal <- ""
  satArgs$model <- c(covstruc[lower.tri(covstruc, diag = TRUE)],
                     paste(varnames, "~ 1"))
  satFit <- do.call(lavaan::lavaan, satArgs)

  ## check for robust estimators
  opts <- lavaan::lavInspect(satFit, "options")
  if (!opts$se %in% c("standard","robust.huber.white"))
    stop(c("Two-Stage estimation requires either se = 'standard' for ",
           "multivariate normal data or se = 'robust.huber.white' to ",
           "correct for non-normality."))

  ## STAGE 2:
  ## fit target model to saturated estimates
  targetArgs <- funArgs
  targetArgs$data <- NULL
  targetArgs$sample.cov <- lavaan::lavInspect(satFit, "cov.ov")
  targetArgs$sample.mean <- lavaan::lavInspect(satFit, "mean.ov")
  targetArgs$sample.nobs <- lavaan::lavInspect(satFit, "nobs")
  targetArgs$se <- "standard"
  targetArgs$sample.cov.rescale <- FALSE
  targetFit <- do.call(fun, targetArgs)

  ## STAGE 0:
  ## fit baseline model (for incremental fit indices)
  baseArgs <- targetArgs
  if (is.null(baseline.model)) {
    basecov <- outer(targetNames, targetNames, function(x, y) paste0(x, " ~~ 0*", y))
    diag(basecov) <- paste(targetNames, "~~", targetNames)
    baseArgs$model <- c(basecov[lower.tri(basecov, diag = TRUE)],
                        paste(targetNames, "~ 1"))
  } else baseArgs$model <- baseline.model
  baseArgs$se <- "standard"
  baseFit <- do.call(lavaan::lavaan, baseArgs)
  if (length(setdiff(lavaan::lavNames(baseFit), targetNames)))
    warning("The baseline model includes variables excluded from the target model.")
  if (length(setdiff(targetNames, lavaan::lavNames(baseFit))))
    warning("The target model includes variables excluded from the baseline model.")

  ## return both models
  out <- new("twostage", saturated = satFit, target = targetFit,
             baseline = baseFit, auxNames = as.character(aux))
  out
}

## methods
setMethod("coef", "twostage", function(object, type = c("free","user")) {
  type <- type[1]
  lavaan::coef(object@target, type = type)
})

setMethod("fitted.values", "twostage",
          function(object, model = c("target","saturated","baseline"),
                   type = "moments", labels = TRUE) {
  model <- model[1]
  lavaan::fitted.values(slot(object, model), type = type, labels = labels)
})
setMethod("fitted", "twostage",
          function(object, model = c("target","saturated","baseline"),
                   type = "moments", labels = TRUE) {
  model <- model[1]
  lavaan::fitted.values(slot(object, model), type = type, labels = labels)
})

setMethod("residuals", "twostage", function(object, type = c("raw","cor","normalized","standardized")) {
  type <- type[1]
  lavaan::residuals(object@target, type = type)
})
setMethod("resid", "twostage", function(object, type = c("raw","cor","normalized","standardized")) {
  type <- type[1]
  lavaan::residuals(object@target, type = type)
})

setMethod("nobs", "twostage",
          function(object, type = c("ntotal","ngroups","n.per.group","norig",
                                    "patterns","coverage")) {
  type <- type[1]
  if (type == "n.per.group") type <- "nobs"
  lavaan::lavInspect(object@saturated, what = type)
})

setMethod("vcov", "twostage", function(object, baseline = FALSE) {
  SLOT <- if (baseline) "baseline" else "target"
  ## calculate model derivatives and complete-data information matrix
  MATS <- twostageMatrices(object, baseline)
  meat <- MATS$H %*% MATS$delta
  bread <- MASS::ginv(t(MATS$delta) %*% meat) # FIXME: why not solve()?
  out <- bread %*% t(meat) %*% MATS$satACOV %*% meat %*% bread
  class(out) <- c("lavaan.matrix.symmetric","matrix")
  if (baseline) {
    rownames(out) <- names(getMethod("coef", "lavaan")(object@baseline))
  } else {
    rownames(out) <- names(getMethod("coef", "twostage")(object))
  }
  colnames(out) <- rownames(out)
  out
})

## chi-squared test results (difference tests not available yet)
setMethod("anova", "twostage", function(object, h1 = NULL, baseline = FALSE) {
  if (is.null(h1)) {
    return(twostageLRT(object, baseline, print = TRUE))
  }
  H0 <- twostageLRT(object, baseline = FALSE)
  H1 <- twostageLRT(h1, baseline = FALSE)
  DF0 <- H0$residual[["df"]]
  DF1 <- H1$residual[["df"]]
  if (DF0 == DF1) stop("Models have the same degrees of freedom.")
  if (min(c(DF0, DF1)) == 0L) return(twostageLRT(object, baseline, print = TRUE))
  parent <- which.min(c(DF0, DF1))
  if (parent == 1L) {
    parent <- H0
    H0 <- H1
    H1 <- parent
    DF0 <- H0$residual[["df"]]
    DF1 <- H1$residual[["df"]]
  }
  DF <- DF0 - DF1
  ## residual-based statistic
  T.res <- H0$residual[["chisq"]] - H1$residual[["chisq"]]
  residual <- c(chisq = T.res, df = DF,
                pvalue = pchisq(T.res, df = DF, lower.tail = FALSE))
  class(residual) <- c("lavaan.vector","numeric")
  ## scaled test statistic
  chisq.naive <- H0$scaled[["chisq.naive"]] - H1$scaled[["chisq.naive"]]
  cc <- (DF0*H0$scaled[["scaling.factor"]] - DF1*H1$scaled[["scaling.factor"]]) / DF
  if (cc < 0) {
    warning("Scaling factor is negative, so it was set to missing.")
    cc <- NA
  }
  scaled <- c(chisq.naive = chisq.naive, scaling.factor = cc,
              chisq.scaled = chisq.naive / cc, DF = DF,
              pvalue = pchisq(chisq.naive / cc, df = DF, lower.tail = FALSE))
  class(scaled) <- c("lavaan.vector","numeric")
  ## return both statistics
  if (lavaan::lavInspect(object@saturated, "options")$se == "standard") {
    cat("Difference test for Browne (1984) residual-based statistics:\n\n")
    print(residual)
  }
  cat("\n\nSatorra-Bentler (2001) scaled difference test:\n\n")
  print(scaled)
  invisible(list(residual = residual, scaled = scaled))
})

setMethod("show", "twostage", function(object) {
  ## show chi-squared test results
  cat("Chi-squared test(s) results, ADJUSTED for missing data:\n\n")
  getMethod("anova", "twostage")(object)
  cat("\n\nChi-squared test results, UNADJUSTED for missing data:\n\n")
  show(object@target)
  invisible(object)
})

setMethod("summary", "twostage", function(object, ...) {
  ## show chi-squared test results AND estimates
  getMethod("show", "twostage")(object)
  cat("\n\nParameter Estimates, with SEs (and tests/CIs) ADJUSTED for missing data:\n\n")
  dots <- list(...)
  if (!"fmi" %in% names(dots)) dots$fmi <- FALSE
  if (!"ci" %in% names(dots)) dots$ci <- TRUE
  if (!"level" %in% names(dots)) dots$level <- .95
  PT <- lavaan::parTable(object@target)
  PT <- PT[PT$group > 0, ]
  PE <- do.call(lavaan::parameterEstimates, c(dots, object = object@target))
  SEs <- sqrt(diag(getMethod("vcov", "twostage")(object)))
  PE$se[PT$free > 0] <- SEs[PT$free]
  PE$z[PT$free > 0] <- PE$est[PT$free > 0] / PE$se[PT$free > 0]
  PE$pvalue[PT$free > 0] <- pnorm(abs(PE$z[PT$free > 0]), lower.tail = FALSE)*2
  if (dots$ci) {
    crit <- qnorm(1 - (1 - dots$level) / 2)
    PE$ci.lower[PT$free > 0] <- PE$est[PT$free > 0] - crit * PE$se[PT$free > 0]
    PE$ci.upper[PT$free > 0] <- PE$est[PT$free > 0] + crit * PE$se[PT$free > 0]
  }
  if (dots$fmi) {
    compVar <- diag(lavaan::vcov(object@target))[PT$free] ## FIXME: need to re-fit model to model-implied moments from Stage 2?
    # compFit <- lavaan::update(object@target, sample.nobs = lavaan::nobs(object@target),
    #                           sample.cov = lavInspect(object@target, "cov.ov"),
    #                           sample.mean = lavInspect(object@target, "mean.ov"))
    # compVar <- diag(lavaan::vcov(compFit))[PT$free]
    missVar <- SEs^2
    PE$fmi[PT$free > 0] <- 1 - compVar / missVar
  }
  PE
})

## (hidden?) function utilized by vcov and anova methods
twostageMatrices <- function(object, baseline) {
  SLOT <- if (baseline) "baseline" else "target"
  ## extract parameter table to isolate estimates by group
  PTsat <- lavaan::parTable(object@saturated)
  nG <- max(PTsat$group)
  isMG <- nG > 1L
  ## model derivatives
  delta <- lavaan::lavInspect(slot(object, SLOT), "delta")
  if (!isMG) delta <- list(delta)
  for (g in 1:nG) {
    covparams <- grep(pattern = "~~", x = rownames(delta[[g]]))
    meanparams <- grep(pattern = "~1", x = rownames(delta[[g]]))
    delta[[g]] <- delta[[g]][c(covparams, meanparams), ]
  }
  ## stack groups' deltas into 1 matrix
  delta <- do.call(rbind, delta)

  ## extract estimated moments from saturated model, and number of moments
  satSigma <- lavaan::lavInspect(object@saturated, "cov.ov")
  satMu <- lavaan::lavInspect(object@saturated, "mean.ov")
  if (!isMG) {
    satSigma <- list(satSigma)
    satMu <- list(satMu)
  }
  if (length(object@auxNames)) {
    an <- object@auxNames
    tn <- lavaan::lavNames(slot(object, SLOT))
    for (g in 1:nG) {
      satSigma[[g]] <- satSigma[[g]][tn, tn]
      satMu[[g]] <- satMu[[g]][tn]
    }
  }
  p <- length(satMu[[1]])
  pStar <- p*(p + 1) / 2
  ## extract model-implied moments
  muHat <- lavaan::lavInspect(slot(object, SLOT), "mean.ov")
  sigmaHat <- lavaan::lavInspect(slot(object, SLOT), "cov.ov")
  if (!isMG) {
    sigmaHat <- list(sigmaHat)
    muHat <- list(muHat)
  }
  shinv <- list()
  for (g in 1:nG) {
    muHat[[g]] <- muHat[[g]][names(satMu[[g]])]
    sigmaHat[[g]] <- sigmaHat[[g]][rownames(satSigma[[g]]), colnames(satSigma[[g]])]
    shinv[[g]] <- solve(sigmaHat[[g]])
  }
  ## assemble complete-data information matrix
  H <- list()
  for (g in 1:nG) H[[g]] <- matrix(0, (pStar + p), (pStar + p))

  if (lavaan::lavInspect(slot(object, SLOT), "options")$estimator == "expected") {
    for (g in 1:nG) {
      H[[g]][1:pStar, 1:pStar] <- .5*lavaan::lav_matrix_duplication_pre_post(shinv[[g]] %x% shinv[[g]])
      H[[g]][(pStar + 1):(pStar + p), (pStar + 1):(pStar + p)] <- shinv[[g]]
    }
  } else {
  	## estimator == "observed"
    dMu <- list()
    for (g in 1:nG) {
      dMu[[g]] <- satMu[[g]] - muHat[[g]]
      H[[g]][1:pStar, 1:pStar] <- lavaan::lav_matrix_duplication_pre_post(shinv[[g]] %x% (shinv[[g]] %*% (satSigma[[g]] + dMu[[g]] %*% t(dMu[[g]])) %*% shinv[[g]] - .5*shinv[[g]]))
      H[[g]][(pStar + 1):(pStar + p), 1:pStar] <- lavaan::lav_matrix_duplication_post(shinv[[g]] %x% (t(dMu[[g]]) %*% shinv[[g]]))
      H[[g]][1:pStar, (pStar + 1):(pStar + p)] <- t(H[[g]][(pStar + 1):(pStar + p), 1:pStar])
      H[[g]][(pStar + 1):(pStar + p), (pStar + 1):(pStar + p)] <- shinv[[g]]
    }
  }
  ## combine into 1 block-diagonal matrix
  H <- do.call(lavaan::lav_matrix_bdiag, H)

  ## asymptotic information and covariance matrices of target model
  satACOV <- lavaan::vcov(object@saturated)
  satInfo <- solve(satACOV * lavaan::nobs(object@saturated))
  ## all(round(acov*N, 8) == round(solve(info), 8))
  ## all(round(acov, 8) == round(solve(info)/N, 8))
  if (length(object@auxNames)) {
  	dimTar <- !(PTsat$lhs %in% an | PTsat$rhs %in% an)
  	dimAux <- PTsat$lhs %in% an | PTsat$rhs %in% an
  	infoTar <- satInfo[dimTar, dimTar]
  	infoAux <- satInfo[dimAux, dimAux]
  	infoAT <- satInfo[dimAux, dimTar]
  	satInfo <- infoTar - t(infoAT) %*% solve(infoAux) %*% infoAT
  	satACOV <- solve(satInfo) / lavaan::nobs(object@saturated)
  }
  list(delta = delta, H = H, satACOV = satACOV, satInfo = satInfo)
}

## (hidden?) function utilized by anova method to test 1 or 2 models
twostageLRT <- function(object, baseline, print = FALSE) {
  SLOT <- if (baseline) "baseline" else "target"
  ## calculate model derivatives and complete-data information matrix
  MATS <- twostageMatrices(object, baseline)
  ## residual-based statistic (Savalei & Bentler, 2009, eq. 8)
  N <- lavaan::nobs(slot(object, SLOT))
  nG <- lavaan::lavInspect(slot(object, SLOT), "ngroups")
  res <- lavaan::residuals(slot(object, SLOT))
  if (nG == 1L) res <- list(res)
  etilde <- do.call(c, lapply(res, function(x) c(lavaan::lav_matrix_vech(x$cov), x$mean)))
  ID <- MATS$satInfo %*% MATS$delta
  T.res <- N*t(etilde) %*% (MATS$satInfo - ID %*% MASS::ginv(t(MATS$delta) %*% ID) %*% t(ID)) %*% etilde # FIXME: why not solve()?
  DF <- lavaan::lavInspect(slot(object, SLOT), "fit")[["df"]]
  pval.res <- pchisq(T.res, df = DF, lower.tail = FALSE)
  residual <- c(chisq = T.res, df = DF, pvalue = pval.res)
  class(residual) <- c("lavaan.vector","numeric")

  ## scaled test statistic (Savalei & Bentler, 2009, eq. 9)
  meat <- MATS$H %*% MATS$delta
  bread <- MASS::ginv(t(MATS$delta) %*% meat) # FIXME: why not solve()?
  cc <- DF / sum(diag(MATS$satACOV %*% (MATS$H - meat %*% bread %*% t(meat))))
  chisq <- lavaan::lavInspect(slot(object, SLOT), "fit")[["chisq"]]
  T.scaled <- cc * chisq
  pval.scaled <- pchisq(T.scaled, df = DF, lower.tail = FALSE)
  scaled <- c(chisq.naive = chisq, scaling.factor = 1 / cc,
              chisq.scaled = T.scaled, df = DF, pvalue = pval.scaled)
  class(scaled) <- c("lavaan.vector","numeric")

  ## return both statistics
  if (print) {
    if (lavaan::lavInspect(object@saturated, "options")$se == "standard") {
      cat("Browne (1984) residual-based test statistic:\n\n")
      print(residual)
    }
    cat("\n\nSatorra-Bentler (2001) scaled test statistic:\n\n")
    print(scaled)
  }
  invisible(list(residual = residual, scaled = scaled))
}


# fitS <- cfa(model = model, data = dat1, missing = "fiml", se = "standard")
# fitR <- cfa(model = model, data = dat1, missing = "fiml", se = "robust.huber.white")
# all(lavInspect(fitS, "information") == lavInspect(fitR, "information"))
# all(vcov(fitS) == vcov(fitR))
