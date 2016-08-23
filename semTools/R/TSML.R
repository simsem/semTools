## Terrence D. Jorgensen
### Last updated: 23 August 2016
### syntax to develop semTools function to implement 2-stage ML

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
  ## check for groups
  if (!is.null(funArgs$group))
    stop(c("Two-Stage estimation is currently only implemented ",
           "for single-group models.  Multiple-group models will be ",
           "implemented in a future version of semTools."))

  ## STAGE 1:
  ## fit saturated model
  targetNames <- lavaan::lavNames(do.call(lavaan::lavaanify, lavaanifyArgs))
  varnames <- c(targetNames, aux)
  covstruc <- outer(varnames, varnames, function(x, y) paste(x, "~~", y))
  satArgs <- funArgs
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
  coef(object@target, type = type)
})

setMethod("fitted.values", "twostage",
          function(object, model = c("target","saturated","baseline"),
                   type = "moments", labels = TRUE) {
  model <- model[1]
  fitted.values(slot(object, model), type = type, labels = labels)
})
setMethod("fitted", "twostage",
          function(object, model = c("target","saturated","baseline"),
                   type = "moments", labels = TRUE) {
  model <- model[1]
  fitted.values(slot(object, model), type = type, labels = labels)
})

setMethod("residuals", "twostage", function(object, type = c("raw","cor","normalized","standardized")) {
  type <- type[1]
  residuals(object@target, type = type)
})
setMethod("resid", "twostage", function(object, type = c("raw","cor","normalized","standardized")) {
  type <- type[1]
  residuals(object@target, type = type)
})

setMethod("nobs", "twostage",
          function(object, type = c("ntotal","ngroups","n.per.group","norig",
                                    "patterns","coverage")) {
  type <- type[1]
  if (type == "n.per.group") type <- "nobs"
  lavaan::lavInspect(object@saturated, what = type)
})

setMethod("vcov", "twostage", function(object, baseline = FALSE) {
  ## calculate model derivatives and complete-data information matrix
  MATS <- twostageMatrices(object, baseline)
  ## return asymptotic covariance matrix of 2-stage estimator
  bread <- solve(t(MATS$delta) %*% MATS$H %*% MATS$delta)
  meat <- MATS$H %*% MATS$delta
  out <- bread %*% t(meat) %*% MATS$satACOV %*% meat %*% bread
  class(out) <- c("lavaan.matrix.symmetric","matrix")
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
  chisq.naïve <- H0$scaled[["chisq.naïve"]] - H1$scaled[["chisq.naïve"]]
  cc <- (DF0*H0$scaled[["scaling.factor"]] - DF1*H1$scaled[["scaling.factor"]]) / DF
  if (cc < 0) {
    warning("Scaling factor is negative, so it was set to missing.")
    cc <- NA
  }
  scaled <- c(chisq.naïve = chisq.naïve, scaling.factor = cc,
              chisq.scaled = chisq.naïve / cc, DF = DF,
              pvalue = pchisq(chisq.naïve / cc, df = DF, lower.tail = FALSE))
  class(scaled) <- c("lavaan.vector","numeric")
  ## return both statistics
  cat("Difference test for Browne (1984) residual-based statistics:\n\n")
  print(residual)
  cat("\n\nSatorra-Bentler (2001) scaled difference test:\n\n")
  print(scaled)
  invisible(list(residual = residual, scaled = scaled))
})

setMethod("show", "twostage", function(object) {
  ## show chi-squared test results
  cat("Chi-squared test(s) results, ADJUSTED for missing data:\n\n")
  anova(object)
  cat("\n\nChi-squared test results, UNADJUSTED for missing data:\n\n")
  show(object@target)
  invisible(object)
})

setMethod("summary", "twostage", function(object, ...) {
  ## show chi-squared test results AND estimates
  show(object)
  cat("\n\nParameter Estimates, with SEs (and tests/CIs) ADJUSTED for missing data:\n\n")
  dots <- list(...)
  if (!"fmi" %in% names(dots)) dots$fmi <- FALSE
  if (!"ci" %in% names(dots)) dots$ci <- TRUE
  if (!"level" %in% names(dots)) dots$level <- .95
  PT <- lavaan::parTable(object@target)
  PT <- PT[PT$group > 0, ]
  PE <- do.call(lavaan::parameterEstimates, c(dots, object = object@target))
  SEs <- sqrt(diag(vcov(object)))
  PE$se[PT$free > 0] <- SEs[PT$free]
  PE$z[PT$free > 0] <- PE$est[PT$free > 0] / PE$se[PT$free > 0]
  PE$pvalue[PT$free > 0] <- pnorm(abs(PE$z[PT$free > 0]), lower.tail = FALSE)*2
  if (dots$ci) {
    crit <- qnorm(1 - (1 - dots$level) / 2)
    PE$ci.lower[PT$free > 0] <- PE$est[PT$free > 0] - crit * PE$se[PT$free > 0]
    PE$ci.upper[PT$free > 0] <- PE$est[PT$free > 0] + crit * PE$se[PT$free > 0]
  }
  if (dots$fmi) {
    compVar <- diag(vcov(object@target))[PT$free] ## FIXME: need to re-fit model to model-implied moments from Stage 2?
    # compFit <- update(object@target, sample.nobs = nobs(object@target),
    #                   sample.cov = lavInspect(object@target, "cov.ov"),
    #                   sample.mean = lavInspect(object@target, "mean.ov"))
    # compVar <- diag(vcov(compFit))[PT$free]
    missVar <- diag(vcov(object))[PT$free]
    PE$fmi[PT$free > 0] <- 1 - compVar / missVar
  }
  PE
})

## (hidden?) function utilized by vcov and anova methods
twostageMatrices <- function(object, baseline) {
  SLOT <- if (baseline) "baseline" else "target"
  ## model derivatives
  delta <- lavaan::lavInspect(slot(object, SLOT), "delta")
  covparams <- grep(pattern = "~~", x = rownames(delta))
  meanparams <- grep(pattern = "~1", x = rownames(delta))
  delta <- delta[c(covparams, meanparams), ]
  ## extract estimated moments from saturated model, and number of moments
  satSigma <- lavaan::lavInspect(object@saturated, "cov.ov")
  satMu <- lavaan::lavInspect(object@saturated, "mean.ov")
  if (length(object@auxNames)) {
    PT <- parTable(object@saturated)
    an <- object@auxNames
    tn <- lavaan::lavNames(slot(object, SLOT))
    satSigma <- satSigma[tn, tn]
    satMu <- satMu[tn]
  }
  p <- length(satMu) # if ngroups > 1, lengths?
  pStar <- p*(p + 1) / 2
  ## extract model-implied moments
  muHat <- lavaan::lavInspect(slot(object, SLOT), "mean.ov")[names(satMu)]
  sigmaHat <- lavaan::lavInspect(slot(object, SLOT), "cov.ov")
  sigmaHat <- sigmaHat[rownames(satSigma), colnames(satSigma)]
  shinv <- solve(sigmaHat)
  ## assemble complete-data information matrix
  H <- matrix(0, (pStar + p), (pStar + p))
  if (lavaan::lavInspect(slot(object, SLOT), "options")$estimator == "observed") {
    H[1:pStar, 1:pStar] <- .5*lavaan::lav_matrix_duplication_pre_post(shinv %x% shinv)
    H[(pStar + 1):(pStar + p), (pStar + 1):(pStar + p)] <- shinv
  } else {
    dMu <- satMu - muHat
    H[1:pStar, 1:pStar] <- lav_matrix_duplication_pre_post(shinv %x% (shinv %*% (satSigma + dMu %*% t(dMu)) %*% shinv - .5*shinv))
    H[(pStar + 1):(pStar + p), 1:pStar] <- lav_matrix_duplication_post(shinv %x% (t(dMu) %*% shinv))
    H[1:pStar, (pStar + 1):(pStar + p)] <- t(H[(pStar + 1):(pStar + p), 1:pStar])
    H[(pStar + 1):(pStar + p), (pStar + 1):(pStar + p)] <- shinv
  }
  ## asymptotic information and covariance matrices of target model
  satACOV <- vcov(object@saturated)
  satInfo <- solve(satACOV * nobs(object@saturated))
  ## all(round(acov*N, 8) == round(solve(info), 8))
  ## all(round(acov, 8) == round(solve(info)/N, 8))
  if (length(object@auxNames)) {
    dimTar <- PT$free[!(PT$lhs %in% an | PT$rhs %in% an)]
    dimAux <- PT$free[PT$lhs %in% an | PT$rhs %in% an]
    infoTar <- satInfo[dimTar, dimTar]
    infoAux <- satInfo[dimAux, dimAux]
    infoAT <- satInfo[dimAux, dimTar]
    satInfo <- infoTar - t(infoAT) %*% solve(infoAux) %*% infoAT
    satACOV <- solve(satInfo) / lavaan::lavInspect(object@saturated, "ntotal")
  }
  list(delta = delta, H = H, satACOV = satACOV, satInfo = satInfo)
}

## (hidden?) function utilized by anova method to test 1 or 2 models
twostageLRT <- function(object, baseline, print = FALSE) {
  SLOT <- if (baseline) "baseline" else "target"
  ## calculate model derivatives and complete-data information matrix
  MATS <- twostageMatrices(object, baseline)

  ## residual-based statistic (Savalei & Bentler, 2009, eq. 8)
  N <- lavaan::lavInspect(slot(object, SLOT), "ntotal")
  res <- residuals(slot(object, SLOT))
  etilde <- c(lavaan::lav_matrix_vech(res$cov), res$mean)
  ID <- MATS$satInfo %*% MATS$delta
  T.res <- N*t(etilde) %*% (MATS$satInfo - ID %*% solve(t(MATS$delta) %*% ID) %*% t(ID)) %*% etilde
  DF <- lavaan::lavInspect(slot(object, SLOT), "fit")[["df"]]
  pval.res <- pchisq(T.res, df = DF, lower.tail = FALSE)
  residual <- c(chisq = T.res, df = DF, pvalue = pval.res)
  class(residual) <- c("lavaan.vector","numeric")

  ## scaled test statistic (Savalei & Bentler, 2009, eq. 9)
  bread <- solve(t(MATS$delta) %*% MATS$H %*% MATS$delta)
  meat <- MATS$H %*% MATS$delta
  cc <- DF / sum(diag(solve(MATS$satInfo) %*% (MATS$H - meat %*% bread %*% t(meat))))
  chisq <- lavaan::lavInspect(slot(object, SLOT), "fit")[["chisq"]]
  T.scaled <- cc * chisq
  pval.scaled <- pchisq(T.scaled, df = DF, lower.tail = FALSE)
  scaled <- c(chisq.naïve = chisq, scaling.factor = 1 / cc,
              chisq.scaled = T.scaled, df = DF, pvalue = pval.scaled)
  class(scaled) <- c("lavaan.vector","numeric")

  ## return both statistics
  if (print) {
    cat("Browne (1984) residual-based test statistic:\n\n")
    print(residual)
    cat("\n\nSatorra-Bentler (2001) scaled test statistic:\n\n")
    print(scaled)
  }
  invisible(list(residual = residual, scaled = scaled))
}


# fitS <- cfa(model = model, data = dat1, missing = "fiml", se = "standard")
# fitR <- cfa(model = model, data = dat1, missing = "fiml", se = "robust.huber.white")
# all(lavInspect(fitS, "information") == lavInspect(fitR, "information"))
# all(vcov(fitS) == vcov(fitR))
