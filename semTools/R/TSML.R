## Terrence D. Jorgensen
### Last updated: 10 January 2021
### semTools function to implement 2-stage ML


## -----------------
## Class and Methods
## -----------------


##' Class for the Results of 2-Stage Maximum Likelihood (TSML) Estimation for
##' Missing Data
##'
##' This class contains the results of 2-Stage Maximum Likelihood (TSML)
##' estimation for missing data.  The \code{summary}, \code{anova}, \code{vcov}
##' methods return corrected \emph{SE}s and test statistics.  Other methods are
##' simply wrappers around the corresponding \code{\linkS4class{lavaan}}
##' methods.
##'
##'
##' @name twostage-class
##' @aliases twostage-class show,twostage-method summary,twostage-method
##' anova,twostage-method vcov,twostage-method coef,twostage-method
##' fitted.values,twostage-method fitted,twostage-method
##' residuals,twostage-method resid,twostage-method nobs,twostage-method
##' @docType class
##'
##' @slot saturated A fitted \code{\linkS4class{lavaan}} object containing the
##'  saturated model results
##' @slot target A fitted \code{\linkS4class{lavaan}} object containing the
##'  target/hypothesized model results
##' @slot baseline A fitted \code{\linkS4class{lavaan}} object containing the
##'  baseline/null model results
##' @slot auxNames A character string (potentially of \code{length == 0}) of any
##'  auxiliary variable names, if used
##'
##' @param object An object of class \code{twostage}.
##' @param ... arguments passed to \code{\link[lavaan]{parameterEstimates}}.
##' @param h1 An object of class \code{twostage} in which \code{object} is
##'        nested, so that their difference in fit can be tested using
##'        \code{anova} (see \bold{Value} section for details).
##' @param baseline \code{logical} indicating whether to return results for the
##'        baseline model, rather than the default target (hypothesized) model.
##' @param type The meaning of this argument varies depending on which method it
##'        it used for. Find detailed descriptions in the \bold{Value} section
##'        under \code{coef}, \code{nobs}, and \code{residuals}.
##' @param model \code{character} naming the slot for which to return the
##'        model-implied sample moments (see \code{fitted.values} description.)
##' @param labels \code{logical} indicating whether the model-implied sample
##'        moments should have (row/column) labels.
##'
##' @return
##'  \item{show}{\code{signature(object = "twostage"):} The \code{show} function
##'   is used to display the results of the \code{anova} method, as well as the
##'   header of the (uncorrected) target model results.}
##'  \item{summary}{\code{signature(object = "twostage", ...):} The summary
##'   function prints the same information from the \code{show} method, but also
##'   provides (and returns) the output of
##'   \code{\link[lavaan]{parameterEstimates}(object@target, ...)} with corrected
##'   \emph{SE}s, test statistics, and confidence intervals.  Additional
##'   arguments can be passed to \code{\link[lavaan]{parameterEstimates}},
##'   including \code{fmi = TRUE} to provide an estimate of the fraction of
##'   missing information.}
##'  \item{anova}{\code{signature(object = "twostage", h1 = NULL, baseline = FALSE):}
##'   The \code{anova} function returns the residual-based \eqn{\chi^2} test
##'   statistic result, as well as the scaled \eqn{\chi^2} test statistic result,
##'   for the model in the \code{target} slot, or for the model in the
##'   \code{baseline} slot if \code{baseline = TRUE}.  The user can also provide
##'   a single additional \code{twostage} object to the \code{h1} argument, in
##'   which case \code{anova} returns residual-based and scaled
##'   (\eqn{\Delta})\eqn{\chi^2} test results, under the assumption that the
##'   models are nested.  The models will be automatically sorted according their
##'   degrees of freedom.}
##'  \item{nobs}{\code{signature(object = "twostage",
##'   type = c("ntotal", "ngroups", "n.per.group", "norig", "patterns", "coverage")):}
##'   The \code{nobs} function will return the total sample sized used in the
##'   analysis by default.  Also available are the number of groups or the sample
##'   size per group, the original sample size (if any rows were deleted because
##'   all variables were missing), the missing data patterns, and the matrix of
##'   coverage (diagonal is the proportion of sample observed on each variable,
##'   and off-diagonal is the proportion observed for both of each pair of
##'   variables).}
##'  \item{coef}{\code{signature(object = "twostage", type = c("free", "user")):}
##'   This is simply a wrapper around the corresponding
##'   \code{\linkS4class{lavaan}} method, providing point estimates from the
##'   \code{target} slot.}
##'  \item{vcov}{\code{signature(object = "twostage", baseline = FALSE):} Returns
##'   the asymptotic covariance matrix of the estimated parameters (corrected for
##'   additional uncertainty due to missing data) for the model in the
##'   \code{target} slot, or for the model in the \code{baseline} slot if
##'   \code{baseline = TRUE}.}
##'  \item{fitted.values, fitted}{\code{signature(object = "twostage",
##'   model = c("target", "saturated", "baseline")):} This is simply a wrapper
##'   around the corresponding \code{\linkS4class{lavaan}} method, providing
##'   model-implied sample moments from the slot specified in the \code{model}
##'   argument.}
##'  \item{residuals, resid}{\code{signature(object = "twostage", type = c("raw",
##'   "cor", "normalized", "standardized")):} This is simply a wrapper around the
##'   corresponding \code{\linkS4class{lavaan}} method, providing residuals of
##'   the specified \code{type} from the \code{target} slot.}
##'
##' @section Objects from the Class: Objects can be created via the
##' \code{\link{twostage}} function.
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##' \email{TJorgensen314@@gmail.com})
##'
##' @seealso \code{\link{twostage}}
##'
##' @examples
##'
##' # See the example from the twostage function
##'
setClass("twostage",
         slots = c(saturated = "lavaan", target = "lavaan", baseline = "lavaan",
                   auxNames = "character"))


##' @rdname twostage-class
##' @aliases show,twostage-method
##' @export
setMethod("show", "twostage", function(object) {
  ## show chi-squared test results
  cat("Chi-squared test(s) results, ADJUSTED for missing data:\n\n")
  getMethod("anova", "twostage")(object)
  cat("\n\nChi-squared test results, UNADJUSTED for missing data:\n\n")
  show(object@target)
  invisible(object)
})


##' @rdname twostage-class
##' @aliases summary,twostage-method
##' @importFrom stats pnorm qnorm
##' @importFrom lavaan parTable
##' @export
setMethod("summary", "twostage", function(object, ...) {
  ## show chi-squared test results AND estimates
  getMethod("show", "twostage")(object)
  cat("\n\nParameter Estimates, with SEs (and tests/CIs) ADJUSTED for missing data:\n\n")
  dots <- list(...)
  if (!"fmi" %in% names(dots)) dots$fmi <- FALSE
  if (!"ci" %in% names(dots)) dots$ci <- TRUE
  if (!"level" %in% names(dots)) dots$level <- .95
  PT <- parTable(object@target)
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


## (hidden) function utilized by vcov and anova methods
##' @importFrom lavaan lavInspect parTable
twostageMatrices <- function(object, baseline) {
  SLOT <- if (baseline) "baseline" else "target"
  ## extract parameter table to isolate estimates by group
  PTsat <- parTable(object@saturated)
  nG <- max(PTsat$group)
  isMG <- nG > 1L
  ## model derivatives
  delta <- lavInspect(slot(object, SLOT), "delta")
  if (!isMG) delta <- list(delta)
  for (g in 1:nG) {
    covparams <- grep(pattern = "~~", x = rownames(delta[[g]]))
    meanparams <- grep(pattern = "~1", x = rownames(delta[[g]]))
    delta[[g]] <- delta[[g]][c(covparams, meanparams), ]
  }
  ## stack groups' deltas into 1 matrix
  delta <- do.call(rbind, delta)

  ## extract estimated moments from saturated model, and number of moments
  satSigma <- lavInspect(object@saturated, "cov.ov")
  satMu <- lavInspect(object@saturated, "mean.ov")
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
  muHat <- lavInspect(slot(object, SLOT), "mean.ov")
  sigmaHat <- lavInspect(slot(object, SLOT), "cov.ov")
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

  if (lavInspect(slot(object, SLOT), "options")$estimator == "expected") {
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
##' @importFrom stats pchisq
##' @importFrom lavaan lavInspect
twostageLRT <- function(object, baseline, print = FALSE) {
  SLOT <- if (baseline) "baseline" else "target"
  ## calculate model derivatives and complete-data information matrix
  MATS <- twostageMatrices(object, baseline)
  ## residual-based statistic (Savalei & Bentler, 2009, eq. 8)
  N <- lavaan::nobs(slot(object, SLOT))
  nG <- lavInspect(slot(object, SLOT), "ngroups")
  res <- lavaan::residuals(slot(object, SLOT))
  if (nG == 1L) res <- list(res)
  etilde <- do.call(c, lapply(res, function(x) c(lavaan::lav_matrix_vech(x$cov), x$mean)))
  ID <- MATS$satInfo %*% MATS$delta
  T.res <- N*t(etilde) %*% (MATS$satInfo - ID %*% MASS::ginv(t(MATS$delta) %*% ID) %*% t(ID)) %*% etilde # FIXME: why not solve()?
  DF <- lavInspect(slot(object, SLOT), "fit")[["df"]]
  pval.res <- pchisq(T.res, df = DF, lower.tail = FALSE)
  residual <- c(chisq = T.res, df = DF, pvalue = pval.res)
  class(residual) <- c("lavaan.vector","numeric")

  ## scaled test statistic (Savalei & Bentler, 2009, eq. 9)
  meat <- MATS$H %*% MATS$delta
  bread <- MASS::ginv(t(MATS$delta) %*% meat) # FIXME: why not solve()?
  cc <- DF / sum(diag(MATS$satACOV %*% (MATS$H - meat %*% bread %*% t(meat))))
  chisq <- lavInspect(slot(object, SLOT), "fit")[["chisq"]]
  T.scaled <- cc * chisq
  pval.scaled <- pchisq(T.scaled, df = DF, lower.tail = FALSE)
  scaled <- c(chisq.naive = chisq, scaling.factor = 1 / cc,
              chisq.scaled = T.scaled, df = DF, pvalue = pval.scaled)
  class(scaled) <- c("lavaan.vector","numeric")

  ## return both statistics
  if (print) {
    if (lavInspect(object@saturated, "options")$se == "standard") {
      cat("Browne (1984) residual-based test statistic:\n\n")
      print(residual)
    }
    cat("\n\nSatorra-Bentler (2001) scaled test statistic:\n\n")
    print(scaled)
  }
  invisible(list(residual = residual, scaled = scaled))
}

##' @rdname twostage-class
##' @aliases anova,twostage-method
##' @importFrom lavaan lavInspect
##' @export
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
  if (lavInspect(object@saturated, "options")$se == "standard") {
    cat("Difference test for Browne (1984) residual-based statistics:\n\n")
    print(residual)
  }
  cat("\n\nSatorra-Bentler (2001) scaled difference test:\n\n")
  print(scaled)
  invisible(list(residual = residual, scaled = scaled))
})


##' @rdname twostage-class
##' @aliases nobs,twostage-method
##' @importFrom lavaan lavInspect
##' @export
setMethod("nobs", "twostage",
function(object, type = c("ntotal","ngroups","n.per.group","norig",
                          "patterns","coverage")) {
  type <- type[1]
  if (type == "n.per.group") type <- "nobs"
  lavInspect(object@saturated, what = type)
})


##' @rdname twostage-class
##' @aliases coef,twostage-method
##' @export
setMethod("coef", "twostage", function(object, type = c("free","user")) {
  type <- type[1]
  lavaan::coef(object@target, type = type)
})


##' @rdname twostage-class
##' @aliases vcov,twostage-method
##' @export
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


##' @rdname twostage-class
##' @aliases fitted.values,twostage-method
##' @export
setMethod("fitted.values", "twostage",
          function(object, model = c("target","saturated","baseline"),
                   type = "moments", labels = TRUE) {
  model <- model[1]
  lavaan::fitted.values(slot(object, model), type = type, labels = labels)
})

##' @rdname twostage-class
##' @aliases fitted,twostage-method
##' @export
setMethod("fitted", "twostage",
          function(object, model = c("target","saturated","baseline"),
                   type = "moments", labels = TRUE) {
  model <- model[1]
  lavaan::fitted.values(slot(object, model), type = type, labels = labels)
})


##' @rdname twostage-class
##' @aliases residuals,twostage-method
##' @export
setMethod("residuals", "twostage",
          function(object, type = c("raw","cor","normalized","standardized")) {
  type <- type[1]
  lavaan::residuals(object@target, type = type)
})

##' @rdname twostage-class
##' @aliases resid,twostage-method
##' @export
setMethod("resid", "twostage",
          function(object, type = c("raw","cor","normalized","standardized")) {
  type <- type[1]
  lavaan::residuals(object@target, type = type)
})



# fitS <- cfa(model = model, data = dat1, missing = "fiml", se = "standard")
# fitR <- cfa(model = model, data = dat1, missing = "fiml", se = "robust.huber.white")
# all(lavInspect(fitS, "information") == lavInspect(fitR, "information"))
# all(vcov(fitS) == vcov(fitR))



## ---------------------
## Constructor Functions
## ---------------------

##' Fit a lavaan model using 2-Stage Maximum Likelihood (TSML) estimation for
##' missing data.
##'
##' This function automates 2-Stage Maximum Likelihood (TSML) estimation,
##' optionally with auxiliary variables.  Step 1 involves fitting a saturated
##' model to the partially observed data set (to variables in the hypothesized
##' model as well as auxiliary variables related to missingness).  Step 2
##' involves fitting the hypothesized model to the model-implied means and
##' covariance matrix (also called the "EM" means and covariance matrix) as if
##' they were complete data.  Step 3 involves correcting the Step-2 standard
##' errors (\emph{SE}s) and chi-squared statistic to account for additional
##' uncertainty due to missing data (using information from Step 1; see
##' References section for sources with formulas).
##'
##' All variables (including auxiliary variables) are treated as endogenous
##' varaibles in the Step-1 saturated model (\code{fixed.x = FALSE}), so data
##' are assumed continuous, although not necessarily multivariate normal
##' (dummy-coded auxiliary variables may be included in Step 1, but categorical
##' endogenous variables in the Step-2 hypothesized model are not allowed).  To
##' avoid assuming multivariate normality, request \code{se =
##' "robust.huber.white"}.  CAUTION: In addition to setting \code{fixed.x =
##' FALSE} and \code{conditional.x = FALSE} in \code{\link[lavaan]{lavaan}},
##' this function will automatically set \code{meanstructure = TRUE},
##' \code{estimator = "ML"}, \code{missing = "fiml"}, and \code{test =
##' "standard"}.  \code{\link[lavaan]{lavaan}}'s \code{se} option can only be
##' set to \code{"standard"} to assume multivariate normality or to
##' \code{"robust.huber.white"} to relax that assumption.
##'
##'
##' @aliases twostage cfa.2stage sem.2stage growth.2stage lavaan.2stage
##' @importFrom lavaan lavInspect
##'
##' @param \dots Arguments passed to the \code{\link[lavaan]{lavaan}} function
##'   specified in the \code{fun} argument.  See also
##'   \code{\link[lavaan]{lavOptions}}.  At a minimum, the user must supply the
##'   first two named arguments to \code{\link[lavaan]{lavaan}} (i.e.,
##'   \code{model} and \code{data}).
##' @param aux An optional character vector naming auxiliary variable(s) in
##'   \code{data}
##' @param fun The character string naming the lavaan function used to fit the
##'   Step-2 hypothesized model (\code{"cfa"}, \code{"sem"}, \code{"growth"}, or
##'   \code{"lavaan"}).
##' @param baseline.model An optional character string, specifying the lavaan
##'   \code{\link[lavaan]{model.syntax}} for a user-specified baseline model.
##'   Interested users can use the fitted baseline model to calculate incremental
##'   fit indices (e.g., CFI and TLI) using the corrected chi-squared values (see
##'   the \code{anova} method in \code{\linkS4class{twostage}}).  If \code{NULL},
##'   the default "independence model" (i.e., freely estimated means and
##'   variances, but all covariances constrained to zero) will be specified
##'   internally.
##'
##' @return The \code{\linkS4class{twostage}} object contains 3 fitted lavaan
##' models (saturated, target/hypothesized, and baseline) as well as the names
##' of auxiliary variables.  None of the individual models provide the correct
##' model results (except the point estimates in the target model are unbiased).
##' Use the methods in \code{\linkS4class{twostage}} to extract corrected
##' \emph{SE}s and test statistics.
##'
##' @author
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @seealso \code{\linkS4class{twostage}}
##'
##' @references
##' Savalei, V., & Bentler, P. M. (2009). A two-stage approach to missing data:
##' Theory and application to auxiliary variables.
##' \emph{Structural Equation Modeling, 16}(3), 477--497.
##' \doi{10.1080/10705510903008238}
##'
##' Savalei, V., & Falk, C. F. (2014). Robust two-stage approach outperforms
##' robust full information maximum likelihood with incomplete nonnormal data.
##' \emph{Structural Equation Modeling, 21}(2), 280--302.
##' \doi{10.1080/10705511.2014.882692}
##'
##' @examples
##'
##' ## impose missing data for example
##' HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
##'                                       "ageyr","agemo","school")]
##' set.seed(12345)
##' HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
##' age <- HSMiss$ageyr + HSMiss$agemo/12
##' HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
##'
##' ## specify CFA model from lavaan's ?cfa help page
##' HS.model <- '
##'   visual  =~ x1 + x2 + x3
##'   textual =~ x4 + x5 + x6
##'   speed   =~ x7 + x8 + x9
##' '
##'
##' ## use ageyr and agemo as auxiliary variables
##' out <- cfa.2stage(model = HS.model, data = HSMiss, aux = c("ageyr","agemo"))
##'
##' ## two versions of a corrected chi-squared test results are shown
##' out
##' ## see Savalei & Bentler (2009) and Savalei & Falk (2014) for details
##'
##' ## the summary additionally provides the parameter estimates with corrected
##' ## standard errors, test statistics, and confidence intervals, along with
##' ## any other options that can be passed to parameterEstimates()
##' summary(out, standardized = TRUE)
##'
##'
##'
##' ## use parameter labels to fit a more constrained model
##' modc <- '
##'   visual  =~ x1 + x2 + x3
##'   textual =~ x4 + x5 + x6
##'   speed   =~ x7 + a*x8 + a*x9
##' '
##' outc <- cfa.2stage(model = modc, data = HSMiss, aux = c("ageyr","agemo"))
##'
##'
##' ## use the anova() method to test this constraint
##' anova(out, outc)
##' ## like for a single model, two corrected statistics are provided
##'
##' @export
twostage <- function(..., aux, fun, baseline.model = NULL) {
  if (all(aux == "")) aux <- NULL
  dots <- list(...)
  if (is.null(dots$model)) stop("lavaan model syntax argument must be named 'model'.")
  ####################### FIXME: also check intersect(names(dots), names(lavOptions()))
  lavaanifyArgs <- dots[intersect(names(dots), names(formals(lavaan::lavaanify)))]
  funArgs <- dots[intersect(names(dots), names(formals(lavaan::lavaan)))] #FIXME: lavOptions too
  ## set some non-optional lavaan arguments
  funArgs$meanstructure <- TRUE
  funArgs$conditional.x <- FALSE
  funArgs$fixed.x <- FALSE
  funArgs$missing <- "fiml"
  funArgs$estimator <- "ML"
  funArgs$test <- "standard"
  if (is.null(funArgs$information)) funArgs$information <- "observed"

  if (funArgs$information[1] == "expected") {
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
  opts <- lavInspect(satFit, "options")
  if (!opts$se %in% c("standard","robust.huber.white"))
    stop(c("Two-Stage estimation requires either se = 'standard' for ",
           "multivariate normal data or se = 'robust.huber.white' to ",
           "correct for non-normality."))

  ## STAGE 2:
  ## fit target model to saturated estimates
  targetArgs <- funArgs
  targetArgs$data <- NULL
  targetArgs$sample.cov <- lavInspect(satFit, "cov.ov")
  targetArgs$sample.mean <- lavInspect(satFit, "mean.ov")
  targetArgs$sample.nobs <- lavInspect(satFit, "nobs")
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


##' @rdname twostage
##' @export
lavaan.2stage <- function(..., aux = NULL, baseline.model = NULL) {
  twostage(..., aux = aux, fun = "lavaan", baseline.model = baseline.model)
}

##' @rdname twostage
##' @export
cfa.2stage <- function(..., aux = NULL, baseline.model = NULL) {
  twostage(..., aux = aux, fun = "cfa", baseline.model = baseline.model)
}

##' @rdname twostage
##' @export
sem.2stage <- function(..., aux = NULL, baseline.model = NULL) {
  twostage(..., aux = aux, fun = "sem", baseline.model = baseline.model)
}

##' @rdname twostage
##' @export
growth.2stage <- function(..., aux = NULL, baseline.model = NULL) {
  twostage(..., aux = aux, fun = "growth", baseline.model = baseline.model)
}


