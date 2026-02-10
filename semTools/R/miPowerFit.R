### Sunthud Pornprasertmanit; with contributions by Terrence D. Jorgensen
### Last updated: 3 February 2026


#' EPC Equivalence Fit Evaluation Using Modification Indices
#'
#' Evaluates model fit from an equivalence-testing perspective by
#' aggregating local EPC-based diagnostics into a global, fit-style
#' assessment. The procedure combines modification indices (MI),
#' expected parameter changes (EPC), statistical power, and confidence
#' intervals relative to a smallest effect size of interest (SESOI).
#'
#' Two complementary local decision rules are implemented:
#'
#' \strong{Method 1 (Power-based; Saris, Satorra, & van der Veld, 2009).}
#' Modification indices, statistical power, and EPC magnitude are jointly
#' evaluated (the J-rule) to classify fixed parameters as misspecified,
#' not misspecified, or inconclusive.
#'
#' \strong{Method 2 (CI-based equivalence testing).}
#' Confidence intervals of EPCs are compared against a trivial
#' misspecification region defined by the SESOI to determine whether
#' fixed parameters are substantially misspecified, trivially misspecified,
#' underpowered, or inconclusive.
#'
#' The resulting local classifications are returned in a single data
#' frame and can be summarized to yield a global equivalence-style fit
#' evaluation.
#'
#' @param lavaanObj A fitted \code{lavaan} object used to evaluate model fit.
#' @param stdLoad Standardized factor loading defining the SESOI for
#'   loading misspecifications. Default is 0.4.
#' @param cor Default standardized correlation defining the SESOI for
#'   covariance misspecifications. Used for both latent and residual
#'   covariances unless overridden.
#' @param corLatent Standardized latent factor correlation defining the
#'   SESOI for latent covariance misspecifications. If \code{NULL},
#'   defaults to \code{cor}.
#' @param corResidual Standardized residual correlation defining the
#'   SESOI for residual covariance misspecifications. If \code{NULL},
#'   defaults to \code{cor}.
#' @param stdBeta Standardized regression coefficient defining the SESOI
#'   for structural misspecifications. Default is 0.1.
#' @param stdIntcept Standardized intercept (Cohen's \emph{d}) defining
#'   the SESOI for intercept misspecifications. Default is 0.2.
#' @param stdSesoi Optional vector of standardized SESOI values. If
#'   provided, overrides operator-specific SESOI definitions.
#' @param sesoi Optional vector of unstandardized SESOI values. If
#'   provided, overrides \code{stdSesoi} and all operator-specific SESOI
#'   arguments.
#' @param cilevel Confidence level for EPC confidence intervals used in
#'   CI-based equivalence testing.
#' @param \dots Additional arguments passed to
#'   \code{\link[lavaan]{modificationIndices}}.
#'
#' @details
#' This function provides a local-to-global equivalence-based alternative
#' to traditional exact-fit evaluation. It is designed to assess whether
#' fixed parameters are substantively misspecified relative to a SESOI,
#' rather than whether a model fits exactly.
#'
#' Models with categorical indicators or unsupported constraints may
#' not be fully supported.
#'
#' @return A data frame with one row per fixed parameter, containing:
#' \enumerate{
#'   \item Parameter identifiers: \code{lhs}, \code{op}, \code{rhs}, and \code{group}.
#'   \item Modification index (\code{mi}) and expected parameter change estimates
#'   (\code{epc}).
#'   \item Unstandardized and standardized smallest effect size of interest values
#'   (\code{sesoi}, \code{std.sesoi}).
#'   \item Power-based decision (\code{decision.pow}) and related diagnostics, including
#'   whether the modification index is statistically significant
#'   (\code{significant.mi}) and whether the misfit at the SESOI has power greater
#'   than 0.80 (\code{high.power}). Decision labels are:
#'   M = Substantially Misspecified,
#'   I = Inconclusive,
#'   NM = Trivially Misspecified,
#'   EPC:M = Substantially Misspecified based on EPC information,
#'   EPC:NM = Trivially Misspecified based on EPC information.
#'   \item EPC-related statistics, including the standard error of the EPC
#'   (\code{se.epc}), confidence interval bounds for the EPC
#'   (\code{lower.epc}, \code{upper.epc}), and confidence interval bounds for the
#'   standardized EPC (\code{lower.std.epc}, \code{upper.std.epc}).
#'   \item Confidence-interval–based equivalence decision (\code{decision.ci}), with
#'   labels:
#'   M = Substantially Misspecified (EPC exceeds the SESOI),
#'   I = Inconclusive,
#'   NM = Trivially Misspecified,
#'   U = Underpowered (CI too wide to evaluate equivalence relative to the SESOI).
#' }
#'
#' @references
#' Saris, W. E., Satorra, A., & van der Veld, W. M. (2009).
#' Testing structural equation models or detection of misspecifications?
#' \emph{Structural Equation Modeling, 16}(4), 561--582.
#'
#' @seealso \code{\link{epcEquivCheck}}
#'
#' @importFrom stats qchisq pchisq qnorm
#' @importFrom lavaan modificationIndices lavNames fitMeasures
#'
#' @aliases epcEquivFit miPowerFit
#'
#' @examples
#'
#' library(lavaan)
#'
#' one.model <- ' onefactor  =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 '
#' fit <- cfa(one.model, data = HolzingerSwineford1939)
#' out <- epcEquivFit(fit)
#' out
#' summary(out)
#'
#' @export
epcEquivFit <- function(lavaanObj,
                       stdLoad = 0.4,
                       cor = 0.1,
                       corLatent = NULL,
                       corResidual = NULL,
                       stdBeta = 0.1,
                       stdIntcept = 0.2,
                       stdSesoi = NULL,
                       sesoi = NULL,
                       cilevel = 0.90, ...) {
  dots <- list(...)

  df_model <- lavaan::fitMeasures(lavaanObj, "df")

  if (is.na(df_model) || df_model <= 0) {
    stop(
      "epcEquivFit() requires a model with positive degrees of freedom (df > 0).\n",
      "The supplied lavaan model has df = ", df_model, ".\n",
      "EPC-based equivalence testing is not defined for just-identified or ",
      "under-identified models.",
      call. = FALSE
    )
  }

  # deprecated aliases
  if (!is.null(dots$intcept)) {
    warning(
      "'intcept' is deprecated; please use 'stdIntcept' instead.",
      call. = FALSE
    )
    stdIntcept <- dots$intcept
  }

  if (!is.null(dots$stdDelta)) {
    warning(
      "'stdDelta' is deprecated; please use 'stdSesoi' instead.",
      call. = FALSE
    )
    stdSesoi <- dots$stdDelta
  }

  if (!is.null(dots$delta)) {
    warning(
      "'delta' is deprecated; please use 'sesoi' instead.",
      call. = FALSE
    )
    sesoi <- dots$delta
  }

  lv_names <- lavaan::lavNames(lavaanObj, type="lv")
  mi <- lavaan::modificationIndices(lavaanObj)
  mi <- mi[mi$op != "==",]
  sigma <- mi[,"epc"] / sqrt(mi[,"mi"])
  sigma[!is.finite(sigma)] <- NA_real_
  if (is.null(corLatent))   corLatent   <- cor
  if (is.null(corResidual)) corResidual <- cor
  if (is.null(sesoi)) {
    if (is.null(stdSesoi))
      stdSesoi <- getTrivialEpc(mi, lv_names=lv_names, stdLoad = stdLoad,
                                corLatent = corLatent, corResidual = corResidual,
                                stdBeta = stdBeta, stdIntcept = stdIntcept)
    if (length(stdSesoi) == 1) stdSesoi <- rep(stdSesoi, nrow(mi))
    sesoi <- unstandardizeEpc(mi, stdSesoi, lavInspectTotalVar(lavaanObj), lavInspectResidualVar(lavaanObj))
  }
  if (length(sesoi) == 1) sesoi <- rep(sesoi, nrow(mi))
  ncp <- (sesoi / sigma)^2
  alpha <- 0.05
  desiredPow <- 0.80
  cutoff <- stats::qchisq(1 - alpha, df = 1)
  pow <- 1 - stats::pchisq(cutoff, df = 1, ncp = ncp)
  sigMI <- mi[,"mi"] > cutoff
  highPow <- pow > desiredPow
  group <- rep(1, nrow(mi))
  if ("group" %in% colnames(mi)) group <- mi[ , "group"]
  decision <- mapply(decisionMIPow, sigMI = sigMI, highPow = highPow,
                     epc = mi[ , "epc"], trivialEpc = sesoi)
  if (is.null(stdSesoi)) stdSesoi <- standardizeEpc(mi, lavInspectTotalVar(lavaanObj),
                                                    lavInspectResidualVar(lavaanObj),
                                                    sesoi = sesoi)
  result <- cbind(mi[ , 1:3], group, as.numeric(mi[ , "mi"]), mi[ , "epc"],
                  sesoi, mi[ , "sepc.all"],
                  stdSesoi, sigMI, highPow, decision)
  # New method
  crit <- abs(stats::qnorm((1 - cilevel)/2))
  seepc <- abs(result[,6]) / sqrt(abs(result[,5]))
  lowerepc <- result[,6] - crit * seepc
  upperepc <- result[,6] + crit * seepc
  stdlowerepc <- standardizeEpc(mi, lavInspectTotalVar(lavaanObj),
                                lavInspectResidualVar(lavaanObj), sesoi = lowerepc)
  stdupperepc <- standardizeEpc(mi, lavInspectTotalVar(lavaanObj),
                                lavInspectResidualVar(lavaanObj), sesoi = upperepc)
  isVar <- mi[,"op"] == "~~" & mi[,"lhs"] == mi[,"rhs"]
  decisionci <- mapply(decisionCIEpc, targetval = as.numeric(stdSesoi),
                       lower = stdlowerepc, upper = stdupperepc,
                       positiveonly = isVar)
  result <- cbind(result, pow, seepc, lowerepc, upperepc, stdlowerepc,
                  stdupperepc, decisionci)
  result <- result[!is.na(decision), ]
  colnames(result) <- c("lhs","op","rhs","group","mi","epc","sesoi",
                        "std.epc","std.sesoi","significant.mi",
                        "high.power","decision.pow","pow","se.epc","lower.epc",
                        "upper.epc","lower.std.epc","upper.std.epc","decision.ci")
  class(result) <- c("epcequivfit.data.frame","lavaan.data.frame","data.frame")

  # backward compatibility alias
  result$target.epc <- result$sesoi
  result$std.target.epc <- result$std.sesoi

  return(result)
}

#FIXME: Change to .Defunct after a few version updates
miPowerFit <- function(...) {
  .Deprecated("epcEquivFit",
              msg = "miPowerFit() has been replaced by epcEquivFit().")
  epcEquivFit(...)
}

#' EPC Equivalence Feasibility Check for Standardized Parameters
#'
#' Performs an EPC-based feasibility check to assess whether a set of
#' standardized population parameters defines a valid population
#' covariance matrix and whether trivially misspecified parameters
#' remain within a user-defined smallest effect size of interest (SESOI).
#' Feasibility is evaluated by constructing implied population models
#' under targeted parameter perturbations and examining EPC behavior
#' using \code{\link{epcEquivFit}}.
#'
#' This function focuses on standardized parameters and supports
#' recursive SEMs with continuous indicators only.
#'
#' @param lavaanObj A fitted \code{lavaan} object representing the target model.
#' @param minRelEffect A scalar in (0, 1) specifying the minimum relative
#'   magnitude of the standardized perturbation to be evaluated. The
#'   default value of 0.75 indicates that perturbations equal to 75\% of
#'   the SESOI are treated as trivial. If EPCs exceed the SESOI under
#'   such perturbations, EPC equivalence testing is not recommended.
#' @param stdLoad Standardized factor loading used to define the SESOI
#'   for loading misspecifications.
#' @param cor Standardized correlation used as a default SESOI for
#'   covariance misspecifications. This value is used for both latent
#'   and residual covariances unless overridden by
#'   \code{corLatent} or \code{corResidual}.
#' @param corLatent Standardized latent factor correlation used to
#'   define the SESOI for latent covariance misspecifications. If
#'   \code{NULL}, defaults to \code{cor}.
#' @param corResidual Standardized residual correlation used to define
#'   the SESOI for indicator residual covariance misspecifications. If
#'   \code{NULL}, defaults to \code{cor}.
#' @param stdBeta Standardized regression coefficient used to define
#'   the SESOI for structural misspecifications.
#'
#' @details
#' The procedure first checks whether the standardized parameters imply
#' a positive definite population covariance matrix. It then evaluates
#' EPC behavior under both positive and negative trivial
#' misspecifications by repeatedly constructing implied population
#' covariance matrices with perturbed parameters
#' (\code{minRelEffect} \eqn{\times} SESOI), refitting the model, and
#' re-evaluating EPCs.
#'
#' Models with categorical indicators, formative indicators, or
#' multiple-group structures are not supported.
#'
#' @return An object of class \code{"epcEquivCheckStd"} containing:
#' \itemize{
#'   \item \code{feasible}: Logical indicator of whether a valid
#'     standardized population model exists.
#'   \item \code{any_M}: Logical indicator of whether any EPC exceeded
#'     the SESOI under the evaluated misspecifications.
#'   \item \code{recommendation}: Character string summarizing feasibility
#'     (e.g., \code{"RECOMMENDED"}, \code{"NOT RECOMMENDED"}).
#'   \item \code{M_table}: Data frame summarizing EPCs exceeding the SESOI,
#'     if any.
#'   \item \code{testeffect}: Data frame reporting the smallest tested
#'     standardized perturbations in each direction.
#' }
#'
#' @importFrom lavaan lavaan
#'
#' @seealso \code{\link{epcEquivFit}}
#'
#' @examples
#'
#' library(lavaan)
#'
#' one.model <- ' onefactor  =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 '
#' fit <- cfa(one.model, data = HolzingerSwineford1939)
#' \donttest{
#' epcEquivCheck(fit)
#' }
#'
#' @export
epcEquivCheck <- function(lavaanObj,
                          minRelEffect = 0.75,
                          stdLoad = 0.4,
                          cor = 0.1,
                          corLatent = NULL,
                          corResidual = NULL,
                          stdBeta = 0.1) {

  tol <- 1e-10

  .new_epcEquivCheckStd_infeasible <- function(reason) {
    out <- list(
      feasible = FALSE,
      any_M = NA,
      recommendation = "NOT APPLICABLE",
      reason = reason,
      M_table = NULL,
      testeffect = NULL
    )
    class(out) <- "epcEquivCheckStd"
    out
  }


  if (minRelEffect <= 0 || minRelEffect >= 1) {
    stop("minRelEffect must be between 0 and 1.")
  }

  if (is.null(corLatent))   corLatent   <- cor
  if (is.null(corResidual)) corResidual <- cor

  # Candidate misspecifications and EPC-based feasibility outputs
  miout <- epcEquivFit(lavaanObj, stdLoad = stdLoad, corLatent = corLatent,
                      corResidual = corResidual, stdBeta = stdBeta)

  # Scope checks: features not supported in this standardized-parameter feasibility check
  if (any(miout$op == "|")) {
    stop("Models with categorical indicators are not supported.")
  }
  if (any(miout$op == "~1")) {
    stop("Models with a mean structure are not supported.")
  }
  if (any(miout$op == "<~")) {
    stop("Models with formative indicators are not supported.")
  }
  if (lavaan::lavInspect(lavaanObj, "ngroups") > 1) {
    stop("Multiple-group models are not supported yet.")
  }

  # Extract standardized parameter matrices
  std <- lavaan::lavInspect(lavaanObj, "std.all")
  lambda <- std$lambda
  corpsi <- cov2cor_safe(std$psi)
  cortheta <- cov2cor_safe(std$theta)
  stdbeta <- std$beta
  if (is.null(stdbeta)) stdbeta <- matrix(0, nrow(corpsi), ncol(corpsi),
                                          dimnames = dimnames(corpsi))

  # ---- Existence check: do standardized parameters define ANY PD population Sigma? ----
  residVarPsi0 <- findFactorResidualVar(
    beta = stdbeta,
    corPsi = corpsi,
    totalVarPsi = rep(1, nrow(stdbeta))
  )
  if (any(!is.finite(residVarPsi0)) || any(residVarPsi0 < -tol)) {
    return(.new_epcEquivCheckStd_infeasible("std_params_not_generative"))
  }
  residVarPsi0[abs(residVarPsi0) < tol] <- 0

  Phi0 <- findFactorTotalCov(
    beta = stdbeta,
    corPsi = corpsi,
    errorVarPsi = residVarPsi0
  )

  Sigma_y0 <- lambda %*% Phi0 %*% t(lambda)

  residVarTheta0 <- findIndResidualVar(
    lambda = lambda,
    totalFactorCov = Phi0,
    totalVarTheta = rep(1, nrow(lambda))
  )
  if (any(!is.finite(residVarTheta0)) || any(residVarTheta0 < -tol)) {
    return(.new_epcEquivCheckStd_infeasible("std_params_not_generative"))
  }
  residVarTheta0[abs(residVarTheta0) < tol] <- 0

  Theta0 <- cor2cov_safe(cortheta, sqrt(residVarTheta0))
  Sigma0 <- Sigma_y0 + Theta0

  if (!isPD(Sigma0)) {
    return(.new_epcEquivCheckStd_infeasible("std_params_not_generative"))
  }

  # ---- Search feasibility for positive perturbations ----
  result <- matrix(NA, nrow(miout), nrow(miout))
  sepc_mat <- matrix(NA, nrow(miout), nrow(miout))
  testeffect <- rep(NA_real_, nrow(miout))
  kseq <- seq(minRelEffect, 0.1, length.out = 10) # decreasing misspecification

  for (i in 1:nrow(miout)) {
    row <- miout[i,]
    found <- FALSE
    k_found <- NA_real_

    for (k in kseq) {

      # Start from base standardized parameter matrices each time
      tlambda <- lambda
      tcorpsi <- corpsi
      tcortheta <- cortheta
      tstdbeta <- stdbeta
      tempRelEffect <- NA_real_

      # Apply a single targeted perturbation (scaled by k)
      if (row$op == "=~") {
        tlambda[row$rhs, row$lhs] <- stdLoad * k
        tempRelEffect <- stdLoad
      } else if (row$op == "~~") {

        # latent (psi) vs observed residual (theta) determined by membership
        if (row$lhs %in% rownames(corpsi)) {
          tcorpsi[row$lhs, row$rhs] <- corLatent * k
          tcorpsi[row$rhs, row$lhs] <- corLatent * k
          tempRelEffect <- corLatent
        } else if (row$lhs %in% rownames(lambda)) {
          tcortheta[row$lhs, row$rhs] <- corResidual * k
          tcortheta[row$rhs, row$lhs] <- corResidual * k
          tempRelEffect <- corResidual
        }

      } else if (row$op == "~") {
        tstdbeta[row$lhs, row$rhs] <- stdBeta * k
        tempRelEffect <- stdBeta
      }

      # Derive implied factor residual variances needed for total factor covariance
      residVarPsi <- findFactorResidualVar(
        beta = tstdbeta,
        corPsi = tcorpsi,
        totalVarPsi = rep(1, nrow(tstdbeta))
      )
      if (any(!is.finite(residVarPsi)) || any(residVarPsi < -tol)) next
      residVarPsi[abs(residVarPsi) < tol] <- 0

      Phi <- tryCatch(
        findFactorTotalCov(
          beta = tstdbeta,
          corPsi = tcorpsi,
          errorVarPsi = residVarPsi
        ),
        error = function(e) NULL
      )
      if (is.null(Phi) || any(!is.finite(Phi))) next

      # Implied indicator covariance (excluding residuals)
      Sigma_y_noTheta <- tlambda %*% Phi %*% t(tlambda)

      # Derive indicator residual variances under standardized total variances (=1)
      residVarTheta <- findIndResidualVar(
        lambda = tlambda,
        totalFactorCov = Phi,
        totalVarTheta = rep(1, nrow(tlambda))
      )
      if (any(!is.finite(residVarTheta)) || any(residVarTheta < -tol)) next
      residVarTheta[abs(residVarTheta) < tol] <- 0

      Theta <- cor2cov_safe(tcortheta, sqrt(residVarTheta))
      Sigma_implied <- Sigma_y_noTheta + Theta

      # Only proceed if Sigma is PD and the model is estimable under that population
      if (isPD(Sigma_implied)) {
        tempout <- tryCatch(
          suppressWarnings(lavaan::lavaan(lavaanObj,
                                  sample.cov  = Sigma_implied,
                                  sample.nobs = 1000000L,
                                  std.lv      = TRUE)),
          error = function(e) NULL
        )
        if (is.null(tempout)) next

        if (lavCheckAdmissibleFit(tempout)) {
          tempEqTest <- epcEquivFit(tempout,
                                   stdLoad = stdLoad,
                                   corLatent = corLatent,
                                   corResidual = corResidual,
                                   stdBeta = stdBeta)
          result[i, ] <- tempEqTest[, "decision.ci"]
          sepc_mat[i, ] <- tempEqTest[,"std.epc"]
          found <- TRUE
          k_found <- k
        }
      }

      if (found) break
    }

    testeffect[i] <- k_found * tempRelEffect
  }

  # ---- Search feasibility for negative perturbations ----
  result2 <- matrix(NA, nrow(miout), nrow(miout))
  sepc_mat2 <- matrix(NA, nrow(miout), nrow(miout))
  testeffect2 <- rep(NA_real_, nrow(miout))
  kseq2 <- seq(minRelEffect, 0.1, length.out = 10) # decreasing misspecification

  for (i in 1:nrow(miout)) {
    row <- miout[i,]
    found2 <- FALSE
    k_found2 <- NA_real_

    for (k2 in kseq2) {

      tlambda2 <- lambda
      tcorpsi2 <- corpsi
      tcortheta2 <- cortheta
      tstdbeta2 <- stdbeta
      tempRelEffect2 <- NA_real_
      if (row$op == "=~") {
        tlambda2[row$rhs, row$lhs] <- -stdLoad * k2
        tempRelEffect2 <- -stdLoad
      } else if (row$op == "~~") {

        if (row$lhs %in% rownames(corpsi)) {
          tcorpsi2[row$lhs, row$rhs] <- -corLatent * k2
          tcorpsi2[row$rhs, row$lhs] <- -corLatent * k2
          tempRelEffect2 <- -corLatent
        } else if (row$lhs %in% rownames(lambda)) {
          tcortheta2[row$lhs, row$rhs] <- -corResidual * k2
          tcortheta2[row$rhs, row$lhs] <- -corResidual * k2
          tempRelEffect2 <- -corResidual
        }

      } else if (row$op == "~") {
        tstdbeta2[row$lhs, row$rhs] <- -stdBeta * k2
        tempRelEffect2 <- -stdBeta
      }

      residVarPsi2 <- findFactorResidualVar(
        beta = tstdbeta2,
        corPsi = tcorpsi2,
        totalVarPsi = rep(1, nrow(tstdbeta2))
      )
      if (any(!is.finite(residVarPsi2)) || any(residVarPsi2 < -tol)) next
      residVarPsi2[abs(residVarPsi2) < tol] <- 0

      Phi2 <- tryCatch(
        findFactorTotalCov(
          beta = tstdbeta2,
          corPsi = tcorpsi2,
          errorVarPsi = residVarPsi2
        ),
        error = function(e) NULL
      )
      if (is.null(Phi2) || any(!is.finite(Phi2))) next

      Sigma_y_noTheta2 <- tlambda2 %*% Phi2 %*% t(tlambda2)

      residVarTheta2 <- findIndResidualVar(
        lambda = tlambda2,
        totalFactorCov = Phi2,
        totalVarTheta = rep(1, nrow(tlambda2))
      )
      if (any(!is.finite(residVarTheta2)) || any(residVarTheta2 < -tol)) next
      residVarTheta2[abs(residVarTheta2) < tol] <- 0

      Theta2 <- cor2cov_safe(tcortheta2, sqrt(residVarTheta2))
      Sigma_implied2 <- Sigma_y_noTheta2 + Theta2

      if (isPD(Sigma_implied2)) {
        tempout2 <- tryCatch(
          suppressWarnings(lavaan(lavaanObj,
                                  sample.cov  = Sigma_implied2,
                                  sample.nobs = 1000000L,
                                  std.lv      = TRUE)),
          error = function(e) NULL
        )
        if (is.null(tempout2)) next

        if (lavCheckAdmissibleFit(tempout2)) {
          tempEqTest2 <- epcEquivFit(tempout2,
                                    stdLoad = stdLoad,
                                    corLatent = corLatent,
                                    corResidual = corResidual,
                                    stdBeta = stdBeta)
          result2[i, ] <- tempEqTest2[, "decision.ci"]
          sepc_mat2[i, ] <- tempEqTest2[,"std.epc"]
          found2 <- TRUE
          k_found2 <- k2
        }
      }

      if (found2) break
    }

    testeffect2[i] <- tempRelEffect2 * k_found2
  }

  resultall <- cbind(result, result2)
  M_pos <- extract_M_table(
    result_mat = result,
    miout      = miout,
    sepc_mat   = sepc_mat,
    direction  = "positive"
  )

  M_neg <- extract_M_table(
    result_mat = result2,
    miout      = miout,
    sepc_mat   = sepc_mat2,
    direction  = "negative"
  )

  M_all <- rbind(M_pos, M_neg)

  T_all <- data.frame(
    lhs = miout$lhs,
    op = miout$op,
    rhs = miout$rhs,
    effect_positive = testeffect,
    effect_negative = testeffect2,
    stringsAsFactors = FALSE
  )
  feasible <- TRUE
  any_M <- any(resultall == "M", na.rm = TRUE)

  recommendation <- if (!feasible) {
    "NOT APPLICABLE"
  } else if (any_M) {
    "NOT RECOMMENDED"
  } else {
    "RECOMMENDED"
  }

  out <- list(
    feasible = feasible,
    any_M = any_M,
    recommendation = recommendation,
    M_table = M_all,
    testeffect = T_all
  )

  class(out) <- "epcEquivCheckStd"
  return(out)
}

## ----------------
## Hidden Functions
## ----------------

# lavInspectTotalVar()
# ------------------------------------------------------------------
# Internal lavaan helper used by epcEquivFit() and EPC-related utilities.
# Extracts total variances of observed and latent variables from a
# fitted lavaan object using the implied covariance matrix. The
# function supports single-group and multi-group models and returns
# group-specific variance vectors for use in EPC standardization and
# variance-based transformations. No EPC estimation or statistical
# inference is performed.

##' @importFrom lavaan lavInspect
lavInspectTotalVar <- function(lavaanObj) {
  result <- list()
  nGroups <- lavInspect(lavaanObj, "ngroups")
  cov.all <- lavInspect(lavaanObj, "cov.all")
  if (nGroups == 1) cov.all <- list(cov.all)
  for (i in 1:nGroups) {
    temp <- diag(cov.all[[i]])
    names(temp) <- rownames(cov.all[[i]])
    result[[i]] <- temp
  }
  return(result)
}

# lavInspectResidualVar()
# ------------------------------------------------------------------
# Internal lavaan helper used by epcEquivFit() and EPC-related utilities.
# Extracts estimated residual variances for indicators (theta) and
# latent variables (psi) from a fitted lavaan object, handling both
# single-group and multi-group models.

##' @importFrom lavaan lavInspect
lavInspectResidualVar <- function(lavaanObj) {
  result <- list()
  nGroups <- lavInspect(lavaanObj, "ngroups")
  est <- lavInspect(lavaanObj, "est")
  if (nGroups == 1) est <- list(est)
  for (i in 1:nGroups) {
    temppsi <- diag(est[[i]]$psi)
    names(temppsi) <- rownames(est[[i]]$psi)
    templambda <- est[[i]]$lambda
    if(ncol(templambda) < nrow(templambda)) {
      temptheta <- diag(est[[i]]$theta)
      names(temptheta) <- rownames(est[[i]]$theta)
      temp <- c(temptheta, temppsi)
    } else {
      temp <- temppsi
    }
    result[[i]] <- temp
  }
  return(result)
}

# getTrivialEpc()
# ------------------------------------------------------------------
# Internal utility used by epcEquivFit() and EPC equivalence diagnostics.
# Assigns operator-specific smallest effect size of interest (SESOI)
# values to each fixed parameter based on its type (e.g., factor
# loadings, structural regressions, latent covariances, residual
# covariances, intercepts). The resulting vector defines the magnitude
# of trivial misspecification used for EPC evaluation.
getTrivialEpc <- function(
    mi,
    lv_names,
    stdLoad = 0.4,
    corLatent = 0.1,
    corResidual = 0.1,
    stdBeta = 0.1,
    stdIntcept = 0.2
) {

  result <- numeric(nrow(mi))

  for (i in seq_len(nrow(mi))) {
    op  <- mi[i, "op"]
    lhs <- mi[i, "lhs"]
    rhs <- mi[i, "rhs"]

    if (op == "=~") {
      result[i] <- stdLoad

    } else if (op == "~~") {
      if (lhs %in% lv_names && rhs %in% lv_names) {
        result[i] <- corLatent
      } else {
        result[i] <- corResidual
      }

    } else if (op == "~1") {
      result[i] <- stdIntcept

    } else if (op == "~") {
      result[i] <- stdBeta

    } else {
      result[i] <- NA_real_
    }
  }

  result
}


# unstandardizeEpc()
# ------------------------------------------------------------------
# Internal utility used by epcEquivFit() and related EPC diagnostics.
# Converts standardized effect-size thresholds (SESOI) back to the
# unstandardized EPC scale using total and residual variances of the
# involved variables. The transformation is operator-specific
# (e.g., loadings, regressions, covariances, intercepts) and provides
# unstandardized quantities required for EPC evaluation.
unstandardizeEpc <- function(mi, sesoi, totalVar, residualVar) {
  name <- names(totalVar[[1]])
  lhsPos <- match(mi[,"lhs"], name)
  rhsPos <- match(mi[,"rhs"], name)
  group <- rep(1, nrow(mi))
  if("group" %in% colnames(mi)) group <- mi[,"group"]
  getVar <- function(pos, group) totalVar[[group]][pos]
  getVarRes <- function(pos, group) residualVar[[group]][pos]
  lhsVar <- mapply(getVar, pos=lhsPos, group=group)
  rhsVar <- mapply(getVar, pos=rhsPos, group=group)
  lhsVarRes <- mapply(getVarRes, pos=lhsPos, group=group)
  rhsVarRes <- mapply(getVarRes, pos=rhsPos, group=group)

  FUN <- function(op, lhsVar, rhsVar, lhsVarRes, rhsVarRes, sesoi) {
    if(op == "|") return(NA)
    lhsSD <- sqrt(lhsVar)
    rhsSD <- sqrt(rhsVar)
    lhsSDRes <- sqrt(lhsVarRes)
    rhsSDRes <- sqrt(rhsVarRes)
    if(!is.numeric(sesoi)) sesoi <- as.numeric(sesoi)
    if(op == "=~") {
      return((rhsSD * sesoi) / lhsSD)
    } else if (op == "~~") {
      return(lhsSDRes * sesoi * rhsSDRes)
    } else if (op == "~1") {
      return(lhsSD * sesoi)
    } else if (op == "~") {
      return((lhsSD * sesoi) / rhsSD)
    } else {
      return(NA)
    }
  }
  sesoi <- mapply(FUN, op=mi[,"op"], lhsVar=lhsVar, rhsVar=rhsVar,
                       lhsVarRes=lhsVarRes, rhsVarRes=rhsVarRes, sesoi=sesoi)
  return(sesoi)
}

# standardizeEpc()
# ------------------------------------------------------------------
# Internal utility used by epcEquivFit() and related EPC diagnostics.
# Transforms unstandardized EPCs or SESOI values into standardized
# effect-size metrics based on total and residual variances of the
# involved variables. The standardization is operator-specific
# (e.g., loadings, regressions, covariances, intercepts) and produces
# standardized quantities suitable for comparison against standardized
# SESOI thresholds.
standardizeEpc <- function(mi, totalVar, residualVar, sesoi = NULL) {
  if(is.null(sesoi)) sesoi <- mi[,"epc"]
  name <- names(totalVar[[1]])
  lhsPos <- match(mi[,"lhs"], name)
  rhsPos <- match(mi[,"rhs"], name)
  group <- rep(1, nrow(mi))
  if("group" %in% colnames(mi)) group <- mi[,"group"]
  getVar <- function(pos, group) totalVar[[group]][pos]
  getVarRes <- function(pos, group) residualVar[[group]][pos]
  lhsVar <- mapply(getVar, pos=lhsPos, group=group)
  rhsVar <- mapply(getVar, pos=rhsPos, group=group)
  lhsVarRes <- mapply(getVarRes, pos=lhsPos, group=group)
  rhsVarRes <- mapply(getVarRes, pos=rhsPos, group=group)
  FUN <- function(op, lhsVar, rhsVar, lhsVarRes, rhsVarRes, sesoi) {
    lhsSD <- sqrt(lhsVar)
    rhsSD <- sqrt(rhsVar)
    lhsSDRes <- sqrt(lhsVarRes)
    rhsSDRes <- sqrt(rhsVarRes)
    if(!is.numeric(sesoi)) sesoi <- as.numeric(sesoi)
    if(op == "=~") {
      #stdload = beta * sdlatent / sdindicator = beta * lhs / rhs
      return((sesoi / rhsSD) * lhsSD)
    } else if (op == "~~") {
      #r = cov / (sd1 * sd2)
      return(sesoi / (lhsSDRes * rhsSDRes))
    } else if (op == "~1") {
      #d = meanDiff/sd
      return(sesoi / lhsSD)
    } else if (op == "~") {
      #beta = b * sdX / sdY = b * rhs / lhs
      return((sesoi / lhsSD) * rhsSD)
    } else {
      return(NA)
    }
  }
  stdSesoi <- mapply(FUN, op=mi[,"op"], lhsVar=lhsVar, rhsVar=rhsVar,
                     lhsVarRes=lhsVarRes, rhsVarRes=rhsVarRes, sesoi=sesoi)
  return(stdSesoi)
}

# decisionMIPow()
# ------------------------------------------------------------------
# Core decision rule for power-based (SSV-style) EPC classification.
# Given the significance of the modification index, an indicator of
# high statistical power, the EPC value, and a smallest effect size
# of interest (SESOI), this function classifies the parameter into
# categories reflecting severity and substantive relevance
# (e.g., M, NM, I, EPC:M, EPC:NM).
decisionMIPow <- function(sigMI, highPow, epc, trivialEpc) {
  if(is.na(sigMI) | is.na(highPow)) return(NA)
  if(sigMI & highPow) {
    if(abs(epc) > abs(trivialEpc)) {
      return("EPC:M")
    } else {
      return("EPC:NM")
    }
  } else if (sigMI & !highPow) {
    return("M")
  } else if (!sigMI & highPow) {
    return("NM")
  } else if (!sigMI & !highPow) {
    return("I")
  } else {
    return(NA)
  }
}

# decisionCIEpc()
# ------------------------------------------------------------------
# Core decision rule for CI-based EPC classification.
# Given an EPC confidence interval and a smallest effect size of
# interest (SESOI), this function classifies the parameter as
# Substantial (M), Not Misspecified (NM), Inconclusive (I), or Underpowered (U).
decisionCIEpc <- function(targetval, lower, upper, positiveonly = FALSE) {
  if(is.na(lower) | is.na(upper)) return(NA)
  if(positiveonly) {
    ciwidth <- upper - lower
    trivialwidth <- targetval
    if (ciwidth > trivialwidth) {
      return("U")
    } else if (lower > targetval) {
      return("M")
    } else if (upper < targetval) {
      return("NM")
    } else {
      return("I")
    }
  } else {
    negtargetval <- -targetval
    ciwidth <- upper - lower
    trivialwidth <- 2*targetval
    if(ciwidth > trivialwidth) {
      return("U")
    } else if(lower > targetval | upper < negtargetval) {
      return("M")
    } else if (upper < targetval & negtargetval < lower) {
      return("NM")
    } else {
      return("I")
    }
  }
}

# isPD()
# ------------------------------------------------------------------
# Internal numerical utility used by epcEquivCheck() and related helpers.
# Checks whether a symmetric matrix is positive definite by verifying
# that all eigenvalues exceed a small tolerance. This function is used
# to screen implied population covariance matrices for admissibility
# before refitting models under perturbed parameters.
isPD <- function(M, tol = 1e-8) {
  if (is.list(M)) {
    return(all(vapply(M, isPD, logical(1), tol = tol)))
  }
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  all(ev > tol)
}

# cov2cor_safe()
# ------------------------------------------------------------------
# Internal numerical utility used by epcEquivCheck() and related helpers.
# Converts a covariance matrix to a correlation matrix while safely
# handling zero or non-positive variances. The function computes
# correlations only for valid sub-blocks, enforces unit diagonals,
# and suppresses numerical artifacts to maintain stability in
# downstream standardized-parameter computations.
cov2cor_safe <- function(S) {
  S <- as.matrix(S)
  p <- nrow(S)
  if (p != ncol(S)) stop("Input must be a square matrix.")

  v <- diag(S)

  # initialize correlation matrix
  R <- matrix(0, p, p)
  dimnames(R) <- dimnames(S)

  # indices with positive variance
  pos <- which(v > 0)

  # normal cov2cor where possible
  if (length(pos) > 0) {
    sd <- sqrt(v[pos])
    R[pos, pos] <- S[pos, pos] / (sd %o% sd)
  }

  # enforce symmetry and unit diagonal
  R[!is.finite(R)] <- 0
  diag(R) <- 1

  R
}

# cor2cov_safe()
# ------------------------------------------------------------------
# Internal numerical utility used by epcEquivCheck() and related helpers.
# Converts a correlation matrix to a covariance matrix given a vector
# of standard deviations, while safely handling zero, negative, or
# non-finite standard deviations. The function operates only on valid
# sub-blocks, enforces symmetry, and removes numerical artifacts to
# ensure stability in subsequent covariance-based computations.
cor2cov_safe <- function(R, sd) {
  R <- as.matrix(R)
  p <- nrow(R)
  if (p != ncol(R)) stop("R must be square.")
  if (length(sd) != p) stop("Length of sd must match dimension of R.")

  # initialize covariance matrix
  S <- matrix(0, p, p)
  dimnames(S) <- dimnames(R)

  # indices with positive sd
  pos <- which(sd > 0 & is.finite(sd))

  if (length(pos) > 0) {
    D <- diag(sd[pos], length(pos))
    S[pos, pos] <- D %*% R[pos, pos, drop = FALSE] %*% D
  }

  # enforce symmetry and clean numerical noise
  S[!is.finite(S)] <- 0
  S
}

# extract_M_table()
# ------------------------------------------------------------------
# Internal helper used by epcEquivCheck().
# Extracts and formats combinations of perturbed and resulting fixed
# parameters that yield Substantial (M) EPC decisions. The function maps
# row–column indices from a decision matrix back to the corresponding
# parameter labels in the \code{epcEquivFit} output and returns a tidy summary
# table for reporting purposes only.
extract_M_table <- function(result_mat, miout, sepc_mat, direction) {

  idx <- which(result_mat == "M", arr.ind = TRUE)
  if (nrow(idx) == 0) return(NULL)

  data.frame(
    perturbed_lhs = miout$lhs[idx[, "row"]],
    perturbed_op  = miout$op[idx[, "row"]],
    perturbed_rhs = miout$rhs[idx[, "row"]],

    resulting_lhs = miout$lhs[idx[, "col"]],
    resulting_op  = miout$op[idx[, "col"]],
    resulting_rhs = miout$rhs[idx[, "col"]],

    resulting_sepc = sepc_mat[cbind(idx[, "row"], idx[, "col"])],

    direction = direction,
    stringsAsFactors = FALSE
  )
}

# print.epcEquivCheckStd()
# ------------------------------------------------------------------
# Internal print method for epcEquivCheckStd objects.
# Formats and displays feasibility and recommendation results from
# standardized-parameter EPC equivalence checks.

#' @export
print.epcEquivCheckStd <- function(x, ...) {
  cat("EPC Equivalence Feasibility (Standardized Parameters)\n")
  cat("----------------------------------------------------\n")

  cat("Feasible standardized population:", x$feasible, "\n")
  cat("Any EPC exceeding SESOI:", x$any_M, "\n")
  cat("Recommendation:", x$recommendation, "\n\n")

  if (x$recommendation == "NOT RECOMMENDED") {
    cat("Non-equivalent EPCs detected (summary):\n")
    print(x$M_table)
  } else if (x$recommendation == "RECOMMENDED") {
    cat("No EPC exceeded the SESOI under tested misspecifications.\n")
  } else {
    cat("Standardized parameters do not define a valid population model.\n")
  }

  invisible(x)
}

# summary.epcequivfit.data.frame()
# ------------------------------------------------------------------
# Internal summary constructor used by epcEquivFit().
# Aggregates and ranks fixed-parameter EPC results from a data-frame
# representation into a structured summaryEpcEquivFit object for
# downstream printing and global decision reporting. This function
# performs data processing only and does not conduct EPC estimation
# or statistical inference.

#' @method summary epcequivfit.data.frame
#' @rdname epcEquivFit
#' @param object An object returned by \code{\link{epcEquivFit}}.
#' @param top Number of top-ranked EPCs to display.
#' @param ssv Logical; whether to include power-based diagnostics.
#' @export
summary.epcequivfit.data.frame <- function(object, ..., top = 5, ssv = FALSE) {
  miout <- object
  stopifnot(is.data.frame(miout))
  num_cols <- c(
    "std.epc", "std.sesoi",
    "lower.std.epc", "upper.std.epc",
    "se.epc", "pow"
  )

  for (nm in num_cols) {
    if (!is.null(miout[[nm]])) {
      miout[[nm]] <- as.numeric(as.character(miout[[nm]]))
    }
  }

  miout$std.sesoi[miout$std.sesoi == 0] <- NA
  miout$se.epc[miout$se.epc == 0] <- NA


  # ---- Derived quantities ----
  miout$severity <- abs(miout$std.epc / miout$std.sesoi)

  miout$ci_gap <- with(miout, pmin(
    abs(lower.std.epc -  std.sesoi),
    abs(upper.std.epc -  std.sesoi),
    abs(lower.std.epc +  std.sesoi),
    abs(upper.std.epc +  std.sesoi),
    na.rm = TRUE
  ))

  # ---- EPC equivalence testing (CI-based; primary) ----
  epc_equiv <- list(
    n_M  = sum(miout$decision.ci == "M",  na.rm = TRUE),
    n_I  = sum(miout$decision.ci == "I",  na.rm = TRUE),
    n_U  = sum(miout$decision.ci == "U",  na.rm = TRUE),
    n_NM = sum(miout$decision.ci == "NM", na.rm = TRUE),
    any_M = any(miout$decision.ci == "M", na.rm = TRUE)
  )

  # Substantially Misspecified EPCs
  top_M <- miout[miout$decision.ci == "M", ]
  top_M <- top_M[is.finite(top_M$severity), ]
  top_M <- top_M[order(top_M$severity, decreasing = TRUE), ]
  if (nrow(top_M) > top) top_M <- top_M[1:top, ]

  # CI-inconclusive EPCs (exclude NM)
  top_I_ci <- miout[miout$decision.ci == "I", ]
  top_I_ci <- top_I_ci[is.finite(top_I_ci$ci_gap), ]
  top_I_ci <- top_I_ci[order(top_I_ci$ci_gap, decreasing = TRUE), ]
  if (nrow(top_I_ci) > top) top_I_ci <- top_I_ci[1:top, ]

  # ---- SSV / power-based diagnostics (secondary) ----
  ssv_out <- list(
    n_M  = sum(miout$decision.pow %in% c("M", "EPC:M"),  na.rm = TRUE),
    n_I  = sum(miout$decision.pow == "I",  na.rm = TRUE),
    n_NM = sum(miout$decision.pow %in% c("NM", "EPC:NM"), na.rm = TRUE),
    any_M = any(miout$decision.pow %in% c("M", "EPC:M"), na.rm = TRUE)
  )

  # Substantially Misspecified SSV
  top_M_pow <- miout[miout$decision.pow == "M", ]
  top_M_pow <- top_M_pow[is.finite(top_M_pow$severity), ]
  top_M_pow <- top_M_pow[order(top_M_pow$severity, decreasing = TRUE), ]
  if (nrow(top_M_pow) > top) top_M_pow <- top_M_pow[1:top, ]

  # Underpowered EPCs (exclude NM)
  top_I_pow <- miout[miout$decision.pow == "I", ]
  top_I_pow <- top_I_pow[is.finite(top_I_pow$pow), ]
  top_I_pow <- top_I_pow[order(top_I_pow$pow), ]  # lowest power first
  if (nrow(top_I_pow) > top) top_I_pow <- top_I_pow[1:top, ]

  out <- list(
    epc_equivalence = epc_equiv,
    ssv = ssv_out,
    top_non_equiv = top_M,
    top_inconclusive_ci = top_I_ci,
    top_non_pow = top_M_pow,
    top_inconclusive_pow = top_I_pow,
    show_ssv = ssv
  )

  class(out) <- "summaryEpcEquivFit"
  out
}

# global_localfit_decision()
# ------------------------------------------------------------------
# Internal decision utility used by print.summaryEpcEquivFit().
# Aggregates local fixed-parameter classifications (e.g., Substantial,
# Inconclusive, Underpowered, Not Misspecified) into a single global
# decision using a conservative priority rule.
global_localfit_decision <- function(n_M, n_I, n_U, n_NM) {

  if (n_U > 0) {
    return("UNDERPOWERED")
  }

  if (n_M > 0) {
    return("SUBSTANTIAL MISSPECIFICATION")
  }

  if (n_I > 0) {
    return("INCONCLUSIVE")
  }

  if (n_NM > 0) {
    return("EQUIVALENT")
  }

  return("NO PARAMETERS EVALUATED")
}

# print.summaryEpcEquivFit()
# ------------------------------------------------------------------
# Internal print method for summaryEpcEquivFit objects.
# Formats and displays results from global EPC equivalence testing
# and the Saris, Satorra, and van der Veld (2009) power-based method.
# This function summarizes aggregate fixed-parameter diagnostics
# for display purposes only and is not part of the analytical
# implementation of epcEquivFit().

#' @export
print.summaryEpcEquivFit <- function(x, ...) {

  cat("Global EPC Evaluation Summary\n")
  cat("--------------------------------\n\n")

  ## ---- EPC equivalence testing (primary) ----
  cat("[1. EPC Equivalence Testing: CI-based]\n")
  cat("Substantially Misspecified (M):",  x$epc_equivalence$n_M,  "\n")
  cat("Inconclusive (I):",    x$epc_equivalence$n_I,  "\n")
  cat("CI-Underpowered (U):",      x$epc_equivalence$n_U,  "\n")
  cat("Trivial / Not Misspecified (NM):", x$epc_equivalence$n_NM, "\n\n")

  global_epc <- global_localfit_decision(
    n_M  = x$epc_equivalence$n_M,
    n_I  = x$epc_equivalence$n_I,
    n_U  = x$epc_equivalence$n_U,
    n_NM = x$epc_equivalence$n_NM
  )

  cat("Global EPC Equivalence Decision:", global_epc, "\n\n")

  if (x$epc_equivalence$n_M > 0) {
    cat("1.1 Top Substantially Misspecified EPCs (ranked by |std.epc / SESOI|):\n")
    print(
      x$top_non_equiv[
        , c("lhs","op","rhs","std.epc","std.sesoi","severity"),
        drop = FALSE
      ]
    )
    cat("\n")
  }

  if (x$epc_equivalence$n_I > 0) {
    cat("1.2 Top CI-inconclusive EPCs\n")
    cat("(ranked by distance to equivalence bounds; larger = needs narrower CI):\n")
    print(
      x$top_inconclusive_ci[
        , c("lhs","op","rhs","lower.std.epc","upper.std.epc","ci_gap"),
        drop = FALSE
      ]
    )
    cat("\n")
  }

  if(isTRUE(x$show_ssv)) {
    ## ---- SSV / power-based diagnostics (secondary) ----
    cat("[2. Saris, Satorra, Van der Veld (2009) / Power-based Diagnostics]\n")
    cat("Substantially Misspecified (M, EPC:M):",   x$ssv$n_M,  "\n")
    cat("Inconclusive (I):", x$ssv$n_I,  "\n")
    cat("Trivial / Not Misspecified (NM, EPC:NM):", x$ssv$n_NM, "\n\n")

    global_ssv <- global_localfit_decision(
      n_M  = x$ssv$n_M,
      n_I  = x$ssv$n_I,
      n_U  = 0,
      n_NM = x$ssv$n_NM
    )

    cat("Global SSV / Power-based Decision:", global_ssv, "\n\n")

    if (x$ssv$n_M > 0) {
      cat("2.1 Top Substantially Misspecified EPCs (ranked by |std.epc / SESOI|):\n")
      print(
        x$top_non_pow[
          , c("lhs","op","rhs","std.epc","std.sesoi","severity"),
          drop = FALSE
        ]
      )
      cat("\n")
    }

    if (x$ssv$n_I > 0) {
      cat("2.2 Top inconclusive EPCs\n")
      cat("(ranked by lowest approximate power):\n")
      print(
        x$top_inconclusive_pow[
          , c("lhs","op","rhs","pow","se.epc"),
          drop = FALSE
        ]
      )
    }
  }

  invisible(x)
}

# augmentBetaWithGamma()
# ------------------------------------------------------------------
# Copied from simsem. Internal utility used by epcEquivCheck().
# Augments the factor regression coefficient matrix (beta) with
# regression paths from exogenous observed variables (gamma) so that
# endogenous and exogenous predictors can be handled within a single
# block-recursive system.
# NOTE: This function is prepared for future support of exogenous covariates
augmentBetaWithGamma <- function(beta, gamma) {
  nf <- nrow(beta)
  nz <- ncol(gamma)
  result <- matrix(0, nf + nz, nf + nz)
  result[(nz + 1):(nz + nf), (nz + 1):(nz + nf)] <- beta
  result[(nz + 1):(nz + nf), 1:nz] <- gamma
  result
}

# augmentPsiWithExogenousCov()
# ------------------------------------------------------------------
# Copied from simsem. Internal utility used by epcEquivCheck().
# Augments the factor residual covariance matrix (psi) with the
# covariance matrix of exogenous observed variables (sigmaxx),
# producing a joint covariance structure compatible with the
# augmented regression system.
# NOTE: This function is prepared for future support of exogenous covariates
augmentPsiWithExogenousCov <- function(psi, sigmaxx) {
  nf <- nrow(psi)
  nz <- ncol(sigmaxx)
  result <- matrix(0, nf + nz, nf + nz)
  result[(nz + 1):(nz + nf), (nz + 1):(nz + nf)] <- psi
  result[1:nz, 1:nz] <- sigmaxx
  result
}

# findRecursiveSet()
# ------------------------------------------------------------------
# Copied from simsem. Internal helper for epcEquivCheck().
# Identifies a recursive (block-triangular) ordering of variables
# implied by the regression coefficient matrix (beta). This ordering
# is required to compute analytical residual variances in recursive
# SEMs.
findRecursiveSet <- function(beta) {
  if (any(diag(beta) != 0))
    stop("Diagonal elements of beta must be zero.")
  result <- list()
  ni <- nrow(beta)
  fix.variable <- rep(FALSE, ni)
  ni.sofar <- 0
  i <- 1
  while (ni.sofar < ni) {
    temp <- findRowZero(beta, fix.variable)
    if (is.null(temp))
      stop("The matrix is not recursive.")
    fix.variable[temp] <- TRUE
    result[[i]] <- temp
    i <- i + 1
    ni.sofar <- ni.sofar + length(temp)
  }
  return(result)
}

# findRowZero()
# ------------------------------------------------------------------
# Internal helper for recursive SEM utilities (copied/adapted from simsem).
# Identifies rows of a square coefficient matrix whose unfixed entries
# are all zero. This is used when constructing a recursive (block-triangular)
# ordering of variables implied by a regression matrix (e.g., beta).
#
# Rows already marked as fixed are ignored, and the function returns the
# indices of rows that contain only zeros among the remaining unfixed columns.
# If no such rows exist, NULL is returned.
#
# This function is used by findRecursiveSet() to iteratively determine
# which variables can be treated as exogenous at each recursion step.
findRowZero <- function(square.matrix, is.row.fixed = FALSE) {
  ni <- nrow(square.matrix)
  if (length(is.row.fixed) == 1) {
    if (is.row.fixed == FALSE)
      is.row.fixed <- rep(FALSE, ni)
  }
  result <- NULL
  desired.zero <- sum(!is.row.fixed)
  for (i in 1:ni) {
    if (is.row.fixed[i] == FALSE) {
      temp <- sum(square.matrix[i, !is.row.fixed] == 0, na.rm = TRUE)
      if (temp == desired.zero)
        result <- c(result, i)
    }
  }
  return(result)
}

# findIndResidualVar()
# ------------------------------------------------------------------
# Copied from simsem. Internal utility used by epcEquivCheck().
# Computes indicator residual variances given factor loadings, total
# factor covariance, and total indicator variances. The resulting
# implied residual variances are used when constructing EPC values
# under trivial misspecification based on a fitted lavaan object.
findIndResidualVar <- function(lambda, totalFactorCov, totalVarTheta = NULL,
                               kappa = NULL, covcov = NULL) {
  ni <- nrow(lambda)
  if (is.null(totalVarTheta)) totalVarTheta <- rep(1, ni)
  factor.part <- lambda %*% totalFactorCov %*% t(lambda)
  error.var <- totalVarTheta - pmax(diag(factor.part), 0)
  if (!is.null(kappa)) error.var <- error.var - diag(kappa %*% covcov %*% t(kappa))
  error.var[(error.var < 0) & (totalVarTheta == 0)] <- 0
  return(as.vector(error.var))
}

# findFactorTotalCov()
# ------------------------------------------------------------------
# Copied from simsem. Internal utility used by epcEquivCheck().
# Computes the total covariance matrix of latent factors implied by
# a recursive regression structure and factor residual covariances.
# The resulting implied factor covariance matrix is used as an
# intermediate quantity for evaluating EPCs under trivial
# misspecification based on a fitted lavaan object.
findFactorTotalCov <- function(beta, psi = NULL, corPsi = NULL,
                               totalVarPsi = NULL, errorVarPsi = NULL,
                               gamma = NULL, covcov = NULL) {
  if (is.null(psi)) {
    if (is.null(errorVarPsi))
      errorVarPsi <- findFactorResidualVar(beta, corPsi, totalVarPsi)
    psi <- suppressWarnings(cor2cov_safe(as.matrix(corPsi), sqrt(errorVarPsi)))
  }
  iden <- diag(nrow(beta))
  inv <- solve(iden - beta)
  facTotalCov <- inv %*% psi %*% t(inv)
  if (!is.null(gamma)) {
    facTotalCov <- facTotalCov + (inv %*% gamma %*% covcov %*% t(gamma) %*% t(inv))
  }
  return(facTotalCov)
}

# findFactorResidualVar()
# ------------------------------------------------------------------
# Copied from simsem. Core analytical routine used by epcEquivCheck().
# Computes factor residual variances implied by a recursive SEM given
# total variances, factor correlations, and regression coefficients.
# These analytically derived residual variances provide the basis
# for evaluating EPCs under trivial misspecification based on a
# fitted lavaan object.
findFactorResidualVar <- function(beta, corPsi, totalVarPsi = NULL,
                                  gamma = NULL, covcov = NULL) {
  stopifnot(nrow(beta) == ncol(beta))
  stopifnot(nrow(corPsi) == ncol(corPsi))
  stopifnot(nrow(beta) == nrow(corPsi))
  if (!is.null(gamma)) {
    if (is.null(totalVarPsi)) totalVarPsi <- rep(1, nrow(beta))
    beta <- augmentBetaWithGamma(beta, gamma)
    corPsi <- augmentPsiWithExogenousCov(corPsi, cov2cor(covcov))
    totalVarPsi <- c(diag(covcov), totalVarPsi)
  }
  if (sum(diag(corPsi)) == 0) diag(corPsi) <- 1
  ni <- nrow(beta)
  set <- findRecursiveSet(beta)
  errorVar <- rep(1, ni)
  if (is.null(totalVarPsi))  totalVarPsi <- rep(1, ni)
  errorVar[set[[1]]] <- totalVarPsi[set[[1]]]
  iv <- NULL
  ivCor <- corPsi[set[[1]], set[[1]]]
  startVar <- totalVarPsi[set[[1]]]
  ivCov <- suppressWarnings(cor2cov_safe(as.matrix(ivCor), sqrt(startVar)))
  for (i in seq_len(length(set) - 1)) {
    iv <- c(iv, set[[i]])
    dv <- set[[i + 1]]
    tempBeta <- matrix(beta[dv, iv], nrow = length(dv), ncol = length(iv))
    var.reg <- (tempBeta %*% ivCov %*% t(tempBeta))
    tempPsi <- corPsi[dv, dv]
    tempPsiSd <- rep(0, length(dv))
    for (j in 1:length(dv)) {
      errorVar[dv[j]] <- totalVarPsi[dv[j]] - var.reg[j, j]
      if (is.na(errorVar[dv[j]]) || errorVar[dv[j]] < 0) {
        tempPsiSd[j] <- NA_real_
      } else {
        tempPsiSd[j] <- sqrt(errorVar[dv[j]])
      }
    }
    if (i < (length(set) - 1)) {
      tempPsi <- suppressWarnings(cor2cov_safe(as.matrix(tempPsi), as.matrix(tempPsiSd))) #### FLAG CHANGES
      real.tempPsi <- matrix(0, length(iv) + length(dv), length(iv) + length(dv))
      real.tempPsi[1:length(iv), 1:length(iv)] <- ivCov
      real.tempPsi[(length(iv) + 1):(length(iv) + length(dv)), (length(iv) +
                                                                  1):(length(iv) + length(dv))] <- tempPsi
      agg <- c(iv, dv)
      blank.path <- matrix(0, nrow = length(iv), ncol = length(agg))
      temp.path2 <- beta[dv, agg]
      temp.path2 <- rbind(blank.path, temp.path2)
      ID <- matrix(0, length(agg), length(agg))
      diag(ID) <- 1
      ivCov <- solve(ID - temp.path2) %*% real.tempPsi %*% t(solve(ID - temp.path2))
    }
  }
  if (!is.null(gamma)) {
    errorVar <- errorVar[(nrow(covcov) + 1):length(errorVar)]
  }
  return(as.vector(errorVar))
}

# lavCheckAdmissibleFit()
# ------------------------------------------------------------------
# Internal admissibility and numerical-sanity check for fitted lavaan
# objects. This function verifies that a fitted model is suitable for
# downstream analytical use (e.g., EPC evaluation, population-based
# refitting) by screening for common estimation pathologies.
#
# The checks include:
#   - model convergence
#   - positive definiteness of observed, implied, residual, and latent
#     covariance matrices
#   - strictly positive residual and latent variances
#   - positive definiteness of the parameter vcov matrix (if available)
#   - positive definiteness of the Hessian matrix (if available)
#
# The function is intentionally conservative: failure of any check
# returns FALSE. It is not intended as a model-fit diagnostic, but as a
# numerical admissibility gate for internal algorithms.
lavCheckAdmissibleFit <- function(fit, tol = 1e-8) {
  cov.ov   <- lavaan::lavInspect(fit, "cov.ov")    # sample cov / polychoric
  sigma    <- lavaan::lavInspect(fit, "sigma.hat") # model-implied cov
  theta    <- lavaan::lavInspect(fit, "theta")     # residual covariances
  psi      <- lavaan::lavInspect(fit, "cov.lv")       # latent covariances
  vc       <- lavaan::lavInspect(fit, "vcov")       # Variance-Covariance Matrix
  # hessian may not exist (e.g., se = "none")
  hess <- try(lavaan::lavInspect(fit, "hessian"), silent = TRUE)
  hess_pd <- if (inherits(hess, "try-error") || is.null(hess)) NA else isPD(hess, tol)

  vccheck <- TRUE
  if(!is.null(vc)) {
    if(all(round(vc,4) == 0)) {
      vccheck <- TRUE
    } else {
      vccheck <- isPD(vc, tol)
    }
  }
  all(
    converged      = lavaan::lavInspect(fit, "converged"),

    cov.ov.pd      = isPD(cov.ov, tol),
    sigma.pd       = isPD(sigma, tol),

    theta.pd       = isPD(theta, tol),
    psi.pd         = isPD(psi, tol),

    vc             = vccheck,

    theta.var.pos  = diag_pos(theta, tol),  # no negative residual variances
    psi.var.pos    = diag_pos(psi, tol),    # no negative latent variances

    hessian.pd     = hess_pd
  )
}

# diag_pos()
# ------------------------------------------------------------------
# Internal numerical utility for model-admissibility checks.
# Tests whether all diagonal elements of a matrix are strictly positive
# up to a numerical tolerance. This is used to screen for negative or
# zero variances in residual (theta) and latent (psi) covariance
# matrices returned by lavaan.
#
# If a list of matrices is supplied (e.g., multigroup models), the
# check is applied recursively to each element.
diag_pos <- function(M, tol = 1e-8) {
  if (is.list(M)) {
    return(all(vapply(M, diag_pos, logical(1), tol = tol)))
  }
  all(diag(M) > tol)
}

