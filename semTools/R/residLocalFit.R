# TO DO: Opdyke distribution in resEquivFit()

#' Residual-Correlation Equivalence Fit Evaluation
#'
#' Evaluates model fit from an equivalence-testing perspective by
#' applying equivalence tests to model residual correlations. The
#' procedure compares each observed residual correlation against a
#' smallest residual correlation of interest (SESOI) and classifies
#' local residual correlations as substantially misspecified,
#' trivially misspecified, inconclusive, or underpowered.
#'
#' The resulting local classifications can be summarized to yield a
#' global residual-correlation equivalence-style fit evaluation.
#'
#' @param lavaanObj A fitted \code{lavaan} object used to evaluate model fit.
#' @param rsesoi Smallest residual correlation of interest. Residual
#'   correlations smaller than this value in absolute magnitude are treated
#'   as trivially misspecified. Default is 0.10.
#' @param method Method used to construct confidence intervals. Available
#'   options are \code{"wald"}, \code{"fisherwald"}, \code{"boot"}, and
#'   \code{"all"}. \code{"wald"} uses residual correlations directly,
#'   \code{"fisherwald"} uses Fisher-z transformed residual correlations,
#'   \code{"boot"} uses bootstrap-based confidence intervals, and
#'   \code{"all"} returns all available methods.
#' @param adjust.method Multiplicity adjustment method applied to confidence
#'   intervals. Available options are \code{"none"}, \code{"bonferroni"},
#'   \code{"max"}, and \code{"all"}. \code{"max"} is available only for
#'   bootstrap-based inference.
#' @param cilevel Confidence level used for equivalence-testing confidence
#'   intervals. Default is 0.90.
#' @param R Number of bootstrap replications used when
#'   \code{method = "boot"} or \code{method = "all"}. Default is 1000.
#' @param seed Random seed used for bootstrap resampling. Default is 123321.
#'
#' @details
#' This function provides a residual-based equivalence-testing alternative
#' to traditional fit-index interpretation. Instead of asking whether the
#' model fits exactly, it evaluates whether each residual correlation is
#' small enough to be considered substantively trivial relative to
#' \code{rsesoi}.
#'
#' For \code{method = "wald"}, equivalence bounds are defined on the
#' residual-correlation scale. For \code{method = "fisherwald"} and
#' \code{method = "boot"}, equivalence bounds are defined on the
#' Fisher-z scale.
#'
#' When \code{adjust.method = "bonferroni"}, confidence intervals are
#' widened using a Bonferroni adjustment across all evaluated residual
#' correlations. When \code{adjust.method = "max"}, bootstrap-max
#' confidence intervals are used to account for multiple residual
#' correlations simultaneously.
#'
#' Bootstrap-based methods require raw data to be available in
#' \code{lavaanObj}. If the model was fitted from a covariance matrix only,
#' \code{method = "boot"} will stop with an error, whereas
#' \code{method = "all"} will skip the bootstrap method with a warning.
#'
#' Models with a mean structure are not currently supported for
#' residual-correlation equivalence testing.
#'
#' @return An object of class \code{resEquivFit}. The object is a list
#' containing:
#' \enumerate{
#'   \item Analysis settings, including \code{rsesoi}, \code{cilevel},
#'   \code{method}, \code{adjust.method}, and sample size information
#'   in \code{N}.
#'   \item \code{wald}, if requested, a data frame containing observed
#'   correlations, model-implied correlations, residual correlations,
#'   standard errors, SESOI bounds, confidence interval bounds, and
#'   equivalence decisions.
#'   \item \code{fisherwald}, if requested, a data frame containing the
#'   same residual-correlation information evaluated on the Fisher-z scale.
#'   \item \code{boot}, if requested and available, a data frame containing
#'   bootstrap-based confidence intervals and equivalence decisions.
#'   \item Bootstrap diagnostics, when applicable, including
#'   \code{bootthetahat}, \code{bootse}, and \code{bootMb}.
#' }
#'
#' Decision labels are:
#' \describe{
#'   \item{M}{Substantially misspecified.}
#'   \item{I}{Inconclusive.}
#'   \item{NM}{Trivially misspecified.}
#'   \item{U}{Underpowered or too imprecise to evaluate equivalence relative
#'   to the SESOI.}
#' }
#'
#' @seealso \code{\link{epcEquivFit}}, \code{\link{resEquivCheck}}
#'
#' @importFrom stats qnorm quantile sd
#' @importFrom lavaan lavInspect lavResiduals
#'
#' @examples
#'
#' library(lavaan)
#'
#' one.model <- ' onefactor =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 '
#' fit <- cfa(one.model, data = HolzingerSwineford1939)
#'
#' out <- resEquivFit(fit)
#' out
#' summary(out)
#'
#' @export
resEquivFit <- function(
    lavaanObj,
    rsesoi = 0.10,
    method = "wald",
    adjust.method = "none",
    cilevel = 0.90,
    R = 1000,
    seed = 123321
) {
  methodargs <- c("wald", "fisherwald", "boot", "all")
  adjustmethodargs <- c("none", "bonferroni", "max", "all")
  method       <- match.arg(method, methodargs)
  adjust.method <- match.arg(adjust.method, adjustmethodargs)

  if (!inherits(lavaanObj, "lavaan")) {
    stop("lavaanObj must be a fitted lavaan object.")
  }

  set.seed(seed)

  ## --- setup ---
  rsesoi <- abs(rsesoi)
  fsesoi <- atanh(rsesoi)
  alpha  <- 1 - cilevel
  crit   <- qnorm(1 - alpha/2) # Symmatric cilevel% confidence interval

  G <- lavInspect(lavaanObj, "nGroups")
  if (lavInspect(lavaanObj, "meanStructure")) {
    warning(
      "Models with a mean structure are not currently supported. ",
      "Residual-correlation equivalence testing is available only for covariance structures."
    )
  }
  N <- lavInspect(lavaanObj, "nObs")

  if(G == 1) {
    obs <- cov2cor(lavInspect(lavaanObj, "observed")$cov)
    imp <- cov2cor(lavInspect(lavaanObj, "implied")$cov)

    res <- lavResiduals(lavaanObj, type = "cor.bollen", se = TRUE)
    resr  <- res$cov
    resse <- res$cov.se

    idx <- lower.tri(obs, diag=FALSE)
    p   <- sum(idx)

    critbon <- qnorm(1 - alpha / (2*p))

    out <- data.frame(
      row   = rownames(obs)[row(obs)[idx]],
      col   = colnames(obs)[col(obs)[idx]],
      obs   = obs[idx],
      imp   = imp[idx],
      diffr = resr[idx],
      se    = resse[idx],
      stringsAsFactors = FALSE
    )
  } else {
    obs <- lapply(lapply(lavInspect(lavaanObj, "observed"), "[[", "cov"), cov2cor)
    imp <- lapply(lapply(lavInspect(lavaanObj, "implied"), "[[", "cov"), cov2cor)

    res <- lavResiduals(lavaanObj, type = "cor.bollen", se = TRUE)
    resr  <- lapply(res, "[[", "cov")
    resse <- lapply(res, "[[", "cov.se")

    idx <- lapply(obs, lower.tri, diag=FALSE)
    p   <- sum(sapply(idx, sum))

    critbon <- qnorm(1 - alpha / (2*p))

    group_labels <- names(obs)
    names(N) <- group_labels
    out <- lapply(seq_along(obs), function(g) {

      idx_g <- idx[[g]]

      data.frame(
        group = group_labels[g],
        row   = rownames(obs[[g]])[row(obs[[g]])[idx_g]],
        col   = colnames(obs[[g]])[col(obs[[g]])[idx_g]],
        obs   = obs[[g]][idx_g],
        imp   = imp[[g]][idx_g],
        diffr = resr[[g]][idx_g],
        se    = resse[[g]][idx_g],
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, out)
  }
  results <- list()
  results$rsesoi <- rsesoi
  results$cilevel <- cilevel
  results$method <- method
  results$adjust.method <- adjust.method
  results$N <- N

  ## ============================================================
  ## Method 1: Wald CI on residual correlation
  ## ============================================================

  if (method %in% c("wald", "all")) {

    fobs <- atanh(out$obs)
    lses <- out$obs - tanh(fobs + fsesoi)
    uses <- out$obs - tanh(fobs - fsesoi)

    out_wald <- out
    out_wald$lowersesoi <- lses
    out_wald$uppersesoi <- uses

    if(adjust.method %in% c("none", "all")) {
      lci <- out$diffr - crit * out$se
      uci <- out$diffr + crit * out$se
      out_wald$cilower <- lci
      out_wald$ciupper <- uci
      out_wald$decision <- mapply(
        decisionCIequivResidual,
        lses, uses, lci, uci
      )
    }


    if(adjust.method %in% c("bonferroni", "all")) {
      lcibon <- out$diffr - critbon * out$se
      ucibon <- out$diffr + critbon * out$se
      out_wald$cilowerbon <- lcibon
      out_wald$ciupperbon <- ucibon
      out_wald$decisionbon <- mapply(
        decisionCIequivResidual,
        lses, uses, lcibon, ucibon
      )
    }
    results$wald <- out_wald
  }

  ## ============================================================
  ## Method 2: Wald CI on Fisher-z residuals
  ## ============================================================

  if (method %in% c("fisherwald", "all")) {

    fobs <- atanh(out$obs)
    lses <- -fsesoi
    uses <- fsesoi

    difff <- atanh(clip(out$obs)) - atanh(clip(out$imp))

    hz <- 1 - ((1 - out$obs^2) / (1 - out$imp^2))
    he <- 1 / (1 - out$imp^2)

    if(G > 1) {
      sed <- sqrt((hz^2 / (N[out$group] - 3)) + (he^2 * out$se^2))
    } else {
      sed <- sqrt((hz^2 / (N - 3)) + (he^2 * out$se^2))
    }

    out_fz <- out
    out_fz$difff   <- difff
    out_fz$lowersesoi <- lses
    out_fz$uppersesoi <- uses

    if(adjust.method %in% c("none", "all")) {
      lci <- difff - crit * sed
      uci <- difff + crit * sed
      out_fz$cilower <- lci
      out_fz$ciupper <- uci
      out_fz$decision <- mapply(
        decisionCIequivResidual,
        lses, uses, lci, uci
      )
    }

    if(adjust.method %in% c("bonferroni", "all")) {
      lcibon <- difff - critbon * sed
      ucibon <- difff + critbon * sed
      out_fz$cilowerbon <- lcibon
      out_fz$ciupperbon <- ucibon
      out_fz$decisionbon <- mapply(
        decisionCIequivResidual,
        lses, uses, lcibon, ucibon
      )
    }

    results$fisherwald <- out_fz
  }

  ## ============================================================
  ## Method 3: Bootstrap (normal or max)
  ## ============================================================

  if (method %in% c("boot", "all")) {
    .has_raw_data <- function(lavaanObj, G) {

      dat <- tryCatch(
        lavInspect(lavaanObj, "data"),
        error = function(e) NULL
      )

      if(G == 1) {
        return(is.matrix(dat))
      } else {
        return(is.list(dat) && is.matrix(dat[[1]]))
      }
    }

    has_data <- .has_raw_data(lavaanObj, G = G)

    if (!has_data) {
      if (method == "boot") {
        stop(
          "Bootstrap-based equivalence testing requires raw data.\n",
          "The supplied lavaan object was fitted without raw data\n",
          "(e.g., using sample.cov =).\n",
          "Refit the model using the `data =` argument."
        )
      }

      if (method == "all") {
        warning(
          "Raw data not found in lavaan object.\n",
          "Bootstrap-based equivalence testing will be skipped."
        )
      }
    }

    if(has_data) {
      dat <- tryCatch(
        lavInspect(lavaanObj, "data"),
        error = function(e) NULL
      )

      if(G > 1){
        FUNG <- function(data, ind, idx, group_var, group_labels) {

          outb <- rep(NA_real_, sum(sapply(idx, sum)))
          datb <- Map(function(d, i) d[i, , drop = FALSE], data, ind)
          datb <- Map(function(d, g) {
            d <- as.data.frame(d)
            d[[group_var]] <- factor(g, levels = group_labels)
            return(d)
            }, datb, group_labels)
          datb <- do.call(rbind, datb)
          fitb <- tryCatch(
            update(lavaanObj, data = datb, warn = FALSE, se="none"),
            error = function(e) NULL
          )
          if (is.null(fitb)) return(outb)
          if (!lavInspect(fitb, "converged")) return(outb)

          Sigma <- lavInspect(fitb, "implied")
          Sigma <- lapply(Sigma, "[[", "cov")
          eigenimp <- sapply(lapply(Sigma, eigen, symmetric=TRUE, only.values=TRUE), "[[", "values")
          if (!all(eigenimp > 0)) return(outb)

          obsb <- lavInspect(fitb, "observed")
          obsb <- lapply(obsb, "[[", "cov")

          sd_o <- lapply(obsb, function(x) sqrt(diag(x)))
          sd_s <- lapply(Sigma, function(x) sqrt(diag(x)))

          r_o <- Map(function(x, s) clip(x / (s %o% s)), obsb, sd_o)
          r_s <- Map(function(x, s) clip(x / (s %o% s)), Sigma, sd_s)

          outlist <- Map(function(o, s, i) atanh(o[i]) - atanh(s[i]), o=r_o, s=r_s, i=idx)
          outvec <- unlist(outlist, use.names = FALSE)
          return(outvec)
        }

        inds <- lapply(seq_len(R), function(x) {lapply(N, function(n) sample(seq_len(n), replace=TRUE))})
        idx <- lapply(obs, lower.tri, diag=FALSE)
        group_var <- lavInspect(lavaanObj, "group")
        theta_hat <- FUNG(data=dat, ind=lapply(N, seq_len), idx=idx, group_var=group_var, group_labels=group_labels)
        Tmat <- sapply(inds, function(ind)
          FUNG(data = dat,
               ind = ind,
               idx = idx,
               group_var = group_var,
               group_labels = group_labels)
        )
        Tmat <- t(Tmat)
      } else {
        FUN <- function(data, ind, idx) {
          outb <- rep(NA_real_, sum(idx))
          fitb <- tryCatch(
            update(lavaanObj, data = data[ind, ], warn = FALSE, se="none"),
            error = function(e) NULL
          )
          if (is.null(fitb)) return(outb)
          if (!lavInspect(fitb, "converged")) return(outb)

          Sigma <- tryCatch(lavInspect(fitb, "implied")$cov, error = function(e) NULL)
          if (is.null(Sigma)) return(outb)
          if (!all(eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values > 0)) return(outb)

          obsb <- tryCatch(lavInspect(fitb, "observed")$cov, error = function(e) NULL)
          if (is.null(obsb)) return(outb)

          sd_o <- sqrt(diag(obsb))
          sd_s <- sqrt(diag(Sigma))

          r_o <- obsb / (sd_o %o% sd_o)
          r_s <- Sigma / (sd_s %o% sd_s)

          r_o <- clip(r_o)
          r_s <- clip(r_s)

          atanh(r_o[idx]) - atanh(r_s[idx])
        }
        inds <- lapply(seq_len(R), function(x) {sample(seq_len(N), replace=TRUE)})
        theta_hat <- FUN(data=dat, ind=seq_len(N), idx=idx)
        Tmat <- sapply(inds, function(ind)
          FUN(data = dat,
              ind = ind,
              idx = idx)
        )
        Tmat <- t(Tmat)
      }


      good <- apply(!is.na(Tmat), 1, all)
      Tmat <- Tmat[good, , drop = FALSE]

      ## normal bootstrap
      fobs <- atanh(out$obs)
      lses <- -fsesoi
      uses <- fsesoi
      seboot <- apply(Tmat, 2, sd, na.rm = TRUE)
      results$bootthetahat <- theta_hat
      results$bootse <- seboot
      out_boot <- out
      out_boot$lowersesoi <- lses
      out_boot$uppersesoi <- uses

      if(adjust.method %in% c("none", "all")) {
        lci <- theta_hat - crit * seboot
        uci <- theta_hat + crit * seboot
        out_boot$cilower <- lci
        out_boot$ciupper <- uci
        out_boot$decision <- mapply(
          decisionCIequivResidual,
          lses, uses, lci, uci
        )
      }

      if(adjust.method %in% c("bonferroni", "all")) {
        lcibon <- theta_hat - critbon * seboot
        ucibon <- theta_hat + critbon * seboot
        out_boot$cilowerbon <- lcibon
        out_boot$ciupperbon <- ucibon
        out_boot$decisionbon <- mapply(
          decisionCIequivResidual,
          lses, uses, lcibon, ucibon
        )
      }

      ## bootstrap-max
      if(adjust.method %in% c("max", "all")) {
        Mb <- apply(abs(sweep(Tmat, 2, theta_hat)), 1, max, na.rm = TRUE)
        results$bootMb <- Mb
        c_alpha <- quantile(Mb, probs = cilevel, na.rm = TRUE)

        lci_m <- theta_hat - c_alpha
        uci_m <- theta_hat + c_alpha
        out_boot$cilowermax <- lci_m
        out_boot$ciuppermax <- uci_m
        out_boot$decisionmax<- mapply(
          decisionCIequivResidual,
          lses, uses, lci_m, uci_m
        )
      }
      results$boot <- out_boot
    }
  }

  class(results) <- "resEquivFit"
  results
}

# update.resEquivFit()
# ------------------------------------------------------------------
# Update the results from resEquivFit with different rsesoi or cilevel.
# This is useful when bootstrap is used and users do not need to rerun
# the bootstrap procedure.

#' @method update resEquivFit
#' @rdname resEquivFit
#' @param object An object returned by \code{\link{resEquivFit}}.
#' @param rsesoi Smallest residual correlation of interest.
#' @param cilevel Confidence level used for equivalence-testing confidence intervals.
#' @export
update.resEquivFit <- function(object,
                               rsesoi = NULL,
                               cilevel = NULL) {

  if (!inherits(object, "resEquivFit")) {
    stop("object must be of class 'resEquivFit'")
  }

  ## -----------------------------
  ## Update rsesoi if requested
  ## -----------------------------

  if (!is.null(rsesoi)) {
    rsesoi <- abs(rsesoi)
    fsesoi <- atanh(rsesoi)
    object$rsesoi <- rsesoi

    for (nm in c("wald", "fisherwald", "boot")) {
      if (!is.null(object[[nm]])) {

        out <- object[[nm]]

        if (nm == "wald") {
          fobs <- atanh(clip(out$obs))
          lses <- out$obs - tanh(fobs + fsesoi)
          uses <- out$obs - tanh(fobs - fsesoi)
        } else {
          lses <- -fsesoi
          uses <-  fsesoi
        }

        out$lowersesoi <- lses
        out$uppersesoi <- uses

        ## recompute decisions only
        if ("decision" %in% names(out)) {
          out$decision <- mapply(
            decisionCIequivResidual,
            lses, uses, out$cilower, out$ciupper
          )
        }

        if ("decisionbon" %in% names(out)) {
          out$decisionbon <- mapply(
            decisionCIequivResidual,
            lses, uses, out$cilowerbon, out$ciupperbon
          )
        }

        if ("decisionmax" %in% names(out)) {
          out$decisionmax <- mapply(
            decisionCIequivResidual,
            lses, uses, out$cilowermax, out$ciuppermax
          )
        }

        object[[nm]] <- out
      }
    }
  }

  ## -----------------------------
  ## Update cilevel if requested
  ## -----------------------------

  if (!is.null(cilevel)) {

    object$cilevel <- cilevel
    alpha  <- 1 - cilevel
    crit   <- qnorm(1 - alpha/2)

    ## Bonferroni correction needs p
    if (!is.null(object$wald)) {
      p <- nrow(object$wald)
      critbon <- qnorm(1 - alpha / (2*p))
    } else if (!is.null(object$fisherwald)) {
      p <- nrow(object$fisherwald)
      critbon <- qnorm(1 - alpha / (2*p))
    } else if (!is.null(object$boot)) {
      p <- nrow(object$boot)
      critbon <- qnorm(1 - alpha / (2*p))
    }

    ## ---- Wald ----
    if (!is.null(object$wald)) {
      out <- object$wald

      lci  <- out$diffr - crit * out$se
      uci  <- out$diffr + crit * out$se

      out$cilower <- lci
      out$ciupper <- uci

      if ("decision" %in% names(out)) {
        out$decision <- mapply(
          decisionCIequivResidual,
          out$lowersesoi, out$uppersesoi,
          lci, uci
        )
      }

      if (object$adjust.method %in% c("bonferroni", "all")) {
        lcibon <- out$diffr - critbon * out$se
        ucibon <- out$diffr + critbon * out$se

        out$cilowerbon <- lcibon
        out$ciupperbon <- ucibon

        if ("decisionbon" %in% names(out)) {
          out$decisionbon <- mapply(
            decisionCIequivResidual,
            out$lowersesoi, out$uppersesoi,
            lcibon, ucibon
          )
        }
      }

      object$wald <- out
    }

    ## ---- Fisher Wald ----
    if (!is.null(object$fisherwald)) {
      out <- object$fisherwald

      difff <- out$difff
      N     <- object$N

      ## need SE recomputation:
      hz <- 1 - ((1 - out$obs^2) / (1 - out$imp^2))
      he <- 1 / (1 - out$imp^2)
      if(length(N) > 1) {
        sed <- sqrt((hz^2 / (N[out$group] - 3)) + (he^2 * out$se^2))
      } else {
        sed <- sqrt((hz^2 / (N - 3)) + (he^2 * out$se^2))
      }

      lci <- difff - crit * sed
      uci <- difff + crit * sed

      out$cilower <- lci
      out$ciupper <- uci

      if ("decision" %in% names(out)) {
        out$decision <- mapply(
          decisionCIequivResidual,
          out$lowersesoi, out$uppersesoi,
          lci, uci
        )
      }

      object$fisherwald <- out
    }

    ## ---- Bootstrap ----
    if (!is.null(object$boot)) {

      out <- object$boot
      theta_hat <- object$bootthetahat
      seboot    <- object$bootse

      ## normal bootstrap
      lci <- theta_hat - crit * seboot
      uci <- theta_hat + crit * seboot

      out$cilower <- lci
      out$ciupper <- uci

      if ("decision" %in% names(out)) {
        out$decision <- mapply(
          decisionCIequivResidual,
          out$lowersesoi, out$uppersesoi,
          lci, uci
        )
      }

      ## max bootstrap
      if (!is.null(object$bootMb)) {
        Mb <- object$bootMb
        c_alpha <- quantile(Mb, probs = cilevel, na.rm = TRUE)

        lci_m <- theta_hat - c_alpha
        uci_m <- theta_hat + c_alpha

        out$cilowermax <- lci_m
        out$ciuppermax <- uci_m

        if ("decisionmax" %in% names(out)) {
          out$decisionmax <- mapply(
            decisionCIequivResidual,
            out$lowersesoi, out$uppersesoi,
            lci_m, uci_m
          )
        }
      }

      object$boot <- out
    }
  }

  object
}

# print.resEquivFit()
# ------------------------------------------------------------------
# Internal print method for resEquivFit objects.

#' @method update resEquivFit
#' @rdname resEquivFit
#' @param object An object returned by \code{\link{resEquivFit}}.
#' @param digits Shown decimal digits.
#' @export
print.resEquivFit <- function(object, digits = 3, ...) {

  cat("Residual Correlation Equivalence Testing\n")
  cat("========================================\n\n")

  cat("SESOI (r):", round(object$rsesoi, digits), "\n")
  cat("CI level  :", round(object$cilevel, digits), "\n")
  cat("Method    :", object$method, "\n")
  cat("Adjustment:", object$adjust.method, "\n\n")

  format_table <- function(df) {
    df_out <- df

    num_cols <- sapply(df_out, is.numeric)
    df_out[num_cols] <- lapply(df_out[num_cols], function(z)
      round(z, digits)
    )

    df_out
  }

  ## -------------------------
  ## Wald
  ## -------------------------
  if (!is.null(object$wald)) {
    cat("---- Wald CI on residual correlation ----\n")
    print(format_table(object$wald))
    cat("\n")
  }

  ## -------------------------
  ## Fisher Wald
  ## -------------------------
  if (!is.null(object$fisherwald)) {
    cat("---- Wald CI on Fisher-z residual ----\n")
    print(format_table(object$fisherwald))
    cat("\n")
  }

  ## -------------------------
  ## Bootstrap
  ## -------------------------
  if (!is.null(object$boot)) {
    cat("---- Bootstrap CI ----\n")
    print(format_table(object$boot))
    cat("\n")
  }

  invisible(object)
}

# summary.resEquivFit()
# ------------------------------------------------------------------
# Internal summary constructor used by resEquivFit().
# Aggregates and ranks residual correlation results from a data-frame
# representation into a structured summaryResEquivFit object for
# downstream printing and global decision reporting. This function
# performs data processing only and does not conduct any estimation
# or statistical inference.

#' @method summary resEquivFit
#' @rdname resEquivFit
#' @param object An object returned by \code{\link{resEquivFit}}.
#' @param method Methods shown in the output.
#' @param adjust.method Adjustment method shown in the output.
#' @param top Shown decimal digits.
#' @export
summary.resEquivFit <- function(
    object,
    method = NULL,
    adjust.method = NULL,
    top = 5,
    ...
) {

  stopifnot(inherits(object, "resEquivFit"))

  countdecision <- function(dec) {
    n_M  <- sum(dec == "M",  na.rm = TRUE)
    n_I  <- sum(dec == "I",  na.rm = TRUE)
    n_U  <- sum(dec == "U",  na.rm = TRUE)
    n_NM <- sum(dec == "NM", na.rm = TRUE)
    c(
      n_M  = n_M,
      n_I  = n_I,
      n_U  = n_U,
      n_NM = n_NM
    )
  }

  object2 <- object[intersect(c("wald", "fisherwald", "boot"), names(object))]
  sumout <- rep(list(NA), length(object2))
  names(sumout) <- names(object2)
  for(i in names(object2)) {
    tempobj <- object2[[i]]
    tempcolnames <- colnames(tempobj)
    if("difff" %in% tempcolnames) {
      misfitmeasures <- "difff"
    } else {
      misfitmeasures <- "diffr"
    }

    namemult <- c("none", "bonferroni", "max")
    namedecmult <- c("decision", "decisionbon", "decisionmax")
    namecilowermult <- c("cilower", "cilowerbon", "cilowermax")
    nameciuppermult <- c("ciupper", "ciupperbon", "ciuppermax")
    tempresult <- list()

    for(j in seq_along(namemult)) {
      if(namedecmult[j] %in% tempcolnames) {
        tempoutput <- list()
        dec <- countdecision(tempobj[[ namedecmult[j] ]])
        tempoutput$n_M <- dec["n_M"]
        tempoutput$n_I <- dec["n_I"]
        tempoutput$n_U <- dec["n_U"]
        tempoutput$n_NM <- dec["n_NM"]
        tempoutput$res_equivalence <- global_localfit_decision(n_M = dec["n_M"],
                                           n_I = dec["n_I"],
                                           n_U = dec["n_U"],
                                           n_NM = dec["n_NM"])
        if(dec["n_M"] > 0) {
          topS <- tempobj[tempobj[,namedecmult[j]] == "M", , drop = FALSE]
          topS <- topS[is.finite(topS[,misfitmeasures]),]
          topS <- topS[order(abs(topS[,misfitmeasures]), decreasing = TRUE),]
          if(nrow(topS) > top) topS <- topS[1:top,]
          tempoutput$top_non_equiv <- topS
        } else {
          tempoutput$top_non_equiv <- NULL
        }

        if(dec["n_I"] > 0) {
          topI <- tempobj[tempobj[,namedecmult[j]] == "I", , drop = FALSE]
          cigap <- pmin(
            abs(topI[,namecilowermult[j]] - topI[,"lowersesoi"]),
            abs(topI[,nameciuppermult[j]] - topI[,"uppersesoi"]),
            abs(topI[,namecilowermult[j]] - topI[,"uppersesoi"]),
            abs(topI[,nameciuppermult[j]] - topI[,"lowersesoi"]),
            na.rm = TRUE
          )
          topI <- topI[is.finite(cigap),]
          topI$ci_gap <- cigap[is.finite(cigap)]
          topI <- topI[order(cigap, decreasing = TRUE),]
          if(nrow(topI) > top) topI <- topI[1:top,]
          tempoutput$top_inconclusive <- topI
        } else {
          tempoutput$top_inconclusive <- NULL
        }
        tempresult[[namemult[j]]] <- tempoutput
      }
    }
    sumout[[i]] <- tempresult
  }

  lastavailmethod <- tail(names(sumout)[lengths(sumout) > 0], 1)
  if(is.null(method)) method <- lastavailmethod
  methodargs <- c("wald", "fisherwald", "boot", "all")
  method       <- match.arg(method, methodargs)

  if(method == "all") {
    availmult <- names(sumout[[lastavailmethod]])
  } else {
    availmult <- names(sumout[[method]])
  }
  if (is.null(adjust.method)) adjust.method <- tail(availmult, 1)
  adjustmethodargs <- c("none", "bonferroni", "max", "all")
  adjust.method <- match.arg(adjust.method, adjustmethodargs)

  out <- list(
    method = method,
    adjustment = adjust.method,
    equivalence = sumout[[method]][[adjust.method]]
  )

  class(out) <- "summaryResEquivFit"
  out
}

# print.summaryResEquivFit()
# ------------------------------------------------------------------
# Internal print method for summaryResEquivFit objects.

#' @export
print.summaryResEquivFit <- function(x, ...) {

  methodargs <- c("wald", "fisherwald", "boot", "all")
  methodnames <- c("Wald CI on residual correlation",
                   "Wald CI on Fisher's residual correlation",
                   "Bootstrap normal CI on Fisher's z residual correlation",
                   "All methods")
  adjustmethodargs <- c("none", "bonferroni", "max", "all")
  adjustmethodnames <- c("No Correction",
                         "Bonferroni Correction",
                         "Bootstrap Max",
                         "All adjustments")

  method_label <- x$method
  method_index <- match(x$method, methodargs)

  if (!is.na(method_index) && methodargs[method_index] != "all") {
    method_label <- methodnames[method_index]
  }

  ## ---- Adjustment label ----
  adjustment_label <- x$adjustment
  adjust_index <- match(x$adjustment, adjustmethodargs)

  if (!is.na(adjust_index) && adjustmethodargs[adjust_index] != "all") {
    adjustment_label <- adjustmethodnames[adjust_index]
  }

  cat("Residual Correlation Equivalence Testing Summary\n")
  cat("------------------------------------------------\n\n")

  cat("Method:", method_label, "\n")
  cat("Multiplicity adjustment:", adjustment_label, "\n\n")

  cat("[1. CI-based Equivalence Testing]\n")
  cat("Severe (M):",  x$equivalence$n_M,  "\n")
  cat("Inconclusive (I):", x$equivalence$n_I, "\n")
  cat("Underpowered (U):", x$equivalence$n_U, "\n")
  cat("Trivial / Equivalent (NM):", x$equivalence$n_NM, "\n\n")

  global_decision <- global_localfit_decision(
    n_M  = x$equivalence$n_M,
    n_I  = x$equivalence$n_I,
    n_U  = x$equivalence$n_U,
    n_NM = x$equivalence$n_NM
  )

  cat("Global Residual Correlation Decision:", global_decision, "\n\n")

  if (x$equivalence$n_M > 0) {
    cat("1.1 Top Severe Residual Correlations\n")
    cat("(ranked by |residual| / SESOI width):\n")
    print(
      x$equivalence$top_non_equiv[
        , c("row","col","obs","imp","diffr"),
        drop = FALSE
      ]
    )
    cat("\n")
  }

  posci <- which(x$adjustment == c("none", "bonferroni", "max"))
  cilowname <- c("cilower", "cilowerbon", "cilowermax")
  ciupname <- c("ciupper", "ciupperbon", "ciuppermax")
  if (x$equivalence$n_I > 0) {
    cat("1.2 Top CI-inconclusive Residual Correlations\n")
    cat("(ranked by distance to equivalence bounds; larger = needs narrower CI):\n")
    print(
      x$equivalence$top_inconclusive[
        , c("row","col",cilowname[posci],ciupname[posci],"lowersesoi","uppersesoi","ci_gap"),
        drop = FALSE
      ]
    )
    cat("\n")
  }

  invisible(x)
}

#' Residual Equivalence Compensatory-Effect Check for Standardized Parameters
#'
#' Performs a residual-correlation compensatory-effect diagnostic to assess
#' whether standardized population parameters define a valid population
#' covariance matrix and whether imposed local misspecifications can produce
#' residual correlations classified as trivial relative to a smallest residual
#' correlation of interest (SESOI).
#'
#' The compensatory effect is evaluated by constructing implied population
#' models under targeted standardized parameter perturbations and examining
#' resulting residual-correlation equivalence decisions. If residual
#' correlations are classified as trivially misspecified under perturbations
#' larger than the SESOI (e.g., 125\% of the SESOI), the compensatory effect
#' is classified as \code{"PRONOUNCED"}.
#'
#' This function operates on standardized parameters and currently supports
#' recursive SEMs with continuous indicators only.
#'
#' @param lavaanObj A fitted \code{lavaan} object representing the target model.
#' @param maxRelEffect A scalar greater than 1 specifying the relative
#'   magnitude of imposed misspecification to be evaluated. The default value
#'   of 1.25 indicates that perturbations producing residual-correlation
#'   changes of 125\% of the SESOI are evaluated.
#' @param rsesoi Smallest residual correlation of interest. Residual
#'   correlations smaller than this value in absolute magnitude are treated
#'   as trivially misspecified. Default is 0.10.
#' @param cilevel Confidence level used in residual-correlation equivalence
#'   testing. Default is 0.90.
#' @param crossLoad Logical. If \code{TRUE}, evaluates imposed cross-loading
#'   misspecifications.
#' @param extraFactor Logical. If \code{TRUE}, evaluates imposed extra-factor
#'   misspecifications using four indicators, with two indicators selected
#'   from each of two factors.
#' @param maxCrossLoadMis Maximum number of cross-loading misspecifications
#'   to evaluate. If the number of candidate misspecifications exceeds this
#'   value, a random subset is selected.
#' @param maxExtraFactorMis Maximum number of extra-factor misspecifications
#'   to evaluate. If the number of candidate misspecifications exceeds this
#'   value, a random subset is selected.
#' @param seed Optional random seed used when subsampling candidate
#'   misspecifications. If \code{NULL}, the random seed is not set.
#' @param verbose Logical. If \code{TRUE}, prints progress messages during
#'   the compensatory-effect search.
#'
#' @details
#' The procedure first verifies whether the standardized parameter values
#' imply a positive definite population covariance matrix. It then evaluates
#' residual-correlation equivalence decisions under imposed misspecifications
#' by constructing population covariance matrices with added cross-loadings
#' or extra-factor loadings, refitting the target model, and applying
#' \code{\link{resEquivFit}}.
#'
#' Cross-loading misspecifications are evaluated in both positive and negative
#' directions. Extra-factor misspecifications are evaluated using the loading
#' pattern that the magnitude of loadings is equal and the signs of loadings
#' are the same as the primary loadings.
#'
#' If at least one imposed misspecification produces a residual-correlation
#' decision of \code{"NM"}, the compensatory effect is labeled
#' \code{"PRONOUNCED"}. Otherwise, it is labeled
#' \code{"NOT PRONOUNCED"}.
#'
#' Models with categorical indicators, formative indicators, mean structures,
#' multiple-group structures, or fewer than two latent variables are not
#' supported. Note that the residual-based equivalence test can evaluate
#' misspecification based on error correlations well. Therefore, the checking
#' function is not designed for one factor because it should perform well.
#'
#' @return An object of class \code{"resEquivCheckStd"} containing:
#' \itemize{
#'   \item \code{feasible}: Logical indicator of whether a valid standardized
#'     population model exists.
#'   \item \code{any_NM}: Logical indicator of whether any tested
#'     misspecification produced a trivially misspecified residual-correlation
#'     decision.
#'   \item \code{compensatory}: Character string summarizing the presence of
#'     the compensatory effect, such as \code{"NOT PRONOUNCED"},
#'     \code{"PRONOUNCED"}, or \code{"NOT APPLICABLE"}.
#'   \item \code{M_crossload}: Data frame summarizing tested cross-loading
#'     perturbations, imposed values, and residual-equivalence decisions.
#'   \item \code{M_extrafactor}: Data frame summarizing tested extra-factor
#'     perturbations, imposed values, and residual-equivalence decisions.
#' }
#'
#' @importFrom lavaan lavaan lavNames lavInspect modificationIndices standardizedSolution
#'
#' @seealso \code{\link{resEquivFit}}
#'
#' @examples
#' library(lavaan)
#'
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#' fit <- cfa(HS.model, data = HolzingerSwineford1939)
#' \donttest{
#' resEquivCheck(fit)
#' }
#'
#' @export
resEquivCheck <- function(lavaanObj,
                          maxRelEffect = 1.25,
                          rsesoi = 0.10,
                          cilevel = 0.90,
                          crossLoad = TRUE,
                          extraFactor = TRUE,
                          maxCrossLoadMis = 100,
                          maxExtraFactorMis = 100,
                          seed = NULL,
                          verbose = FALSE) {

  if(!is.null(seed)) set.seed(seed)
  # Function to check from decision results from each residual correlation
  global_localfit_decision_vec <- function(vec) {
    if (any(vec == "U", na.rm = TRUE)) {
      return("U")
    }
    if (any(vec == "M", na.rm = TRUE)) {
      return("M")
    }
    if (any(vec == "I", na.rm = TRUE)) {
      return("I")
    }
    if (any(vec == "NM", na.rm = TRUE)) {
      return("NM")
    }
    return(NA_character_)
  }

  # Tolerance for optim
  tol <- 1e-10

  # Return errors obtained from nonpositive definite model-implied covariance
  .new_resEquivCheckStd_infeasible <- function(reason, feasible = FALSE) {
    out <- list(
      feasible = feasible,
      any_NM = NA,
      compensatory = "NOT APPLICABLE",
      reason = reason,
      M_crossload = NULL,
      M_extrafactor = NULL
    )
    class(out) <- "resEquivCheckStd"
    out
  }

  if(!crossLoad & !extraFactor) {
    stop("Either crossLoad or extraFactor must be TRUE.")
  }
  if (maxRelEffect <= 1) {
    stop("maxRelEffect must be greater than 1.")
  }
  lvnames <- lavaan::lavNames(lavaanObj, type="lv")
  if(length(lvnames) < 2) {
    stop("This function supports a model with at least two latent variables.")
  }

  miout <- lavaan::modificationIndices(lavaanObj)

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
  if (is.null(stdbeta)) {
    stdbeta <- matrix(
      0,
      nrow(corpsi),
      ncol(corpsi),
      dimnames = dimnames(corpsi)
    )
  }

  # ---- Existence check: do standardized parameters define ANY PD population Sigma? ----
  residVarPsi0 <- findFactorResidualVar(
    beta = stdbeta,
    corPsi = corpsi,
    totalVarPsi = rep(1, nrow(stdbeta))
  )
  if (any(!is.finite(residVarPsi0)) || any(residVarPsi0 < -tol)) {
    return(.new_resEquivCheckStd_infeasible("std_params_not_generative"))
  }
  residVarPsi0[abs(residVarPsi0) < tol] <- 0

  Phi <- findFactorTotalCov(
    beta = stdbeta,
    corPsi = corpsi,
    errorVarPsi = residVarPsi0
  )

  Sigma_y_noTheta <- lambda %*% Phi %*% t(lambda)

  residVarTheta0 <- findIndResidualVar(
    lambda = lambda,
    totalFactorCov = Phi,
    totalVarTheta = rep(1, nrow(lambda))
  )
  if (any(!is.finite(residVarTheta0)) || any(residVarTheta0 < -tol)) {
    return(.new_resEquivCheckStd_infeasible("std_params_not_generative"))
  }
  residVarTheta0[abs(residVarTheta0) < tol] <- 0

  Theta <- cor2cov_safe(cortheta, sqrt(residVarTheta0))
  Sigma_implied <- Sigma_y_noTheta + Theta

  if (!isPD(Sigma_implied)) {
    return(.new_resEquivCheckStd_infeasible("std_params_not_generative"))
  }

  # Get the position of lower triangular elements
  idx <- lower.tri(Sigma_implied, diag=FALSE)

  # Perturbation in cross-loadings
  mioutcrossload <- miout[miout$op == "=~",]
  if(crossLoad && nrow(mioutcrossload) > 0) {
    if(nrow(mioutcrossload) > maxCrossLoadMis) {
      mioutcrossload <- mioutcrossload[sort(sample(1:nrow(mioutcrossload), maxExtraFactorMis)),]
    }
    resultcrossload <- matrix(NA, nrow(mioutcrossload), 2)
    testeffectcrossload <- matrix(NA_real_, nrow(mioutcrossload), 2)

    for (i in 1:nrow(mioutcrossload)) {
      if(verbose) cat("Check cross-loading misspecification", i, "/", nrow(mioutcrossload), "\n")
      row <- mioutcrossload[i,]
      k <- maxRelEffect

      # Positive Side
      FUN <- function(x) {
        lambda2 <- lambda
        lambda2[row$rhs, row$lhs] <- x
        Sigma_y_noTheta2 <- lambda2 %*% Phi %*% t(lambda2)
        implied <- cov2cor(Sigma_implied)
        obs <- cov2cor(Sigma_y_noTheta2 + Theta)
        fimplied <- atanh(implied[idx])
        fobs <- atanh(obs[idx])
        (max(abs(fimplied - fobs)) - k*rsesoi)^2
      }
      opt <- optimize(FUN, interval = c(0, 1))
      lambda2 <- lambda
      lambda2[row$rhs, row$lhs] <- opt$minimum
      testeffectcrossload[i, 1] <- opt$minimum
      Sigma_y_noTheta2 <- lambda2 %*% Phi %*% t(lambda2)
      Sigma_Implied2 <- Sigma_y_noTheta2 + Theta

      # Only proceed if Sigma is PD and the model is estimable under that population
      if (isPD(Sigma_Implied2)) {
        tempout <- tryCatch(
          suppressWarnings(lavaan::lavaan(lavaanObj,
                                          sample.cov  = Sigma_Implied2,
                                          sample.nobs = 1000000L,
                                          std.lv      = TRUE)),
          error = function(e) NULL
        )
        if (is.null(tempout)) next

        if (suppressWarnings(lavCheckAdmissibleFit(tempout))) {
          tempResTest <- resEquivFit(tempout,
                                    rsesoi = rsesoi,
                                    method = "wald",
                                    adjust.method = "none",
                                    cilevel = cilevel)
          resultcrossload[i, 1] <- global_localfit_decision_vec(tempResTest$wald$decision)
        }
      }

      # Negative Side
      optneg <- optimize(FUN, interval = c(-1, 0))
      lambda2 <- lambda
      lambda2[row$rhs, row$lhs] <- optneg$minimum
      testeffectcrossload[i, 2] <- optneg$minimum
      Sigma_y_noTheta2 <- lambda2 %*% Phi %*% t(lambda2)
      Sigma_Implied2 <- Sigma_y_noTheta2 + Theta

      # Only proceed if Sigma is PD and the model is estimable under that population
      if (isPD(Sigma_Implied2)) {
        tempout <- tryCatch(
          suppressWarnings(lavaan::lavaan(lavaanObj,
                                          sample.cov  = Sigma_Implied2,
                                          sample.nobs = 1000000L,
                                          std.lv      = TRUE)),
          error = function(e) NULL
        )
        if (is.null(tempout)) next

        if (lavCheckAdmissibleFit(tempout)) {
          tempResTest <- resEquivFit(tempout,
                                     rsesoi = rsesoi,
                                     method = "wald",
                                     adjust.method = "none",
                                     cilevel = cilevel)
          resultcrossload[i, 2] <- global_localfit_decision_vec(tempResTest$wald$decision)
        }
      }
    }
  } else {
    resultcrossload <- NULL
    testeffectcrossload <- NULL
  }

  # Add extra factor
  if(extraFactor) {
    stdtable <- lavaan::standardizedSolution(lavaanObj)
    stdtable <- stdtable[stdtable[,"op"] == "=~", ]

    # 1. remove items with dual loadings
    item_count <- table(stdtable$rhs)
    single_items <- names(item_count[item_count == 1])
    stdtable <- stdtable[stdtable$rhs %in% single_items, ]

    # 2. remove factors with less than three indicators
    fac_count <- table(stdtable$lhs)
    valid_factors <- names(fac_count[fac_count >= 3])
    stdtable <- stdtable[stdtable$lhs %in% valid_factors, ]

    # 3. Proceed if remaining factors are at least two
    if (length(unique(stdtable$lhs)) > 1) {

      # split indicators by factor
      items_by_factor <- split(stdtable$rhs, stdtable$lhs)

      # choose two factors, then choose two items from each factor
      factor_pairs <- combn(names(items_by_factor), 2, simplify = FALSE)

      extra_factor_table <- do.call(rbind, lapply(factor_pairs, function(fp) {

        f1 <- fp[1]
        f2 <- fp[2]

        comb1 <- combn(items_by_factor[[f1]], 2, simplify = FALSE)
        comb2 <- combn(items_by_factor[[f2]], 2, simplify = FALSE)

        out <- expand.grid(
          i = seq_along(comb1),
          j = seq_along(comb2)
        )

        out2 <- matrix(NA_character_, nrow(out), 4)
        for(i in seq_len(nrow(out))){
          out2[i, 1:2] <- comb1[[out[i,1]]]
          out2[i, 3:4] <- comb2[[out[i,2]]]
        }
        out2
      }))

      # table with four possible indicators
      extra_factor_table <- as.data.frame(extra_factor_table)
      colnames(extra_factor_table) <- paste0("indicator", 1:4)
      rownames(extra_factor_table) <- NULL

      # if the number of rows is over maxExtraFactorMis, randomly select rows
      if(nrow(extra_factor_table) > maxExtraFactorMis) {
        extra_factor_table <- extra_factor_table[sort(sample(1:nrow(extra_factor_table), maxExtraFactorMis)),]
      }

      resultextrafactor <- matrix(NA, nrow(extra_factor_table), 1)
      testeffectextrafactor <- matrix(NA, nrow(extra_factor_table), 1)

      for (i in 1:nrow(extra_factor_table)) {
        if(verbose) cat("Check extra-factor misspecification", i, "/", nrow(extra_factor_table), "\n")
        row <- unlist(extra_factor_table[i,])
        k <- maxRelEffect
        posextrafac <- ncol(lambda) + 1
        signvec <- ifelse(rowSums(lambda[row, , drop = FALSE]) >= 0, 1, -1)

        FUN <- function(x) {
          lambda2 <- cbind(lambda, 0)
          lambda2[row, posextrafac] <- signvec*x
          Phi2 <- diag(posextrafac)
          Phi2[1:ncol(lambda), 1:ncol(lambda)] <- Phi
          Sigma_y_noTheta2 <- lambda2 %*% Phi2 %*% t(lambda2)
          implied <- cov2cor(Sigma_implied)
          obs <- cov2cor(Sigma_y_noTheta2 + Theta)
          fimplied <- atanh(implied[idx])
          fobs <- atanh(obs[idx])
          (max(abs(fimplied - fobs)) - k*rsesoi)^2
        }
        opt <- optimize(FUN, interval = c(0, 1))
        lambda2 <- cbind(lambda, 0)
        lambda2[row, posextrafac] <- signvec*opt$minimum
        testeffectextrafactor[i, 1] <- opt$minimum
        Phi2 <- diag(posextrafac)
        Phi2[1:ncol(lambda), 1:ncol(lambda)] <- Phi
        Sigma_y_noTheta2 <- lambda2 %*% Phi2 %*% t(lambda2)
        Sigma_Implied2 <- Sigma_y_noTheta2 + Theta

        # Only proceed if Sigma is PD and the model is estimable under that population
        if (isPD(Sigma_Implied2)) {
          tempout <- tryCatch(
            suppressWarnings(lavaan::lavaan(lavaanObj,
                                            sample.cov  = Sigma_Implied2,
                                            sample.nobs = 1000000L,
                                            std.lv      = TRUE)),
            error = function(e) NULL
          )
          if (is.null(tempout)) next

          if (suppressWarnings(lavCheckAdmissibleFit(tempout))) {
            tempResTest <- resEquivFit(tempout,
                                       rsesoi = rsesoi,
                                       method = "wald",
                                       adjust.method = "none",
                                       cilevel = cilevel)
            resultextrafactor[i, 1] <- global_localfit_decision_vec(tempResTest$wald$decision)
          }
        }
      }
    } else {
      resultextrafactor <- NULL
      testeffectextrafactor <- NULL
    }
  } else {
    resultextrafactor <- NULL
    testeffectextrafactor <- NULL
  }

  # Check for any NM. Note that "I" and "U" are ignored because the tested model is
  # based on sample size of 1,000,000. Thus, expect all perturbations return "M"
  all_results <- c(as.vector(resultcrossload), as.vector(resultextrafactor))
  all_results <- all_results[!is.na(all_results)]
  any_NM <- if (length(all_results) == 0) {
    return(.new_resEquivCheckStd_infeasible("no_illegal_imposed_misspecifications",
                                            feasible=TRUE))
  } else {
    any(all_results == "NM", na.rm = TRUE)
  }
  if(!is.null(resultcrossload)) {
    colnames(resultcrossload) <- c("decision_pos", "decision_neg")
    colnames(testeffectcrossload) <- c("imposed_val_pos", "imposed_val_neg")
    M_crossload <- data.frame(mioutcrossload[,1:3], testeffectcrossload, resultcrossload)
  } else {
    M_crossload <- NULL
  }
  if(!is.null(resultextrafactor)) {
    colnames(resultextrafactor) <- "decision"
    colnames(testeffectextrafactor) <- "imposed"
    M_extrafactor <- data.frame(extra_factor_table, testeffectextrafactor, resultextrafactor)
  } else {
    M_extrafactor <- NULL
  }
  feasible <- TRUE
  compensatory <- if (!feasible) {
    "NOT APPLICABLE"
  } else if (any_NM) {
    "PRONOUNCED"
  } else {
    "NOT PRONOUNCED"
  }
  out <- list(
    feasible = feasible,
    any_NM = any_NM,
    compensatory = compensatory,
    M_crossload = M_crossload,
    M_extrafactor = M_extrafactor
  )
  class(out) <- "resEquivCheckStd"
  return(out)
}

# print.resEquivCheckStd()
# ------------------------------------------------------------------
# Internal print method for resEquivCheckStd objects.

#' @method print resEquivCheckStd
#' @rdname resEquivCheckStd
#' @param object An object returned by \code{\link{resEquivCheck}}.
#' @param max.print The number of imposed misspecification shown in the output
#' @export
print.resEquivCheckStd <- function(x, ..., max.print = 6L) {
  cat("Residual Equivalence Compensatory Effect (Standardized Parameters)\n")
  cat("----------------------------------------------------\n")

  cat("Feasible standardized population:", x$feasible, "\n")

  if(x$feasible) {
    cat("Any residual correlation less than SESOI:", x$any_NM, "\n")
    cat("Compensatory Effect:", x$compensatory, "\n\n")
  }

  if (x$compensatory == "PRONOUNCED") {
    cat("Residual correlation below SESOI under tested perturbations (summary):\n")
    printmaxprint <- FALSE
    tempcl <- x$M_crossload
    if(!is.null(tempcl)) {
      tempclpos <- tempcl[!is.na(tempcl[,"decision_pos"]) & tempcl[,"decision_pos"] == "NM",]
      tempclneg <- tempcl[!is.na(tempcl[,"decision_neg"]) & tempcl[,"decision_neg"] == "NM",]
      if(nrow(tempclpos) > 0) {
        cat("-- Positive cross-loading perturbation:\n")
        tempclpos_print <- tempclpos[
          , !colnames(tempclpos) %in% c("decision_neg", "imposed_val_neg"),
          drop = FALSE
        ]
        print(head(tempclpos_print, n = max.print))
        if(nrow(tempclpos) > max.print) printmaxprint <- TRUE
      }
      if(nrow(tempclneg) > 0) {
        cat("-- Negative cross-loading perturbation:\n")
        tempclneg_print <- tempclneg[
          , !colnames(tempclneg) %in% c("decision_pos", "imposed_val_pos"),
          drop = FALSE
        ]
        print(head(tempclneg_print, n = max.print))
        if(nrow(tempclneg) > max.print) printmaxprint <- TRUE
      }
    }
    tempef <- x$M_extrafactor
    if(!is.null(tempef) && nrow(tempef) > 0) {

      cat("-- Extra-factor-loading perturbation:\n")
      print(head(tempef, n = max.print))
      if(nrow(tempef) > max.print) printmaxprint <- TRUE
    }
    if(printmaxprint) {
      cat("Please adjust max.print in print(..., max.print=6L) to see more NM results.")
    }
  } else if (x$compensatory == "NOT PRONOUNCED") {
    cat("No residual correlation below the SESOI under tested perturbations.\n")
  } else {
    if (x$reason == "no_illegal_imposed_misspecifications") {
      cat("No imposed misspecification provides valid tested results.\n")
    } else {
      cat("Standardized parameters do not define a valid population model.\n")
    }
  }

  invisible(x)
}


### Helper for Residual Equivalence Test

decisionCIequivResidual <- function(lsesoi, usesoi, lci, uci) {
  if (is.na(lci) || is.na(uci)) return(NA_character_)
  ciwidth    <- uci - lci
  sesoiwidth <- usesoi - lsesoi
  if (ciwidth > sesoiwidth) {
    "U"
  } else if (lci > usesoi || uci < lsesoi) {
    "M"
  } else if (lci > lsesoi && uci < usesoi) {
    "NM"
  } else {
    "I"
  }
}

clip <- function(x, eps = 1e-12) pmin(pmax(x, -1 + eps), 1 - eps)


