## Title: Compute more fit indices
## Authors: Terrence Jorgensen <TJorgensen314@gmail.com>
##          Sunthud Pornprasertmanit <psunthud@ku.edu>,
##          Aaron Boulton <aboulton@ku.edu>,
##          Ruben Arslan <rubenarslan@gmail.com>
## Last updated: 2 June 2018
## Description: Calculations for promising alternative fit indices
##----------------------------------------------------------------------------



#' Calculate more fit indices
#'
#' Calculate more fit indices that are not already provided in lavaan.
#'
#' Gamma Hat (gammaHat; West, Taylor, & Wu, 2012) is a global fit index which
#' can be computed (assuming equal number of indicators across groups) by
#'
#' \deqn{ gammaHat =\frac{p}{p + 2 \times \frac{\chi^{2}_{k} - df_{k}}{N}} ,}
#'
#' where \eqn{p} is the number of variables in the model, \eqn{\chi^{2}_{k}} is
#' the \eqn{\chi^2} test statistic value of the target model, \eqn{df_{k}} is
#' the degree of freedom when fitting the target model, and \eqn{N} is the
#' sample size (or sample size minus the number of groups if \code{mimic} is
#' set to \code{"EQS"}).
#'
#' Adjusted Gamma Hat (adjGammaHat; West, Taylor, & Wu, 2012) is a global fit
#' index which can be computed by
#'
#' \deqn{ adjGammaHat = \left(1 - \frac{K \times p \times (p + 1)}{2 \times
#' df_{k}} \right) \times \left( 1 - gammaHat \right) ,}
#'
#' where \eqn{K} is the number of groups (please refer to Dudgeon, 2004 for the
#' multiple-group adjustment for agfi*).
#'
#' Corrected Akaike Information Criterion (aic.smallN; Burnham & Anderson,
#' 2003) is a corrected version of AIC for small sample size, often abbreviated
#' AICc:
#'
#' \deqn{ aic.smallN = AIC + \frac{2k(k + 1)}{N - k - 1},}
#'
#' where \eqn{AIC} is the original AIC: \eqn{-2 \times LL + 2k} (where \eqn{k}
#' = the number of estimated parameters in the target model). Note that AICc is
#' a small-sample correction derived for univariate regression models, so it is
#' probably \emph{not} appropriate for comparing SEMs.
#'
#' Corrected Bayesian Information Criterion (bic.priorN; Kuha, 2004) is similar
#' to BIC but explicitly specifying the sample size on which the prior is based
#' (\eqn{N_{prior}}).
#'
#' \deqn{ bic.priorN = f + k\log{(1 + N/N_{prior})},}
#'
#' Stochastic information criterion (SIC; Preacher, 2006) is similar to AIC or
#' BIC. This index will account for model complexity in the model's function
#' form, in addition to the number of free parameters. This index will be
#' provided only when the \eqn{\chi^2} value is not scaled. The SIC can be
#' computed by
#'
#' \deqn{ sic = \frac{1}{2}\left(f - \log{\det{I(\hat{\theta})}}\right),}
#'
#' where \eqn{I(\hat{\theta})} is the information matrix of the parameters.
#'
#' Hannan-Quinn Information Criterion (hqc; Hannan & Quinn, 1979) is used for
#' model selection similar to AIC or BIC.
#'
#' \deqn{ hqc = f + 2k\log{(\log{N})},}
#'
#' Note that if Satorra--Bentler or Yuan--Bentler's method is used, the fit
#' indices using the scaled \eqn{\chi^2} values are also provided.
#'
#' See \code{\link{nullRMSEA}} for the further details of the computation of
#' RMSEA of the null model.
#'
#'
#' @importFrom lavaan lavInspect
#'
#' @param object The lavaan model object provided after running the \code{cfa},
#' \code{sem}, \code{growth}, or \code{lavaan} functions.
#' @param fit.measures Additional fit measures to be calculated. All additional
#' fit measures are calculated by default
#' @param nPrior The sample size on which prior is based. This argument is used
#' to compute BIC*.
#' @return \enumerate{
#'  \item \code{gammaHat}: Gamma Hat
#'  \item \code{adjGammaHat}: Adjusted Gamma Hat
#'  \item \code{baseline.rmsea}: RMSEA of the Baseline (Null) Model
#'  \item \code{aic.smallN}: Corrected (for small sample size) Akaike Information Criterion
#'  \item \code{bic.priorN}: Bayesian Information Criterion with specified prior sample size
#'  \item \code{sic}: Stochastic Information Criterion
#'  \item \code{hqc}: Hannan-Quinn Information Criterion
#'  \item \code{gammaHat.scaled}: Gamma Hat using scaled \eqn{\chi^2}
#'  \item \code{adjGammaHat.scaled}: Adjusted Gamma Hat using scaled \eqn{\chi^2}
#'  \item \code{baseline.rmsea.scaled}: RMSEA of the Baseline (Null) Model using scaled \eqn{\chi^2}
#' }
#' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#'
#' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
#'
#' Aaron Boulton (University of North Carolina, Chapel Hill; \email{aboulton@@email.unc.edu})
#'
#' Ruben Arslan (Humboldt-University of Berlin, \email{rubenarslan@@gmail.com})
#'
#' Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
#'
#' @seealso \itemize{ \item \code{\link{miPowerFit}} For the modification
#' indices and their power approach for model fit evaluation \item
#' \code{\link{nullRMSEA}} For RMSEA of the null model }
#'
#' @references Burnham, K., & Anderson, D. (2003). \emph{Model selection and
#' multimodel inference: A practical--theoretic approach}. New York, NY:
#' Springer--Verlag.
#'
#' Dudgeon, P. (2004). A note on extending Steiger's (1998) multiple sample
#' RMSEA adjustment to other noncentrality parameter-based statistic.
#' \emph{Structural Equation Modeling, 11}(3), 305--319.
#' doi:10.1207/s15328007sem1103_1
#'
#' Kuha, J. (2004). AIC and BIC: Comparisons of assumptions and performance.
#' \emph{Sociological Methods Research, 33}(2), 188--229.
#' doi:10.1177/0049124103262065
#'
#' Preacher, K. J. (2006). Quantifying parsimony in structural equation
#' modeling. \emph{Multivariate Behavioral Research, 43}(3), 227-259.
#' doi:10.1207/s15327906mbr4103_1
#'
#' West, S. G., Taylor, A. B., & Wu, W. (2012). Model fit and model selection
#' in structural equation modeling. In R. H. Hoyle (Ed.), \emph{Handbook of
#' Structural Equation Modeling} (pp. 209--231). New York, NY: Guilford.
#' @examples
#'
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#'
#' fit <- cfa(HS.model, data = HolzingerSwineford1939)
#' moreFitIndices(fit)
#'
#' fit2 <- cfa(HS.model, data = HolzingerSwineford1939, estimator = "mlr")
#' moreFitIndices(fit2)
#'
#' @export
moreFitIndices <- function(object, fit.measures = "all", nPrior = 1) {
  ## check for validity of user-specified "fit.measures" argument
  fit.choices <- c("gammaHat","adjGammaHat","baseline.rmsea",
                   "gammaHat.scaled","adjGammaHat.scaled","baseline.rmsea.scaled",
                   "aic.smallN","bic.priorN","hqc","sic")
  flags <- setdiff(fit.measures, c("all", fit.choices))
  if (length(flags)) stop(paste("Argument 'fit.measures' includes invalid options:",
                                paste(flags, collapse = ", "),
                                "Please choose 'all' or among the following:",
                                paste(fit.choices, collapse = ", "), sep = "\n"))
  if ("all" %in% fit.measures) fit.measures <- fit.choices

  # Extract fit indices information from lavaan object
  fit <- lavInspect(object, "fit")
  # Get the number of variable
  p <- length(lavaan::lavNames(object, type = "ov", group = 1))
  # Get the number of parameters
  nParam <- fit["npar"]

  # Find the number of groups
  ngroup <- lavInspect(object, "ngroups")
  # Get number of observations
  N <- n <- lavInspect(object, "ntotal")
  if (lavInspect(object, "options")$mimic == "EQS") n <- n - ngroup

  # Calculate -2*log(likelihood)
  f <- -2 * fit["logl"]

  # Compute fit indices
  result <- list()
  if (length(grep("gamma", fit.measures, ignore.case = TRUE))) {
    gammaHat <- p / (p + 2 * ((fit["chisq"] - fit["df"]) / n))
    adjGammaHat <- 1 - (((ngroup * p * (p + 1)) / 2) / fit["df"]) * (1 - gammaHat)
    result["gammaHat"] <- gammaHat
    result["adjGammaHat"] <- adjGammaHat
    if (!lavInspect(object, "options")$test %in% c("standard","bollen.stine")) {
      gammaHatScaled <- p / (p + 2 * ((fit["chisq.scaled"] - fit["df.scaled"]) / n))
      adjGammaHatScaled <- 1 - (((ngroup * p * (p + 1)) / 2) / fit["df.scaled"]) * (1 - gammaHatScaled)
      result["gammaHat.scaled"] <- gammaHatScaled
      result["adjGammaHat.scaled"] <- adjGammaHatScaled
    }
  }
  if (length(grep("rmsea", fit.measures))) {
    result["baseline.rmsea"] <- nullRMSEA(object, silent = TRUE)
    if (!lavInspect(object, "options")$test %in% c("standard","bollen.stine")) {
      result["baseline.rmsea.scaled"] <- nullRMSEA(object, scaled = TRUE, silent = TRUE)
    }
  }
  if (!is.na(f)) {
    if ("aic.smallN" %in% fit.measures) {
      warning('AICc (aic.smallN) was developed for univariate linear models.',
              ' It is probably not appropriate to use AICc to compare SEMs.')
      result["aic.smallN"] <- fit[["aic"]] + (2 * nParam * (nParam + 1)) / (N - nParam - 1)
    }
    if ("bic.priorN" %in% fit.measures) {
      result["bic.priorN"] <- f + log(1 + N/nPrior) * nParam
    }
    if ("hqc" %in% fit.measures) result["hqc"] <- f + 2 * log(log(N)) * nParam
    if ("sic" %in% fit.measures) result["sic"] <- sic(f, object)
  }
  class(result) <- c("lavaan.vector","numeric")
  unlist(result[fit.measures])
}



#' Calculate the RMSEA of the null model
#'
#' Calculate the RMSEA of the null (baseline) model
#'
#' RMSEA of the null model is calculated similar to the formula provided in the
#' \code{lavaan} package. The standard formula of RMSEA is
#'
#' \deqn{ RMSEA =\sqrt{\frac{\chi^2}{N \times df} - \frac{1}{N}} \times
#' \sqrt{G} }
#'
#' where \eqn{\chi^2} is the chi-square test statistic value of the target
#' model, \eqn{N} is the total sample size, \eqn{df} is the degree of freedom
#' of the hypothesized model, \eqn{G} is the number of groups. Kenny proposed
#' in his website that
#'
#' "A reasonable rule of thumb is to examine the RMSEA for the null model and
#' make sure that is no smaller than 0.158. An RMSEA for the model of 0.05 and
#' a TLI of .90, implies that the RMSEA of the null model is 0.158.  If the
#' RMSEA for the null model is less than 0.158, an incremental measure of fit
#' may not be that informative."
#'
#' See also \url{http://davidakenny.net/cm/fit.htm}
#'
#'
#' @importFrom lavaan lavInspect
#'
#' @param object The lavaan model object provided after running the \code{cfa},
#' \code{sem}, \code{growth}, or \code{lavaan} functions.
#' @param scaled If \code{TRUE}, the scaled (or robust, if available) RMSEA
#'   is returned. Ignored if a robust test statistic was not requested.
#' @param silent If \code{TRUE}, do not print anything on the screen.
#'
#' @return A value of RMSEA of the null model (a \code{numeric} vector)
#'   returned invisibly.
#'
#' @author
#' Ruben Arslan (Humboldt-University of Berlin, \email{rubenarslan@@gmail.com})
#'
#' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{miPowerFit}} For the modification indices and their
#'      power approach for model fit evaluation
#'   \item \code{\link{moreFitIndices}} For other fit indices
#' }
#'
#' @references Kenny, D. A., Kaniskan, B., & McCoach, D. B. (2015). The
#' performance of RMSEA in models with small degrees of freedom.
#' \emph{Sociological Methods Research, 44}(3), 486--507.
#' doi:10.1177/0049124114543236
#'
#' @examples
#'
#' HS.model <- ' visual  =~ x1 + x2 + x3
#'               textual =~ x4 + x5 + x6
#'               speed   =~ x7 + x8 + x9 '
#'
#' fit <- cfa(HS.model, data = HolzingerSwineford1939)
#' nullRMSEA(fit)
#'
#' @export
nullRMSEA <- function(object, scaled = FALSE, silent = FALSE) {
  fit <- lavaan::update(object, model = lavaan::lav_partable_independence(object))
	fits <- lavaan::fitMeasures(fit, fit.measures = c("rmsea","rmsea.scaled",
	                                                  "rmsea.robust"))
	if (scaled) {
	  RMSEA <- fits["rmsea.robust"]
	  if (is.na(RMSEA)) RMSEA <- fits["rmsea.scaled"]
	  if (is.na(RMSEA)) RMSEA <- fits["rmsea"]
	} else RMSEA <- fits["rmsea"]

	if (!silent) {
	  cat("The baseline model's RMSEA =", RMSEA, "\n\n")
		if (RMSEA < 0.158 ) {
			cat("CFI, TLI, and other incremental fit indices may not be very",
			    "informative because the baseline model's RMSEA < 0.158",
			    "(Kenny, Kaniskan, & McCoach, 2015). \n")
		}
	}
	invisible(RMSEA)
}



## Stochastic Information Criterion
## f = minimized discrepancy function
## lresults = lavaan sem output object
sic <- function(f, lresults = NULL) {

  E.inv <- lavaan::lavTech(lresults, "inverted.information")
  if (inherits(E.inv, "try-error")) {
    return(as.numeric(NA))
  }
  E <- MASS::ginv(E.inv) * lavaan::nobs(lresults)

  eigvals <- eigen(E, symmetric = TRUE, only.values = TRUE)$values
  # only positive ones
  eigvals <- eigvals[ eigvals > sqrt(.Machine$double.eps)]

  DET <- prod(eigvals)

  ## check singular
  if (DET <= 0) return(NA)

  ## return SIC
  f + log(DET)
}



#' \emph{k}-factor correction for \eqn{chi^2} test statistic
#'
#' Calculate \emph{k}-factor correction for \eqn{chi^2} model-fit test
#' statistic to adjust for small sample size.
#'
#' The \emph{k}-factor correction (Nevitt & Hancock, 2004) is a global fit
#' index which can be computed by:
#'
#' \deqn{ kc = 1 - \frac{2 \times P + 4 \times K + 5}{6 \times N}}
#'
#' where \eqn{N} is the sample size when using normal likelihood, or \eqn{N -
#' 1} when using \code{likelihood = 'wishart'}.
#'
#'
#' @importFrom lavaan lavInspect
#' @importFrom stats pchisq
#'
#' @param fit0 The lavaan model object provided after running the \code{cfa},
#' \code{sem}, \code{growth}, or \code{lavaan} functions.
#' @param fit1 Optional additional \linkS4class{lavaan} model, in which
#' \code{fit0} is nested.  If \code{fit0} has fewer \emph{df} than \code{fit1},
#' the models will be swapped, still on the assumption that they are nested.
#' @param \dots Additional arguments to the \code{\link[lavaan]{lavTestLRT}}
#' function.
#' @return A numeric vector including the unadjusted (naive) chi-squared test
#' statistic, the \emph{k}-factor correction, the corrected test statistic, the
#' \emph{df} for the test, and the \emph{p} value for the test under the null
#' hypothesis that the model fits perfectly (or that the 2 models have
#' equivalent fit).
#' @author Terrence D. Jorgensen (University of Amsterdam;
#' \email{TJorgensen314@@gmail.com})
#' @references Nevitt, J., & Hancock, G. R. (2004). Evaluating small sample
#' approaches for model test statistics in structural equation modeling.
#' \emph{Multivariate Behavioral Research, 39}(3), 439--478.
#' doi:10.1207/S15327906MBR3903_3
#' @examples
#'
#' HS.model <- '
#'     visual  =~ x1 + b1*x2 + x3
#'     textual =~ x4 + b2*x5 + x6
#'     speed   =~ x7 + b3*x8 + x9
#' '
#' fit1 <- cfa(HS.model, data = HolzingerSwineford1939)
#' ## test a single model (implicitly compared to a saturated model)
#' chisqSmallN(fit1)
#'
#' ## fit a more constrained model
#' fit0 <- cfa(HS.model, data = HolzingerSwineford1939, orthogonal = TRUE)
#' ## compare 2 models
#' chisqSmallN(fit1, fit0)
#'
#' @export
chisqSmallN <- function(fit0, fit1 = NULL, ...) {
  ## if there are 2 models, order them by DF
  if (!is.null(fit1)) {
    DF0 <- lavaan::fitMeasures(fit0, "df")
    DF1 <- lavaan::fitMeasures(fit1, "df")
    if (DF0 == DF1) stop("Models have the same degrees of freedom.")
    parent <- which.min(c(DF0, DF1))
    if (parent == 1L) {
      parent <- fit0
      fit0 <- fit1
      fit1 <- parent
    }
    #if (min(c(DF0, DF1)) == 0L) fit1 <- NULL
  }
  ## calculate k-factor correction
  N <- lavInspect(fit0, "ntotal")
  Ng <- lavInspect(fit0, "ngroups")
  if (!lavInspect(fit0, "options")$sample.cov.rescale) N <- N - Ng
  P <- length(lavaan::lavNames(fit0))
  K <- length(lavaan::lavNames(fit0, type = "lv")) # count latent factors
  if (!is.null(fit1)) {
    N1 <- lavInspect(fit1, "ntotal")
    Ng1 <- lavInspect(fit1, "ngroups")
    if (!lavInspect(fit1, "options")$sample.cov.rescale) N1 <- N1 - Ng1
    if (N != N1) stop("Unequal sample sizes")
    if (P != length(lavaan::lavNames(fit1))) stop("Unequal number of variables")
    K <- max(K, length(lavaan::lavNames(fit1, type = "lv")))
  }
  kc <- 1 - ((2*P + 4*K + 5) / (6*N))
  if (is.null(fit1)) {
    scaled <- lavInspect(fit0, "options")$test %in%
      c("satorra.bentler","yuan.bentler","mean.var.adjusted","scaled.shifted")
    chi <- lavaan::fitMeasures(fit0)[[if (scaled) "chisq.scaled" else "chisq"]]
    DF <- lavaan::fitMeasures(fit0)[["df"]]
  } else {
    AOV <- lavaan::lavTestLRT(fit0, fit1, ...)
    chi <- AOV[["Chisq diff"]][2]
    DF <- AOV[["Df diff"]][2]
  }
  out <- c(naive.chisq = chi, `k-factor` = kc, adj.chisq = chi*kc,
           df = DF, pvalue = pchisq(chi*kc, DF, lower.tail = FALSE))
  class(out) <- c("lavaan.vector","numeric")
  out
}


