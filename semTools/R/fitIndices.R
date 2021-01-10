### Title: Compute more fit indices
### Authors: Terrence D. Jorgensen <TJorgensen314@gmail.com>
###          Sunthud Pornprasertmanit <psunthud@ku.edu>,
###          Aaron Boulton <aboulton@ku.edu>,
###          Ruben Arslan <rubenarslan@gmail.com>
### Last updated: 10 January 2021
### Description: Calculations for promising alternative fit indices



##' Calculate more fit indices
##'
##' Calculate more fit indices that are not already provided in lavaan.
##'
##' Gamma Hat (gammaHat; West, Taylor, & Wu, 2012) is a global fit index which
##' can be computed (assuming equal number of indicators across groups) by
##'
##' \deqn{ gammaHat =\frac{p}{p + 2 \times \frac{\chi^{2}_{k} - df_{k}}{N}} ,}
##'
##' where \eqn{p} is the number of variables in the model, \eqn{\chi^{2}_{k}} is
##' the \eqn{\chi^2} test statistic value of the target model, \eqn{df_{k}} is
##' the degree of freedom when fitting the target model, and \eqn{N} is the
##' sample size (or sample size minus the number of groups if \code{mimic} is
##' set to \code{"EQS"}).
##'
##' Adjusted Gamma Hat (adjGammaHat; West, Taylor, & Wu, 2012) is a global fit
##' index which can be computed by
##'
##' \deqn{ adjGammaHat = \left(1 - \frac{K \times p \times (p + 1)}{2 \times
##' df_{k}} \right) \times \left( 1 - gammaHat \right) ,}
##'
##' where \eqn{K} is the number of groups (please refer to Dudgeon, 2004 for the
##' multiple-group adjustment for agfi*).
##'
##' Corrected Akaike Information Criterion (aic.smallN; Burnham & Anderson,
##' 2003) is a corrected version of AIC for small sample size, often abbreviated
##' AICc:
##'
##' \deqn{ aic.smallN = AIC + \frac{2k(k + 1)}{N - k - 1},}
##'
##' where \eqn{AIC} is the original AIC: \eqn{-2 \times LL + 2k} (where \eqn{k}
##' = the number of estimated parameters in the target model). Note that AICc is
##' a small-sample correction derived for univariate regression models, so it is
##' probably \emph{not} appropriate for comparing SEMs.
##'
##' Corrected Bayesian Information Criterion (bic.priorN; Kuha, 2004) is similar
##' to BIC but explicitly specifying the sample size on which the prior is based
##' (\eqn{N_{prior}}).
##'
##' \deqn{ bic.priorN = f + k\log{(1 + N/N_{prior})},}
##'
##' Stochastic information criterion (SIC; Preacher, 2006) is similar to AIC or
##' BIC. This index will account for model complexity in the model's function
##' form, in addition to the number of free parameters. This index will be
##' provided only when the \eqn{\chi^2} value is not scaled. The SIC can be
##' computed by
##'
##' \deqn{ sic = \frac{1}{2}\left(f - \log{\det{I(\hat{\theta})}}\right),}
##'
##' where \eqn{I(\hat{\theta})} is the information matrix of the parameters.
##'
##' Hannan-Quinn Information Criterion (hqc; Hannan & Quinn, 1979) is used for
##' model selection similar to AIC or BIC.
##'
##' \deqn{ hqc = f + 2k\log{(\log{N})},}
##'
##' Note that if Satorra--Bentler or Yuan--Bentler's method is used, the fit
##' indices using the scaled \eqn{\chi^2} values are also provided.
##'
##' See \code{\link{nullRMSEA}} for the further details of the computation of
##' RMSEA of the null model.
##'
##'
##' @importFrom lavaan lavInspect
##'
##' @param object The lavaan model object provided after running the \code{cfa},
##' \code{sem}, \code{growth}, or \code{lavaan} functions.
##' @param fit.measures Additional fit measures to be calculated. All additional
##' fit measures are calculated by default
##' @param nPrior The sample size on which prior is based. This argument is used
##' to compute BIC*.
##' @return \enumerate{
##'  \item \code{gammaHat}: Gamma Hat
##'  \item \code{adjGammaHat}: Adjusted Gamma Hat
##'  \item \code{baseline.rmsea}: RMSEA of the Baseline (Null) Model
##'  \item \code{aic.smallN}: Corrected (for small sample size) Akaike Information Criterion
##'  \item \code{bic.priorN}: Bayesian Information Criterion with specified prior sample size
##'  \item \code{sic}: Stochastic Information Criterion
##'  \item \code{hqc}: Hannan-Quinn Information Criterion
##'  \item \code{gammaHat.scaled}: Gamma Hat using scaled \eqn{\chi^2}
##'  \item \code{adjGammaHat.scaled}: Adjusted Gamma Hat using scaled \eqn{\chi^2}
##'  \item \code{baseline.rmsea.scaled}: RMSEA of the Baseline (Null) Model using scaled \eqn{\chi^2}
##' }
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' Aaron Boulton (University of North Carolina, Chapel Hill; \email{aboulton@@email.unc.edu})
##'
##' Ruben Arslan (Humboldt-University of Berlin, \email{rubenarslan@@gmail.com})
##'
##' Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
##'
##' @seealso \itemize{ \item \code{\link{miPowerFit}} For the modification
##' indices and their power approach for model fit evaluation \item
##' \code{\link{nullRMSEA}} For RMSEA of the null model }
##'
##' @references Burnham, K., & Anderson, D. (2003). \emph{Model selection and
##' multimodel inference: A practical--theoretic approach}. New York, NY:
##' Springer--Verlag.
##'
##' Dudgeon, P. (2004). A note on extending Steiger's (1998) multiple sample
##' RMSEA adjustment to other noncentrality parameter-based statistic.
##' \emph{Structural Equation Modeling, 11}(3), 305--319.
##' \doi{10.1207/s15328007sem1103_1}
##'
##' Kuha, J. (2004). AIC and BIC: Comparisons of assumptions and performance.
##' \emph{Sociological Methods Research, 33}(2), 188--229.
##' \doi{10.1177/0049124103262065}
##'
##' Preacher, K. J. (2006). Quantifying parsimony in structural equation
##' modeling. \emph{Multivariate Behavioral Research, 43}(3), 227-259.
##' \doi{10.1207/s15327906mbr4103_1}
##'
##' West, S. G., Taylor, A. B., & Wu, W. (2012). Model fit and model selection
##' in structural equation modeling. In R. H. Hoyle (Ed.), \emph{Handbook of
##' Structural Equation Modeling} (pp. 209--231). New York, NY: Guilford.
##' @examples
##'
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ x7 + x8 + x9 '
##'
##' fit <- cfa(HS.model, data = HolzingerSwineford1939)
##' moreFitIndices(fit)
##'
##' fit2 <- cfa(HS.model, data = HolzingerSwineford1939, estimator = "mlr")
##' moreFitIndices(fit2)
##'
##' @export
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
    if (any(grepl(pattern = "scaled", x = names(fit)))) {
      gammaHatScaled <- p / (p + 2 * ((fit["chisq.scaled"] - fit["df.scaled"]) / n))
      adjGammaHatScaled <- 1 - (((ngroup * p * (p + 1)) / 2) / fit["df.scaled"]) * (1 - gammaHatScaled)
      result["gammaHat.scaled"] <- gammaHatScaled
      result["adjGammaHat.scaled"] <- adjGammaHatScaled
    }
  }
  if (length(grep("rmsea", fit.measures))) {
    result["baseline.rmsea"] <- nullRMSEA(object, silent = TRUE)
    if (any(grepl(pattern = "scaled", x = names(fit)))) {
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



##' Calculate the RMSEA of the null model
##'
##' Calculate the RMSEA of the null (baseline) model
##'
##' RMSEA of the null model is calculated similar to the formula provided in the
##' \code{lavaan} package. The standard formula of RMSEA is
##'
##' \deqn{ RMSEA =\sqrt{\frac{\chi^2}{N \times df} - \frac{1}{N}} \times
##' \sqrt{G} }
##'
##' where \eqn{\chi^2} is the chi-square test statistic value of the target
##' model, \eqn{N} is the total sample size, \eqn{df} is the degree of freedom
##' of the hypothesized model, \eqn{G} is the number of groups. Kenny proposed
##' in his website that
##'
##' "A reasonable rule of thumb is to examine the RMSEA for the null model and
##' make sure that is no smaller than 0.158. An RMSEA for the model of 0.05 and
##' a TLI of .90, implies that the RMSEA of the null model is 0.158.  If the
##' RMSEA for the null model is less than 0.158, an incremental measure of fit
##' may not be that informative."
##'
##' See also \url{http://davidakenny.net/cm/fit.htm}
##'
##'
##' @importFrom lavaan lavInspect
##'
##' @param object The lavaan model object provided after running the \code{cfa},
##' \code{sem}, \code{growth}, or \code{lavaan} functions.
##' @param scaled If \code{TRUE}, the scaled (or robust, if available) RMSEA
##'   is returned. Ignored if a robust test statistic was not requested.
##' @param silent If \code{TRUE}, do not print anything on the screen.
##'
##' @return A value of RMSEA of the null model (a \code{numeric} vector)
##'   returned invisibly.
##'
##' @author
##' Ruben Arslan (Humboldt-University of Berlin, \email{rubenarslan@@gmail.com})
##'
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @seealso
##' \itemize{
##'   \item \code{\link{miPowerFit}} For the modification indices and their
##'      power approach for model fit evaluation
##'   \item \code{\link{moreFitIndices}} For other fit indices
##' }
##'
##' @references Kenny, D. A., Kaniskan, B., & McCoach, D. B. (2015). The
##' performance of RMSEA in models with small degrees of freedom.
##' \emph{Sociological Methods Research, 44}(3), 486--507.
##' \doi{10.1177/0049124114543236}
##'
##' @examples
##'
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ x7 + x8 + x9 '
##'
##' fit <- cfa(HS.model, data = HolzingerSwineford1939)
##' nullRMSEA(fit)
##'
##' @export
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
#TODO: update to extract f from lresults. Make public?
sic <- function(f, lresults = NULL) {
  ## p. 596 of doi:10.1007/s10519-004-5587-0 says to use observed Fisher information
  E.inv <- lavaan::lavTech(lresults, "inverted.information.observed")
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



##' Small-\emph{N} correction for \eqn{chi^2} test statistic
##'
##' Calculate small-\emph{N} corrections for \eqn{chi^2} model-fit test
##' statistic to adjust for small sample size (relative to model size).
##'
##' Four finite-sample adjustments to the chi-squared statistic are currently
##' available, all of which are described in Shi et al. (2018). These all
##' assume normally distributed data, and may not work well with severely
##' nonnormal data. Deng et al. (2018, section 4) review proposed small-\emph{N}
##' adjustments that do not assume normality, which rarely show promise, so
##' they are not implemented here. This function currently will apply
##' small-\emph{N} adjustments to scaled test statistics with a warning that
##' they do not perform well (Deng et al., 2018).
##'
##' @importFrom lavaan lavInspect lavNames
##' @importFrom stats pchisq
##' @importFrom methods getMethod
##'
##' @param fit0,fit1 \linkS4class{lavaan} object(s) provided after running the
##'   \code{cfa}, \code{sem}, \code{growth}, or \code{lavaan} functions.
##'   \linkS4class{lavaan.mi} object(s) also accepted.
##' @param smallN.method \code{character} indicating the small-\emph{N}
##'   correction method to use. Multiple may be chosen (all of which assume
##'   normality), as described in Shi et al. (2018):
##'   \code{c("swain","yuan.2015","yuan.2005","bartlett")}. Users may also
##'   simply select \code{"all"}.
##' @param \dots Additional arguments to the \code{\link[lavaan]{lavTestLRT}} or
##'   \code{\link{lavTestLRT.mi}} functions. Ignored when \code{is.null(fit1)}.
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'   imputations from pooled results. Ignored unless \code{fit0} (and
##'   optionally \code{fit1}) is a \linkS4class{lavaan.mi} object. See
##'   \code{\link{lavTestLRT.mi}} for a description of options and defaults.
##'
##' @return A \code{list} of \code{numeric} vectors: one for the originally
##'   requested statistic(s), along with one per requested \code{smallN.method}.
##'   All include the the (un)adjusted test statistic, its \emph{df}, and the
##'   \emph{p} value for the test under the null hypothesis that the model fits
##'   perfectly (or that the 2 models have equivalent fit).
##'   The adjusted chi-squared statistic(s) also include(s) the scaling factor
##'   for the small-\emph{N} adjustment.
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Deng, L., Yang, M., & Marcoulides, K. M. (2018). Structural equation
##'   modeling with many variables: A systematic review of issues and
##'   developments. \emph{Frontiers in Psychology, 9}, 580.
##'   \doi{10.3389/fpsyg.2018.00580}
##'
##'   Shi, D., Lee, T., & Terry, R. A. (2018). Revisiting the model
##'   size effect in structural equation modeling.
##'   \emph{Structural Equation Modeling, 25}(1), 21--40.
##'   \doi{10.1080/10705511.2017.1369088}
##'
##' @examples
##'
##' HS.model <- '
##'     visual  =~ x1 + b1*x2 + x3
##'     textual =~ x4 + b2*x5 + x6
##'     speed   =~ x7 + b3*x8 + x9
##' '
##' fit1 <- cfa(HS.model, data = HolzingerSwineford1939[1:50,])
##' ## test a single model (implicitly compared to a saturated model)
##' chisqSmallN(fit1)
##'
##' ## fit a more constrained model
##' fit0 <- cfa(HS.model, data = HolzingerSwineford1939[1:50,],
##'             orthogonal = TRUE)
##' ## compare 2 models
##' chisqSmallN(fit1, fit0)
##'
##' @export
chisqSmallN <- function(fit0, fit1 = NULL,
                        smallN.method = if (is.null(fit1)) c("swain","yuan.2015") else "yuan.2005",
                        ..., omit.imps = c("no.conv","no.se")) {
  if ("all" %in% smallN.method) smallN.method <- c("swain","yuan.2015",
                                                   "yuan.2005","bartlett")
  smallN.method <- intersect(tolower(smallN.method),
                             c("swain","yuan.2015","yuan.2005","bartlett"))
  if (!any(smallN.method %in% c("swain","yuan.2015","yuan.2005","bartlett")))
    stop('No recognized options for "smallN.method" argument')

  ## check class
  if (!inherits(fit0, what = c("lavaan","lavaanList")))
    stop("this function is only applicable to fitted lavaan models")

  ## if there are 2 models...
  if (!is.null(fit1)) {

    ## check classes
    if (!inherits(fit1, what = c("lavaan","lavaanList")))
      stop("this function is only applicable to fitted lavaan models")
    modClass <- unique(sapply(list(fit0, fit1), class))
    if (length(modClass) > 1L) stop('All models must be of the same class (e.g.,',
                                    ' cannot compare lavaan objects to lavaan.mi)')

    ## check order of DF
    suppressMessages(DF0 <- getMethod("fitMeasures", class(fit0))(fit0, fit.measures = "df",
                                                                  omit.imps = omit.imps)[1])
    suppressMessages(DF1 <- getMethod("fitMeasures", class(fit1))(fit1, fit.measures = "df",
                                                                  omit.imps = omit.imps)[1])
    if (DF0 == DF1) stop("Models have the same degrees of freedom.")
    parent <- which.min(c(DF0, DF1))
    if (parent == 1L) {
      parent <- fit0
      fit0 <- fit1
      fit1 <- parent

      parent <- DF0
      DF0 <- DF1
      DF1 <- parent
    }
    if (min(c(DF0, DF1)) == 0L) {
      message('Less restricted model has df=0, so chi-squared difference ',
              'not needed to compare models. Using only the restricted ',
              "model's chi-squared statistic.")
      fit1 <- NULL
    }
  }

  ## check whether methods can be used
  if (!is.null(fit1)) {

    if (any(smallN.method %in% c("yuan.2015","swain"))) {
      message('Swain(1975) and Yuan (2015) corrections depend on the number ',
              'of free parameters, so it is unavailable for model comparison.')
      smallN.method <- smallN.method[-which(smallN.method %in% c("yuan.2015","swain"))]
    }

    if (!length(smallN.method)) {
      stop('No valid options for "smallN.method" argument')
    } else warning('Small-N corrections developed for single models, not for ',
                   'model comparison. Experimentally applying correction to ',
                   'chi-squared difference statistic, which might be invalid.')
  }

  ## save quantities relevant across correction methods
  N <- lavInspect(fit0, "ntotal")
  Ng <- lavInspect(fit0, "ngroups")
  if (!lavInspect(fit0, "options")$sample.cov.rescale) N <- N - Ng
  P <- length(lavNames(fit0))
  K <- length(lavNames(fit0, type = "lv")) # count latent factors

  if (is.null(fit1)) {
    FIT <- getMethod("fitMeasures", class(fit0))(fit0,
                                                 ## lavaan.mi arguments ignored
                                                 ## for lavaan objects
                                                 omit.imps = omit.imps,
                                                 asymptotic = TRUE,
                                                 fit.measures = c("npar","chisq",
                                                                  "df","pvalue"))
    scaled <- any(grepl(pattern = "scaled", x = names(FIT)))
    if (scaled) warning('Small-N corrections developed assuming normality, but',
                        ' a scaled test was requested. Applying correction(s) ',
                        'to the scaled test statistic, but this has not ',
                        'performed well in past simulations.')
    NPAR <- FIT[["npar"]]
    chi <- FIT[[if (scaled) "chisq.scaled" else "chisq"]]
    DF  <- FIT[[if (scaled) "df.scaled" else "df"]]
    PV  <- FIT[[if (scaled) "pvalue.scaled" else "pvalue"]]

  } else {
    ## Compare to a second model. Check matching stats.
    N1 <- lavInspect(fit1, "ntotal")
    Ng1 <- lavInspect(fit1, "ngroups")
    if (!lavInspect(fit1, "options")$sample.cov.rescale) N1 <- N1 - Ng1
    if (N != N1) stop("Unequal sample sizes")
    if (P != length(lavNames(fit1))) stop("Unequal number of observed variables")
    K1 <- length(lavNames(fit1, type = "lv"))
    if (K != K1 && any(smallN.method %in% c("yuan.2005","bartlett"))) {
      warning("Unequal number of latent variables (k). Unclear how to apply ",
              "Yuan (2005) or Bartlett (2015) corrections when comparing ",
              "models with different k. Experimentally using the larger ",
              "model's k, but there is no evidence this is valid.")
      K <- max(K, K1)
    }

    try(AOV <- compareFit(fit0, fit1, argsLRT = list(...), indices = FALSE),
        silent = TRUE)
    if (inherits(AOV, "try-error")) stop('Model comparison failed. Try using ',
                                         'lavTestLRT() to investigate why.')

    if (inherits(fit0, "lavaan")) {
      if (grepl("scaled", attr(AOV@nested, "heading"), ignore.case = TRUE))
        warning('Small-N corrections developed assuming normality, but scaled ',
                'tests were requested. Applying correction(s) to the scaled test',
                ' statistic, but this has not performed well in past simulations.')
      chi <- AOV@nested[["Chisq diff"]][2]
      DF  <- AOV@nested[["Df diff"]][2]
      PV  <- AOV@nested[["Pr(>Chisq)"]][2]

    } else if (inherits(fit0, "lavaan.mi")) {
      if (any(grepl("scaled", colnames(AOV@nested), ignore.case = TRUE)))
        warning('Small-N corrections developed assuming normality, but scaled ',
                'tests were requested. Applying correction(s) to the scaled test',
                ' statistic, but this has not performed well in past simulations.')
      chi <- AOV@nested[1, 1]
      DF  <- AOV@nested[1, 2]
      PV  <- AOV@nested[1, 3]
    }

  }


  ## empty list to store correction(s)
  out <- list()
  out[[ lavInspect(fit0, "options")$test ]] <- c(chisq = chi, df = DF,
                                                 pvalue = PV)
  class(out[[1]]) <- c("lavaan.vector","numeric")

  ## calculate Swain's (1975) correction
  ## (undefined for model comparison)
  if ("swain" %in% smallN.method) {
    Q <- (sqrt(1 + 8*NPAR) - 1) / 2
    num <- P*(2*P^2 + 3*P - 1) - Q*(2*Q^2 + 3*Q - 1)
    SC <- 1 - num / (12*DF*N)
    out[["swain"]] <- c(chisq = chi*SC, df = DF,
                        pvalue = pchisq(chi*SC, DF, lower.tail = FALSE),
                        smallN.factor = SC)
    class(out[["swain"]]) <- c("lavaan.vector","numeric")
  }

  ## calculate Yuan's (2015) correction
  ## (undefined for model comparison)
  if ("yuan.2015" %in% smallN.method) {
    ## numerator uses actual N regardless of sample.cov.rescale
    SC <- (lavInspect(fit0, "ntotal") - (2.381 + .361*P + .006*NPAR)) / N
    out[["yuan.2015"]] <- c(chisq = chi*SC, df = DF,
                            pvalue = pchisq(chi*SC, DF, lower.tail = FALSE),
                            smallN.factor = SC)
    class(out[["yuan.2015"]]) <- c("lavaan.vector","numeric")
  }

  ## calculate Yuan's (2005) correction
  if ("yuan.2005" %in% smallN.method) {
    SC <- 1 - ((2*P + 2*K + 7) / (6*N))
    out[["yuan.2005"]] <- c(chisq = chi*SC, df = DF,
                            pvalue = pchisq(chi*SC, DF, lower.tail = FALSE),
                            smallN.factor = SC)
    class(out[["yuan.2005"]]) <- c("lavaan.vector","numeric")
  }

  ## calculate Bartlett's (1950) k-factor correction (ONLY appropriate for EFA)
  if ("bartlett" %in% smallN.method) {
    message('Bartlett\'s k-factor correction was developed for EFA models, ',
            'not for general SEMs.')
    SC <- 1 - ((2*P + 4*K + 5) / (6*N))
    out[["bartlett"]] <- c(chisq = chi*SC, df = DF,
                           pvalue = pchisq(chi*SC, DF, lower.tail = FALSE),
                           smallN.factor = SC)
    class(out[["bartlett"]]) <- c("lavaan.vector","numeric")
  }

  out[c(lavInspect(fit0, "options")$test, smallN.method)]
}


