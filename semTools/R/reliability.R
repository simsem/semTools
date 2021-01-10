### Sunthud Pornprasertmanit, Terrence D. Jorgensen, Yves Rosseel
### Last updated: 10 January 2021


## -------------
## reliability()
## -------------


##' Calculate reliability values of factors
##'
##' Calculate reliability values of factors by coefficient omega
##'
##' The coefficient alpha (Cronbach, 1951) can be calculated by
##'
##' \deqn{ \alpha = \frac{k}{k - 1}\left[ 1 - \frac{\sum^{k}_{i = 1}
##' \sigma_{ii}}{\sum^{k}_{i = 1} \sigma_{ii} + 2\sum_{i < j} \sigma_{ij}}
##' \right],}
##'
##' where \eqn{k} is the number of items in a factor, \eqn{\sigma_{ii}} is the
##' item \emph{i} observed variances, \eqn{\sigma_{ij}} is the observed
##' covariance of items \emph{i} and \emph{j}.
##'
##' The coefficient omega (Bollen, 1980; see also Raykov, 2001) can be
##' calculated by
##'
##' \deqn{ \omega_1 =\frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
##' Var\left( \psi \right)}{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
##' Var\left( \psi \right) + \sum^{k}_{i = 1} \theta_{ii} + 2\sum_{i < j}
##' \theta_{ij} }, }
##'
##' where \eqn{\lambda_i} is the factor loading of item \emph{i}, \eqn{\psi} is
##' the factor variance, \eqn{\theta_{ii}} is the variance of measurement errors
##' of item \emph{i}, and \eqn{\theta_{ij}} is the covariance of measurement
##' errors from item \emph{i} and \emph{j}.
##'
##' The second coefficient omega (Bentler, 1972, 2009) can be calculated by
##'
##' \deqn{ \omega_2 = \frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
##' Var\left( \psi \right)}{\bold{1}^\prime \hat{\Sigma} \bold{1}}, }
##'
##' where \eqn{\hat{\Sigma}} is the model-implied covariance matrix, and
##' \eqn{\bold{1}} is the \eqn{k}-dimensional vector of 1. The first and the
##' second coefficients omega will have the same value when the model has simple
##' structure, but different values when there are (for example) cross-loadings
##' or method factors. The first coefficient omega can be viewed as the
##' reliability controlling for the other factors (like \eqn{\eta^2_{partial}} in
##' ANOVA). The second coefficient omega can be viewed as the unconditional
##' reliability (like \eqn{\eta^2} in ANOVA).
##'
##' The third coefficient omega (McDonald, 1999), which is sometimes referred to
##' hierarchical omega, can be calculated by
##'
##' \deqn{ \omega_3 =\frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
##' Var\left( \psi \right)}{\bold{1}^\prime \Sigma \bold{1}}, }
##'
##' where \eqn{\Sigma} is the observed covariance matrix. If the model fits the
##' data well, the third coefficient omega will be similar to the
##' \eqn{\omega_2}. Note that if there is a directional effect in the model, all
##' coefficients omega will use the total factor variances, which is calculated
##' by \code{\link[lavaan]{lavInspect}(object, "cov.lv")}.
##'
##' In conclusion, \eqn{\omega_1}, \eqn{\omega_2}, and \eqn{\omega_3} are
##' different in the denominator. The denominator of the first formula assumes
##' that a model is congeneric factor model where measurement errors are not
##' correlated. The second formula accounts for correlated measurement errors.
##' However, these two formulas assume that the model-implied covariance matrix
##' explains item relationships perfectly. The residuals are subject to sampling
##' error. The third formula use observed covariance matrix instead of
##' model-implied covariance matrix to calculate the observed total variance.
##' This formula is the most conservative method in calculating coefficient
##' omega.
##'
##' The average variance extracted (AVE) can be calculated by
##'
##' \deqn{ AVE = \frac{\bold{1}^\prime
##' \textrm{diag}\left(\Lambda\Psi\Lambda^\prime\right)\bold{1}}{\bold{1}^\prime
##' \textrm{diag}\left(\hat{\Sigma}\right) \bold{1}}, }
##'
##' Note that this formula is modified from Fornell & Larcker (1981) in the case
##' that factor variances are not 1. The proposed formula from Fornell & Larcker
##' (1981) assumes that the factor variances are 1. Note that AVE will not be
##' provided for factors consisting of items with dual loadings. AVE is the
##' property of items but not the property of factors.
##'
##' Regarding categorical indicators, coefficient alpha and AVE are calculated
##' based on polychoric correlations. The coefficient alpha from this function
##' differs from the standard alpha calculation, which does not assume items are
##' continuous, so numerically weighted categories can be treated as numeric.
##' Researchers may check the \code{alpha} function in the \code{psych} package
##' for the standard coefficient alpha calculation.
##'
##' Item thresholds are not accounted for. Coefficient omega for categorical
##' items, however, is calculated by accounting for both item covariances and
##' item thresholds using Green and Yang's (2009, formula 21) approach. Three
##' types of coefficient omega indicate different methods to calculate item
##' total variances. The original formula from Green and Yang is equivalent to
##' \eqn{\omega_3} in this function. Green and Yang did not propose a method for
##' calculating reliability with a mixture of categorical and continuous
##' indicators, and we are currently unaware of an appropriate method.
##' Therefore, when \code{reliability} detects both categorical and continuous
##' indicators in the model, an error is returned. If the categorical indicators
##' load on a different factor(s) than continuous indicators, then reliability
##' can be calculated separately for those scales by fitting separate models and
##' submitting each to the \code{reliability} function.
##'
##'
##' @importFrom lavaan lavInspect lavNames
##' @importFrom methods getMethod
##'
##' @param object A \code{\linkS4class{lavaan}} or
##'   \code{\linkS4class{lavaan.mi}} object, expected to contain only
##'   exogenous common factors (i.e., a CFA model).
##' @param return.total \code{logical} indicating whether to return a final
##'   column containing the reliability of a composite of all items. Ignored
##'   in 1-factor models, and should only be set \code{TRUE} if all factors
##'   represent scale dimensions that could nonetheless be collapsed to a
##'   single scale composite (scale sum or scale mean).
##' @param dropSingle \code{logical} indicating whether to exclude factors
##'   defined by a single indicator from the returned results. If \code{TRUE}
##'   (default), single indicators will still be included in the \code{total}
##'   column when \code{return.total = TRUE}.
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'   imputations from pooled results.  Can include any of
##'   \code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
##'   default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (\code{"no.npd"}) would exclude any imputations which
##'   yielded a nonpositive definite covariance matrix for observed or
##'   latent variables, which would include any "improper solutions" such
##'   as Heywood cases.  NPD solutions are not excluded by default because
##'   they are likely to occur due to sampling error, especially in small
##'   samples.  However, gross model misspecification could also cause
##'   NPD solutions, users can compare pooled results with and without
##'   this setting as a sensitivity analysis to see whether some
##'   imputations warrant further investigation.
##'
##' @return Reliability values (coefficient alpha, coefficients omega, average
##'   variance extracted) of each factor in each group. If there are multiple
##'   factors, a \code{total} column can optionally be included.
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##'   Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
##'
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @seealso \code{\link{reliabilityL2}} for reliability value of a desired
##' second-order factor, \code{\link{maximalRelia}} for the maximal reliability
##' of weighted composite
##'
##' @references
##' Bollen, K. A. (1980). Issues in the comparative measurement of
##' political democracy. \emph{American Sociological Review, 45}(3), 370--390.
##' \doi{10.2307/2095172}
##'
##' Bentler, P. M. (1972). A lower-bound method for the dimension-free
##' measurement of internal consistency. \emph{Social Science Research, 1}(4),
##' 343--357. \doi{10.1016/0049-089X(72)90082-8}
##'
##' Bentler, P. M. (2009). Alpha, dimension-free, and model-based internal
##' consistency reliability. \emph{Psychometrika, 74}(1), 137--143.
##' \doi{10.1007/s11336-008-9100-1}
##'
##' Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
##' tests. \emph{Psychometrika, 16}(3), 297--334. \doi{10.1007/BF02310555}
##'
##' Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models
##' with unobservable variables and measurement errors. \emph{Journal of
##' Marketing Research, 18}(1), 39--50. \doi{10.2307/3151312}
##'
##' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
##' structural equation modeling: An alternative to coefficient alpha.
##' \emph{Psychometrika, 74}(1), 155--167. \doi{10.1007/s11336-008-9099-3}
##'
##' McDonald, R. P. (1999). \emph{Test theory: A unified treatment}. Mahwah, NJ:
##' Erlbaum.
##'
##' Raykov, T. (2001). Estimation of congeneric scale reliability using
##' covariance structure analysis with nonlinear constraints \emph{British
##' Journal of Mathematical and Statistical Psychology, 54}(2), 315--323.
##' \doi{10.1348/000711001159582}
##'
##' @examples
##'
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ x7 + x8 + x9 '
##'
##' fit <- cfa(HS.model, data = HolzingerSwineford1939)
##' reliability(fit)
##' reliability(fit, return.total = TRUE)
##'
##' @export
reliability <- function(object, return.total = FALSE, dropSingle = TRUE,
                        omit.imps = c("no.conv","no.se")) {

  ngroups <- lavInspect(object, "ngroups") #TODO: adapt to multiple levels
  nlevels <- lavInspect(object, "nlevels")
  nblocks <- ngroups*nlevels #FIXME: always true?
  group.label <- if (ngroups > 1L) lavInspect(object, "group.label") else NULL
  #FIXME? lavInspect(object, "level.labels")
  clus.label <- if (nlevels > 1L) c("within", lavInspect(object, "cluster")) else NULL
  if (nblocks > 1L) {
    block.label <- paste(rep(group.label, each = nlevels), clus.label,
                         sep = if (ngroups > 1L && nlevels > 1L) "_" else "")
  }

  ## parameters in GLIST format (not flat, need block-level list)
  if (inherits(object, "lavaan")) {
    param <- lavInspect(object, "est")
    ve <- lavInspect(object, "cov.lv") # model-implied latent covariance matrix
    S <- object@h1$implied$cov # observed sample covariance matrix (already a list)

    if (nblocks == 1L) {
      param <- list(param)
      ve <- list(ve)
    }

  } else if (inherits(object, "lavaan.mi")) {
    useImps <- rep(TRUE, length(object@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

    param <- object@coefList[[ useImps[1] ]] # first admissible as template
    coefList <- object@coefList[useImps]
    phiList <- object@phiList[useImps]
    ## add block-level list per imputation?
    if (nblocks == 1L) {
      param <- list(param)
      for (i in 1:m) {
        coefList[[i]] <- list(coefList[[i]])
        phiList[[i]] <- list(phiList[[i]])
      }
    }
    S <- vector("list", nblocks) # pooled observed covariance matrix
    ve <- vector("list", nblocks)
    ## loop over blocks
    for (b in 1:nblocks) {

      ## param:  loop over GLIST elements
      for (mat in names(param[[b]])) {
        matList <- lapply(coefList, function(i) i[[b]][[mat]])
        param[[b]][[mat]] <- Reduce("+", matList) / length(matList)
      } # mat

      ## pooled observed covariance matrix
      covList <- lapply(object@h1List[useImps], function(i) i$implied$cov[[b]])
      S[[b]] <- Reduce("+", covList) / m

      ## pooled model-implied latent covariance matrix
      ve[[b]] <- Reduce("+", lapply(phiList, "[[", i = b) ) / m

    } # b

  }

  if (nblocks == 1L) {
    SigmaHat <- getMethod("fitted", class(object))(object)["cov"] # retain list format
  } else {
    SigmaHat <- sapply(getMethod("fitted", class(object))(object),
                       "[[", "cov", simplify = FALSE)
  }

  ly <- lapply(param, "[[", "lambda")
  te <- lapply(param, "[[", "theta")

	anyCategorical <- lavInspect(object, "categorical")
	threshold <- if (anyCategorical) getThreshold(object) else NULL
	latScales <- if (anyCategorical) getScales(object) else NULL

	result <- list()
	warnHigher <- FALSE
	## loop over i blocks (groups/levels)
	for (i in 1:nblocks) {
	  ## extract factor and indicator names
	  indNames <- rownames(ly[[i]])
	  facNames <- colnames(ly[[i]])
	  ## distinguish between categorical, continuous, and latent indicators
	  nameArgs <- list(object = object)
	  if (nblocks > 1L) nameArgs$block <- i
	  ordNames <- do.call(lavNames, c(nameArgs, list(type = "ov.ord")))
	  latInds  <- do.call(lavNames, c(nameArgs, list(type = "lv.ind")))
	  higher   <- setdiff(facNames, latInds)
	  ## keep track of factor indices, to drop higher-order factors,
	  ## and optionally single-indicator factors
	  idx.drop <- numeric(0)

	  ## relevant quantities
		common <- (apply(ly[[i]], 2, sum)^2) * diag(ve[[i]])
		truevar <- ly[[i]] %*% ve[[i]] %*% t(ly[[i]])
		## vectors to store results for each factor
		error <- rep(NA, length(common))
		alpha <- rep(NA, length(common))
		total <- rep(NA, length(common))
		omega1 <- omega2 <- omega3 <- rep(NA, length(common))
		impliedTotal <- rep(NA, length(common))
		avevar <- rep(NA, length(common))

		## loop over j factors
		for (j in 1:length(common)) {
			index <- which(ly[[i]][,j] != 0)
			## check for categorical (or mixed) indicators
			categorical <-      any(indNames[index] %in% ordNames)
			if (categorical && !all(indNames[index] %in% ordNames)) {
			  stop('Reliability cannot be computed for factors with combinations ',
			       'of categorical and continuous (including latent) indicators.')
			}
			## check for latent indicators
			if (length(latInds) && facNames[j] %in% higher) {
			  warnHigher <- TRUE
			  idx.drop <- c(idx.drop, j)
			  next
			}
			if (dropSingle && length(index) == 1L) {
			  idx.drop <- c(idx.drop, j)
			  next
			}

			sigma <- S[[i]][index, index, drop = FALSE]
			alpha[j] <- computeAlpha(sigma)
			faccontrib <- ly[[i]][,j, drop = FALSE] %*% ve[[i]][j,j, drop = FALSE] %*% t(ly[[i]][,j, drop = FALSE])
			truefac <- diag(faccontrib[index, index, drop = FALSE])
			trueitem <- diag(truevar[index, index, drop = FALSE])
			erritem <- diag(te[[i]][index, index, drop = FALSE])
			if (sum(abs(trueitem - truefac)) < 0.00001) {
				avevar[j] <- sum(trueitem) / sum(trueitem + erritem)
			} else {
				avevar[j] <- NA
			}
			if (categorical) {
				omega1[j] <- omegaCat(truevar = faccontrib[index, index, drop = FALSE],
				                      threshold = threshold[[i]][index],
				                      scales = latScales[[i]][index],
				                      denom = faccontrib[index, index, drop = FALSE] + te[[i]][index, index, drop = FALSE])
				omega2[j] <- omegaCat(truevar = faccontrib[index, index, drop = FALSE],
				                      threshold = threshold[[i]][index],
				                      scales = latScales[[i]][index],
				                      denom = SigmaHat[[i]][index, index, drop = FALSE])
				omega3[j] <- omegaCat(truevar = faccontrib[index, index, drop = FALSE],
				                      threshold = threshold[[i]][index],
				                      scales = latScales[[i]][index],
				                      denom = sigma)
			} else {
			  commonfac <- sum(faccontrib[index, index, drop = FALSE])
			  error[j] <- sum(te[[i]][index, index, drop = FALSE])
			  impliedTotal[j] <- sum(SigmaHat[[i]][index, index, drop = FALSE])
			  total[j] <- sum(sigma)

			  omega1[j] <- commonfac / (commonfac + error[j])
				omega2[j] <- commonfac / impliedTotal[j]
				omega3[j] <- commonfac / total[j]
			}
			## end loop over j factors
		}

		if (return.total & length(facNames) > 1L) {
		  alpha <- c(alpha, total = computeAlpha(S[[i]]))
		  #FIXME: necessary?    names(alpha) <- c(names(common), "total")
		  if (categorical) {
		    omega1 <- c(omega1, total = omegaCat(truevar = truevar,
		                                         threshold = threshold[[i]],
		                                         scales = latScales[[i]],
		                                         denom = truevar + te[[i]]))
		    omega2 <- c(omega2, total = omegaCat(truevar = truevar,
		                                         threshold = threshold[[i]],
		                                         scales = latScales[[i]],
		                                         denom = SigmaHat[[i]]))
		    omega3 <- c(omega3, total = omegaCat(truevar = truevar,
		                                         threshold = threshold[[i]],
		                                         scales = latScales[[i]],
		                                         denom = S[[i]]))
		  } else {
		    omega1 <- c(omega1, total = sum(truevar) / (sum(truevar) + sum(te[[i]])))
		    omega2 <- c(omega2, total = sum(truevar) / (sum(SigmaHat[[i]])))
		    omega3 <- c(omega3, total = sum(truevar) / (sum(S[[i]])))
		  }
		  avevar <- c(avevar,
		              total = sum(diag(truevar)) / sum((diag(truevar) + diag(te[[i]]))))
		}

		result[[i]] <- rbind(alpha = alpha, omega = omega1, omega2 = omega2,
		                     omega3 = omega3, avevar = avevar)
		colnames(result[[i]])[1:length(facNames)] <- facNames
		if (return.total & length(facNames) > 1L) {
		  colnames(result[[i]])[ ncol(result[[i]]) ] <- "total"
		}
		if (length(idx.drop)) {
		  result[[i]] <- result[[i]][ , -idx.drop]
		  ## reset indices for next block (could have different model/variables)
		  idx.drop <- numeric(0)
		}
		## end loop over blocks
	}

	if (anyCategorical) message("For constructs with categorical indicators, ",
	                            "the alpha and the average variance extracted ",
	                            "are calculated from polychoric (polyserial) ",
	                            "correlations, not from Pearson correlations.\n")
	if (warnHigher) message('Higher-order factors were ignored.\n')

	## drop list structure?
	if (nblocks == 1L) {
		result <- result[[1]]
	} else names(result) <- block.label

	result
}



## ---------------
## reliabilityL2()
## ---------------

##' Calculate the reliability values of a second-order factor
##'
##' Calculate the reliability values (coefficient omega) of a second-order
##' factor
##'
##' The first formula of the coefficient omega (in the
##' \code{\link{reliability}}) will be mainly used in the calculation. The
##' model-implied covariance matrix of a second-order factor model can be
##' separated into three sources: the second-order common-factor variance,
##' the residual variance of the first-order common factors (i.e., not
##' accounted for by the second-order factor), and the measurement error of
##' observed indicators:
##'
##' \deqn{ \hat{\Sigma} = \Lambda \bold{B} \Phi_2 \bold{B}^{\prime}
##' \Lambda^{\prime} + \Lambda \Psi_{u} \Lambda^{\prime} + \Theta, }
##'
##' where \eqn{\hat{\Sigma}} is the model-implied covariance matrix,
##' \eqn{\Lambda} contains first-order factor loadings, \eqn{\bold{B}} contains
##' second-order factor loadings, \eqn{\Phi_2} is the covariance matrix of the
##' second-order factor(s), \eqn{\Psi_{u}} is the covariance matrix of residuals
##' from first-order factors, and \eqn{\Theta} is the covariance matrix of the
##' measurement errors from observed indicators. Thus, we can calculate the
##' proportion of variance of a composite score calculated from the observed
##' indicators (e.g., a total score or scale mean) that is attributable to the
##' second-order factor, i.e. coefficient omega at Level 1:
##'
##' \deqn{ \omega_{L1} = \frac{\bold{1}^{\prime} \Lambda \bold{B} \Phi_2
##' \bold{B}^{\prime} \Lambda^{\prime} \bold{1}}{\bold{1}^{\prime} \Lambda
##' \bold{B} \Phi_2 \bold{B} ^{\prime} \Lambda^{\prime} \bold{1} +
##' \bold{1}^{\prime} \Lambda \Psi_{u} \Lambda^{\prime} \bold{1} +
##' \bold{1}^{\prime} \Theta \bold{1}}, }
##'
##' where \eqn{\bold{1}} is the \emph{k}-dimensional vector of 1 and \emph{k} is
##' the number of observed variables.
##'
##' The model-implied covariance matrix among first-order factors (\eqn{\Phi_1})
##' can be calculated as:
##'
##' \deqn{ \Phi_1 = \bold{B} \Phi_2 \bold{B}^{\prime} + \Psi_{u}, }
##'
##' Thus, the proportion of variance among first-order common factors that is
##' attributable to the second-order factor (i.e., coefficient omega at Level 2)
##' can be calculated as:
##'
##' \deqn{ \omega_{L2} = \frac{\bold{1_F}^{\prime} \bold{B} \Phi_2
##' \bold{B}^{\prime} \bold{1_F}}{\bold{1_F}^{\prime} \bold{B} \Phi_2
##' \bold{B}^{\prime} \bold{1_F} + \bold{1_F}^{\prime} \Psi_{u} \bold{1_F}}, }
##'
##' where \eqn{\bold{1_F}} is the \emph{F}-dimensional vector of 1 and \emph{F}
##' is the number of first-order factors. This Level-2 omega can be interpreted
##' as an estimate of the reliability of a hypothetical composite calculated
##' from error-free observable variables representing the first-order common
##' factors. This might only be meaningful as a thought experiment.
##'
##' An additional thought experiment is possible: If the observed indicators
##' contained only the second-order common-factor variance and unsystematic
##' measurement error, then there would be no first-order common factors because
##' their unique variances would be excluded from the observed measures. An
##' estimate of this hypothetical composite reliability can be calculated as the
##' partial coefficient omega at Level 1, or the proportion of observed
##' variance explained by the second-order factor after partialling out the
##' uniqueness from the first-order factors:
##'
##' \deqn{ \omega_{L1} = \frac{\bold{1}^{\prime} \Lambda \bold{B} \Phi_2
##' \bold{B}^{\prime} \Lambda^{\prime} \bold{1}}{\bold{1}^{\prime} \Lambda
##' \bold{B} \Phi_2 \bold{B}^{\prime} \Lambda^{\prime} \bold{1} +
##' \bold{1}^{\prime} \Theta \bold{1}}, }
##'
##' Note that if the second-order factor has a direct factor loading on some
##' observed variables, the observed variables will be counted as first-order
##' factors, which might not be desirable.
##'
##'
##' @importFrom lavaan lavInspect
##'
##' @param object A \code{\linkS4class{lavaan}} or
##'   \code{\linkS4class{lavaan.mi}} object, expected to contain a least one
##'   exogenous higher-order common factor.
##' @param secondFactor The name of a single second-order factor in the
##'   model fitted in \code{object}. The function must be called multiple
##'   times to estimate reliability for each higher-order factor.
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'        imputations from pooled results.  Can include any of
##'        \code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
##'        default setting, which excludes any imputations that did not
##'        converge or for which standard errors could not be computed.  The
##'        last option (\code{"no.npd"}) would exclude any imputations which
##'        yielded a nonpositive definite covariance matrix for observed or
##'        latent variables, which would include any "improper solutions" such
##'        as Heywood cases.  NPD solutions are not excluded by default because
##'        they are likely to occur due to sampling error, especially in small
##'        samples.  However, gross model misspecification could also cause
##'        NPD solutions, users can compare pooled results with and without
##'        this setting as a sensitivity analysis to see whether some
##'        imputations warrant further investigation.
##'
##' @return Reliability values at Levels 1 and 2 of the second-order factor, as
##'   well as the partial reliability value at Level 1
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' @seealso
##'   \code{\link{reliability}} for the reliability of the first-order factors.
##'
##' @examples
##'
##' HS.model3 <- ' visual  =~ x1 + x2 + x3
##'                textual =~ x4 + x5 + x6
##'                speed   =~ x7 + x8 + x9
##' 			         higher =~ visual + textual + speed'
##'
##' fit6 <- cfa(HS.model3, data = HolzingerSwineford1939)
##' reliability(fit6) # Should provide a warning for the endogenous variables
##' reliabilityL2(fit6, "higher")
##'
##' @export
reliabilityL2 <- function(object, secondFactor,
                          omit.imps = c("no.conv","no.se")) {
  secondFactor <- as.character(secondFactor)[1] # only one at a time

  ngroups <- lavInspect(object, "ngroups") #TODO: adapt to multiple levels
  nlevels <- lavInspect(object, "nlevels")
  nblocks <- ngroups*nlevels #FIXME: always true?
  group.label <- if (ngroups > 1L) lavInspect(object, "group.label") else NULL
  #FIXME? lavInspect(object, "level.labels")
  clus.label <- if (nlevels > 1L) c("within", lavInspect(object, "cluster")) else NULL
  if (nblocks > 1L) {
    block.label <- paste(rep(group.label, each = nlevels), clus.label,
                         sep = if (ngroups > 1L && nlevels > 1L) "_" else "")
  }

  ## parameters in GLIST format (not flat, need block-level list)
  if (inherits(object, "lavaan")) {
    param <- lavInspect(object, "est")
    ve <- lavInspect(object, "cov.lv") # model-implied latent covariance matrix
    S <- object@h1$implied$cov # observed sample covariance matrix (already a list)

    if (nblocks == 1L) {
      param <- list(param)
      ve <- list(ve)
    }

  } else if (inherits(object, "lavaan.mi")) {
    useImps <- rep(TRUE, length(object@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

    param <- object@coefList[[ useImps[1] ]] # first admissible as template
    coefList <- object@coefList[useImps]
    phiList <- object@phiList[useImps]
    ## add block-level list per imputation?
    if (nblocks == 1L) {
      param <- list(param)
      for (i in 1:m) {
        coefList[[i]] <- list(coefList[[i]])
        phiList[[i]] <- list(phiList[[i]])
      }
    }

    S <- vector("list", nblocks) # pooled observed covariance matrix
    ve <- vector("list", nblocks)
    ## loop over blocks
    for (b in 1:nblocks) {

      ## param:  loop over GLIST elements
      for (mat in names(param[[b]])) {
        matList <- lapply(coefList, function(i) i[[b]][[mat]])
        param[[b]][[mat]] <- Reduce("+", matList) / length(matList)
      } # mat

      ## pooled observed covariance matrix
      covList <- lapply(object@h1List[useImps], function(i) i$implied$cov[[b]])
      S[[b]] <- Reduce("+", covList) / m

      ## pooled model-implied latent covariance matrix
      ve[[b]] <- Reduce("+", lapply(phiList, "[[", i = b) ) / m

    } # b

  }

  if (nblocks == 1L) {
    SigmaHat <- getMethod("fitted", class(object))(object)["cov"] # retain list format
  } else {
    SigmaHat <- sapply(getMethod("fitted", class(object))(object),
                       "[[", "cov", simplify = FALSE)
  }

  ly <- lapply(param, "[[", "lambda")
  te <- lapply(param, "[[", "theta")
  ps <- lapply(param, "[[", "psi")
	be <- lapply(param, "[[", "beta")

	result <- list()
	for (i in 1:nblocks) {

		# Prepare for higher-order reliability
		l2var <- ve[[i]][secondFactor, secondFactor, drop = FALSE]
		l2load <- be[[1]][,secondFactor]
		indexl2 <- which(l2load != 0)
		commonl2 <- (sum(l2load)^2) * l2var
		errorl2 <- sum(ps[[i]][indexl2, indexl2, drop = FALSE])

		# Prepare for lower-order reliability
		indexl1 <- which(apply(ly[[i]][,indexl2, drop = FALSE], 1, function(x) sum(x != 0)) > 0)
		l1load <- ly[[i]][,indexl2] %*% as.matrix(be[[1]][indexl2, secondFactor, drop = FALSE])
		commonl1 <- (sum(l1load)^2) * l2var
		errorl1 <- sum(te[[i]][indexl1, indexl1, drop = FALSE])
		uniquel1 <- 0
		for (j in seq_along(indexl2)) {
			uniquel1 <- uniquel1 + (sum(ly[[i]][,indexl2[j]])^2) * ps[[i]][indexl2[j], indexl2[j], drop = FALSE]
		}

		# Adjustment for direct loading from L2 to observed variables
		if (any(ly[[i]][,secondFactor] != 0)) {
			indexind <- which(ly[[i]][,secondFactor] != 0)
			if (length(intersect(indexind, indexl1)) > 0)
			  stop("Direct and indirect loadings of higher-order factor to observed",
			       " variables are specified at the same time.")
			commonl2 <- sum(c(ly[[i]][,secondFactor], l2load))^2 * l2var
			errorl2 <- errorl2 + sum(te[[i]][indexind, indexind, drop = FALSE])
			commonl1 <- sum(c(ly[[i]][,secondFactor], l1load))^2 * l2var
			errorl1 <- errorl1 + sum(te[[i]][indexind, indexind, drop = FALSE])
		}

		# Calculate Reliability
		omegaL1 <- commonl1 / (commonl1 + uniquel1 + errorl1)
		omegaL2 <- commonl2 / (commonl2 + errorl2)
		partialOmegaL1 <- commonl1 / (commonl1 + errorl1)
		result[[i]] <- c(omegaL1 = omegaL1, omegaL2 = omegaL2, partialOmegaL1 = partialOmegaL1)
	}

	if (nblocks == 1L) {
	  result <- result[[1]]
	} else names(result) <- block.label

	result
}



## --------------
## maximalRelia()
## --------------

##' Calculate maximal reliability
##'
##' Calculate maximal reliability of a scale
##'
##' Given that a composite score (\eqn{W}) is a weighted sum of item scores:
##'
##' \deqn{ W = \bold{w}^\prime \bold{x} ,}
##'
##' where \eqn{\bold{x}} is a \eqn{k \times 1} vector of the scores of each
##' item, \eqn{\bold{w}} is a \eqn{k \times 1} weight vector of each item, and
##' \eqn{k} represents the number of items. Then, maximal reliability is
##' obtained by finding \eqn{\bold{w}} such that reliability attains its maximum
##' (Li, 1997; Raykov, 2012). Note that the reliability can be obtained by
##'
##' \deqn{ \rho = \frac{\bold{w}^\prime \bold{S}_T \bold{w}}{\bold{w}^\prime
##' \bold{S}_X \bold{w}}}
##'
##' where \eqn{\bold{S}_T} is the covariance matrix explained by true scores and
##' \eqn{\bold{S}_X} is the observed covariance matrix. Numerical method is used
##' to find \eqn{\bold{w}} in this function.
##'
##' For continuous items, \eqn{\bold{S}_T} can be calculated by
##'
##' \deqn{ \bold{S}_T = \Lambda \Psi \Lambda^\prime,}
##'
##' where \eqn{\Lambda} is the factor loading matrix and \eqn{\Psi} is the
##' covariance matrix among factors. \eqn{\bold{S}_X} is directly obtained by
##' covariance among items.
##'
##' For categorical items, Green and Yang's (2009) method is used for
##' calculating \eqn{\bold{S}_T} and \eqn{\bold{S}_X}. The element \eqn{i} and
##' \eqn{j} of \eqn{\bold{S}_T} can be calculated by
##'
##' \deqn{ \left[\bold{S}_T\right]_{ij} = \sum^{C_i - 1}_{c_i = 1} \sum^{C_j -
##' 1}_{c_j - 1} \Phi_2\left( \tau_{x_{c_i}}, \tau_{x_{c_j}}, \left[ \Lambda
##' \Psi \Lambda^\prime \right]_{ij} \right) - \sum^{C_i - 1}_{c_i = 1}
##' \Phi_1(\tau_{x_{c_i}}) \sum^{C_j - 1}_{c_j - 1} \Phi_1(\tau_{x_{c_j}}),}
##'
##' where \eqn{C_i} and \eqn{C_j} represents the number of thresholds in Items
##' \eqn{i} and \eqn{j}, \eqn{\tau_{x_{c_i}}} represents the threshold \eqn{c_i}
##' of Item \eqn{i}, \eqn{\tau_{x_{c_j}}} represents the threshold \eqn{c_i} of
##' Item \eqn{j}, \eqn{ \Phi_1(\tau_{x_{c_i}})} is the cumulative probability of
##' \eqn{\tau_{x_{c_i}}} given a univariate standard normal cumulative
##' distribution and \eqn{\Phi_2\left( \tau_{x_{c_i}}, \tau_{x_{c_j}}, \rho
##' \right)} is the joint cumulative probability of \eqn{\tau_{x_{c_i}}} and
##' \eqn{\tau_{x_{c_j}}} given a bivariate standard normal cumulative
##' distribution with a correlation of \eqn{\rho}
##'
##' Each element of \eqn{\bold{S}_X} can be calculated by
##'
##' \deqn{ \left[\bold{S}_T\right]_{ij} = \sum^{C_i - 1}_{c_i = 1} \sum^{C_j -
##' 1}_{c_j - 1} \Phi_2\left( \tau_{V_{c_i}}, \tau_{V_{c_j}}, \rho^*_{ij}
##' \right) - \sum^{C_i - 1}_{c_i = 1} \Phi_1(\tau_{V_{c_i}}) \sum^{C_j -
##' 1}_{c_j - 1} \Phi_1(\tau_{V_{c_j}}),}
##'
##' where \eqn{\rho^*_{ij}} is a polychoric correlation between Items \eqn{i}
##' and \eqn{j}.
##'
##'
##' @importFrom lavaan lavInspect lavNames
##'
##' @param object A \code{\linkS4class{lavaan}} or
##'   \code{\linkS4class{lavaan.mi}} object, expected to contain only
##'   exogenous common factors (i.e., a CFA model).
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'        imputations from pooled results.  Can include any of
##'        \code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
##'        default setting, which excludes any imputations that did not
##'        converge or for which standard errors could not be computed.  The
##'        last option (\code{"no.npd"}) would exclude any imputations which
##'        yielded a nonpositive definite covariance matrix for observed or
##'        latent variables, which would include any "improper solutions" such
##'        as Heywood cases.  NPD solutions are not excluded by default because
##'        they are likely to occur due to sampling error, especially in small
##'        samples.  However, gross model misspecification could also cause
##'        NPD solutions, users can compare pooled results with and without
##'        this setting as a sensitivity analysis to see whether some
##'        imputations warrant further investigation.
##'
##' @return Maximal reliability values of each group. The maximal-reliability
##'   weights are also provided. Users may extracted the weighted by the
##'   \code{attr} function (see example below).
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' @seealso \code{\link{reliability}} for reliability of an unweighted
##'   composite score
##'
##' @references
##' Li, H. (1997). A unifying expression for the maximal reliability of a linear
##' composite. \emph{Psychometrika, 62}(2), 245--249. \doi{10.1007/BF02295278}
##'
##' Raykov, T. (2012). Scale construction and development using structural
##' equation modeling. In R. H. Hoyle (Ed.), \emph{Handbook of structural
##' equation modeling} (pp. 472--494). New York, NY: Guilford.
##'
##' @examples
##'
##' total <- 'f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 '
##' fit <- cfa(total, data = HolzingerSwineford1939)
##' maximalRelia(fit)
##'
##' # Extract the weight
##' mr <- maximalRelia(fit)
##' attr(mr, "weight")
##'
##' @export
maximalRelia <- function(object, omit.imps = c("no.conv","no.se")) {
  ngroups <- lavInspect(object, "ngroups") #TODO: adapt to multiple levels
  nlevels <- lavInspect(object, "nlevels")
  nblocks <- ngroups*nlevels #FIXME: always true?
  group.label <- if (ngroups > 1L) lavInspect(object, "group.label") else NULL
  #FIXME? lavInspect(object, "level.labels")
  clus.label <- if (nlevels > 1L) c("within", lavInspect(object, "cluster")) else NULL
  if (nblocks > 1L) {
    block.label <- paste(rep(group.label, each = nlevels), clus.label,
                         sep = if (ngroups > 1L && nlevels > 1L) "_" else "")
  }

  ## parameters in GLIST format (not flat, need block-level list)
  if (inherits(object, "lavaan")) {
    param <- lavInspect(object, "est")
    ve <- lavInspect(object, "cov.lv") # model-implied latent covariance matrix
    S <- object@h1$implied$cov # observed sample covariance matrix (already a list)

    if (nblocks == 1L) {
      param <- list(param)
      ve <- list(ve)
    }

  } else if (inherits(object, "lavaan.mi")) {
    useImps <- rep(TRUE, length(object@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

    param <- object@coefList[[ useImps[1] ]] # first admissible as template
    coefList <- object@coefList[useImps]
    phiList <- object@phiList[useImps]
    ## add block-level list per imputation?
    if (nblocks == 1L) {
      param <- list(param)
      for (i in 1:m) {
        coefList[[i]] <- list(coefList[[i]])
        phiList[[i]] <- list(phiList[[i]])
      }
    }
    S <- vector("list", nblocks) # pooled observed covariance matrix
    ve <- vector("list", nblocks)
    ## loop over blocks
    for (b in 1:nblocks) {

      ## param:  loop over GLIST elements
      for (mat in names(param[[b]])) {
        matList <- lapply(coefList, function(i) i[[b]][[mat]])
        param[[b]][[mat]] <- Reduce("+", matList) / length(matList)
      } # mat

      ## pooled observed covariance matrix
      covList <- lapply(object@h1List[useImps], function(i) i$implied$cov[[b]])
      S[[b]] <- Reduce("+", covList) / m

      ## pooled model-implied latent covariance matrix
      ve[[b]] <- Reduce("+", lapply(phiList, "[[", i = b) ) / m

    } # b

  }

  if (nblocks == 1L) {
    SigmaHat <- getMethod("fitted", class(object))(object)["cov"] # retain list format
  } else {
    SigmaHat <- sapply(getMethod("fitted", class(object))(object),
                       "[[", "cov", simplify = FALSE)
  }

  ly <- lapply(param, "[[", "lambda")
  te <- lapply(param, "[[", "theta")

  categorical <- lavInspect(object, "categorical")
  threshold <- if (categorical) getThreshold(object) else NULL

  result <- list()
  for (i in 1:nblocks) {
    truevar <- ly[[i]] %*% ve[[i]] %*% t(ly[[i]])
    varnames <- colnames(truevar)
    if (categorical) {
      invstdvar <- 1 / sqrt(diag(SigmaHat[[i]]))
      polyr <- diag(invstdvar) %*% truevar %*% diag(invstdvar)
      nitem <- ncol(SigmaHat[[i]])
      result[[i]] <- calcMaximalReliaCat(polyr, threshold[[i]], S[[i]], nitem, varnames)
    } else {
      result[[i]] <- calcMaximalRelia(truevar, S[[i]], varnames)
    }
  }

  if (nblocks == 1L) {
    result <- result[[1]]
  } else names(result) <- block.label

  result
}



## ----------------
## Hidden Functions
## ----------------

computeAlpha <- function(S) {
  k <- nrow(S)
  k/(k - 1) * (1.0 - sum(diag(S)) / sum(S))
}

#' @importFrom stats cov2cor pnorm
omegaCat <- function(truevar, threshold, scales, denom) {
  ## must be in standardized latent scale
  R <- diag(scales) %*% truevar %*% diag(scales)

	## denom could be model-implied polychoric correlation assuming diagonal theta,
	##       model-implied polychoric correlation accounting for error covariances,
	##       or "observed" polychoric correlation matrix.
  ## If parameterization="theta", standardize the polychoric coVARIANCE matrix
	denom <- cov2cor(denom)

	nitem <- ncol(denom)
	## initialize sums of cumulative probabilities
	sumnum <- 0 # numerator
	addden <- 0 # denominator
	## loop over all pairs of items
	for (j in 1:nitem) {
  	for (jp in 1:nitem) {
  	  ## initialize sums of cumulative probabilities *per item*
  		sumprobn2 <- 0
  		addprobn2 <- 0
  		## for each pair of items, loop over all their thresholds
  		t1 <- threshold[[j]]  * scales[j] # on standardized latent scale
  		t2 <- threshold[[jp]] * scales[jp]
  		for (c in 1:length(t1)) {
    		for (cp in 1:length(t2)) {
    			sumprobn2 <- sumprobn2 + p2(t1[c], t2[cp], R[j, jp])
    			addprobn2 <- addprobn2 + p2(t1[c], t2[cp], denom[j, jp])
    		}
  		}
  		sumprobn1 <- sum(pnorm(t1))
  		sumprobn1p <- sum(pnorm(t2))
  		sumnum <- sumnum + (sumprobn2 - sumprobn1 * sumprobn1p)
  		addden <- addden + (addprobn2 - sumprobn1 * sumprobn1p)
  	}
	}
	reliab <- sumnum / addden
	reliab
}


p2 <- function(t1, t2, r) {
	mnormt::pmnorm(c(t1, t2), c(0,0), matrix(c(1, r, r, 1), 2, 2))
}


# polycorLavaan <- function(object) {
# 	ngroups <- lavInspect(object, "ngroups")
# 	coef <- lavInspect(object, "est")
# 	targettaunames <- NULL
# 	if (ngroups == 1L) {
# 		targettaunames <- rownames(coef$tau)
# 	} else {
# 		targettaunames <- rownames(coef[[1]]$tau)
# 	}
# 	barpos <- sapply(strsplit(targettaunames, ""), function(x) which(x == "|"))
# 	varnames <- unique(apply(data.frame(targettaunames, barpos - 1), MARGIN = 1,
# 	                         FUN = function(x) substr(x[1], 1, x[2])))
# 	if (length(varnames))
# 	script <- ""
# 	for (i in 2:length(varnames)) {
# 		temp <- paste0(varnames[1:(i - 1)], collapse = " + ")
# 		temp <- paste0(varnames[i], "~~", temp, "\n")
# 		script <- paste(script, temp)
# 	}
# 	newobject <- refit(script, object)
# 	if (ngroups == 1L) {
# 		return(lavInspect(newobject, "est")$theta)
# 	}
# 	lapply(lavInspect(newobject, "est"), "[[", "theta")
# }

##' @importFrom lavaan lavInspect lavNames
getThreshold <- function(object) {
	ngroups <- lavInspect(object, "ngroups") #TODO: add nlevels when capable
	ordnames <- lavNames(object, "ov.ord")
	EST <- lavInspect(object, "est")

	if (ngroups == 1L) {
	  thresholds <- EST$tau[,"threshold"]
	  result <- lapply(ordnames,
	                   function(nn) thresholds[grepl(nn, names(thresholds))])
	  names(result) <- ordnames
	  ## needs to be within a list when called above within block-loops
	  result <- list(result)

	} else {
	  allThr <- EST[which(names(EST) == "tau")]
	  ## convert 1-column matrices to vectors, preserving rownames
	  thresholds <- sapply(allThr, "[", j = "threshold", simplify = FALSE)
	  result <- list()
		group.label <- lavInspect(object, "group.label")

		for (g in 1:ngroups) {
		  result[[ group.label[g] ]] <- lapply(ordnames, function(nn) {
		    thresholds[[g]][ grepl(nn, names(thresholds[[g]])) ]
		  })
		  names(result[[ group.label[g] ]]) <- ordnames
		}

	}

	return(result)
}

##' @importFrom lavaan lavInspect lavNames
getScales <- function(object) {
  ngroups <- lavInspect(object, "ngroups") #TODO: add nlevels when capable
  ordnames <- lavNames(object, "ov.ord") #TODO: use to allow mix of cat/con vars
  EST <- lavInspect(object, "est")

  if (ngroups == 1L) {
    result <- list(EST$delta[,"scales"])
  } else {
    result <- lapply(EST[which(names(EST) == "delta")],
                     function(x) x[,"scales"])
    names(result) <- lavInspect(object, "group.label")
  }

  return(result)
}

invGeneralRelia <- function(w, truevar, totalvar) {
	1 - (t(w) %*% truevar %*% w) / (t(w) %*% totalvar %*% w)
}

#' @importFrom stats pnorm
invGeneralReliaCat <- function(w, polyr, threshold, denom, nitem) {
	# denom could be polychoric correlation, model-implied correlation, or model-implied without error correlation
	upper <- matrix(NA, nitem, nitem)
	lower <- matrix(NA, nitem, nitem)
	for (j in 1:nitem) {
  	for (jp in 1:nitem) {
  		sumprobn2 <- 0
  		addprobn2 <- 0
  		t1 <- threshold[[j]]
  		t2 <- threshold[[jp]]
  		for (c in 1:length(t1)) {
  		for (cp in 1:length(t2)) {
  			sumprobn2 <- sumprobn2 + p2(t1[c], t2[cp], polyr[j, jp])
  			addprobn2 <- addprobn2 + p2(t1[c], t2[cp], denom[j, jp])
  		}
  		}
  		sumprobn1 <- sum(pnorm(t1))
  		sumprobn1p <- sum(pnorm(t2))
  		upper[j, jp] <- (sumprobn2 - sumprobn1 * sumprobn1p)
  		lower[j, jp] <- (addprobn2 - sumprobn1 * sumprobn1p)
  	}
	}
	1 - (t(w) %*% upper %*% w) / (t(w) %*% lower %*% w)
}

#' @importFrom stats nlminb
calcMaximalRelia <- function(truevar, totalvar, varnames) {
	start <- rep(1, nrow(truevar))
	out <- nlminb(start, invGeneralRelia, truevar = truevar, totalvar = totalvar)
	if (out$convergence != 0) stop("The numerical method for finding the maximal",
	                               " reliability did not converge.")
	result <- 1 - out$objective
	weight <- out$par / mean(out$par)
	names(weight) <- varnames
	attr(result, "weight") <- weight
	result
}

#' @importFrom stats nlminb
calcMaximalReliaCat <- function(polyr, threshold, denom, nitem, varnames) {
	start <- rep(1, nrow(polyr))
	out <- nlminb(start, invGeneralReliaCat, polyr = polyr, threshold = threshold,
	              denom = denom, nitem = nitem)
	if (out$convergence != 0) stop("The numerical method for finding the maximal",
	                               " reliability did not converge.")
	result <- 1 - out$objective
	weight <- out$par / mean(out$par)
	names(weight) <- varnames
	attr(result, "weight") <- weight
	result
}


