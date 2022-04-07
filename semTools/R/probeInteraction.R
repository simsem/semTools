### Sunthud Pornprasertmanit & Terrence D. Jorgensen
### Last updated: 7 April 2022



## --------
## 2-way MC
## --------

##' Probing two-way interaction on the no-centered or mean-centered latent
##' interaction
##'
##' Probing interaction for simple intercept and simple slope for the
##' no-centered or mean-centered latent two-way interaction
##'
##' Before using this function, researchers need to make the products of the
##' indicators between the first-order factors using mean centering (Marsh, Wen,
##' & Hau, 2004). Note that the double-mean centering may not be appropriate for
##' probing interaction if researchers are interested in simple intercepts. The
##' mean or double-mean centering can be done by the \code{\link{indProd}}
##' function. The indicator products can be made for all possible combination or
##' matched-pair approach (Marsh et al., 2004). Next, the hypothesized model
##' with the regression with latent interaction will be used to fit all original
##' indicators and the product terms. See the example for how to fit the product
##' term below. Once the lavaan result is obtained, this function will be used
##' to probe the interaction.
##'
##' Let that the latent interaction model regressing the dependent variable
##' (\eqn{Y}) on the independent varaible (\eqn{X}) and the moderator (\eqn{Z})
##' be \deqn{ Y = b_0 + b_1X + b_2Z + b_3XZ + r, } where \eqn{b_0} is the
##' estimated intercept or the expected value of \eqn{Y} when both \eqn{X} and
##' \eqn{Z} are 0, \eqn{b_1} is the effect of \eqn{X} when \eqn{Z} is 0,
##' \eqn{b_2} is the effect of \eqn{Z} when \eqn{X} is 0, \eqn{b_3} is the
##' interaction effect between \eqn{X} and \eqn{Z}, and \eqn{r} is the residual
##' term.
##'
##' For probing two-way interaction, the simple intercept of the independent
##' variable at each value of the moderator (Aiken & West, 1991; Cohen, Cohen,
##' West, & Aiken, 2003; Preacher, Curran, & Bauer, 2006) can be obtained by
##' \deqn{ b_{0|X = 0, Z} = b_0 + b_2Z. }
##'
##' The simple slope of the independent varaible at each value of the moderator
##' can be obtained by \deqn{ b_{X|Z} = b_1 + b_3Z. }
##'
##' The variance of the simple intercept formula is \deqn{ Var\left(b_{0|X = 0,
##' Z}\right) = Var\left(b_0\right) + 2ZCov\left(b_0, b_2\right) +
##' Z^2Var\left(b_2\right) } where \eqn{Var} denotes the variance of a parameter
##' estimate and \eqn{Cov} denotes the covariance of two parameter estimates.
##'
##' The variance of the simple slope formula is \deqn{ Var\left(b_{X|Z}\right) =
##' Var\left(b_1\right) + 2ZCov\left(b_1, b_3\right) + Z^2Var\left(b_3\right) }
##'
##' Wald \emph{z} statistic is used for test statistic (even for objects of
##' class \code{\linkS4class{lavaan.mi}}).
##'
##'
##' @importFrom lavaan lavInspect parTable
##' @importFrom stats pnorm
##' @importFrom methods getMethod
##'
##' @param fit A fitted \code{\linkS4class{lavaan}} or
##'   \code{\linkS4class{lavaan.mi}} object with a latent 2-way interaction.
##' @param nameX \code{character} vector of all 3 factor names used as the
##'   predictors. The lower-order factors must be listed first, and the final
##'   name must be the latent interaction factor.
##' @param nameY The name of factor that is used as the dependent variable.
##' @param modVar The name of factor that is used as a moderator. The effect of
##'   the other independent factor will be probed at each value of the
##'   moderator variable listed in \code{valProbe}.
##' @param valProbe The values of the moderator that will be used to probe the
##'   effect of the focal predictor.
##' @param group In multigroup models, the label of the group for which the
##'   results will be returned. Must correspond to one of
##'   \code{\link[lavaan]{lavInspect}(fit, "group.label")}, or an integer
##'   corresponding to which of those group labels.
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'        imputations from pooled results. Ignored unless \code{fit} is of
##'        class \code{\linkS4class{lavaan.mi}}. Can include any of
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
##' @return A list with two elements:
##' \enumerate{
##'  \item \code{SimpleIntercept}: The intercepts given each value of the
##'   moderator. This element will be \code{NULL} unless the factor intercept is
##'   estimated (e.g., not fixed at 0).
##'  \item \code{SimpleSlope}: The slopes given each value of the moderator.
##' }
##' In each element, the first column represents the values of the moderators
##' specified in the \code{valProbe} argument. The second column is the simple
##' intercept or simple slope. The third column is the \emph{SE} of the simple
##' intercept or simple slope. The fourth column is the Wald (\emph{z})
##' statistic. The fifth column is the \emph{p} value testing whether the simple
##' intercepts or slopes are different from 0.
##'
##' @author
##' Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @seealso \itemize{
##'  \item \code{\link{indProd}} For creating the indicator products with no
##'   centering, mean centering, double-mean centering, or residual centering.
##'  \item \code{\link{probe3WayMC}} For probing the three-way latent interaction
##'   when the results are obtained from mean-centering, or double-mean centering
##'  \item \code{\link{probe2WayRC}} For probing the two-way latent interaction
##'   when the results are obtained from residual-centering approach.
##'  \item \code{\link{probe3WayRC}} For probing the two-way latent interaction
##'   when the results are obtained from residual-centering approach.
##'  \item \code{\link{plotProbe}} Plot the simple intercepts and slopes of the
##'   latent interaction.
##' }
##'
##' @references
##'
##' Tutorial:
##'
##' Schoemann, A. M., & Jorgensen, T. D. (2021). Testing and interpreting
##' latent variable interactions using the \code{semTools} package.
##' \emph{Psych, 3}(3), 322--335. \doi{10.3390/psych3030024}
##'
##' Background literature:
##'
##' Aiken, L. S., & West, S. G. (1991). \emph{Multiple regression: Testing
##' and interpreting interactions}. Newbury Park, CA: Sage.
##'
##' Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). \emph{Applied
##' multiple regression/correlation analysis for the behavioral sciences}
##' (3rd ed.). New York, NY: Routledge.
##'
##' Marsh, H. W., Wen, Z., & Hau, K. T. (2004). Structural equation models of
##' latent interactions: Evaluation of alternative estimation strategies and
##' indicator construction. \emph{Psychological Methods, 9}(3), 275--300.
##' \doi{10.1037/1082-989X.9.3.275}
##'
##' Preacher, K. J., Curran, P. J., & Bauer, D. J. (2006). Computational tools
##' for probing interactions in multiple linear regression, multilevel modeling,
##' and latent curve analysis. \emph{Journal of Educational and Behavioral
##' Statistics, 31}(4), 437--448. \doi{10.3102/10769986031004437}
##'
##' @examples
##'
##' dat2wayMC <- indProd(dat2way, 1:3, 4:6)
##'
##' model1 <- "
##' f1 =~ x1 + x2 + x3
##' f2 =~ x4 + x5 + x6
##' f12 =~ x1.x4 + x2.x5 + x3.x6
##' f3 =~ x7 + x8 + x9
##' f3 ~ f1 + f2 + f12
##' f12 ~~ 0*f1 + 0*f2
##' x1 + x4 + x1.x4 + x7 ~ 0*1 # identify latent means
##' f1 + f2 + f12 + f3 ~ NA*1
##' "
##'
##' fitMC2way <- sem(model1, data = dat2wayMC, meanstructure = TRUE)
##' summary(fitMC2way)
##'
##' probe2WayMC(fitMC2way, nameX = c("f1", "f2", "f12"), nameY = "f3",
##'             modVar = "f2", valProbe = c(-1, 0, 1))
##'
##'
##' ## can probe multigroup models, one group at a time
##' dat2wayMC$g <- 1:2
##'
##' model2 <- "
##' f1 =~ x1 + x2 + x3
##' f2 =~ x4 + x5 + x6
##' f12 =~ x1.x4 + x2.x5 + x3.x6
##' f3 =~ x7 + x8 + x9
##' f3 ~ c(b1.g1, b1.g2)*f1 + c(b2.g1, b2.g2)*f2 + c(b12.g1, b12.g2)*f12
##' f12 ~~ 0*f1 + 0*f2
##' x1 + x4 + x1.x4 + x7 ~ 0*1 # identify latent means
##' f1 + f2 + f12 ~ NA*1
##' f3 ~ NA*1 + c(b0.g1, b0.g2)*1
##' "
##' fit2 <- sem(model2, data = dat2wayMC, group = "g")
##' probe2WayMC(fit2, nameX = c("f1", "f2", "f12"), nameY = "f3",
##'             modVar = "f2", valProbe = c(-1, 0, 1)) # group = 1 by default
##' probe2WayMC(fit2, nameX = c("f1", "f2", "f12"), nameY = "f3",
##'             modVar = "f2", valProbe = c(-1, 0, 1), group = 2)
##'
##' @export
probe2WayMC <- function(fit, nameX, nameY, modVar, valProbe, group = 1L,
                        omit.imps = c("no.conv","no.se")) {
  ## TDJ: verify class
  if (inherits(fit, "lavaan")) {
    est <- lavInspect(fit, "est")[[group]]

  } else if (inherits(fit, "lavaan.mi")) {
    useImps <- rep(TRUE, length(fit@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(fit@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(fit@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(fit@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(fit@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    ## custom removal by imputation number
    rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
    if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
    ## whatever is left
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

  } else stop('"fit" must inherit from lavaan or lavaan.mi class', call. = FALSE)


	# Check whether modVar is correct
	if (is.character(modVar)) modVar <- match(modVar, nameX)
	if (is.na(modVar) || !(modVar %in% 1:2))
	  stop("The moderator name is not in the name of independent factors or not 1 or 2.")

	## TDJ: If multigroup, check group %in% group.label
	nG <- lavInspect(fit, "ngroups")
	if (nG > 1L) {
	  group.label <- lavInspect(fit, "group.label")
	  ## assign numeric to character
	  if (is.numeric(group)) {
	    if (group %in% 1:nG) {
	      group <- group.label[group]
	    } else group <- as.character(group)
	  } else group <- as.character(group)
	  ## check that character is a group
	  if (!as.character(group) %in% group.label)
	    stop('"group" must be a character string naming a group of interest, or ',
	         'an integer corresponding to a group in  lavInspect(fit, "group.label")')
	  group.number <- which(group.label == group)
	} else group.number <- 1L

	# Extract all varEst
	if (inherits(fit, "lavaan")) {
	  varEst <- lavaan::vcov(fit)
	} else if (inherits(fit, "lavaan.mi")) {
	  varEst <- getMethod("vcov", "lavaan.mi")(fit, omit.imps = omit.imps)
	}

	## Check whether the outcome's intercept is estimated
	PT <- parTable(fit)
	if (lavInspect(fit, "options")$meanstructure) {
	  targetcol <- PT$label[PT$lhs == nameY & PT$op == "~1" & PT$group == group.number]
	  if (targetcol == "") {
	    ## no custom label, use default
	    targetcol <- paste0(nameY, "~1")
	    if (nG > 1L && group.number > 1L) {
	      targetcol <- paste0(targetcol, ".g", group.number)
	    }
	  }
	  ## check it is actually estimated (thus, has sampling variance)
	  estimateIntcept <- targetcol %in% rownames(varEst)

	} else estimateIntcept <- FALSE


	## Get the parameter estimates for that group
	if (nG > 1L) {

	  if (inherits(fit, "lavaan")) {
	    est <- lavInspect(fit, "est")[[group]]
	  } else if (inherits(fit, "lavaan.mi")) {
	    est <- list()
	    GLIST <- fit@coefList[useImps]
	    est$beta <- Reduce("+", lapply(GLIST, function(i) i[[group]]$beta)) / m
	    if (estimateIntcept) {
	      est$alpha <- Reduce("+", lapply(GLIST, function(i) i[[group]]$alpha)) / m
	    }
	  }

	} else {
	  ## single-group model

	  if (inherits(fit, "lavaan")) {
	    est <- lavInspect(fit, "est")
	  } else if (inherits(fit, "lavaan.mi")) {
	    est <- list()
	    est$beta <- Reduce("+", lapply(fit@coefList[useImps], "[[", i = "beta")) / m
	    if (estimateIntcept) {
	      est$alpha <- Reduce("+", lapply(fit@coefList[useImps], "[[", i = "alpha")) / m
	    }
	  }

	}

	# Compute the intercept of no-centering
	betaNC <- matrix(est$beta[nameY, nameX], ncol = 1,
	                 dimnames = list(nameX, nameY))

	pvalue <- function(x) (1 - pnorm(abs(x))) * 2

	resultIntcept <- NULL
	resultSlope <- NULL
	if (estimateIntcept) {
		# Extract SE from centered result
	  newRows <- which(PT$lhs == nameY & PT$op == "~" & PT$rhs %in% nameX & PT$group == group.number)
	  newLabels <- PT$label[newRows]
	  if (any(newLabels == "")) for (i in which(newLabels == "")) {
	    newLabels[i] <- paste0(nameY, "~", nameX[i],
	                           ifelse(nG > 1L && group.number > 1L, no = "",
	                                  yes = paste0(".g", group.number)))
	  }
		targetcol <- c(targetcol, newLabels)

		# Transform it to non-centering SE
		usedVar <-  varEst[targetcol, targetcol]
		usedBeta <- rbind(est$alpha[nameY,], betaNC)

		# Change the order of usedVar and usedBeta if the moderator variable is listed first
		if (modVar == 1) {
			usedVar <- usedVar[c(1, 3, 2, 4), c(1, 3, 2, 4)]
			usedBeta <- usedBeta[c(1, 3, 2, 4)]
		}

		# Find simple intercept
		simpleIntcept <- usedBeta[1] + usedBeta[3] * valProbe
		varIntcept <- usedVar[1, 1] + 2 * valProbe * usedVar[1, 3] + (valProbe^2) * usedVar[3, 3]
		zIntcept <- simpleIntcept/sqrt(varIntcept)
		pIntcept <- pvalue(zIntcept)
		resultIntcept <- data.frame(valProbe, simpleIntcept, sqrt(varIntcept), zIntcept, pIntcept)
		colnames(resultIntcept) <- c(nameX[modVar], "est", "se", "z", "pvalue")
		class(resultIntcept) <- c("lavaan.data.frame","data.frame")

		# Find simple slope
		simpleSlope <- usedBeta[2] + usedBeta[4] * valProbe
		varSlope <- usedVar[2, 2] + 2 * valProbe * usedVar[2, 4] + (valProbe^2) * usedVar[4, 4]
		zSlope <- simpleSlope / sqrt(varSlope)
		pSlope <- pvalue(zSlope)
		resultSlope <- data.frame(valProbe, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "est", "se", "z", "pvalue")
		class(resultSlope) <- c("lavaan.data.frame","data.frame")

	} else {
	  newRows <- which(PT$lhs == nameY & PT$op == "~" & PT$rhs %in% nameX & PT$group == group.number)
	  targetcol <- PT$label[newRows]
	  if (any(targetcol == "")) for (i in which(targetcol == "")) {
      targetcol[i] <- paste0(nameY, "~", PT$rhs[ newRows[i] ],
                             ifelse(nG > 1L && group.number > 1L, no = "",
                                    yes = paste0(".g", group.number)))
	  }

		# Transform it to non-centering SE
		usedVar <-  varEst[targetcol, targetcol]
		usedBeta <- betaNC

		# Change the order of usedVar and usedBeta if the moderator variable is listed first
		if(modVar == 1) {
			usedVar <- usedVar[c(2, 1, 3), c(2, 1, 3)]
			usedBeta <- usedBeta[c(2, 1, 3)]
		}

		# Find simple slope
		simpleSlope <- usedBeta[1] + usedBeta[3] * valProbe
		varSlope <- usedVar[1, 1] + 2 * valProbe * usedVar[1, 3] + (valProbe^2) * usedVar[3, 3]
		zSlope <- simpleSlope/sqrt(varSlope)
		pSlope <- pvalue(zSlope)
		resultSlope <- data.frame(valProbe, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "est", "se", "z", "pvalue")
		class(resultSlope) <- c("lavaan.data.frame","data.frame")
	}

	list(SimpleIntcept = resultIntcept, SimpleSlope = resultSlope)
}



## --------
## 2-way RC
## --------

##' Probing two-way interaction on the residual-centered latent interaction
##'
##' Probing interaction for simple intercept and simple slope for the
##' residual-centered latent two-way interaction (Geldhof et al., 2013)
##'
##' Before using this function, researchers need to make the products of the
##' indicators between the first-order factors and residualize the products by
##' the original indicators (Lance, 1988; Little, Bovaird, & Widaman, 2006). The
##' process can be automated by the \code{\link{indProd}} function. Note that
##' the indicator products can be made for all possible combination or
##' matched-pair approach (Marsh et al., 2004). Next, the hypothesized model
##' with the regression with latent interaction will be used to fit all original
##' indicators and the product terms. To use this function the model must be fit
##' with a mean structure. See the example for how to fit the product term
##' below. Once the lavaan result is obtained, this function will be used to
##' probe the interaction.
##'
##' The probing process on residual-centered latent interaction is based on
##' transforming the residual-centered result into the no-centered result. See
##' Geldhof et al. (2013) for further details. Note that this approach based on
##' a strong assumption that the first-order latent variables are normally
##' distributed. The probing process is applied after the no-centered result
##' (parameter estimates and their covariance matrix among parameter estimates)
##' has been computed. See the \code{\link{probe2WayMC}} for further details.
##'
##'
##' @importFrom lavaan lavInspect parTable
##' @importFrom stats pnorm
##' @importFrom methods getMethod
##'
##' @param fit A fitted \code{\linkS4class{lavaan}} or
##'   \code{\linkS4class{lavaan.mi}} object with a latent 2-way interaction.
##' @param nameX \code{character} vector of all 3 factor names used as the
##'   predictors. The lower-order factors must be listed first, and the final
##'   name must be the latent interaction factor.
##' @param nameY The name of factor that is used as the dependent variable.
##' @param modVar The name of factor that is used as a moderator. The effect of
##'   the other independent factor will be probed at each value of the
##'   moderator variable listed in \code{valProbe}.
##' @param valProbe The values of the moderator that will be used to probe the
##'   effect of the focal predictor.
##' @param group In multigroup models, the label of the group for which the
##'   results will be returned. Must correspond to one of
##'   \code{\link[lavaan]{lavInspect}(fit, "group.label")}, or an integer
##'   corresponding to which of those group labels.
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'        imputations from pooled results. Ignored unless \code{fit} is of
##'        class \code{\linkS4class{lavaan.mi}}. Can include any of
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
##' @return A list with two elements:
##' \enumerate{
##'  \item \code{SimpleIntercept}: The intercepts given each value of the
##'   moderator. This element will be \code{NULL} unless the factor intercept is
##'   estimated (e.g., not fixed at 0).
##'  \item \code{SimpleSlope}: The slopes given each value of the moderator.
##' }
##' In each element, the first column represents the values of the moderators
##' specified in the \code{valProbe} argument. The second column is the simple
##' intercept or simple slope. The third column is the standard error of the
##' simple intercept or simple slope. The fourth column is the Wald (\emph{z})
##' statistic. The fifth column is the \emph{p} value testing whether the simple
##' intercepts or slopes are different from 0.
##'
##' @author
##' Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @seealso \itemize{
##'  \item \code{\link{indProd}} For creating the indicator products with no
##'   centering, mean centering, double-mean centering, or residual centering.
##'  \item \code{\link{probe2WayMC}} For probing the two-way latent interaction
##'   when the results are obtained from mean-centering, or double-mean centering
##'  \item \code{\link{probe3WayMC}} For probing the three-way latent interaction
##'   when the results are obtained from mean-centering, or double-mean centering
##'  \item \code{\link{probe3WayRC}} For probing the two-way latent interaction
##'   when the results are obtained from residual-centering approach.
##'  \item \code{\link{plotProbe}} Plot the simple intercepts and slopes of the
##'   latent interaction.
##' }
##' @references
##'
##' Tutorial:
##'
##' Schoemann, A. M., & Jorgensen, T. D. (2021). Testing and interpreting
##' latent variable interactions using the \code{semTools} package.
##' \emph{Psych, 3}(3), 322--335. \doi{10.3390/psych3030024}
##'
##' Background literature:
##'
##' Lance, C. E. (1988). Residual centering, exploratory and confirmatory
##' moderator analysis, and decomposition of effects in path models containing
##' interactions. \emph{Applied Psychological Measurement, 12}(2), 163--175.
##' \doi{10.1177/014662168801200205}
##'
##' Little, T. D., Bovaird, J. A., & Widaman, K. F. (2006). On the merits of
##' orthogonalizing powered and product terms: Implications for modeling
##' interactions. \emph{Structural Equation Modeling, 13}(4), 497--519.
##' \doi{10.1207/s15328007sem1304_1}
##'
##' Marsh, H. W., Wen, Z., & Hau, K. T. (2004). Structural equation models of
##' latent interactions: Evaluation of alternative estimation strategies and
##' indicator construction. \emph{Psychological Methods, 9}(3), 275--300.
##' \doi{10.1037/1082-989X.9.3.275}
##'
##' Geldhof, G. J., Pornprasertmanit, S., Schoemann, A. M., & Little, T. D.
##' (2013). Orthogonalizing through residual centering: Extended applications
##' and caveats. \emph{Educational and Psychological Measurement, 73}(1), 27--46.
##' \doi{10.1177/0013164412445473}
##'
##' @examples
##'
##' dat2wayRC <- orthogonalize(dat2way, 1:3, 4:6)
##'
##' model1 <- "
##' f1 =~ x1 + x2 + x3
##' f2 =~ x4 + x5 + x6
##' f12 =~ x1.x4 + x2.x5 + x3.x6
##' f3 =~ x7 + x8 + x9
##' f3 ~ f1 + f2 + f12
##' f12 ~~ 0*f1 + 0*f2
##' x1 + x4 + x1.x4 + x7 ~ 0*1 # identify latent means
##' f1 + f2 + f12 + f3 ~ NA*1
##' "
##'
##' fitRC2way <- sem(model1, data = dat2wayRC, meanstructure = TRUE)
##' summary(fitRC2way)
##'
##' probe2WayRC(fitRC2way, nameX = c("f1", "f2", "f12"), nameY = "f3",
##'             modVar = "f2", valProbe = c(-1, 0, 1))
##'
##'
##' ## can probe multigroup models, one group at a time
##' dat2wayRC$g <- 1:2
##'
##' model2 <- "
##' f1  =~ x1 + x2 + x3
##' f2  =~ x4 + x5 + x6
##' f12 =~ x1.x4 + x2.x5 + x3.x6
##' f3  =~ x7 + x8 + x9
##' f3 ~ c(b1.g1, b1.g2)*f1 + c(b2.g1, b2.g2)*f2 + c(b12.g1, b12.g2)*f12
##' f12 ~~ 0*f1 + 0*f2
##' x1 + x4 + x1.x4 + x7 ~ 0*1 # identify latent means
##' f1 + f2 + f12 ~ NA*1
##' f3 ~ NA*1 + c(b0.g1, b0.g2)*1
##' "
##' fit2 <- sem(model2, data = dat2wayRC, group = "g")
##' probe2WayRC(fit2, nameX = c("f1", "f2", "f12"), nameY = "f3",
##'             modVar = "f2", valProbe = c(-1, 0, 1)) # group = 1 by default
##' probe2WayRC(fit2, nameX = c("f1", "f2", "f12"), nameY = "f3",
##'             modVar = "f2", valProbe = c(-1, 0, 1), group = 2)
##'
##' @export
probe2WayRC <- function(fit, nameX, nameY, modVar, valProbe, group = 1L,
                        omit.imps = c("no.conv","no.se")) {
  ## TDJ: verify class
  if (inherits(fit, "lavaan")) {
    est <- lavInspect(fit, "est")[[group]]

  } else if (inherits(fit, "lavaan.mi")) {
    useImps <- rep(TRUE, length(fit@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(fit@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(fit@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(fit@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(fit@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    ## custom removal by imputation number
    rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
    if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
    ## whatever is left
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

  } else stop('"fit" must inherit from lavaan or lavaan.mi class', call. = FALSE)

  if (!lavInspect(fit, "options")$meanstructure)
    stop('This function requires the model to be fit with a mean structure.',
         call. = FALSE)


  # Check whether modVar is correct
	if (is.character(modVar))	modVar <- match(modVar, nameX)
	if (is.na(modVar) || !(modVar %in% 1:2))
	  stop("The moderator name is not in the name of independent factors or not 1 or 2.")

	## TDJ: If multigroup, check group %in% group.label
	nG <- lavInspect(fit, "ngroups")
	if (nG > 1L) {
	  group.label <- lavInspect(fit, "group.label")
	  ## assign numeric to character
	  if (is.numeric(group)) {
	    if (group %in% 1:nG) {
	      group <- group.label[group]
	    } else group <- as.character(group)
	  } else group <- as.character(group)
	  ## check that character is a group
	  if (!as.character(group) %in% group.label)
	    stop('"group" must be a character string naming a group of interest, or ',
	         'an integer corresponding to a group in  lavInspect(fit, "group.label")')
	  group.number <- which(group.label == group)
	} else group.number <- 1L

	# Extract all varEst
	if (inherits(fit, "lavaan")) {
	  varEst <- lavaan::vcov(fit)
	} else if (inherits(fit, "lavaan.mi")) {
	  varEst <- getMethod("vcov", "lavaan.mi")(fit, omit.imps = omit.imps)
	}

	## Check whether the outcome's intercept is estimated
	PT <- parTable(fit)
	if (lavInspect(fit, "options")$meanstructure) {
	  targetcol <- PT$label[PT$lhs == nameY & PT$op == "~1" & PT$group == group.number]
	  if (targetcol == "") {
	    ## no custom label, use default
	    targetcol <- paste0(nameY, "~1")
	    if (nG > 1L && group.number > 1L) {
	      targetcol <- paste0(targetcol, ".g", group.number)
	    }
	  }
	  ## check it is actually estimated (thus, has sampling variance)
	  estimateIntcept <- targetcol %in% rownames(varEst)

	} else estimateIntcept <- FALSE


	## Get the parameter estimates for that group
	if (nG > 1L) {

	  if (inherits(fit, "lavaan")) {
	    est <- lavInspect(fit, "est")[[group]]
	  } else if (inherits(fit, "lavaan.mi")) {
	    est <- list()
	    GLIST <- fit@coefList[useImps]
	    est$beta <- Reduce("+", lapply(GLIST, function(i) i[[group]]$beta)) / m
	    est$alpha <- Reduce("+", lapply(GLIST, function(i) i[[group]]$alpha)) / m
	    est$psi <- Reduce("+", lapply(GLIST, function(i) i[[group]]$psi)) / m
	  }

	} else {
	  ## single-group model

	  if (inherits(fit, "lavaan")) {
	    est <- lavInspect(fit, "est")
	  } else if (inherits(fit, "lavaan.mi")) {
	    est <- list()
	    est$beta <- Reduce("+", lapply(fit@coefList[useImps], "[[", i = "beta")) / m
	    est$alpha <- Reduce("+", lapply(fit@coefList[useImps], "[[", i = "alpha")) / m
	    est$psi <- Reduce("+", lapply(fit@coefList[useImps], "[[", i = "psi")) / m
	  }

	}

	# Find the mean and covariance matrix of independent factors
	varX <- est$psi[nameX, nameX]
	meanX <- matrix(est$alpha[nameX,], ncol = 1, dimnames = list(NULL, "intcept"))

	# Find the intercept, regression coefficients, and residual variance of residual-centered regression
	intceptRC <- est$alpha[nameY,]
	resVarRC <- est$psi[nameY, nameY]
	betaRC <- matrix(est$beta[nameY, nameX], ncol = 1,
	                 dimnames = list(nameX, nameY))

	# Find the number of observations
	numobs <- lavInspect(fit, "nobs")[group.number]

	# Compute SSRC
	meanXwith1 <- rbind(1, meanX)
	varXwith0 <- cbind(0, rbind(0, varX))
	SSRC <- numobs * (varXwith0 + (meanXwith1 %*% t(meanXwith1)))

	# Compute Mean(Y) and Var(Y)
	betaRCWithIntcept <- rbind(intceptRC, betaRC)
	meanY <- t(meanXwith1) %*% betaRCWithIntcept
	varY <- (t(betaRCWithIntcept) %*% SSRC %*% betaRCWithIntcept)/numobs - meanY^2 + resVarRC

	# Compute Cov(Y, X)
	covY <- as.matrix((varX %*% betaRC)[1:2,])

	# Compute E(XZ)
	meanX[3] <- meanX[1] * meanX[2] + varX[1, 2]

	# Compute Var(XZ)
	varX[3, 3] <- meanX[1]^2 * varX[2, 2] + meanX[2]^2 * varX[1, 1] + 2 * meanX[1] * meanX[2] * varX[1, 2] + varX[1, 1] * varX[2, 2] + varX[1, 2]^2

	# Compute Cov(X, XZ), Cov(Z, XZ)
	varX[1, 3] <- varX[3, 1] <- meanX[1] * varX[1, 2] + meanX[2] * varX[1, 1]
	varX[2, 3] <- varX[3, 2] <- meanX[1] * varX[2, 2] + meanX[2] * varX[1, 2]

	# Compute Cov(Y, XZ) and regression coefficients of no-centering
	betaNC <- solve(varX[1:2,1:2], covY - rbind(varX[1,3] * betaRC[3,1], varX[2, 3] * betaRC[3,1]))
	betaNC <- rbind(betaNC, betaRC[3, 1])
	covY <- rbind(covY, (varX %*% betaNC)[3, 1])

	# Aggregate the non-centering sufficient statistics (Just show how to do but not necessary)
	fullCov <- rbind(cbind(varX, covY), c(covY, varY))
	fullMean <- rbind(meanX, meanY)

	# Compute the intercept of no-centering
	intceptNC <- meanY - t(betaNC) %*% meanX

	# Compute SSNC
	betaNCWithIntcept <- rbind(intceptNC, betaNC)
	meanXwith1 <- rbind(1, meanX)
	varXwith0 <- rbind(0, cbind(0, varX))
	SSNC <- numobs * (varXwith0 + (meanXwith1 %*% t(meanXwith1)))

	# Compute residual variance on non-centering
	resVarNC <- varY - (t(betaNCWithIntcept) %*% SSNC %*% betaNCWithIntcept)/numobs + meanY^2

	pvalue <- function(x) (1 - pnorm(abs(x))) * 2

	resultIntcept <- NULL
	resultSlope <- NULL
	if (estimateIntcept) {
		# Extract SE from residual centering
	  newRows <- which(PT$lhs == nameY & PT$op == "~" & PT$rhs %in% nameX & PT$group == group.number)
	  newLabels <- PT$label[newRows]
	  if (any(newLabels == "")) for (i in which(newLabels == "")) {
	    newLabels[i] <- paste0(nameY, "~", nameX[i],
	                           ifelse(nG > 1L && group.number > 1L, no = "",
	                                  yes = paste0(".g", group.number)))
	  }
	  targetcol <- c(targetcol, newLabels)
		varEstSlopeRC <- varEst[targetcol, targetcol]

		# Transform it to non-centering SE
		usedVar <-  as.numeric(resVarNC/resVarRC) * (varEstSlopeRC %*% SSRC %*% solve(SSNC))
		usedBeta <- betaNCWithIntcept

		# Change the order of usedVar and usedBeta if the moderator variable is listed first
		if (modVar == 1) {
			usedVar <- usedVar[c(1, 3, 2, 4), c(1, 3, 2, 4)]
			usedBeta <- usedBeta[c(1, 3, 2, 4)]
		}

		# Find simple intercept
		simpleIntcept <- usedBeta[1] + usedBeta[3] * valProbe
		varIntcept <- usedVar[1, 1] + 2 * valProbe * usedVar[1, 3] + (valProbe^2) * usedVar[3, 3]
		zIntcept <- simpleIntcept/sqrt(varIntcept)
		pIntcept <- pvalue(zIntcept)
		resultIntcept <- data.frame(valProbe, simpleIntcept, sqrt(varIntcept), zIntcept, pIntcept)
		colnames(resultIntcept) <- c(nameX[modVar], "est", "se", "z", "pvalue")
		class(resultIntcept) <- c("lavaan.data.frame","data.frame")

		# Find simple slope
		simpleSlope <- usedBeta[2] + usedBeta[4] * valProbe
		varSlope <- usedVar[2, 2] + 2 * valProbe * usedVar[2, 4] + (valProbe^2) * usedVar[4, 4]
		zSlope <- simpleSlope/sqrt(varSlope)
		pSlope <- pvalue(zSlope)
		resultSlope <- data.frame(valProbe, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "est", "se", "z", "pvalue")
		class(resultSlope) <- c("lavaan.data.frame","data.frame")

	} else {
	  newRows <- which(PT$lhs == nameY & PT$op == "~" & PT$rhs %in% nameX & PT$group == group.number)
	  targetcol <- PT$label[newRows]
	  if (any(targetcol == "")) for (i in which(targetcol == "")) {
	    targetcol[i] <- paste0(nameY, "~", PT$rhs[ newRows[i] ],
	                           ifelse(nG > 1L && group.number > 1L, no = "",
	                                  yes = paste0(".g", group.number)))
	  }
	  varEstSlopeRC <- varEst[targetcol, targetcol]

		# Transform it to non-centering SE
		usedVar <-  as.numeric(resVarNC/resVarRC) * (varEstSlopeRC %*% SSRC[2:4, 2:4] %*% solve(SSNC[2:4, 2:4]))
		usedBeta <- betaNC

		# Change the order of usedVar and usedBeta if the moderator variable is listed first
		if (modVar == 1) {
			usedVar <- usedVar[c(2, 1, 3), c(2, 1, 3)]
			usedBeta <- usedBeta[c(2, 1, 3)]
		}

		# Find simple slope
		simpleSlope <- usedBeta[1] + usedBeta[3] * valProbe
		varSlope <- usedVar[1, 1] + 2 * valProbe * usedVar[1, 3] + (valProbe^2) * usedVar[3, 3]
		zSlope <- simpleSlope/sqrt(varSlope)
		pSlope <- pvalue(zSlope)
		resultSlope <- data.frame(valProbe, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "est", "se", "z", "pvalue")
		class(resultSlope) <- c("lavaan.data.frame","data.frame")
	}

	list(SimpleIntcept = resultIntcept, SimpleSlope = resultSlope)
}



## --------
## 3-way MC
## --------

##' Probing three-way interaction on the no-centered or mean-centered latent
##' interaction
##'
##' Probing interaction for simple intercept and simple slope for the
##' no-centered or mean-centered latent two-way interaction
##'
##' Before using this function, researchers need to make the products of the
##' indicators between the first-order factors using mean centering (Marsh, Wen,
##' & Hau, 2004). Note that the double-mean centering may not be appropriate for
##' probing interaction if researchers are interested in simple intercepts. The
##' mean or double-mean centering can be done by the \code{\link{indProd}}
##' function. The indicator products can be made for all possible combination or
##' matched-pair approach (Marsh et al., 2004). Next, the hypothesized model
##' with the regression with latent interaction will be used to fit all original
##' indicators and the product terms. See the example for how to fit the product
##' term below. Once the lavaan result is obtained, this function will be used
##' to probe the interaction.
##'
##' Let that the latent interaction model regressing the dependent variable
##' (\eqn{Y}) on the independent varaible (\eqn{X}) and two moderators (\eqn{Z}
##' and \eqn{W}) be \deqn{ Y = b_0 + b_1X + b_2Z + b_3W + b_4XZ + b_5XW + b_6ZW
##' + b_7XZW + r, } where \eqn{b_0} is the estimated intercept or the expected
##' value of \eqn{Y} when \eqn{X}, \eqn{Z}, and \eqn{W} are 0, \eqn{b_1} is the
##' effect of \eqn{X} when \eqn{Z} and \eqn{W} are 0, \eqn{b_2} is the effect of
##' \eqn{Z} when \eqn{X} and \eqn{W} is 0, \eqn{b_3} is the effect of \eqn{W}
##' when \eqn{X} and \eqn{Z} are 0, \eqn{b_4} is the interaction effect between
##' \eqn{X} and \eqn{Z} when \eqn{W} is 0, \eqn{b_5} is the interaction effect
##' between \eqn{X} and \eqn{W} when \eqn{Z} is 0, \eqn{b_6} is the interaction
##' effect between \eqn{Z} and \eqn{W} when \eqn{X} is 0, \eqn{b_7} is the
##' three-way interaction effect between \eqn{X}, \eqn{Z}, and \eqn{W}, and
##' \eqn{r} is the residual term.
##'
##' For probing three-way interaction, the simple intercept of the independent
##' variable at the specific values of the moderators (Aiken & West, 1991) can
##' be obtained by \deqn{ b_{0|X = 0, Z, W} = b_0 + b_2Z + b_3W + b_6ZW. }
##'
##' The simple slope of the independent varaible at the specific values of the
##' moderators can be obtained by \deqn{ b_{X|Z, W} = b_1 + b_3Z + b_4W + b_7ZW.
##' }
##'
##' The variance of the simple intercept formula is \deqn{ Var\left(b_{0|X = 0,
##' Z, W}\right) = Var\left(b_0\right) + Z^2Var\left(b_2\right) +
##' W^2Var\left(b_3\right) + Z^2W^2Var\left(b_6\right) + 2ZCov\left(b_0,
##' b_2\right) + 2WCov\left(b_0, b_3\right) + 2ZWCov\left(b_0, b_6\right) +
##' 2ZWCov\left(b_2, b_3\right) + 2Z^2WCov\left(b_2, b_6\right) +
##' 2ZW^2Cov\left(b_3, b_6\right) } where \eqn{Var} denotes the variance of a
##' parameter estimate and \eqn{Cov} denotes the covariance of two parameter
##' estimates.
##'
##' The variance of the simple slope formula is \deqn{ Var\left(b_{X|Z,
##' W}\right) = Var\left(b_1\right) + Z^2Var\left(b_4\right) +
##' W^2Var\left(b_5\right) + Z^2W^2Var\left(b_7\right) + 2ZCov\left(b_1,
##' b_4\right) + 2WCov\left(b_1, b_5\right) + 2ZWCov\left(b_1, b_7\right) +
##' 2ZWCov\left(b_4, b_5\right) + 2Z^2WCov\left(b_4, b_7\right) +
##' 2ZW^2Cov\left(b_5, b_7\right) }
##'
##' Wald \emph{z} statistic is used for test statistic (even for objects of
##' class \code{\linkS4class{lavaan.mi}}).
##'
##'
##' @importFrom lavaan lavInspect parTable
##' @importFrom stats pnorm
##' @importFrom methods getMethod
##'
##' @param fit A fitted \code{\linkS4class{lavaan}} or
##'   \code{\linkS4class{lavaan.mi}} object with a latent 2-way interaction.
##' @param nameX \code{character} vector of all 7 factor names used as the
##'   predictors. The 3 lower-order factors must be listed first, followed by
##'   the 3 second-order factors (specifically, the 4th element must be the
##'   interaction between the factors listed first and second, the 5th element
##'   must be the interaction between the factors listed first and third, and
##'   the 6th element must be the interaction between the factors listed second
##'   and third). The final name will be the factor representing the 3-way
##'   interaction.
##' @param nameY The name of factor that is used as the dependent variable.
##' @param modVar The name of two factors that are used as the moderators. The
##'   effect of the independent factor on each combination of the moderator
##'   variable values will be probed.
##' @param valProbe1 The values of the first moderator that will be used to
##'   probe the effect of the independent factor.
##' @param valProbe2 The values of the second moderator that will be used to
##'   probe the effect of the independent factor.
##' @param group In multigroup models, the label of the group for which the
##'   results will be returned. Must correspond to one of
##'   \code{\link[lavaan]{lavInspect}(fit, "group.label")}.
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'        imputations from pooled results. Ignored unless \code{fit} is of
##'        class \code{\linkS4class{lavaan.mi}}. Can include any of
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
##' @return A list with two elements:
##' \enumerate{
##'  \item \code{SimpleIntercept}: The intercepts given each combination of
##'   moderator values. This element will be shown only if the factor intercept
##'   is estimated (e.g., not fixed at 0).
##'  \item \code{SimpleSlope}: The slopes given each combination of moderator
##'   values.
##' }
##' In each element, the first column represents values of the first moderator
##' specified in the \code{valProbe1} argument. The second column represents
##' values of the second moderator specified in the \code{valProbe2} argument.
##' The third column is the simple intercept or simple slope. The fourth column
##' is the standard error of the simple intercept or simple slope. The fifth
##' column is the Wald (\emph{z}) statistic. The sixth column is the \emph{p}
##' value testing whether the simple intercepts or slopes are different from 0.
##'
##' @author
##' Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @seealso \itemize{
##'  \item \code{\link{indProd}} For creating the indicator products with no
##'   centering, mean centering, double-mean centering, or residual centering.
##'  \item \code{\link{probe2WayMC}} For probing the two-way latent interaction
##'   when the results are obtained from mean-centering, or double-mean centering
##'  \item \code{\link{probe2WayRC}} For probing the two-way latent interaction
##'   when the results are obtained from residual-centering approach.
##'  \item \code{\link{probe3WayRC}} For probing the two-way latent interaction
##'   when the results are obtained from residual-centering approach.
##'  \item \code{\link{plotProbe}} Plot the simple intercepts and slopes of the
##'   latent interaction.
##' }
##'
##' @references
##' Tutorial:
##'
##' Schoemann, A. M., & Jorgensen, T. D. (2021). Testing and interpreting
##' latent variable interactions using the \code{semTools} package.
##' \emph{Psych, 3}(3), 322--335. \doi{10.3390/psych3030024}
##'
##' Background literature:
##'
##' Aiken, L. S., & West, S. G. (1991). \emph{Multiple regression: Testing
##' and interpreting interactions}. Newbury Park, CA: Sage.
##'
##' Marsh, H. W., Wen, Z., & Hau, K. T. (2004). Structural equation models of
##' latent interactions: Evaluation of alternative estimation strategies and
##' indicator construction. \emph{Psychological Methods, 9}(3), 275--300.
##' \doi{10.1037/1082-989X.9.3.275}
##'
##' @examples
##'
##' dat3wayMC <- indProd(dat3way, 1:3, 4:6, 7:9)
##'
##' model3 <- " ## define latent variables
##' f1 =~ x1 + x2 + x3
##' f2 =~ x4 + x5 + x6
##' f3 =~ x7 + x8 + x9
##' ## 2-way interactions
##' f12 =~ x1.x4 + x2.x5 + x3.x6
##' f13 =~ x1.x7 + x2.x8 + x3.x9
##' f23 =~ x4.x7 + x5.x8 + x6.x9
##' ## 3-way interaction
##' f123 =~ x1.x4.x7 + x2.x5.x8 + x3.x6.x9
##' ## outcome variable
##' f4 =~ x10 + x11 + x12
##'
##' ## latent regression model
##' f4 ~ b1*f1 + b2*f2 + b3*f3 + b12*f12 + b13*f13 + b23*f23 + b123*f123
##'
##' ## orthogonal terms among predictors
##' f1 ~~ 0*f12 + 0*f13 + 0*f123
##' f2 ~~ 0*f12 + 0*f23 + 0*f123
##' f3 ~~ 0*f13 + 0*f23 + 0*f123
##' f12 + f13 + f23 ~~ 0*f123
##'
##' ## identify latent means
##' x1 + x4 + x7 + x1.x4 + x1.x7 + x4.x7 + x1.x4.x7 + x10 ~ 0*1
##' f1 + f2 + f3 + f12 + f13 + f23 + f123 + f4 ~ NA*1
##' "
##'
##' fitMC3way <- sem(model3, data = dat3wayMC, meanstructure = TRUE)
##' summary(fitMC3way)
##'
##' probe3WayMC(fitMC3way, nameX = c("f1" ,"f2" ,"f3",
##'                                  "f12","f13","f23", # the order matters!
##'                                  "f123"),           # 3-way interaction
##'             nameY = "f4", modVar = c("f1", "f2"),
##'             valProbe1 = c(-1, 0, 1), valProbe2 = c(-1, 0, 1))
##'
##' @export
probe3WayMC <- function(fit, nameX, nameY, modVar, valProbe1, valProbe2,
                        group = 1L, omit.imps = c("no.conv","no.se")) {
  ## TDJ: verify class
  if (inherits(fit, "lavaan")) {
    est <- lavInspect(fit, "est")[[group]]

  } else if (inherits(fit, "lavaan.mi")) {
    useImps <- rep(TRUE, length(fit@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(fit@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(fit@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(fit@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(fit@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    ## custom removal by imputation number
    rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
    if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
    ## whatever is left
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

  } else stop('"fit" must inherit from lavaan or lavaan.mi class', call. = FALSE)


  # Check whether modVar is correct
	if (is.character(modVar)) modVar <- match(modVar, nameX)
	if ((NA %in% modVar) || !(do.call("&", as.list(modVar %in% 1:3))))
	  stop("The moderator name is not in the list of independent factors and is not 1, 2 or 3.")

	## TDJ: If multigroup, check group %in% group.label
	nG <- lavInspect(fit, "ngroups")
	if (nG > 1L) {
	  group.label <- lavInspect(fit, "group.label")
	  ## assign numeric to character
	  if (is.numeric(group)) {
	    if (group %in% 1:nG) {
	      group <- group.label[group]
	    } else group <- as.character(group)
	  } else group <- as.character(group)
	  ## check that character is a group
	  if (!as.character(group) %in% group.label)
	    stop('"group" must be a character string naming a group of interest, or ',
	         'an integer corresponding to a group in  lavInspect(fit, "group.label")')
	  group.number <- which(group.label == group)
	} else group.number <- 1L

	# Extract all varEst
	if (inherits(fit, "lavaan")) {
	  varEst <- lavaan::vcov(fit)
	} else if (inherits(fit, "lavaan.mi")) {
	  varEst <- getMethod("vcov", "lavaan.mi")(fit, omit.imps = omit.imps)
	}

	## Check whether the outcome's intercept is estimated
	PT <- parTable(fit)
	if (lavInspect(fit, "options")$meanstructure) {
	  targetcol <- PT$label[PT$lhs == nameY & PT$op == "~1" & PT$group == group.number]
	  if (targetcol == "") {
	    ## no custom label, use default
	    targetcol <- paste0(nameY, "~1")
	    if (nG > 1L && group.number > 1L) {
	      targetcol <- paste0(targetcol, ".g", group.number)
	    }
	  }
	  ## check it is actually estimated (thus, has sampling variance)
	  estimateIntcept <- targetcol %in% rownames(varEst)

	} else estimateIntcept <- FALSE


	## Get the parameter estimates for that group
	if (nG > 1L) {

	  if (inherits(fit, "lavaan")) {
	    est <- lavInspect(fit, "est")[[group]]
	  } else if (inherits(fit, "lavaan.mi")) {
	    est <- list()
	    GLIST <- fit@coefList[useImps]
	    est$beta <- Reduce("+", lapply(GLIST, function(i) i[[group]]$beta)) / m
	    if (estimateIntcept) {
	      est$alpha <- Reduce("+", lapply(GLIST, function(i) i[[group]]$alpha)) / m
	    }
	  }

	} else {
	  ## single-group model

	  if (inherits(fit, "lavaan")) {
	    est <- lavInspect(fit, "est")
	  } else if (inherits(fit, "lavaan.mi")) {
	    est <- list()
	    est$beta <- Reduce("+", lapply(fit@coefList[useImps], "[[", i = "beta")) / m
	    if (estimateIntcept) {
	      est$alpha <- Reduce("+", lapply(fit@coefList[useImps], "[[", i = "alpha")) / m
	    }
	  }

	}


	# Compute the intercept of no-centering
	betaNC <- matrix(est$beta[nameY, nameX], ncol = 1,
	                 dimnames = list(nameX, nameY))

	pvalue <- function(x) (1 - pnorm(abs(x))) * 2

	# Find the order to rearrange
	ord <- c(setdiff(1:3, modVar), modVar)
	ord <- c(ord, 7 - rev(ord))

	resultIntcept <- NULL
	resultSlope <- NULL
	if(estimateIntcept) {
		# Extract SE from centered result
	  newRows <- which(PT$lhs == nameY & PT$op == "~" & PT$rhs %in% nameX & PT$group == group.number)
	  newLabels <- PT$label[newRows]
	  if (any(newLabels == "")) for (i in which(newLabels == "")) {
	    newLabels[i] <- paste0(nameY, "~", nameX[i],
	                           ifelse(nG > 1L && group.number > 1L, no = "",
	                                  yes = paste0(".g", group.number)))
	  }
	  targetcol <- c(targetcol, newLabels)

		# Transform it to non-centering SE
		usedVar <-  varEst[targetcol, targetcol]
		usedBeta <- rbind(est$alpha[nameY,], betaNC)
		if (sum(diag(usedVar) < 0) > 0)
		  stop("This method does not work. The resulting calculation provided negative standard errors.") # JG: edited this error

		# Change the order of usedVar and usedBeta if the moderator variable is listed first
		usedVar <- usedVar[c(1, ord+1, 8), c(1, ord+1, 8)]
		usedBeta <- usedBeta[c(1, ord+1, 8)]

		# Find probe value
		val <- expand.grid(valProbe1, valProbe2)

		# Find simple intercept
		simpleIntcept <- usedBeta[1] + usedBeta[3] * val[,1] + usedBeta[4] * val[,2] + usedBeta[7] * val[,1] * val[,2]
		varIntcept <- usedVar[1, 1] + val[,1]^2 * usedVar[3, 3] + val[,2]^2 * usedVar[4, 4] + val[,1]^2 * val[,2]^2 * usedVar[7, 7] + 2 * val[,1] * usedVar[1, 3] + 2 * val[,2] * usedVar[1, 4] + 2 * val[,1] * val[,2] * usedVar[1, 7] + 2 * val[,1] * val[,2] * usedVar[3, 4] + 2 * val[,1]^2 * val[,2] * usedVar[3, 7] + 2* val[,1] * val[,2]^2 * usedVar[4, 7]
		zIntcept <- simpleIntcept / sqrt(varIntcept)
		pIntcept <- pvalue(zIntcept)
		resultIntcept <- data.frame(val, simpleIntcept, sqrt(varIntcept), zIntcept, pIntcept)
		colnames(resultIntcept) <- c(nameX[modVar], "est", "se", "z", "pvalue")
		class(resultIntcept) <- c("lavaan.data.frame","data.frame")

		# Find simple slope
		simpleSlope <- usedBeta[2] + usedBeta[5] * val[,1] + usedBeta[6] * val[,2] + usedBeta[8] * val[,1] * val[,2]
		varSlope <- usedVar[2, 2] + val[,1]^2 * usedVar[5, 5] + val[,2]^2 * usedVar[6, 6] + val[,1]^2 * val[,2]^2 * usedVar[8, 8] + 2 * val[,1] * usedVar[2, 5] + 2 * val[,2] * usedVar[2, 6] + 2 * val[,1] * val[,2] * usedVar[2, 8] + 2 * val[,1] * val[,2] * usedVar[5, 6] + 2 * val[,1]^2 * val[,2] * usedVar[5, 8] + 2 * val[,1] * val[,2]^2 * usedVar[6, 8]
		zSlope <- simpleSlope / sqrt(varSlope)
		pSlope <- pvalue(zSlope)
		resultSlope <- data.frame(val, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "est", "se", "z", "pvalue")
		class(resultSlope) <- c("lavaan.data.frame","data.frame")

	} else {
	  newRows <- which(PT$lhs == nameY & PT$op == "~" & PT$rhs %in% nameX & PT$group == group.number)
	  targetcol <- PT$label[newRows]
	  if (any(targetcol == "")) for (i in which(targetcol == "")) {
	    targetcol[i] <- paste0(nameY, "~", PT$rhs[ newRows[i] ],
	                           ifelse(nG > 1L && group.number > 1L, no = "",
	                                  yes = paste0(".g", group.number)))
	  }

		# Transform it to non-centering SE
		usedVar <-  varEst[targetcol, targetcol]
		usedBeta <- betaNC
		if (sum(diag(usedVar) < 0) > 0) stop("This method does not work. The resulting calculation provided negative standard errors.") # JG: edited this error

		# Change the order of usedVar and usedBeta if the moderator variable is listed first
		usedVar <- usedVar[c(ord, 7), c(ord, 7)]
		usedBeta <- usedBeta[c(ord, 7)]

		# Find probe value
		val <- expand.grid(valProbe1, valProbe2)

		# Find simple slope
		simpleSlope <- usedBeta[1] + usedBeta[4] * val[,1] + usedBeta[5] * val[,2] + usedBeta[7] * val[,1] * val[,2]
		varSlope <- usedVar[1, 1] + val[,1]^2 * usedVar[4, 4] + val[,2]^2 * usedVar[5, 5] + val[,1]^2 * val[,2]^2 * usedVar[7, 7] + 2 * val[,1] * usedVar[1, 4] + 2 * val[,2] * usedVar[1, 5] + 2 * val[,1] * val[,2] * usedVar[1, 7] + 2 * val[,1] * val[,2] * usedVar[4, 5] + 2 * val[,1]^2 * val[,2] * usedVar[4, 7] + 2 * val[,1] * val[,2]^2 * usedVar[5, 7]
		zSlope <- simpleSlope / sqrt(varSlope)
		pSlope <- pvalue(zSlope)
		resultSlope <- data.frame(val, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "est", "se", "z", "pvalue")
		class(resultSlope) <- c("lavaan.data.frame","data.frame")
	}

	list(SimpleIntcept = resultIntcept, SimpleSlope = resultSlope)
}



## --------
## 3-way RC
## --------

##' Probing three-way interaction on the residual-centered latent interaction
##'
##' Probing interaction for simple intercept and simple slope for the
##' residual-centered latent three-way interaction (Geldhof et al., 2013)
##'
##' Before using this function, researchers need to make the products of the
##' indicators between the first-order factors and residualize the products by
##' the original indicators (Lance, 1988; Little, Bovaird, & Widaman, 2006). The
##' process can be automated by the \code{\link{indProd}} function. Note that
##' the indicator products can be made for all possible combination or
##' matched-pair approach (Marsh et al., 2004). Next, the hypothesized model
##' with the regression with latent interaction will be used to fit all original
##' indicators and the product terms (Geldhof et al., 2013). To use this
##' function the model must be fit with a mean structure. See the example for
##' how to fit the product term below. Once the lavaan result is obtained, this
##' function will be used to probe the interaction.
##'
##' The probing process on residual-centered latent interaction is based on
##' transforming the residual-centered result into the no-centered result. See
##' Geldhof et al. (2013) for further details. Note that this approach based on
##' a strong assumption that the first-order latent variables are normally
##' distributed. The probing process is applied after the no-centered result
##' (parameter estimates and their covariance matrix among parameter estimates)
##' has been computed. See the \code{\link{probe3WayMC}} for further details.
##'
##'
##' @importFrom lavaan lavInspect parTable
##' @importFrom stats pnorm
##' @importFrom methods getMethod
##'
##' @param fit A fitted \code{\linkS4class{lavaan}} or
##'   \code{\linkS4class{lavaan.mi}} object with a latent 2-way interaction.
##' @param nameX \code{character} vector of all 7 factor names used as the
##'   predictors. The 3 lower-order factors must be listed first, followed by
##'   the 3 second-order factors (specifically, the 4th element must be the
##'   interaction between the factors listed first and second, the 5th element
##'   must be the interaction between the factors listed first and third, and
##'   the 6th element must be the interaction between the factors listed second
##'   and third). The final name will be the factor representing the 3-way
##'   interaction.
##' @param nameY The name of factor that is used as the dependent variable.
##' @param modVar The name of two factors that are used as the moderators. The
##'   effect of the independent factor on each combination of the moderator
##'   variable values will be probed.
##' @param valProbe1 The values of the first moderator that will be used to
##'   probe the effect of the independent factor.
##' @param valProbe2 The values of the second moderator that will be used to
##'   probe the effect of the independent factor.
##' @param group In multigroup models, the label of the group for which the
##'   results will be returned. Must correspond to one of
##'   \code{\link[lavaan]{lavInspect}(fit, "group.label")}.
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'        imputations from pooled results. Ignored unless \code{fit} is of
##'        class \code{\linkS4class{lavaan.mi}}. Can include any of
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
##' @return A list with two elements:
##' \enumerate{
##'  \item \code{SimpleIntercept}: The intercepts given each value of the moderator.
##'    This element will be shown only if the factor intercept is estimated
##'    (e.g., not fixed as 0).
##'  \item \code{SimpleSlope}: The slopes given each value of the moderator.
##' }
##' In each element, the first column represents values of the first moderator
##' specified in the \code{valProbe1} argument. The second column represents
##' values of the second moderator specified in the \code{valProbe2} argument.
##' The third column is the simple intercept or simple slope. The fourth column
##' is the \emph{SE} of the simple intercept or simple slope. The fifth column
##' is the Wald (\emph{z}) statistic. The sixth column is the \emph{p} value
##' testing whether the simple intercepts or slopes are different from 0.
##'
##' @author
##' Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @seealso \itemize{
##'  \item \code{\link{indProd}} For creating the indicator products with no
##'   centering, mean centering, double-mean centering, or residual centering.
##'  \item \code{\link{probe2WayMC}} For probing the two-way latent interaction
##'   when the results are obtained from mean-centering, or double-mean centering
##'  \item \code{\link{probe3WayMC}} For probing the three-way latent interaction
##'   when the results are obtained from mean-centering, or double-mean centering
##'  \item \code{\link{probe2WayRC}} For probing the two-way latent interaction
##'   when the results are obtained from residual-centering approach.
##'  \item \code{\link{plotProbe}} Plot the simple intercepts and slopes of the
##'   latent interaction.
##' }
##'
##' @references
##' Tutorial:
##'
##' Schoemann, A. M., & Jorgensen, T. D. (2021). Testing and interpreting
##' latent variable interactions using the \code{semTools} package.
##' \emph{Psych, 3}(3), 322--335. \doi{10.3390/psych3030024}
##'
##' Background literature:
##'
##' Geldhof, G. J., Pornprasertmanit, S., Schoemann, A., & Little,
##' T. D. (2013). Orthogonalizing through residual centering: Extended
##' applications and caveats. \emph{Educational and Psychological Measurement,
##' 73}(1), 27--46. \doi{10.1177/0013164412445473}
##'
##' Lance, C. E. (1988). Residual centering, exploratory and confirmatory
##' moderator analysis, and decomposition of effects in path models containing
##' interactions. \emph{Applied Psychological Measurement, 12}(2), 163--175.
##' \doi{10.1177/014662168801200205}
##'
##' Little, T. D., Bovaird, J. A., & Widaman, K. F. (2006). On the merits of
##' orthogonalizing powered and product terms: Implications for modeling
##' interactions. \emph{Structural Equation Modeling, 13}(4), 497--519.
##' \doi{10.1207/s15328007sem1304_1}
##'
##' Marsh, H. W., Wen, Z., & Hau, K. T. (2004). Structural equation models of
##' latent interactions: Evaluation of alternative estimation strategies and
##' indicator construction. \emph{Psychological Methods, 9}(3), 275--300.
##' \doi{10.1037/1082-989X.9.3.275}
##'
##' Pornprasertmanit, S., Schoemann, A. M., Geldhof, G. J., & Little, T. D.
##' (submitted). \emph{Probing latent interaction estimated with a residual
##' centering approach.}
##'
##' @examples
##'
##' dat3wayRC <- orthogonalize(dat3way, 1:3, 4:6, 7:9)
##'
##' model3 <- " ## define latent variables
##' f1 =~ x1 + x2 + x3
##' f2 =~ x4 + x5 + x6
##' f3 =~ x7 + x8 + x9
##' ## 2-way interactions
##' f12 =~ x1.x4 + x2.x5 + x3.x6
##' f13 =~ x1.x7 + x2.x8 + x3.x9
##' f23 =~ x4.x7 + x5.x8 + x6.x9
##' ## 3-way interaction
##' f123 =~ x1.x4.x7 + x2.x5.x8 + x3.x6.x9
##' ## outcome variable
##' f4 =~ x10 + x11 + x12
##'
##' ## latent regression model
##' f4 ~ b1*f1 + b2*f2 + b3*f3 + b12*f12 + b13*f13 + b23*f23 + b123*f123
##'
##' ## orthogonal terms among predictors
##' f1 ~~ 0*f12 + 0*f13 + 0*f123
##' f2 ~~ 0*f12 + 0*f23 + 0*f123
##' f3 ~~ 0*f13 + 0*f23 + 0*f123
##' f12 + f13 + f23 ~~ 0*f123
##'
##' ## identify latent means
##' x1 + x4 + x7 + x1.x4 + x1.x7 + x4.x7 + x1.x4.x7 + x10 ~ 0*1
##' f1 + f2 + f3 + f12 + f13 + f23 + f123 + f4 ~ NA*1
##' "
##'
##' fitRC3way <- sem(model3, data = dat3wayRC, meanstructure = TRUE)
##' summary(fitRC3way)
##'
##' probe3WayMC(fitRC3way, nameX = c("f1" ,"f2" ,"f3",
##'                                  "f12","f13","f23", # the order matters!
##'                                  "f123"),           # 3-way interaction
##'             nameY = "f4", modVar = c("f1", "f2"),
##'             valProbe1 = c(-1, 0, 1), valProbe2 = c(-1, 0, 1))
##'
##' @export
probe3WayRC <- function(fit, nameX, nameY, modVar, valProbe1, valProbe2,
                        group = 1L, omit.imps = c("no.conv","no.se")) {
  ## TDJ: verify class
  if (inherits(fit, "lavaan")) {
    est <- lavInspect(fit, "est")[[group]]

  } else if (inherits(fit, "lavaan.mi")) {
    useImps <- rep(TRUE, length(fit@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(fit@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(fit@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(fit@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(fit@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    ## custom removal by imputation number
    rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
    if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
    ## whatever is left
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

  } else stop('"fit" must inherit from lavaan or lavaan.mi class', call. = FALSE)

  if (!lavInspect(fit, "options")$meanstructure)
    stop('This function requires the model to be fit with a mean structure.',
         call. = FALSE)

  # Check whether modVar is correct
	if (is.character(modVar)) modVar <- match(modVar, nameX)
	if ((NA %in% modVar) || !(do.call("&", as.list(modVar %in% 1:3))))
	  stop("The moderator name is not in the list of independent factors and is ",
	       "not 1, 2 or 3.") # JG: Changed error

	## TDJ: If multigroup, check group %in% group.label
	nG <- lavInspect(fit, "ngroups")
	if (nG > 1L) {
	  group.label <- lavInspect(fit, "group.label")
	  ## assign numeric to character
	  if (is.numeric(group)) {
	    if (group %in% 1:nG) {
	      group <- group.label[group]
	    } else group <- as.character(group)
	  } else group <- as.character(group)
	  ## check that character is a group
	  if (!as.character(group) %in% group.label)
	    stop('"group" must be a character string naming a group of interest, or ',
	         'an integer corresponding to a group in  lavInspect(fit, "group.label")')
	  group.number <- which(group.label == group)
	} else group.number <- 1L

	# Extract all varEst
	if (inherits(fit, "lavaan")) {
	  varEst <- lavaan::vcov(fit)
	} else if (inherits(fit, "lavaan.mi")) {
	  varEst <- getMethod("vcov", "lavaan.mi")(fit, omit.imps = omit.imps)
	}

	## Check whether the outcome's intercept is estimated
	PT <- parTable(fit)
	if (lavInspect(fit, "options")$meanstructure) {
	  targetcol <- PT$label[PT$lhs == nameY & PT$op == "~1" & PT$group == group.number]
	  if (targetcol == "") {
	    ## no custom label, use default
	    targetcol <- paste0(nameY, "~1")
	    if (nG > 1L && group.number > 1L) {
	      targetcol <- paste0(targetcol, ".g", group.number)
	    }
	  }
	  ## check it is actually estimated (thus, has sampling variance)
	  estimateIntcept <- targetcol %in% rownames(varEst)

	} else estimateIntcept <- FALSE


	## Get the parameter estimates for that group
	if (nG > 1L) {

	  if (inherits(fit, "lavaan")) {
	    est <- lavInspect(fit, "est")[[group]]
	  } else if (inherits(fit, "lavaan.mi")) {
	    est <- list()
	    GLIST <- fit@coefList[useImps]
	    est$beta <- Reduce("+", lapply(GLIST, function(i) i[[group]]$beta)) / m
	    est$alpha <- Reduce("+", lapply(GLIST, function(i) i[[group]]$alpha)) / m
	    est$psi <- Reduce("+", lapply(GLIST, function(i) i[[group]]$psi)) / m
	  }

	} else {
	  ## single-group model

	  if (inherits(fit, "lavaan")) {
	    est <- lavInspect(fit, "est")
	  } else if (inherits(fit, "lavaan.mi")) {
	    est <- list()
	    est$beta <- Reduce("+", lapply(fit@coefList[useImps], "[[", i = "beta")) / m
	    est$alpha <- Reduce("+", lapply(fit@coefList[useImps], "[[", i = "alpha")) / m
	    est$psi <- Reduce("+", lapply(fit@coefList[useImps], "[[", i = "psi")) / m
	  }

	}


	# Find the mean and covariance matrix of independent factors
	varX <- est$psi[nameX, nameX]
	meanX <- matrix(est$alpha[nameX,], ncol = 1, dimnames = list(NULL, "intcept"))

	# Find the intercept, regression coefficients, and residual variance of residual-centered regression
	intceptRC <- est$alpha[nameY,]
	resVarRC <- est$psi[nameY, nameY]
	if (resVarRC < 0) stop("The residual variance is negative. The model did not converge!") # JG: Changed error
	betaRC <- as.matrix(est$beta[nameY, nameX]); colnames(betaRC) <- nameY

	# Find the number of observations
	numobs <- lavInspect(fit, "nobs")[group.number]

	# Compute SSRC
	meanXwith1 <- rbind(1, meanX)
	varXwith0 <- cbind(0, rbind(0, varX))
	SSRC <- numobs * (varXwith0 + (meanXwith1 %*% t(meanXwith1)))

	# Compute Mean(Y) and Var(Y)
	betaRCWithIntcept <- rbind(intceptRC, betaRC)
	meanY <- t(meanXwith1) %*% betaRCWithIntcept
	varY <- (t(betaRCWithIntcept) %*% SSRC %*% betaRCWithIntcept)/numobs - meanY^2 + resVarRC

	# Compute Cov(Y, X)
	covY <- as.matrix((varX %*% betaRC)[1:3,])

	# Compute E(XZ), E(XW), E(ZW), E(XZW)
	meanX[4] <- expect2NormProd(meanX[c(1,2)], varX[c(1,2), c(1,2)])
	meanX[5] <- expect2NormProd(meanX[c(1,3)], varX[c(1,3), c(1,3)])
	meanX[6] <- expect2NormProd(meanX[c(2,3)], varX[c(2,3), c(2,3)])
	meanX[7] <- expect3NormProd(meanX[1:3], varX[1:3, 1:3])

	# Compute Var(XZ), Var(XW), Var(ZW), Var(XZW)
	varX[4, 4] <- var2NormProd(meanX[c(1,2)], varX[c(1,2), c(1,2)])
	varX[5, 5] <- var2NormProd(meanX[c(1,3)], varX[c(1,3), c(1,3)])
	varX[6, 6] <- var2NormProd(meanX[c(2,3)], varX[c(2,3), c(2,3)])
	varX[7, 7] <- var3NormProd(meanX[1:3], varX[1:3, 1:3])

	# Compute All covariances
	varX[4, 1] <- varX[1, 4] <- expect3NormProd(meanX[c(1, 2, 1)], varX[c(1, 2, 1),c(1, 2, 1)]) - expect2NormProd(meanX[c(1,2)], varX[c(1,2),c(1,2)]) * meanX[1]
	varX[5, 1] <- varX[1, 5] <- expect3NormProd(meanX[c(1, 3, 1)], varX[c(1, 3, 1),c(1, 3, 1)]) - expect2NormProd(meanX[c(1,3)], varX[c(1,3),c(1,3)]) * meanX[1]
	varX[6, 1] <- varX[1, 6] <- expect3NormProd(meanX[c(2, 3, 1)], varX[c(2, 3, 1),c(2, 3, 1)]) - expect2NormProd(meanX[c(2,3)], varX[c(2,3),c(2,3)]) * meanX[1]
	varX[7, 1] <- varX[1, 7] <- expect4NormProd(meanX[c(1,2,3,1)], varX[c(1,2,3,1),c(1,2,3,1)]) - expect3NormProd(meanX[c(1,2,3)], varX[c(1,2,3),c(1,2,3)]) * meanX[1]

	varX[4, 2] <- varX[2, 4] <- expect3NormProd(meanX[c(1, 2, 2)], varX[c(1, 2, 2),c(1, 2, 2)]) - expect2NormProd(meanX[c(1,2)], varX[c(1,2),c(1,2)]) * meanX[2]
	varX[5, 2] <- varX[2, 5] <- expect3NormProd(meanX[c(1, 3, 2)], varX[c(1, 3, 2),c(1, 3, 2)]) - expect2NormProd(meanX[c(1,3)], varX[c(1,3),c(1,3)]) * meanX[2]
	varX[6, 2] <- varX[2, 6] <- expect3NormProd(meanX[c(2, 3, 2)], varX[c(2, 3, 2),c(2, 3, 2)]) - expect2NormProd(meanX[c(2,3)], varX[c(2,3),c(2,3)]) * meanX[2]
	varX[7, 2] <- varX[2, 7] <- expect4NormProd(meanX[c(1,2,3,2)], varX[c(1,2,3,2),c(1,2,3,2)]) - expect3NormProd(meanX[c(1,2,3)], varX[c(1,2,3),c(1,2,3)]) * meanX[2]

	varX[4, 3] <- varX[3, 4] <- expect3NormProd(meanX[c(1, 2, 3)], varX[c(1, 2, 3),c(1, 2, 3)]) - expect2NormProd(meanX[c(1,2)], varX[c(1,2),c(1,2)]) * meanX[3]
	varX[5, 3] <- varX[3, 5] <- expect3NormProd(meanX[c(1, 3, 3)], varX[c(1, 3, 3),c(1, 3, 3)]) - expect2NormProd(meanX[c(1,3)], varX[c(1,3),c(1,3)]) * meanX[3]
	varX[6, 3] <- varX[3, 6] <- expect3NormProd(meanX[c(2, 3, 3)], varX[c(2, 3, 3),c(2, 3, 3)]) - expect2NormProd(meanX[c(2,3)], varX[c(2,3),c(2,3)]) * meanX[3]
	varX[7, 3] <- varX[3, 7] <- expect4NormProd(meanX[c(1,2,3,3)], varX[c(1,2,3,3),c(1,2,3,3)]) - expect3NormProd(meanX[c(1,2,3)], varX[c(1,2,3),c(1,2,3)]) * meanX[3]

	varX[5, 4] <- varX[4, 5] <- expect4NormProd(meanX[c(1,3,1,2)], varX[c(1,3,1,2),c(1,3,1,2)]) - expect2NormProd(meanX[c(1,3)], varX[c(1,3),c(1,3)]) * expect2NormProd(meanX[c(1,2)], varX[c(1,2),c(1,2)])
	varX[6, 4] <- varX[4, 6] <- expect4NormProd(meanX[c(2,3,1,2)], varX[c(2,3,1,2),c(2,3,1,2)]) - expect2NormProd(meanX[c(2,3)], varX[c(2,3),c(2,3)]) * expect2NormProd(meanX[c(1,2)], varX[c(1,2),c(1,2)])
	varX[7, 4] <- varX[4, 7] <- expect5NormProd(meanX[c(1,2,3,1,2)], varX[c(1,2,3,1,2),c(1,2,3,1,2)]) - expect3NormProd(meanX[c(1, 2, 3)], varX[c(1, 2, 3),c(1, 2, 3)]) * expect2NormProd(meanX[c(1,2)], varX[c(1,2),c(1,2)])

	varX[6, 5] <- varX[5, 6] <- expect4NormProd(meanX[c(2,3,1,3)], varX[c(2,3,1,3),c(2,3,1,3)]) - expect2NormProd(meanX[c(2,3)], varX[c(2,3),c(2,3)]) * expect2NormProd(meanX[c(1,3)], varX[c(1,3),c(1,3)])
	varX[7, 5] <- varX[5, 7] <- expect5NormProd(meanX[c(1,2,3,1,3)], varX[c(1,2,3,1,3),c(1,2,3,1,3)]) - expect3NormProd(meanX[c(1, 2, 3)], varX[c(1, 2, 3),c(1, 2, 3)]) * expect2NormProd(meanX[c(1,3)], varX[c(1,3),c(1,3)])
	varX[7, 6] <- varX[6, 7] <- expect5NormProd(meanX[c(1,2,3,2,3)], varX[c(1,2,3,2,3),c(1,2,3,2,3)]) - expect3NormProd(meanX[c(1, 2, 3)], varX[c(1, 2, 3),c(1, 2, 3)]) * expect2NormProd(meanX[c(2,3)], varX[c(2,3),c(2,3)])

	# Find the meanX and varX without XZW
	meanXReducedWith1 <- rbind(1, as.matrix(meanX[1:6]))
	varXReducedWith0 <- cbind(0, rbind(0, varX[1:6, 1:6]))
	SSMCReduced <- numobs * (varXReducedWith0 + (meanXReducedWith1 %*% t(meanXReducedWith1)))

	# Find product of main and two-way onto three-way
	covXZWwith0 <- rbind(0, as.matrix(varX[7, 1:6]))
	meanXZWwith1 <- meanX[7] * meanXReducedWith1
	SSXZW <- numobs * (covXZWwith0 + meanXZWwith1) # should the mean vector be squared (postmultiplied by its transpose)?

	# Compute a vector and b4, b5, b6
	a <- solve(SSMCReduced) %*% as.matrix(SSXZW)
	betaTemp <- betaRC[4:6] - (as.numeric(betaRC[7]) * a[5:7])
	betaTemp <- c(betaTemp, betaRC[7])

	# Compute Cov(Y, XZ) and regression coefficients of no-centering
	betaNC <- solve(varX[1:3,1:3], as.matrix(covY) - (t(varX[4:7, 1:3]) %*% as.matrix(betaTemp)))
	betaNC <- rbind(as.matrix(betaNC), as.matrix(betaTemp))
	covY <- rbind(covY, as.matrix((varX %*% betaNC)[4:7, 1]))

	# Aggregate the non-centering sufficient statistics (Just show how to do but not necessary)
	fullCov <- rbind(cbind(varX, covY), c(covY, varY))
	fullMean <- rbind(meanX, meanY)

	# Compute the intercept of no-centering
	intceptNC <- meanY - t(betaNC) %*% meanX

	# Compute SSNC
	betaNCWithIntcept <- rbind(intceptNC, betaNC)
	meanXwith1 <- rbind(1, meanX) #JG: redundant
	varXwith0 <- rbind(0, cbind(0, varX)) #JG: redundant
	SSNC <- numobs * (varXwith0 + (meanXwith1 %*% t(meanXwith1)))

	# Compute residual variance on non-centering
	resVarNC <- varY - (t(betaNCWithIntcept) %*% SSNC %*% betaNCWithIntcept)/numobs + meanY^2


	pvalue <- function(x) (1 - pnorm(abs(x))) * 2

	# Find the order to rearrange
	ord <- c(setdiff(1:3, modVar), modVar)
	ord <- c(ord, 7 - rev(ord))

	resultIntcept <- NULL
	resultSlope <- NULL
	if (estimateIntcept) {
		# Extract SE from residual centering
	  newRows <- which(PT$lhs == nameY & PT$op == "~" & PT$rhs %in% nameX & PT$group == group.number)
	  newLabels <- PT$label[newRows]
	  if (any(newLabels == "")) for (i in which(newLabels == "")) {
	    newLabels[i] <- paste0(nameY, "~", nameX[i],
	                           ifelse(nG > 1L && group.number > 1L, no = "",
	                                  yes = paste0(".g", group.number)))
	  }
	  targetcol <- c(targetcol, newLabels)
		varEstSlopeRC <- varEst[targetcol, targetcol]

		# Transform it to non-centering SE
		usedVar <-  as.numeric(resVarNC/resVarRC) * (varEstSlopeRC %*% SSRC %*% solve(SSNC))
		usedBeta <- betaNCWithIntcept
		if (sum(diag(usedVar) < 0) > 0) stop("This method does not work. The resulting calculation provided negative standard errors.") # JG: edited this error

		# Change the order of usedVar and usedBeta if the moderator variable is listed first
		usedVar <- usedVar[c(1, ord+1, 8), c(1, ord+1, 8)]
		usedBeta <- usedBeta[c(1, ord+1, 8)]

		# Find probe value
		val <- expand.grid(valProbe1, valProbe2)

		# Find simple intercept
		simpleIntcept <- usedBeta[1] + usedBeta[3] * val[,1] + usedBeta[4] * val[,2] + usedBeta[7] * val[,1] * val[,2]
		varIntcept <- usedVar[1, 1] + val[,1]^2 * usedVar[3, 3] + val[,2]^2 * usedVar[4, 4] + val[,1]^2 * val[,2]^2 * usedVar[7, 7] + 2 * val[,1] * usedVar[1, 3] + 2 * val[,2] * usedVar[1, 4] + 2 * val[,1] * val[,2] * usedVar[1, 7] + 2 * val[,1] * val[,2] * usedVar[3, 4] + 2 * val[,1]^2 * val[,2] * usedVar[3, 7] + 2* val[,1] * val[,2]^2 * usedVar[4, 7]
		zIntcept <- simpleIntcept / sqrt(varIntcept)
		pIntcept <- pvalue(zIntcept)
		resultIntcept <- data.frame(val, simpleIntcept, sqrt(varIntcept), zIntcept, pIntcept)
		colnames(resultIntcept) <- c(nameX[modVar], "est", "se", "z", "pvalue")
		class(resultIntcept) <- c("lavaan.data.frame","data.frame")

		# Find simple slope
		simpleSlope <- usedBeta[2] + usedBeta[5] * val[,1] + usedBeta[6] * val[,2] + usedBeta[8] * val[,1] * val[,2]
		varSlope <- usedVar[2, 2] + val[,1]^2 * usedVar[5, 5] + val[,2]^2 * usedVar[6, 6] + val[,1]^2 * val[,2]^2 * usedVar[8, 8] + 2 * val[,1] * usedVar[2, 5] + 2 * val[,2] * usedVar[2, 6] + 2 * val[,1] * val[,2] * usedVar[2, 8] + 2 * val[,1] * val[,2] * usedVar[5, 6] + 2 * val[,1]^2 * val[,2] * usedVar[5, 8] + 2 * val[,1] * val[,2]^2 * usedVar[6, 8]
		zSlope <- simpleSlope / sqrt(varSlope)
		pSlope <- pvalue(zSlope)
		resultSlope <- data.frame(val, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "est", "se", "z", "pvalue")
		class(resultSlope) <- c("lavaan.data.frame","data.frame")

	} else {
	  newRows <- which(PT$lhs == nameY & PT$op == "~" & PT$rhs %in% nameX & PT$group == group.number)
	  targetcol <- PT$label[newRows]
	  if (any(targetcol == "")) for (i in which(targetcol == "")) {
	    targetcol[i] <- paste0(nameY, "~", PT$rhs[ newRows[i] ],
	                           ifelse(nG > 1L && group.number > 1L, no = "",
	                                  yes = paste0(".g", group.number)))
	  }
		varEstSlopeRC <- varEst[targetcol, targetcol]

		# Transform it to non-centering SE
		usedVar <-  as.numeric(resVarNC/resVarRC) * (varEstSlopeRC %*% SSRC[2:8, 2:8] %*% solve(SSNC[2:8, 2:8]))
		usedBeta <- betaNC
		if(sum(diag(usedVar) < 0) > 0) stop("This method does not work. The resulting calculation provided negative standard errors.") # JG: edited this error

		# Change the order of usedVar and usedBeta if the moderator variable is listed first
		usedVar <- usedVar[c(ord, 7), c(ord, 7)]
		usedBeta <- usedBeta[c(ord, 7)]

		# Find probe value
		val <- expand.grid(valProbe1, valProbe2)

		# Find simple slope
		simpleSlope <- usedBeta[1] + usedBeta[4] * val[,1] + usedBeta[5] * val[,2] + usedBeta[7] * val[,1] * val[,2]
		varSlope <- usedVar[1, 1] + val[,1]^2 * usedVar[4, 4] + val[,2]^2 * usedVar[5, 5] + val[,1]^2 * val[,2]^2 * usedVar[7, 7] + 2 * val[,1] * usedVar[1, 4] + 2 * val[,2] * usedVar[1, 5] + 2 * val[,1] * val[,2] * usedVar[1, 7] + 2 * val[,1] * val[,2] * usedVar[4, 5] + 2 * val[,1]^2 * val[,2] * usedVar[4, 7] + 2 * val[,1] * val[,2]^2 * usedVar[5, 7]
		zSlope <- simpleSlope / sqrt(varSlope)
		pSlope <- pvalue(zSlope)
		resultSlope <- data.frame(val, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "est", "se", "z", "pvalue")
		class(resultSlope) <- c("lavaan.data.frame","data.frame")
	}

	list(SimpleIntcept = resultIntcept, SimpleSlope = resultSlope)
}



## -----------------
## Plotting Function
## -----------------

##' Plot a latent interaction
##'
##' This function will plot the line graphs representing the simple effect of
##' the independent variable given the values of the moderator. For multigroup
##' models, it will only generate a plot for 1 group, as specified in the
##' function used to obtain the first argument.
##'
##'
##' @param object The result of probing latent interaction obtained from
##'   \code{\link{probe2WayMC}}, \code{\link{probe2WayRC}},
##'   \code{\link{probe3WayMC}}, or \code{\link{probe3WayRC}} function.
##' @param xlim The vector of two numbers: the minimum and maximum values of the
##'   independent variable
##' @param xlab The label of the x-axis
##' @param ylab The label of the y-axis
##' @param legend \code{logical}. If \code{TRUE} (default), a legend is printed.
##' @param legendArgs \code{list} of arguments passed to \code{\link{legend}}
##'   function if \code{legend=TRUE}.
##' @param \dots Any addition argument for the \code{\link{plot}} function
##'
##' @return None. This function will plot the simple main effect only.
##'
##' @author
##' Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @seealso \itemize{
##'  \item \code{\link{indProd}} For creating the indicator products with no
##'   centering, mean centering, double-mean centering, or residual centering.
##'  \item \code{\link{probe2WayMC}} For probing the two-way latent interaction
##'   when the results are obtained from mean-centering, or double-mean centering
##'  \item \code{\link{probe3WayMC}} For probing the three-way latent interaction
##'   when the results are obtained from mean-centering, or double-mean centering
##'  \item \code{\link{probe2WayRC}} For probing the two-way latent interaction
##'   when the results are obtained from residual-centering approach.
##'  \item \code{\link{probe3WayRC}} For probing the two-way latent interaction
##'   when the results are obtained from residual-centering approach.
##' }
##'
##' @references
##'
##' Schoemann, A. M., & Jorgensen, T. D. (2021). Testing and interpreting
##' latent variable interactions using the \code{semTools} package.
##' \emph{Psych, 3}(3), 322--335. \doi{10.3390/psych3030024}
##'
##' @examples
##'
##' library(lavaan)
##'
##' dat2wayMC <- indProd(dat2way, 1:3, 4:6)
##'
##' model1 <- "
##' f1 =~ x1 + x2 + x3
##' f2 =~ x4 + x5 + x6
##' f12 =~ x1.x4 + x2.x5 + x3.x6
##' f3 =~ x7 + x8 + x9
##' f3 ~ f1 + f2 + f12
##' f12 ~~ 0*f1
##' f12 ~~ 0*f2
##' x1 ~ 0*1
##' x4 ~ 0*1
##' x1.x4 ~ 0*1
##' x7 ~ 0*1
##' f1 ~ NA*1
##' f2 ~ NA*1
##' f12 ~ NA*1
##' f3 ~ NA*1
##' "
##'
##' fitMC2way <- sem(model1, data = dat2wayMC, meanstructure = TRUE)
##' result2wayMC <- probe2WayMC(fitMC2way, nameX = c("f1", "f2", "f12"),
##'                             nameY = "f3", modVar = "f2", valProbe = c(-1, 0, 1))
##' plotProbe(result2wayMC, xlim = c(-2, 2))
##'
##'
##' dat3wayMC <- indProd(dat3way, 1:3, 4:6, 7:9)
##'
##' model3 <- "
##' f1 =~ x1 + x2 + x3
##' f2 =~ x4 + x5 + x6
##' f3 =~ x7 + x8 + x9
##' f12 =~ x1.x4 + x2.x5 + x3.x6
##' f13 =~ x1.x7 + x2.x8 + x3.x9
##' f23 =~ x4.x7 + x5.x8 + x6.x9
##' f123 =~ x1.x4.x7 + x2.x5.x8 + x3.x6.x9
##' f4 =~ x10 + x11 + x12
##' f4 ~ f1 + f2 + f3 + f12 + f13 + f23 + f123
##' f1 ~~ 0*f12
##' f1 ~~ 0*f13
##' f1 ~~ 0*f123
##' f2 ~~ 0*f12
##' f2 ~~ 0*f23
##' f2 ~~ 0*f123
##' f3 ~~ 0*f13
##' f3 ~~ 0*f23
##' f3 ~~ 0*f123
##' f12 ~~ 0*f123
##' f13 ~~ 0*f123
##' f23 ~~ 0*f123
##' x1 ~ 0*1
##' x4 ~ 0*1
##' x7 ~ 0*1
##' x10 ~ 0*1
##' x1.x4 ~ 0*1
##' x1.x7 ~ 0*1
##' x4.x7 ~ 0*1
##' x1.x4.x7 ~ 0*1
##' f1 ~ NA*1
##' f2 ~ NA*1
##' f3 ~ NA*1
##' f12 ~ NA*1
##' f13 ~ NA*1
##' f23 ~ NA*1
##' f123 ~ NA*1
##' f4 ~ NA*1
##' "
##'
##' fitMC3way <- sem(model3, data = dat3wayMC, std.lv = FALSE,
##'                  meanstructure = TRUE)
##' result3wayMC <- probe3WayMC(fitMC3way, nameX = c("f1", "f2", "f3", "f12",
##'                                                  "f13", "f23", "f123"),
##'                             nameY = "f4", modVar = c("f1", "f2"),
##'                             valProbe1 = c(-1, 0, 1), valProbe2 = c(-1, 0, 1))
##' plotProbe(result3wayMC, xlim = c(-2, 2))
##'
##' @export
plotProbe <- function(object, xlim, xlab = "Indepedent Variable",
                      ylab = "Dependent Variable", legend = TRUE,
                      legendArgs = list(), ...) {
  if (length(xlim) != 2)  stop("The x-limit should be specified as a numeric",
                               " vector with the length of 2.")

  # Extract simple slope
  slope <- object$SimpleSlope

  # Check whether the object is the two-way or three-way interaction result
  numInt <- 2
  if (ncol(slope) == 6) numInt <- 3
  estSlope <- slope[, ncol(slope) - 3]

  # Get whether the simple slope is significant. If so, the resulting lines will be
  # shown as red. If not, the line will be black.
  estSlopeSig <- (slope[, ncol(slope)] < 0.05) + 1

  # Extract simple intercept. If the simple intercept is not provided, the intercept
  # will be fixed as 0.
  estIntercept <- NULL
  if (!is.null(object$SimpleIntcept))
    estIntercept <- object$SimpleIntcept[, ncol(slope) - 3]
  if (numInt == 2) {
    if (is.null(legendArgs$title)) legendArgs$title <- colnames(slope)[1]
    if (is.null(legendArgs$legend)) legendArgs$legend <- slope[, 1]
    plotSingleProbe(estSlope, estIntercept, xlim = xlim, xlab = xlab, ylab = ylab,
                    colLine = estSlopeSig, legend = legend,
                    legendArgs = legendArgs, ...)
  } else if (numInt == 3) {
    # Three-way interaction; separate lines for the first moderator, separate graphs
    # for the second moderator
    mod2 <- unique(slope[, 2])
    mod1 <- unique(slope[, 1])

    # Use multiple graphs in a figure
    if (length(mod2) == 2) {
      obj <- par(mfrow = c(1, 2))
    } else if (length(mod2) == 3) {
      obj <- par(mfrow = c(1, 3))
    } else if (length(mod2) > 3) {
      obj <- par(mfrow = c(2, ceiling(length(mod2)/2)))
    } else if (length(mod2) == 1) {
      # Intentionally leaving as blank
    } else stop("Some errors occur")

    for (i in 1:length(mod2)) {
      select <- slope[, 2] == mod2[i]
      if (is.null(legendArgs$title)) legendArgs$title <- colnames(slope)[1]
      if (is.null(legendArgs$legend)) legendArgs$legend <- mod1
      plotSingleProbe(estSlope[select], estIntercept[select], xlim = xlim,
                      xlab = xlab, ylab = ylab, colLine = estSlopeSig[select],
                      main = paste(colnames(slope)[2], "=", mod2[i]),
                      legend = legend, legendArgs = legendArgs, ...)
    }
    if (length(mod2) > 1)
      par(obj)
  } else {
    stop("Please make sure that the object argument is obtained from",
         " 'probe2wayMC', 'probe2wayRC', 'probe3wayMC', or 'probe3wayRC'.")
  }
}



## ----------------
## Hidden Functions
## ----------------

## Find the expected value of the product of two normal variates
## m = the mean of each normal variate
## s = the covariance matrix of all variates
expect2NormProd <- function(m, s) return(prod(m) + s[1, 2])



## Find the expected value of the product of three normal variates
## m = the mean of each normal variate
## s = the covariance matrix of all variates
expect3NormProd <- function(m, s) {
	return(prod(m) + m[3] * s[1, 2] + m[2] * s[1, 3] + m[1] * s[2, 3])
}



## Find the expected value of the product of four normal variates
## m = the mean of each normal variate
## s = the covariance matrix of all variates
expect4NormProd <- function(m, s) {
	first <- prod(m)
	com <- utils::combn(1:4, 2)
	forSecond <- function(draw, meanval, covval, index) {
		draw2 <- setdiff(index, draw)
		prod(meanval[draw2]) * covval[draw[1], draw[2]]
	}
	second <- sum(apply(com, 2, forSecond, meanval=m, covval=s, index=1:4))

	com2 <- com[,1:3] #select only first three terms containing the first element only
	forThird <- function(draw, covval, index) {
		draw2 <- setdiff(index, draw)
		covval[draw[1], draw[2]] * covval[draw2[1], draw2[2]]
	}
	third <- sum(apply(com2, 2, forThird, covval=s, index=1:4))
	return(first + second + third)
}



## Find the expected value of the product of five normal variates
## m = the mean of each normal variate
## s = the covariance matrix of all variates
expect5NormProd <- function(m, s) {
	first <- prod(m)
	com <- utils::combn(1:5, 2)
	forSecond <- function(draw, meanval, covval, index) {
		draw2 <- setdiff(index, draw)
		prod(meanval[draw2]) * covval[draw[1], draw[2]]
	}
	second <- sum(apply(com, 2, forSecond, meanval=m, covval=s, index=1:5))

	com2 <- utils::combn(1:5, 4)
	forThirdOuter <- function(index, m, s, indexall) {
		targetMean <- m[setdiff(indexall, index)]
		cominner <- utils::combn(index, 2)[,1:3] #select only first three terms containing the first element only
		forThirdInner <- function(draw, covval, index) {
			draw2 <- setdiff(index, draw)
			covval[draw[1], draw[2]] * covval[draw2[1], draw2[2]]
		}
		thirdInner <- targetMean * sum(apply(cominner, 2, forThirdInner, covval=s, index=index))
		return(thirdInner)
	}
	third <- sum(apply(com2, 2, forThirdOuter, m=m, s=s, indexall=1:5))
	return(first + second + third)
}



## Find the variance of the product of two normal variates
## m = the mean of each normal variate
## s = the covariance matrix of all variates
var2NormProd <- function(m, s) {
	first <- m[2]^2 * s[1, 1] + m[1]^2 * s[2, 2]
	second <- 2 * m[1] * m[2] * s[1, 2]
	third <- s[1, 1] * s[2, 2]
	fourth <- s[1, 2]^2
	return(first + second + third + fourth)
}



## Find the variance of the product of three normal variates
## m = the mean of each normal variate
## s = the covariance matrix of all variates
var3NormProd <- function(m, s) {
	com <- utils::combn(1:3, 2)
	forFirst <- function(draw, meanval, covval, index) {
		# draw = 2, 3; draw2 = 1
		draw2 <- setdiff(index, draw)
		term1 <- meanval[draw[1]]^2 * meanval[draw[2]]^2 * covval[draw2, draw2]
		term2 <- 2 * meanval[draw2]^2 * meanval[draw[1]] * meanval[draw[2]] * covval[draw[1], draw[2]]
		term3 <- (meanval[draw2]^2 * covval[draw[1], draw[1]] * covval[draw[2], draw[2]]) + (meanval[draw2]^2 * covval[draw[1], draw[2]]^2)
		term4 <- 4 * meanval[draw[1]] * meanval[draw[2]] * covval[draw2, draw2] * covval[draw[1], draw[2]]
		term5 <- 6 * meanval[draw[1]] * meanval[draw[2]] * covval[draw2, draw[1]] * covval[draw2, draw[2]]
		term1 + term2 + term3 + term4 + term5
	}
	first <- sum(apply(com, 2, forFirst, meanval=m, covval=s, index=1:3))
	second <- prod(diag(s))
	third <- 2 * s[3, 3] * s[1, 2]^2 + 2 * s[2, 2] * s[1, 3]^2 + 2 * s[1, 1] * s[2, 3]^2
	fourth <- 8 * s[1, 2] * s[1, 3] * s[2, 3]
	return(first + second + third + fourth)
}



## plotSingleProbe : plot the probing interaction result specific for only one moderator
## estSlope = slope of each line
## estIntercept = intercept of each line
## xlim = the minimum and maximum values of the independent variable (x-axis)
## xlab = the label for the independent variable
## ylab = the lable for the dependent variable
## main = the title of the graph
## colLine = the color of each line
## legend = whether to print a legend
## legendArgs = arguments to pass to legend() function
plotSingleProbe <- function(estSlope, estIntercept = NULL, xlim,
                            xlab = "Indepedent Variable",
                            ylab = "Dependent Variable", main = NULL,
                            colLine = "black", legend = TRUE,
                            legendArgs = list(), ...) {
	if (is.null(estIntercept)) estIntercept <- rep(0, length(estSlope))
	if (length(colLine) == 1) colLine <- rep(colLine, length(estSlope))
	lower <- estIntercept + (xlim[1] * estSlope)
	upper <- estIntercept + (xlim[2] * estSlope)
	ylim <- c(min(c(lower, upper)), max(c(lower, upper)))
	plot(cbind(xlim, ylim), xlim = xlim, ylim = ylim, type = "n",
	     xlab = xlab, ylab = ylab, main = main, ...)
	for (i in 1:length(estSlope)) {
		lines(cbind(xlim, c(lower[i], upper[i])),
		      col = colLine[i], lwd = 1.5, lty = i)
	}
	if (legend) {
		positionX <- 0.25
		if (all(estSlope > 0)) positionX <- 0.01
		if (all(estSlope < 0)) positionX <- 0.50
		if (is.null(legendArgs$x)) legendArgs$x <- positionX * (xlim[2] - xlim[1]) + xlim[1]
		if (is.null(legendArgs$y)) legendArgs$y <- 0.99 * (ylim[2] - ylim[1]) + ylim[1]
		if (is.null(legendArgs$col)) legendArgs$col <- colLine
		if (is.null(legendArgs$lty)) legendArgs$lty <- 1:length(estSlope)
		do.call(graphics::legend, legendArgs)
	}
}


