### Sunthud Pornprasertmanit
### Last updated: 10 January 2021


##' Modification indices and their power approach for model fit evaluation
##'
##' The model fit evaluation approach using modification indices and expected
##' parameter changes.
##'
##' In the lavaan object, one can inspect the modification indices and expected
##' parameter changes. Those values can be used to evaluate model fit by two
##' methods.
##'
##' First, Saris, Satorra, and van der Veld (2009, pp. 570-573) used the power
##' to detect modification indices and expected parameter changes to evaluate
##' model fit. First, one should evaluate whether the modification index of each
##' parameter is significant. Second, one should evaluate whether the power to
##' detect a target expected parameter change is high enough. If the
##' modification index is not significant and the power is high, there is no
##' misspecification. If the modification index is significant and the power is
##' low, the fixed parameter is misspecified. If the modification index is
##' significant and the power is high, the expected parameter change is
##' investigated. If the expected parameter change is large (greater than the
##' the target expected parameter change), the parameter is misspecified. If the
##' expected parameter change is low (lower than the target expected parameter
##' change), the parameter is not misspecificied. If the modification index is
##' not significant and the power is low, the decision is inconclusive.
##'
##' Second, the confidence intervals of the expected parameter changes are
##' formed. These confidence intervals are compared with the range of trivial
##' misspecification, which could be (-\code{delta}, \code{delta}) or (0,
##' \code{delta}) for nonnegative parameters. If the confidence intervals are
##' outside of the range of trivial misspecification, the fixed parameters are
##' severely misspecified. If the confidence intervals are inside the range of
##' trivial misspecification, the fixed parameters are trivially misspecified.
##' If confidence intervals are overlapped the range of trivial
##' misspecification, the decision is inconclusive.
##'
##' @aliases miPowerFit miPowerFit
##' @importFrom lavaan lavInspect
##' @importFrom stats qnorm qchisq pchisq
##'
##' @param lavaanObj The lavaan model object used to evaluate model fit
##' @param stdLoad The amount of standardized factor loading that one would like
##' to be detected (rejected). The default value is 0.4, which is suggested by
##' Saris and colleagues (2009, p. 571).
##' @param cor The amount of factor or error correlations that one would like to
##' be detected (rejected). The default value is 0.1, which is suggested by
##' Saris and colleagues (2009, p. 571).
##' @param stdBeta The amount of standardized regression coefficients that one
##' would like to be detected (rejected). The default value is 0.1, which is
##' suggested by Saris and colleagues (2009, p. 571).
##' @param intcept The amount of standardized intercept (similar to Cohen's
##' \emph{d} that one would like to be detected (rejected). The default value is
##' 0.2, which is equivalent to a low effect size proposed by Cohen (1988,
##' 1992).
##' @param stdDelta The vector of the standardized parameters that one would
##' like to be detected (rejected). If this argument is specified, the value
##' here will overwrite the other arguments above. The order of the vector must
##' be the same as the row order from modification indices from the
##' \code{lavaan} object. If a single value is specified, the value will be
##' applied to all parameters.
##' @param delta The vector of the unstandardized parameters that one would like
##' to be detected (rejected). If this argument is specified, the value here
##' will overwrite the other arguments above. The order of the vector must be
##' the same as the row order from modification indices from the \code{lavaan}
##' object. If a single value is specified, the value will be applied to all
##' parameters.
##' @param cilevel The confidence level of the confidence interval of expected
##' parameter changes. The confidence intervals are used in the equivalence
##' testing.
##' @return A data frame with these variables:
##'  \enumerate{
##'   \item lhs: The left-hand side variable, with respect to the operator in
##'    in the lavaan \code{\link[lavaan]{model.syntax}}
##'   \item op: The lavaan syntax operator: "~~" represents covariance,
##'     "=~" represents factor loading, "~" represents regression, and
##'     "~1" represents intercept.
##'   \item rhs: The right-hand side variable
##'   \item group: The level of the group variable for the parameter in question
##'   \item mi: The modification index of the fixed parameter
##'   \item epc: The expected parameter change if the parameter is freely
##'    estimated
##'   \item target.epc: The target expected parameter change that represents
##'    the minimum size of misspecification that one would like to be detected
##'    by the test with a high power
##'   \item std.epc: The standardized expected parameter change if the parameter
##'    is freely estimated
##'   \item std.target.epc: The standardized target expected parameter change
##'   \item significant.mi: Represents whether the modification index value is
##'     significant
##'   \item high.power: Represents whether the power is enough to detect the
##'    target expected parameter change
##'   \item decision.pow: The decision whether the parameter is misspecified
##'    or not based on Saris et al's method: \code{"M"} represents the parameter
##'    is misspecified, \code{"NM"} represents the parameter is not misspecified,
##'    \code{"EPC:M"} represents the parameter is misspecified decided by
##'    checking the expected parameter change value, \code{"EPC:NM"} represents
##'    the parameter is not misspecified decided by checking the expected
##'    parameter change value, and \code{"I"} represents the decision is
##'    inconclusive.
##'   \item se.epc: The standard errors of the expected parameter changes.
##'   \item lower.epc: The lower bound of the confidence interval of expected
##'    parameter changes.
##'   \item upper.epc: The upper bound of the confidence interval of expected
##'    parameter changes.
##'   \item lower.std.epc: The lower bound of the confidence interval of
##'    standardized expected parameter changes.
##'   \item upper.std.epc: The upper bound of the confidence interval of
##'    standardized expected parameter changes.
##'   \item decision.ci: The decision whether the parameter is misspecified or
##'    not based on the confidence interval method: \code{"M"} represents the
##'    parameter is misspecified, \code{"NM"} represents the parameter is not
##'    misspecified, and \code{"I"} represents the decision is inconclusive.
##' }
##'
##'  The row numbers matches with the results obtained from the
##'  \code{inspect(object, "mi")} function.
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##' @seealso \code{\link{moreFitIndices}} For the additional fit indices
##' information
##' @references Cohen, J. (1988). \emph{Statistical power analysis for the
##' behavioral sciences} (2nd ed.). Hillsdale, NJ: Erlbaum.
##'
##' Cohen, J. (1992). A power primer. \emph{Psychological Bulletin, 112}(1),
##' 155--159. \doi{10.1037/0033-2909.112.1.155}
##'
##' Saris, W. E., Satorra, A., & van der Veld, W. M. (2009). Testing structural
##' equation models or detection of misspecifications? \emph{Structural Equation
##' Modeling, 16}(4), 561--582. \doi{10.1080/10705510903203433}
##' @examples
##'
##' library(lavaan)
##'
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ x7 + x8 + x9 '
##'
##' fit <- cfa(HS.model, data = HolzingerSwineford1939,
##'            group = "sex", meanstructure = TRUE)
##' miPowerFit(fit)
##'
##' model <- '
##'   # latent variable definitions
##'      ind60 =~ x1 + x2 + x3
##'      dem60 =~ y1 + a*y2 + b*y3 + c*y4
##'      dem65 =~ y5 + a*y6 + b*y7 + c*y8
##'
##'   # regressions
##'     dem60 ~ ind60
##'     dem65 ~ ind60 + dem60
##'
##'   # residual correlations
##'     y1 ~~ y5
##'     y2 ~~ y4 + y6
##'     y3 ~~ y7
##'     y4 ~~ y8
##'     y6 ~~ y8
##' '
##' fit2 <- sem(model, data = PoliticalDemocracy, meanstructure = TRUE)
##' miPowerFit(fit2, stdLoad = 0.3, cor = 0.2, stdBeta = 0.2, intcept = 0.5)
##'
##' @export
miPowerFit <- function(lavaanObj, stdLoad = 0.4, cor = 0.1, stdBeta = 0.1,
                       intcept = 0.2, stdDelta = NULL, delta = NULL,
                       cilevel = 0.90) {
	mi <- lavInspect(lavaanObj, "mi")
	mi <- mi[mi$op != "==",]
	sigma <- mi[,"epc"] / sqrt(mi[,"mi"])
	if (is.null(delta)) {
		if (is.null(stdDelta))
		  stdDelta <- getTrivialEpc(mi, stdLoad = stdLoad, cor = cor,
		                            stdBeta = stdBeta, intcept = intcept)
		if (length(stdDelta) == 1) stdDelta <- rep(stdDelta, nrow(mi))
		delta <- unstandardizeEpc(mi, stdDelta, findTotalVar(lavaanObj))
	}
	if (length(delta) == 1) delta <- rep(delta, nrow(mi))
	ncp <- (delta / sigma)^2
	alpha <- 0.05
	desiredPow <- 0.80
	cutoff <- qchisq(1 - alpha, df = 1)
	pow <- 1 - pchisq(cutoff, df = 1, ncp = ncp)
	sigMI <- mi[,"mi"] > cutoff
	highPow <- pow > desiredPow
	group <- rep(1, nrow(mi))
	if ("group" %in% colnames(mi)) group <- mi[ , "group"]
	decision <- mapply(decisionMIPow, sigMI = sigMI, highPow = highPow,
	                   epc = mi[ , "epc"], trivialEpc = delta)
	if (is.null(stdDelta)) stdDelta <- standardizeEpc(mi, findTotalVar(lavaanObj),
	                                                  delta = delta)
	result <- cbind(mi[ , 1:3], group, as.numeric(mi[ , "mi"]), mi[ , "epc"],
	                delta, standardizeEpc(mi, findTotalVar(lavaanObj)),
	                stdDelta, sigMI, highPow, decision)
	# New method
	crit <- abs(qnorm((1 - cilevel)/2))
	seepc <- abs(result[,6]) / sqrt(abs(result[,5]))
	lowerepc <- result[,6] - crit * seepc
	upperepc <- result[,6] + crit * seepc
	stdlowerepc <- standardizeEpc(mi, findTotalVar(lavaanObj), delta = lowerepc)
	stdupperepc <- standardizeEpc(mi, findTotalVar(lavaanObj), delta = upperepc)
	isVar <- mi[,"op"] == "~~" & mi[,"lhs"] == mi[,"rhs"]
	decisionci <- mapply(decisionCIEpc, targetval = as.numeric(stdDelta),
	                     lower = stdlowerepc, upper = stdupperepc,
	                     positiveonly = isVar)
	result <- cbind(result, seepc, lowerepc, upperepc, stdlowerepc,
	                stdupperepc, decisionci)
	result <- result[!is.na(decision), ]
	colnames(result) <- c("lhs","op","rhs","group","mi","epc","target.epc",
	                      "std.epc","std.target.epc","significant.mi",
	                      "high.power","decision.pow","se.epc","lower.epc",
	                      "upper.epc","lower.std.epc","upper.std.epc","decision.ci")
	result <- format(result, scientific = FALSE, digits = 4)
	return(result)
}



## ----------------
## Hidden Functions
## ----------------

## totalFacVar: Find total factor variances when regression coeffient matrix
## and factor residual covariance matrix are specified

totalFacVar <- function(beta, psi) {
	ID <- diag(nrow(psi))
  total <- solve(ID - beta) %*% psi %*% t(solve(ID - beta))
  return(diag(total))
}

## findTotalVar: find the total indicator and factor variances

##' @importFrom lavaan lavInspect
findTotalVar <- function(lavaanObj) {
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

## getTrivialEpc: find the trivial misspecified expected parameter changes
## given the type of parameters in each row of modification indices

getTrivialEpc <- function(mi, stdLoad=0.4, cor=0.1, stdBeta=0.1, intcept=0.2) {
	op <- mi[,"op"]
	result <- gsub("=~", stdLoad, op)
	result <- gsub("~~", cor, result)
	result <- gsub("~1", intcept, result)
	result <- gsub("~", stdBeta, result)
	return(result)
}

## unstandardizeEpc: Transform from standardized EPC to unstandardized EPC

unstandardizeEpc <- function(mi, delta, totalVar) {
	name <- names(totalVar[[1]])
	lhsPos <- match(mi[,"lhs"], name)
	rhsPos <- match(mi[,"rhs"], name)
	group <- rep(1, nrow(mi))
	if("group" %in% colnames(mi)) group <- mi[,"group"]
	getVar <- function(pos, group) totalVar[[group]][pos]
	lhsVar <- mapply(getVar, pos=lhsPos, group=group)
	rhsVar <- mapply(getVar, pos=rhsPos, group=group)
	FUN <- function(op, lhsVar, rhsVar, delta) {
		if(op == "|") return(NA)
		lhsSD <- sqrt(lhsVar)
		rhsSD <- sqrt(rhsVar)
		if(!is.numeric(delta)) delta <- as.numeric(delta)
		if(op == "=~") {
			return((rhsSD * delta) / lhsSD)
		} else if (op == "~~") {
			return(lhsSD * delta * rhsSD)
		} else if (op == "~1") {
			return(lhsSD * delta)
		} else if (op == "~") {
			return((lhsSD * delta) / rhsSD)
		} else {
			return(NA)
		}
	}
	unstdDelta <- mapply(FUN, op=mi[,"op"], lhsVar=lhsVar, rhsVar=rhsVar, delta=delta)
	return(unstdDelta)
}

## unstandardizeEpc: Transform from unstandardized EPC to standardized EPC.
## If delta is null, the unstandardized epc from the modification indices
## data.frame are used

standardizeEpc <- function(mi, totalVar, delta = NULL) {
	if(is.null(delta)) delta <- mi[,"epc"]
	name <- names(totalVar[[1]])
	lhsPos <- match(mi[,"lhs"], name)
	rhsPos <- match(mi[,"rhs"], name)
	group <- rep(1, nrow(mi))
	if("group" %in% colnames(mi)) group <- mi[,"group"]
	getVar <- function(pos, group) totalVar[[group]][pos]
	lhsVar <- mapply(getVar, pos=lhsPos, group=group)
	rhsVar <- mapply(getVar, pos=rhsPos, group=group)
	FUN <- function(op, lhsVar, rhsVar, delta) {
		lhsSD <- sqrt(lhsVar)
		rhsSD <- sqrt(rhsVar)
		if(!is.numeric(delta)) delta <- as.numeric(delta)
		if(op == "=~") {
			#stdload = beta * sdlatent / sdindicator = beta * lhs / rhs
			return((delta / rhsSD) * lhsSD)
		} else if (op == "~~") {
			#r = cov / (sd1 * sd2)
			return(delta / (lhsSD * rhsSD))
		} else if (op == "~1") {
			#d = meanDiff/sd
			return(delta / lhsSD)
		} else if (op == "~") {
			#beta = b * sdX / sdY = b * rhs / lhs
			return((delta / lhsSD) * rhsSD)
		} else {
			return(NA)
		}
	}
	stdDelta <- mapply(FUN, op=mi[,"op"], lhsVar=lhsVar, rhsVar=rhsVar, delta=delta)
	return(stdDelta)
}

## decisionMIPow: provide the decision given the significance of modification
## indices and power to detect trivial misspecification

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

decisionCIEpc <- function(targetval, lower, upper, positiveonly = FALSE) {
	if(is.na(lower) | is.na(upper)) return(NA)
	if(positiveonly) {
		if(lower > targetval) {
			return("M")
		} else if (upper < targetval) {
			return("NM")
		} else {
			return("I")
		}
	} else {
		negtargetval <- -targetval
		if(lower > targetval | upper < negtargetval) {
			return("M")
		} else if (upper < targetval & negtargetval < lower) {
			return("NM")
		} else {
			return("I")
		}
	}
}


