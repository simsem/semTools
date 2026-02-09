### Sunthud Pornprasertmanit
### Last updated: 9 February 2026


##' Partial Measurement Invariance Testing Across Groups
##'
##' This test will provide partial invariance testing by (a) freeing a parameter
##' one-by-one from nested model and compare with the original nested model or
##' (b) fixing (or constraining) a parameter one-by-one from the parent model
##' and compare with the original parent model. This function only works with
##' congeneric models. The `partialInvariance` is used for continuous
##' variable. The `partialInvarianceCat` is used for categorical variables.
##'
##' There are four types of partial invariance testing:
##'
##' \itemize{
##'  \item Partial weak invariance. The model named `fit.configural`
##' from the list of models is compared with the model named `fit.loadings`.
##' Each loading will be freed or fixed from the metric and configural
##' invariance models respectively. The modified models are compared with the
##' original model. Note that the objects in the list of models must have the
##' names of `"fit.configural"` and `"fit.loadings"`. Users may use "metric",
##' "weak", "loading", or "loadings" in the `type` argument. Note that, for
##' testing invariance on marker variables, other variables will be assigned as
##' marker variables automatically.
##'
##'  \item Partial strong invariance. The model
##' named `fit.loadings` from the list of models is compared with the model
##' named either `fit.intercepts` or 'fit.thresholds'. Each intercept will be
##' freed or fixed from the scalar and metric invariance models respectively.
##' The modified models are compared with the original model. Note that the
##' objects in the list of models must have the names of "fit.loadings" and
##' either "fit.intercepts" or "fit.thresholds". Users may use "scalar",
##' "strong", "intercept", "intercepts", "threshold", or "thresholds" in the
##' `type` argument. Note that, for testing invariance on marker variables,
##' other variables will be assigned as marker variables automatically. Note
##' that if all variables are dichotomous, scalar invariance testing is not
##' available.
##'
##'  \item Partial strict invariance. The model named either
##' 'fit.intercepts' or 'fit.thresholds' (or 'fit.loadings') from the list of
##' models is compared with the model named 'fit.residuals'. Each residual
##' variance will be freed or fixed from the strict and scalar (or metric)
##' invariance models respectively. The modified models are compared with the
##' original model. Note that the objects in the list of models must have the
##' names of "fit.residuals" and either "fit.intercepts", "fit.thresholds", or
##' "fit.loadings". Users may use "strict", "residual", "residuals", "error", or
##' "errors" in the `type` argument.
##'
##'  \item Partial mean invariance. The
##' model named either 'fit.intercepts' or 'fit.thresholds' (or 'fit.residuals'
##' or 'fit.loadings') from the list of models is compared with the model named
##' 'fit.means'. Each factor mean will be freed or fixed from the means and
##' scalar (or strict or metric) invariance models respectively. The modified
##' models are compared with the original model. Note that the objects in the
##' list of models must have the names of "fit.means" and either
##' "fit.residuals", "fit.intercepts", "fit.thresholds", or "fit.loadings".
##' Users may use "means" or "mean" in the `type` argument. }
##'
##' Two types of comparisons are used in this function:
##' \enumerate{
##' \item `free`: The nested model is used as a template. Then, one
##' parameter indicating the differences between two models is free. The new
##' model is compared with the nested model. This process is repeated for all
##' differences between two models. The likelihood-ratio test and the difference
##' in CFI are provided.
##'
##' \item `fix`: The parent model is used as a template. Then, one parameter
##' indicating the differences between two models is fixed or constrained to be
##' equal to other parameters. The new model is then compared with the parent
##' model. This process is repeated for all differences between two models. The
##' likelihood-ratio test and the difference in CFI are provided.
##'
##' \item `wald`: This method is similar to the `fix` method. However,
##' instead of building a new model and compare them with likelihood-ratio test,
##' multivariate wald test is used to compare equality between parameter
##' estimates. See [lavaan::lavTestWald()] for further details. Note
##' that if any rows of the contrast cannot be summed to 0, the Wald test is not
##' provided, such as comparing two means where one of the means is fixed as 0.
##' This test statistic is not as accurate as likelihood-ratio test provided in
##' `fix`. I provide it here in case that likelihood-ratio test fails to
##' converge.
##' }
##'
##' Note that this function does not adjust for the inflated Type I error rate
##' from multiple tests. The degree of freedom of all tests would be the number
##' of groups minus 1.
##'
##' The details of standardized estimates and the effect size used for each
##' parameters are provided in the vignettes by running
##' `vignette("partialInvariance")`.
##'
##' @importFrom lavaan lavInspect parTable
##' @aliases partialInvariance partialInvarianceCat
##'
##' @param fit A list of models for invariance testing. Each model should be
##'   assigned by appropriate names (see details).
##' @param type The types of invariance testing: "metric", "scalar", "strict",
##'   or "means"
##' @param free A vector of variable names that are free across groups in
##'   advance. If partial mean invariance is tested, this argument represents a
##'   vector of factor names that are free across groups.
##' @param fix A vector of variable names that are constrained to be equal
##'   across groups in advance. If partial mean invariance is tested, this
##'   argument represents a vector of factor names that are fixed across groups.
##' @param refgroup The reference group used to make the effect size comparison
##'   with the other groups.
##' @param poolvar If `TRUE`, the variances are pooled across group for
##'   standardization. Otherwise, the variances of the reference group are used
##'   for standardization.
##' @param p.adjust The method used to adjust p values. See
##'   [stats::p.adjust()] for the options for adjusting p values. The
##'   default is to not use any corrections.
##' @param fbound The z-scores of factor that is used to calculate the effect
##'   size of the loading difference proposed by Millsap and Olivera-Aguilar
##'   (2012).
##' @param return.fit Return the submodels fitted by this function
##' @param method The method used to calculate likelihood ratio test. See
##'   [lavaan::lavTestLRT()] for available options
##'
##' @return A list of results are provided. The list will consists of at least
##' two elements:
##' \enumerate{
##'  \item `estimates`: The results of parameter estimates including pooled
##'   estimates (`poolest`), the estimates for each group, standardized
##'   estimates for each group (`std`), the difference in standardized
##'   values, and the effect size statistic (*q* for factor loading
##'   difference and *h* for error variance difference). See the details of
##'   this effect size statistic by running `vignette("partialInvariance")`.
##'   In the `partialInvariance` function, the additional effect statistics
##'   proposed by Millsap and Olivera-Aguilar (2012) are provided. For factor
##'   loading, the additional outputs are the observed mean difference
##'   (`diff_mean`), the mean difference if factor scores are low
##'   (`low_fscore`), and the mean difference if factor scores are high
##'   (`high_fscore`). The low factor score is calculated by (a) finding the
##'   factor scores that its *z* score equals -`bound` (the default is
##'   \eqn{-2}) from all groups and (b) picking the minimum value among the
##'   factor scores. The high factor score is calculated by (a) finding the
##'   factor scores that its *z* score equals `bound` (default = 2)
##'   from all groups and (b) picking the maximum value among the factor scores.
##'   For measurement intercepts, the additional outputs are the observed means
##'   difference (`diff_mean`) and the proportion of the differences in the
##'   intercepts over the observed means differences (`propdiff`). For error
##'   variances, the additional outputs are the proportion of the difference in
##'   error variances over the difference in observed variances (`propdiff`).
##'
##'  \item `results`: Statistical tests as well as the change in CFI are
##'   provided. \eqn{\chi^2} and *p* value are provided for all methods.
##'
##'  \item `models`: The submodels used in the `free` and `fix`
##'   methods, as well as the nested and parent models. The nested and parent
##'   models will be changed from the original models if `free` or
##'   `fit` arguments are specified.
##' }
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' @references Millsap, R. E., & Olivera-Aguilar, M. (2012). Investigating
##' measurement invariance using confirmatory factor analysis. In R. H. Hoyle
##' (Ed.), *Handbook of structural equation modeling* (pp. 380--392). New
##' York, NY: Guilford.
##'
##' @examples
##'
##' ## Conduct weak invariance testing manually by using fixed-factor
##' ## method of scale identification
##'
##' library(lavaan)
##'
##' conf <- "
##' f1 =~ NA*x1 + x2 + x3
##' f2 =~ NA*x4 + x5 + x6
##' f1 ~~ c(1, 1)*f1
##' f2 ~~ c(1, 1)*f2
##' "
##'
##' weak <- "
##' f1 =~ NA*x1 + x2 + x3
##' f2 =~ NA*x4 + x5 + x6
##' f1 ~~ c(1, NA)*f1
##' f2 ~~ c(1, NA)*f2
##' "
##'
##' configural <- cfa(conf, data = HolzingerSwineford1939, std.lv = TRUE, group="school")
##' weak <- cfa(weak, data = HolzingerSwineford1939, group="school", group.equal="loadings")
##' models <- list(fit.configural = configural, fit.loadings = weak)
##' partialInvariance(models, "metric")
##'
# \donttest{
# partialInvariance(models, "metric", free = "x5") # "x5" is free across groups in advance
# partialInvariance(models, "metric", fix = "x4") # "x4" is fixed across groups in advance
#
# ## Use the result from the measurementInvariance function
# HW.model <- ' visual =~ x1 + x2 + x3
#               textual =~ x4 + x5 + x6
#               speed =~ x7 + x8 + x9 '
#
# models2 <- measurementInvariance(model = HW.model, data=HolzingerSwineford1939,
#                                  group="school")
# partialInvariance(models2, "scalar")
#
# ## Conduct weak invariance testing manually by using fixed-factor
# ## method of scale identification for dichotomous variables
#
# f <- rnorm(1000, 0, 1)
# u1 <- 0.9*f + rnorm(1000, 1, sqrt(0.19))
# u2 <- 0.8*f + rnorm(1000, 1, sqrt(0.36))
# u3 <- 0.6*f + rnorm(1000, 1, sqrt(0.64))
# u4 <- 0.7*f + rnorm(1000, 1, sqrt(0.51))
# u1 <- as.numeric(cut(u1, breaks = c(-Inf, 0, Inf)))
# u2 <- as.numeric(cut(u2, breaks = c(-Inf, 0.5, Inf)))
# u3 <- as.numeric(cut(u3, breaks = c(-Inf, 0, Inf)))
# u4 <- as.numeric(cut(u4, breaks = c(-Inf, -0.5, Inf)))
# g <- rep(c(1, 2), 500)
# dat2 <- data.frame(u1, u2, u3, u4, g)
#
# configural2 <- "
# f1 =~ NA*u1 + u2 + u3 + u4
# u1 | c(t11, t11)*t1
# u2 | c(t21, t21)*t1
# u3 | c(t31, t31)*t1
# u4 | c(t41, t41)*t1
# f1 ~~ c(1, 1)*f1
# f1 ~ c(0, NA)*1
# u1 ~~ c(1, 1)*u1
# u2 ~~ c(1, NA)*u2
# u3 ~~ c(1, NA)*u3
# u4 ~~ c(1, NA)*u4
# "
#
# outConfigural2 <- cfa(configural2, data = dat2, group = "g",
#                       parameterization = "theta", estimator = "wlsmv",
#                       ordered = c("u1", "u2", "u3", "u4"))
#
# weak2 <- "
# f1 =~ NA*u1 + c(f11, f11)*u1 + c(f21, f21)*u2 + c(f31, f31)*u3 + c(f41, f41)*u4
# u1 | c(t11, t11)*t1
# u2 | c(t21, t21)*t1
# u3 | c(t31, t31)*t1
# u4 | c(t41, t41)*t1
# f1 ~~ c(1, NA)*f1
# f1 ~ c(0, NA)*1
# u1 ~~ c(1, 1)*u1
# u2 ~~ c(1, NA)*u2
# u3 ~~ c(1, NA)*u3
# u4 ~~ c(1, NA)*u4
# "
#
# outWeak2 <- cfa(weak2, data = dat2, group = "g", parameterization = "theta",
#                 estimator = "wlsmv", ordered = c("u1", "u2", "u3", "u4"))
# modelsCat <- list(fit.configural = outConfigural2, fit.loadings = outWeak2)
#
# partialInvarianceCat(modelsCat, type = "metric")
#
# partialInvarianceCat(modelsCat, type = "metric", free = "u2")
# partialInvarianceCat(modelsCat, type = "metric", fix = "u3")
#
# ## Use the result from the measurementInvarianceCat function
#
# model <- ' f1 =~ u1 + u2 + u3 + u4
#            f2 =~ u5 + u6 + u7 + u8'
#
# modelsCat2 <- measurementInvarianceCat(model = model, data = datCat, group = "g",
# 	                                      parameterization = "theta",
# 	                                      estimator = "wlsmv", strict = TRUE)
#
# partialInvarianceCat(modelsCat2, type = "scalar")
# }
##'
##' @export
partialInvariance <- function(fit, type, free = NULL, fix = NULL, refgroup = 1,
                              poolvar = TRUE, p.adjust = "none", fbound = 2,
                              return.fit = FALSE, method = "satorra.bentler.2001") {

	type <- tolower(type)
	numType <- 0
	fit1 <- fit0 <- NULL
	# fit0 = Nested model, fit1 = Parent model
	if(type %in% c("metric", "weak", "loading", "loadings")) {
		numType <- 1
		if(all(c("fit.configural", "fit.loadings") %in% names(fit))) {
			fit1 <- fit$fit.configural
			fit0 <- fit$fit.loadings
		} else {
			stop("The elements named 'fit.configural' and 'fit.loadings' are needed in the 'fit' argument")
		}
	} else if (type %in% c("scalar", "strong", "intercept", "intercepts", "threshold", "thresholds")) {
		numType <- 2
		if(all(c("fit.loadings", "fit.intercepts") %in% names(fit))) {
			fit1 <- fit$fit.loadings
			fit0 <- fit$fit.intercepts
		} else {
			stop("The elements named 'fit.loadings' and 'fit.intercepts' are needed in the 'fit' argument")
		}
	} else if (type %in% c("strict", "residual", "residuals", "error", "errors")) {
		numType <- 3
		if(all(c("fit.intercepts", "fit.residuals") %in% names(fit))) {
			fit1 <- fit$fit.intercepts
			fit0 <- fit$fit.residuals
		} else {
			stop("The elements named 'fit.intercepts' and 'fit.residuals' are needed in the 'fit' argument")
		}
	} else if (type %in% c("means", "mean")) {
		numType <- 4
		if("fit.means" %in% names(fit)) {
			fit0 <- fit$fit.means
			if("fit.residuals" %in% names(fit)) {
				fit1 <- fit$fit.residuals
			} else if ("fit.intercepts" %in% names(fit)) {
				fit1 <- fit$fit.intercepts
			} else {
				stop("The elements named either 'fit.residuals' or 'fit.intercepts	' is needed in the 'fit' argument")
			}
		} else {
			stop("The elements named 'fit.means' is needed in the 'fit' argument")
		}
	} else {
		stop("Please specify the correct type of measurement invariance. See the help page.")
	}
	pt1 <- parTable(fit1)
	pt0 <- parTable(fit0)
	pt0$start <- pt0$est <- pt0$se <- NULL
	pt1$start <- pt1$est <- pt1$se <- NULL
	pt1$label[substr(pt1$label, 1, 1) == "." & substr(pt1$label, nchar(pt1$label), nchar(pt1$label)) == "."] <- ""
	pt0$label[substr(pt0$label, 1, 1) == "." & substr(pt0$label, nchar(pt0$label), nchar(pt0$label)) == "."] <- ""
	namept1 <- paramNameFromPt(pt1)
	namept0 <- paramNameFromPt(pt0)
	if(length(table(table(pt0$rhs[pt0$op == "=~"]))) != 1) stop("The model is not congeneric. This function does not support non-congeneric model.")
	varfree <- varnames <- unique(pt0$rhs[pt0$op == "=~"])
	facnames <- unique(pt0$lhs[(pt0$op == "=~") & (pt0$rhs %in% varnames)])
	facrepresent <- table(pt0$lhs[(pt0$op == "=~") & (pt0$rhs %in% varnames)], pt0$rhs[(pt0$op == "=~") & (pt0$rhs %in% varnames)])
	if(any(apply(facrepresent, 2, function(x) sum(x != 0)) > 1)) stop("The model is not congeneric. This function does not support non-congeneric model.")
	facList <- list()
	for(i in 1:nrow(facrepresent)) {
		facList[[i]] <- colnames(facrepresent)[facrepresent[i,] > 0]
	}
	names(facList) <- rownames(facrepresent)
	facList <- facList[match(names(facList), facnames)]
	fixLoadingFac <- list()
	for(i in seq_along(facList)) {
		select <- pt1$lhs == names(facList)[i] & pt1$op == "=~" & pt1$rhs %in% facList[[i]] & pt1$group == 1 & pt1$free == 0 & (!is.na(pt1$ustart) & pt1$ustart > 0)
		fixLoadingFac[[i]] <- pt1$rhs[select]
	}
	names(fixLoadingFac) <- names(facList)
	fixIntceptFac <- list()
	for(i in seq_along(facList)) {
		select <- pt1$op == "~1" & pt1$rhs %in% facList[[i]] & pt1$group == 1 & pt1$free == 0
		fixIntceptFac[[i]] <- pt1$rhs[select]
	}
	names(fixIntceptFac) <- names(facList)

	ngroups <- max(pt0$group)
	neach <- lavInspect(fit0, "nobs")
	groupvar <- lavInspect(fit0, "group")
	grouplab <- lavInspect(fit0, "group.label")
	if(!is.numeric(refgroup)) refgroup <- which(refgroup == grouplab)
	grouporder <- 1:ngroups
	grouporder <- c(refgroup, setdiff(grouporder, refgroup))
	grouplaborder <- grouplab[grouporder]
	complab <- paste(grouplaborder[2:ngroups], "vs.", grouplaborder[1])
	if(ngroups <= 1) stop("Well, the number of groups is 1. Measurement invariance across 'groups' cannot be done.")

	if(numType == 4) {
		if(!all(c(free, fix) %in% facnames)) stop("'free' and 'fix' arguments should consist of factor names because mean invariance is tested.")
	} else {
		if(!all(c(free, fix) %in% varnames)) stop("'free' and 'fix' arguments should consist of variable names.")
	}
	result <- fixCon <- freeCon <- NULL
	estimates <- NULL
	listFreeCon <- listFixCon <- list()
	beta <- lavaan::coef(fit1)
	beta0 <- lavaan::coef(fit0)
	waldMat <- matrix(0, ngroups - 1, length(beta))
	if(numType == 1) {
		if(!is.null(free) | !is.null(fix)) {
			if(!is.null(fix)) {
				facinfix <- findFactor(fix, facList)
				dup <- duplicated(facinfix)
				for(i in seq_along(fix)) {
					if(dup[i]) {
						pt0 <- constrainParTable(pt0, facinfix[i], "=~", fix[i], 1:ngroups)
						pt1 <- constrainParTable(pt1, facinfix[i], "=~", fix[i], 1:ngroups)
					} else {
						oldmarker <- fixLoadingFac[[facinfix[i]]]
						if(length(oldmarker) > 0) {
							oldmarkerval <- pt1$ustart[pt1$lhs == facinfix[i] & pt1$op == "=~" & pt1$rhs == oldmarker & pt1$group == 1]
							if(oldmarker == fix[i]) {
								pt0 <- fixParTable(pt0, facinfix[i], "=~", fix[i], 1:ngroups, oldmarkerval)
								pt1 <- fixParTable(pt1, facinfix[i], "=~", fix[i], 1:ngroups, oldmarkerval)
 							} else {
								pt0 <- freeParTable(pt0, facinfix[i], "=~", oldmarker, 1:ngroups)
								pt0 <- constrainParTable(pt0, facinfix[i], "=~", oldmarker, 1:ngroups)
								pt1 <- freeParTable(pt1, facinfix[i], "=~", oldmarker, 1:ngroups)
								pt0 <- fixParTable(pt0, facinfix[i], "=~", fix[i], 1:ngroups, oldmarkerval)
								pt1 <- fixParTable(pt1, facinfix[i], "=~", fix[i], 1:ngroups, oldmarkerval)
								fixLoadingFac[[facinfix[i]]] <- fix[i]
							}
						} else {
							pt0 <- constrainParTable(pt0, facinfix[i], "=~", fix[i], 1:ngroups)
							pt1 <- constrainParTable(pt1, facinfix[i], "=~", fix[i], 1:ngroups)
						}
					}
				}
			}
			if(!is.null(free)) {
				facinfree <- findFactor(free, facList)
				for(i in seq_along(free)) {
					# Need to change marker variable if fixed
					oldmarker <- fixLoadingFac[[facinfree[i]]]
					if(length(oldmarker) > 0 && oldmarker == free[i]) {
						oldmarkerval <- pt1$ustart[pt1$lhs == facinfix[i] & pt1$op == "=~" & pt1$rhs == oldmarker & pt1$group == 1]
						candidatemarker <- setdiff(facList[[facinfree[i]]], free[i])[1]
						pt0 <- freeParTable(pt0, facinfree[i], "=~", free[i], 1:ngroups)
						pt1 <- freeParTable(pt1, facinfree[i], "=~", free[i], 1:ngroups)
						pt0 <- fixParTable(pt0, facinfix[i], "=~", candidatemarker, 1:ngroups, oldmarkerval)
						pt1 <- fixParTable(pt1, facinfix[i], "=~", candidatemarker, 1:ngroups, oldmarkerval)
						fixLoadingFac[[facinfix[i]]] <- candidatemarker
					} else {
						pt0 <- freeParTable(pt0, facinfree[i], "=~", free[i], 1:ngroups)
						pt1 <- freeParTable(pt1, facinfree[i], "=~", free[i], 1:ngroups)
					}
				}
			}
			namept1 <- paramNameFromPt(pt1)
			namept0 <- paramNameFromPt(pt0)
			fit0 <- refit(pt0, fit0)
			fit1 <- refit(pt1, fit1)
			beta <- lavaan::coef(fit1)
			beta0 <- lavaan::coef(fit0)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		obsmean <- sapply(lavInspect(fit0, "sampstat"), "[[", "mean") #FIXME: there might not be a mean structure
		obsmean <- obsmean[,grouporder]
		obsdiff <- obsmean[,2:ngroups, drop = FALSE] - matrix(obsmean[,1], nrow(obsmean), ngroups - 1)
		obsdiff <- obsdiff[varfree, , drop = FALSE]
		colnames(obsdiff) <- paste0("diff_mean:", complab)

		estimates <- matrix(NA, length(varfree), ngroups + 1)
		stdestimates <- matrix(NA, length(varfree), ngroups)
		colnames(estimates) <- c("poolest", paste0("load:", grouplab))
		colnames(stdestimates) <- paste0("std:", grouplab)
		esstd <- esz <- matrix(NA, length(varfree), ngroups - 1)
		colnames(esstd) <- paste0("diff_std:", complab)
		colnames(esz) <- paste0("q:", complab)
		esdiff <- matrix(NA, length(varfree), ngroups - 1)

		# Extract facmean, facsd, load, tau -> lowdiff, highdiff
		lowdiff <- matrix(NA, length(varfree), ngroups - 1)
		highdiff <- matrix(NA, length(varfree), ngroups - 1)
		colnames(lowdiff) <- paste0("low_fscore:", complab)
		colnames(highdiff) <- paste0("high_fscore:", complab)

		fixCon <- freeCon <- matrix(NA, length(varfree), 4)
		waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
		colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")
		index <- which((pt1$rhs %in% varfree) & (pt1$op == "=~") & (pt1$group == 1))
		facinfix <- findFactor(fix, facList)
		varinfixvar <- unlist(facList[facinfix])
		varinfixvar <- setdiff(varinfixvar, setdiff(varinfixvar, varfree))
		indexfixvar <- which((pt1$rhs %in% varinfixvar) & (pt1$op == "=~") & (pt1$group == 1))
		varnonfixvar <- setdiff(varfree, varinfixvar)
		indexnonfixvar <- setdiff(index, indexfixvar)

		pos <- 1
		for(i in seq_along(indexfixvar)) {
			runnum <- indexfixvar[i]
			temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
			tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
			if(!is(tryresult, "try-error")) {
				compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
			}
			listFixCon <- c(listFixCon, tryresult)
			temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			estimates[pos, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				loadVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[pos, 2:ncol(estimates)] <- loadVal
				facVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~~", pt0$lhs[runnum], 1:ngroups)
				totalVal <- sapply(lavaan::fitted.values(tempfit0), function(x, v) x$cov[v, v], v = pt0$rhs[runnum])
				names(facVal) <- names(totalVal) <- grouplab
				ifelse(poolvar, refFacVal <- poolVariance(facVal, neach), refFacVal <- facVal[refgroup])
				ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
				stdLoadVal <- loadVal * sqrt(refFacVal) / sqrt(refTotalVal)
				stdestimates[pos,] <- stdLoadVal
				stdLoadVal <- stdLoadVal[grouporder]
				esstd[pos,] <- stdLoadVal[2:ngroups] - stdLoadVal[1]
				if(any(abs(stdLoadVal) > 0.9999)) warning(paste("Standardized Loadings of", pt0$rhs[runnum], "in some groups are less than -1 or over 1. The standardized loadings used in Fisher z transformation are changed to -0.9999 or 0.9999."))
				stdLoadVal[stdLoadVal > 0.9999] <- 0.9999
				stdLoadVal[stdLoadVal < -0.9999] <- -0.9999
				zLoadVal <- atanh(stdLoadVal)
				esz[pos,] <- zLoadVal[2:ngroups] - zLoadVal[1]

				facMean <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~1", "", 1:ngroups)
				wlow <- min(facMean - fbound * sqrt(facVal))
				whigh <- max(facMean + fbound * sqrt(facVal))
				intVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$rhs[runnum], "~1", "", 1:ngroups)
				loadVal <- loadVal[grouporder]
				intVal <- intVal[grouporder]
				loaddiff <- loadVal[2:ngroups] - loadVal[1]
				intdiff <- intVal[2:ngroups] - intVal[1]
				lowdiff[pos,] <- intdiff + wlow * loaddiff
				highdiff[pos,] <- intdiff + whigh * loaddiff
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			waldCon[pos,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
			pos <- pos + 1
		}

		facinvarfree <- findFactor(varnonfixvar, facList)
		for(i in seq_along(indexnonfixvar)) {
			runnum <- indexnonfixvar[i]
			# Need to change marker variable if fixed
			oldmarker <- fixLoadingFac[[facinvarfree[i]]]
			if(length(oldmarker) > 0 && oldmarker == varnonfixvar[i]) {
				candidatemarker <- setdiff(facList[[facinvarfree[i]]], varnonfixvar[i])[1]
				temp <- freeParTable(pt1, facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups)
				temp <- fixParTable(temp, facinvarfree[i], "=~", candidatemarker, 1:ngroups, ustart = 1)
				temp <- constrainParTable(temp, facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups)
				newparent <- freeParTable(pt1, facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups)
				newparent <- fixParTable(newparent, facinvarfree[i], "=~", candidatemarker, 1:ngroups, ustart = 1)
				newparentresult <- try(newparentfit <- refit(newparent, fit1), silent = TRUE)
				if(!is(newparentresult, "try-error")) {
					tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
					if(!is(tryresult, "try-error")) {
						compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, newparentfit, method = method), silent = TRUE)
						if(!is(compresult, "try-error")) fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(newparentfit, tempfit))
					}
					waldCon[pos,] <- waldConstraint(newparentfit, newparent, waldMat, cbind(facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups))
				}
			} else {
				temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
				tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
				if(!is(tryresult, "try-error")) {
					compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
					if(!is(compresult, "try-error"))  fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
				}
				waldCon[pos,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
			}
			listFixCon <- c(listFixCon, tryresult)
			if(length(oldmarker) > 0 && oldmarker == varnonfixvar[i]) {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
			} else {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			}
			estimates[pos, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error")) freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				loadVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[pos, 2:ncol(estimates)] <- loadVal
				facVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~~", pt0$lhs[runnum], 1:ngroups)
				totalVal <- sapply(lavaan::fitted.values(tempfit0), function(x, v) x$cov[v, v], v = pt0$rhs[runnum])
				names(facVal) <- names(totalVal) <- grouplab
				ifelse(poolvar, refFacVal <- poolVariance(facVal, neach), refFacVal <- facVal[refgroup])
				ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
				stdLoadVal <- loadVal * sqrt(refFacVal) / sqrt(refTotalVal)
				stdestimates[pos,] <- stdLoadVal
				stdLoadVal <- stdLoadVal[grouporder]
				esstd[pos,] <- stdLoadVal[2:ngroups] - stdLoadVal[1]
				if(any(abs(stdLoadVal) > 0.9999)) warning(paste("Standardized Loadings of", pt0$rhs[runnum], "in some groups are less than -1 or over 1. The standardized loadings used in Fisher z transformation are changed to -0.9999 or 0.9999."))
				stdLoadVal[stdLoadVal > 0.9999] <- 0.9999
				stdLoadVal[stdLoadVal < -0.9999] <- -0.9999
				zLoadVal <- atanh(stdLoadVal)
				esz[pos,] <- zLoadVal[2:ngroups] - zLoadVal[1]

				facMean <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~1", "", 1:ngroups)
				wlow <- min(facMean - fbound * sqrt(facVal))
				whigh <- max(facMean + fbound * sqrt(facVal))
				intVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$rhs[runnum], "~1", "", 1:ngroups)
				loadVal <- loadVal[grouporder]
				intVal <- intVal[grouporder]
				loaddiff <- loadVal[2:ngroups] - loadVal[1]
				intdiff <- intVal[2:ngroups] - intVal[1]
				lowdiff[pos,] <- intdiff + wlow * loaddiff
				highdiff[pos,] <- intdiff + whigh * loaddiff
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			pos <- pos + 1
		}
		freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
		fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
		waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)
		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[c(indexfixvar, indexnonfixvar)]
		estimates <- cbind(estimates, stdestimates, esstd, esz, obsdiff, lowdiff, highdiff)
		result <- cbind(freeCon, fixCon, waldCon)
	} else if (numType == 2) {
		if(!is.null(free) | !is.null(fix)) {
			if(!is.null(fix)) {
				facinfix <- findFactor(fix, facList)
				dup <- duplicated(facinfix)
				for(i in seq_along(fix)) {
					if(dup[i]) {
						pt0 <- constrainParTable(pt0, fix[i], "~1", "", 1:ngroups)
						pt1 <- constrainParTable(pt1, fix[i], "~1", "", 1:ngroups)
					} else {
						oldmarker <- fixIntceptFac[[facinfix[i]]]
						if(length(oldmarker) > 0) {
							oldmarkerval <- pt1$ustart[pt1$lhs == fix[i] & pt1$op == "~1" & pt1$rhs == "" & pt1$group == 1]
							if(oldmarker == fix[i]) {
								pt0 <- fixParTable(pt0, fix[i], "~1", "", 1:ngroups, oldmarkerval)
								pt1 <- fixParTable(pt1, fix[i], "~1", "", 1:ngroups, oldmarkerval)
							} else {
								pt0 <- freeParTable(pt0, oldmarker, "~1", "", 1:ngroups)
								pt0 <- constrainParTable(pt0, oldmarker, "~1", "", 1:ngroups)
								pt1 <- freeParTable(pt1, oldmarker, "~1", "", 1:ngroups)
								pt0 <- fixParTable(pt0, fix[i], "~1", "", 1:ngroups, oldmarkerval)
								pt1 <- fixParTable(pt1, fix[i], "~1", "", 1:ngroups, oldmarkerval)
								fixIntceptFac[[facinfix[i]]] <- fix[i]
							}
						} else {
							pt0 <- constrainParTable(pt0, fix[i], "~1", "", 1:ngroups)
							pt1 <- constrainParTable(pt1, fix[i], "~1", "", 1:ngroups)
						}
					}
				}
			}
			if(!is.null(free)) {
				facinfree <- findFactor(free, facList)
				for(i in seq_along(free)) {
					# Need to change marker variable if fixed
					oldmarker <- fixIntceptFac[[facinfree[i]]]
					if(length(oldmarker) > 0 && oldmarker == free[i]) {
						oldmarkerval <- pt1$ustart[pt1$lhs == oldmarker & pt1$op == "~1" & pt1$rhs == "" & pt1$group == 1]
						candidatemarker <- setdiff(facList[[facinfree[i]]], free[i])[1]
						pt0 <- freeParTable(pt0, free[i], "~1", "", 1:ngroups)
						pt1 <- freeParTable(pt1, free[i], "~1", "", 1:ngroups)
						pt0 <- fixParTable(pt0, candidatemarker, "~1", "", 1:ngroups, oldmarkerval)
						pt1 <- fixParTable(pt1, candidatemarker, "~1", "", 1:ngroups, oldmarkerval)
						fixIntceptFac[[facinfix[i]]] <- candidatemarker
					} else {
						pt0 <- freeParTable(pt0, free[i], "~1", "", 1:ngroups)
						pt1 <- freeParTable(pt1, free[i], "~1", "", 1:ngroups)
					}
				}
			}
			namept1 <- paramNameFromPt(pt1)
			namept0 <- paramNameFromPt(pt0)
			fit0 <- refit(pt0, fit0)
			fit1 <- refit(pt1, fit1)
			beta <- lavaan::coef(fit1)
			beta0 <- lavaan::coef(fit0)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		obsmean <- sapply(lavInspect(fit0, "sampstat"), "[[", "mean") #FIXME: there might not be a mean structure
		obsmean <- obsmean[,grouporder]
		obsdiff <- obsmean[,2:ngroups, drop = FALSE] - matrix(obsmean[,1], nrow(obsmean), ngroups - 1)
		obsdiff <- obsdiff[varfree, , drop = FALSE]
		colnames(obsdiff) <- paste0("diff_mean:", complab)

		# Prop diff
		propdiff <- matrix(NA, length(varfree), ngroups - 1)
		colnames(propdiff) <- paste0("propdiff:", complab)

		estimates <- matrix(NA, length(varfree), ngroups + 1)
		stdestimates <- matrix(NA, length(varfree), ngroups)
		colnames(estimates) <- c("poolest", paste0("int:", grouplab))
		colnames(stdestimates) <- paste0("std:", grouplab)
		esstd <- matrix(NA, length(varfree), ngroups - 1)
		colnames(esstd) <- paste0("diff_std:", complab)
		fixCon <- freeCon <- matrix(NA, length(varfree), 4)
		waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
		colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")
		index <- which((pt1$lhs %in% varfree) & (pt1$op == "~1") & (pt1$group == 1))
		facinfix <- findFactor(fix, facList)
		varinfixvar <- unlist(facList[facinfix])
		varinfixvar <- setdiff(varinfixvar, setdiff(varinfixvar, varfree))
		indexfixvar <- which((pt1$lhs %in% varinfixvar) & (pt1$op == "~1") & (pt1$group == 1))
		varnonfixvar <- setdiff(varfree, varinfixvar)
		indexnonfixvar <- setdiff(index, indexfixvar)

		pos <- 1
		for(i in seq_along(varinfixvar)) {
			runnum <- indexfixvar[i]
			temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
			tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
			if(!is(tryresult, "try-error")) {
				compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
			}
			listFixCon <- c(listFixCon, tryresult)
			temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			estimates[pos, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				intVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[pos, 2:ncol(estimates)] <- intVal
				totalVal <- sapply(lavaan::fitted.values(tempfit0), function(x, v) x$cov[v, v], v = pt0$lhs[runnum])
				ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
				stdIntVal <- intVal / sqrt(refTotalVal)
				stdestimates[pos,] <- stdIntVal
				stdIntVal <- stdIntVal[grouporder]
				esstd[pos,] <- stdIntVal[2:ngroups] - stdIntVal[1]

				intVal <- intVal[grouporder]
				propdiff[pos,] <- (intVal[2:ngroups] - intVal[1]) / obsdiff[pos,]
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			waldCon[pos,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
			pos <- pos + 1
		}

		facinvarfree <- findFactor(varfree, facList)
		for(i in seq_along(varnonfixvar)) {
			runnum <- indexnonfixvar[i]
			# Need to change marker variable if fixed
			oldmarker <- fixIntceptFac[[facinvarfree[i]]]
			if(length(oldmarker) > 0 && oldmarker == varfree[i]) {
				candidatemarker <- setdiff(facList[[facinvarfree[i]]], varfree[i])[1]
				temp <- freeParTable(pt1, varfree[i], "~1", "", 1:ngroups)
				temp <- constrainParTable(temp, varfree[i], "~1", "", 1:ngroups)
				temp <- fixParTable(temp, candidatemarker, "~1", "", 1:ngroups)
				newparent <- freeParTable(pt1, varfree[i], "~1", "", 1:ngroups)
				newparent <- fixParTable(newparent, candidatemarker, "~1", "", 1:ngroups)
				newparentresult <- try(newparentfit <- refit(newparent, fit1), silent = TRUE)
				if(!is(newparentresult, "try-error")) {
					tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
					if(!is(tryresult, "try-error")) {
						compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, newparentfit, method = method), silent = TRUE)
						if(!is(compresult, "try-error")) fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(newparentfit, tempfit))
					}
					waldCon[pos,] <- waldConstraint(newparentfit, newparent, waldMat, cbind(varfree[i], "~1", "", 1:ngroups))
				}
			} else {
				temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
				tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
				if(!is(tryresult, "try-error")) {
					compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
					if(!is(compresult, "try-error"))  fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
				}
				waldCon[pos,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
			}
			listFixCon <- c(listFixCon, tryresult)
			if(length(oldmarker) > 0 && oldmarker == varfree[i]) {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
			} else {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			}
			estimates[pos, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				intVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[pos, 2:ncol(estimates)] <- intVal
				totalVal <- sapply(lavaan::fitted.values(tempfit0), function(x, v) x$cov[v, v], v = pt0$lhs[runnum])
				ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
				stdIntVal <- intVal / sqrt(refTotalVal)
				stdestimates[pos,] <- stdIntVal
				stdIntVal <- stdIntVal[grouporder]
				esstd[pos,] <- stdIntVal[2:ngroups] - stdIntVal[1]

				intVal <- intVal[grouporder]
				propdiff[pos,] <- (intVal[2:ngroups] - intVal[1]) / obsdiff[pos,]
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			pos <- pos + 1
		}
		freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
		fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
		waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)

		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[c(indexfixvar, indexnonfixvar)]
		estimates <- cbind(estimates, stdestimates, esstd, obsdiff, propdiff)
		result <- cbind(freeCon, fixCon, waldCon)
	} else if (numType == 3) {
		if(!is.null(free) | !is.null(fix)) {
			if(!is.null(fix)) {
				for(i in seq_along(fix)) {
					pt0 <- constrainParTable(pt0, fix[i], "~~", fix[i], 1:ngroups)
					pt1 <- constrainParTable(pt1, fix[i], "~~", fix[i], 1:ngroups)
				}
			}
			if(!is.null(free)) {
				for(i in seq_along(free)) {
					pt0 <- freeParTable(pt0, free[i], "~~", free[i], 1:ngroups)
					pt1 <- freeParTable(pt1, free[i], "~~", free[i], 1:ngroups)
				}
			}
			namept1 <- paramNameFromPt(pt1)
			namept0 <- paramNameFromPt(pt0)
			fit0 <- refit(pt0, fit0)
			fit1 <- refit(pt1, fit1)
			beta <- lavaan::coef(fit1)
			beta0 <- lavaan::coef(fit0)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		# Prop diff
		propdiff <- matrix(NA, length(varfree), ngroups - 1)
		colnames(propdiff) <- paste0("propdiff:", complab)

		estimates <- matrix(NA, length(varfree), ngroups + 1)
		stdestimates <- matrix(NA, length(varfree), ngroups)
		colnames(estimates) <- c("poolest", paste0("errvar:", grouplab))
		colnames(stdestimates) <- paste0("std:", grouplab)
		esstd <- esz <- matrix(NA, length(varfree), ngroups - 1)
		colnames(esstd) <- paste0("diff_std:", complab)
		colnames(esz) <- paste0("h:", complab)
		fixCon <- freeCon <- matrix(NA, length(varfree), 4)
		waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
		colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")
		index <- which((pt1$lhs %in% varfree) & (pt1$op == "~~") & (pt1$lhs == pt1$rhs) & (pt1$group == 1))
		for(i in seq_along(index)) {
			runnum <- index[i]
			temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
			tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
			if(!is(tryresult, "try-error")) {
				compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[i,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
			}
			listFixCon <- c(listFixCon, tryresult)
			temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			estimates[i, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[i,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				errVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[i, 2:ncol(estimates)] <- errVal
				totalVal <- sapply(lavaan::fitted.values(tempfit0), function(x, v) x$cov[v, v], v = pt0$rhs[runnum])
				ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
				stdErrVal <- errVal / sqrt(refTotalVal)
				stdestimates[i,] <- stdErrVal
				stdErrVal <- stdErrVal[grouporder]
				esstd[i,] <- stdErrVal[2:ngroups] - stdErrVal[1]
				if(any(abs(stdErrVal) > 0.9999)) warning(paste("The uniqueness of", pt0$rhs[runnum], "in some groups are over 1. The uniqueness used in arctan transformation are changed to 0.9999."))
				stdErrVal[stdErrVal > 0.9999] <- 0.9999
				zErrVal <- asin(sqrt(stdErrVal))
				esz[i,] <- zErrVal[2:ngroups] - zErrVal[1]

				errVal <- errVal[grouporder]
				totalVal <- totalVal[grouporder]
				errdiff <- errVal[2:ngroups] - errVal[1]
				totaldiff <- totalVal[2:ngroups] - totalVal[1]
				propdiff[i,] <- errdiff / totaldiff
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			waldCon[i,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
		}

		freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
		fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
		waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)
		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[index]
		estimates <- cbind(estimates, stdestimates, esstd, esz, propdiff)
		result <- cbind(freeCon, fixCon, waldCon)
	} else if (numType == 4) {
		varfree <- facnames
		if(!is.null(free) | !is.null(fix)) {
			if(!is.null(fix)) {
				for(i in seq_along(fix)) {
					pt0 <- constrainParTable(pt0, fix[i], "~1", "", 1:ngroups)
					pt1 <- constrainParTable(pt1, fix[i], "~1", "", 1:ngroups)
				}
			}
			if(!is.null(free)) {
				for(i in seq_along(free)) {
					pt0 <- freeParTable(pt0, free[i], "~1", "", 1:ngroups)
					pt1 <- freeParTable(pt1, free[i], "~1", "", 1:ngroups)
				}
			}
			namept1 <- paramNameFromPt(pt1)
			namept0 <- paramNameFromPt(pt0)
			fit0 <- refit(pt0, fit0)
			fit1 <- refit(pt1, fit1)
			beta <- lavaan::coef(fit1)
			beta0 <- lavaan::coef(fit0)
			waldMat <- matrix(0, ngroups - 1, length(beta))
			varfree <- setdiff(varfree, c(free, fix))
		}

		estimates <- matrix(NA, length(varfree), ngroups + 1)
		stdestimates <- matrix(NA, length(varfree), ngroups)
		colnames(estimates) <- c("poolest", paste0("mean:", grouplab))
		colnames(stdestimates) <- paste0("std:", grouplab)
		esstd <- matrix(NA, length(varfree), ngroups - 1)
		colnames(esstd) <- paste0("diff_std:", complab)
		fixCon <- freeCon <- matrix(NA, length(varfree), 4)
		waldCon <- matrix(NA, length(varfree), 3)
		colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
		colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
		colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")
		index <- which((pt1$lhs %in% varfree) & (pt1$op == "~1") & (pt1$group == 1))
		for(i in seq_along(index)) {
			runnum <- index[i]
			isfree <- pt1$free[runnum] != 0
			if(isfree) {
				temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
			} else {
				temp <- fixParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 2:ngroups, ustart = pt1$ustart[runnum])
			}
			tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
			if(!is(tryresult, "try-error")) {
				compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
				if(!is(compresult, "try-error"))  fixCon[i,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
			}
			listFixCon <- c(listFixCon, tryresult)
			isfree0 <- pt0$free[runnum] != 0
			if(isfree0) {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
			} else {
				temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
			}
			estimates[i, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
			tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
			if(!is(tryresult0, "try-error")) {
				compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
				if(!is(compresult0, "try-error"))  freeCon[i,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
				meanVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
				estimates[i, 2:ncol(estimates)] <- meanVal
				facVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~~", pt0$lhs[runnum], 1:ngroups)
				ifelse(poolvar, refFacVal <- poolVariance(facVal, neach), refFacVal <- facVal[refgroup])
				stdMeanVal <- meanVal / sqrt(refFacVal)
				stdestimates[i,] <- stdMeanVal
				stdMeanVal <- stdMeanVal[grouporder]
				esstd[i,] <- stdMeanVal[2:ngroups] - stdMeanVal[1]
			}
			listFreeCon <- c(listFreeCon, tryresult0)
			waldCon[i,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
		}
		freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
		fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
		waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)

		rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[index]
		estimates <- cbind(estimates, stdestimates, esstd)
		result <- cbind(freeCon, fixCon, waldCon)
	}
	if(return.fit) {
		return(invisible(list(estimates = estimates, results = result, models = list(free = listFreeCon, fix = listFixCon, nested = fit0, parent = fit1))))
	} else {
		return(list(estimates = estimates, results = result))
	}
}


##' @importFrom lavaan lavInspect parTable
##' @rdname partialInvariance
##' @export
partialInvarianceCat <- function(fit, type, free = NULL, fix = NULL,
                                 refgroup = 1, poolvar = TRUE,
                                 p.adjust = "none", return.fit = FALSE,
                                 method = "satorra.bentler.2001") {

  type <- tolower(type)
  numType <- 1
  fit1 <- fit0 <- NULL
  # fit0 = Nested model, fit1 = Parent model
  if (type %in% c("metric", "weak", "loading", "loadings")) {
    numType <- 1
    if (all(c("fit.configural", "fit.loadings") %in% names(fit))) {
      fit1 <- fit$fit.configural
      fit0 <- fit$fit.loadings
    } else {
      stop("The elements named 'fit.configural' and 'fit.loadings' are needed",
           " in the 'fit' argument")
    }
  } else if (type %in% c("scalar", "strong", "intercept", "intercepts",
                         "threshold", "thresholds")) {
    numType <- 2
    if (all(c("fit.loadings", "fit.thresholds") %in% names(fit))) {
      fit1 <- fit$fit.loadings
      fit0 <- fit$fit.thresholds
    } else {
      stop("The elements named 'fit.loadings' and 'fit.thresholds' are needed",
           " in the 'fit' argument")
    }
  } else if (type %in% c("strict", "residual", "residuals", "error", "errors")) {
    numType <- 3
    if ("fit.residuals" %in% names(fit)) {
      fit0 <- fit$fit.residuals
      if ("fit.thresholds" %in% names(fit)) {
        fit1 <- fit$fit.thresholds
      } else if ("fit.loadings" %in% names(fit)) {
        fit1 <- fit$fit.loadings
      } else {
        stop("The element named either 'fit.thresholds' or 'fit.loadings' is",
             " needed in the 'fit' argument")
      }
    } else {
      stop("The element named 'fit.residuals' is needed in the 'fit' argument")
    }
  } else if (type %in% c("means", "mean")) {
    numType <- 4
    if ("fit.means" %in% names(fit)) {
      fit0 <- fit$fit.means
      if("fit.residuals" %in% names(fit)) {
        fit1 <- fit$fit.residuals
      } else if ("fit.thresholds" %in% names(fit)) {
        fit1 <- fit$fit.thresholds
      } else if ("fit.loadings" %in% names(fit)) {
        fit1 <- fit$fit.loadings
      } else {
        stop("The element named either 'fit.residuals', 'fit.thresholds',",
             " or 'fit.loadings' is needed in the 'fit' argument")
      }
    } else {
      stop("The element named 'fit.means' is needed in the 'fit' argument")
    }
  } else {
    stop("Please specify the correct type of measurement invariance. See the help page.")
  }
  pt1 <- parTable(fit1)
  pt0 <- parTable(fit0)
  pt0$start <- pt0$est <- pt0$se <- NULL
  pt1$start <- pt1$est <- pt1$se <- NULL

  pt1$label[substr(pt1$label, 1, 1) == "." & substr(pt1$label, nchar(pt1$label),
                                                    nchar(pt1$label)) == "."] <- ""
  pt0$label[substr(pt0$label, 1, 1) == "." & substr(pt0$label, nchar(pt0$label),
                                                    nchar(pt0$label)) == "."] <- ""
  namept1 <- paramNameFromPt(pt1)
  namept0 <- paramNameFromPt(pt0)
  if (length(table(table(pt0$rhs[pt0$op == "=~"]))) != 1)
    stop("The model is not congeneric. This function does not support non-congeneric model.")
  varfree <- varnames <- unique(pt0$rhs[pt0$op == "=~"])
  facnames <- unique(pt0$lhs[(pt0$op == "=~") & (pt0$rhs %in% varnames)])
  facrepresent <- table(pt0$lhs[(pt0$op == "=~") & (pt0$rhs %in% varnames)],
                        pt0$rhs[(pt0$op == "=~") & (pt0$rhs %in% varnames)])
  if (any(apply(facrepresent, 2, function(x) sum(x != 0)) > 1))
    stop("The model is not congeneric. This function does not support non-congeneric model.")
  facList <- list()
  for (i in 1:nrow(facrepresent)) {
    facList[[i]] <- colnames(facrepresent)[facrepresent[i,] > 0]
  }
  names(facList) <- rownames(facrepresent)
  facList <- facList[match(names(facList), facnames)]
  fixLoadingFac <- list()
  for (i in seq_along(facList)) {
    select <- pt1$lhs == names(facList)[i] & pt1$op == "=~" & pt1$rhs %in% facList[[i]] & pt1$group == 1 & pt1$free == 0 & (!is.na(pt1$ustart) & pt1$ustart > 0)
    fixLoadingFac[[i]] <- pt1$rhs[select]
  }
  names(fixLoadingFac) <- names(facList)

  # Find the number of thresholds
  # Check whether the factor configuration is the same across gorups

  conParTable <- lapply(pt1, "[", pt1$op == "==")
  group1pt <- lapply(pt1, "[", pt1$group != 1)

  numThreshold <- table(sapply(group1pt, "[", group1pt$op == "|")[,"lhs"])
  plabelthres <- split(group1pt$plabel[group1pt$op == "|"], group1pt$lhs[group1pt$op == "|"])
  numFixedThreshold <- sapply(lapply(plabelthres, function(vec) !is.na(match(vec, conParTable$lhs)) | !is.na(match(vec, conParTable$rhs))), sum)[names(numThreshold)]

  #numFixedThreshold <- table(sapply(group1pt, "[", group1pt$op == "|" & group1pt$eq.id != 0)[,"lhs"])
  fixIntceptFac <- list()
  for (i in seq_along(facList)) {
    tmp <- numFixedThreshold[facList[[i]]]
    if (all(tmp > 1)) {
      fixIntceptFac[[i]] <- integer(0)
    } else {
      fixIntceptFac[[i]] <- names(which.max(tmp))[1]
    }
  }
  names(fixIntceptFac) <- names(facList)

  ngroups <- max(pt0$group)
  neach <- lavInspect(fit0, "nobs")
  groupvar <- lavInspect(fit0, "group")
  grouplab <- lavInspect(fit0, "group.label")
  if (!is.numeric(refgroup)) refgroup <- which(refgroup == grouplab)
  grouporder <- 1:ngroups
  grouporder <- c(refgroup, setdiff(grouporder, refgroup))
  grouplaborder <- grouplab[grouporder]
  complab <- paste(grouplaborder[2:ngroups], "vs.", grouplaborder[1])
  if (ngroups <= 1) stop("Well, the number of groups is 1. Measurement",
                         " invariance across 'groups' cannot be done.")

  if (numType == 4) {
    if (!all(c(free, fix) %in% facnames))
      stop("'free' and 'fix' arguments should consist of factor names because",
           " mean invariance is tested.")
  } else {
    if (!all(c(free, fix) %in% varnames))
      stop("'free' and 'fix' arguments should consist of variable names.")
  }
  result <- fixCon <- freeCon <- NULL
  estimates <- NULL
  listFreeCon <- listFixCon <- list()
  beta <- lavaan::coef(fit1)
  beta0 <- lavaan::coef(fit0)
  waldMat <- matrix(0, ngroups - 1, length(beta))
  if (numType == 1) {
    if (!is.null(free) | !is.null(fix)) {
      if (!is.null(fix)) {
        facinfix <- findFactor(fix, facList)
        dup <- duplicated(facinfix)
        for (i in seq_along(fix)) {
          if (dup[i]) {
            pt0 <- constrainParTable(pt0, facinfix[i], "=~", fix[i], 1:ngroups)
            pt1 <- constrainParTable(pt1, facinfix[i], "=~", fix[i], 1:ngroups)
          } else {
            oldmarker <- fixLoadingFac[[facinfix[i]]]
            if (length(oldmarker) > 0) {
              oldmarkerval <- pt1$ustart[pt1$lhs == facinfix[i] & pt1$op == "=~" & pt1$rhs == oldmarker & pt1$group == 1]
              if (oldmarker == fix[i]) {
                pt0 <- fixParTable(pt0, facinfix[i], "=~", fix[i], 1:ngroups, oldmarkerval)
                pt1 <- fixParTable(pt1, facinfix[i], "=~", fix[i], 1:ngroups, oldmarkerval)
              } else {
                pt0 <- freeParTable(pt0, facinfix[i], "=~", oldmarker, 1:ngroups)
                pt0 <- constrainParTable(pt0, facinfix[i], "=~", oldmarker, 1:ngroups)
                pt1 <- freeParTable(pt1, facinfix[i], "=~", oldmarker, 1:ngroups)
                pt0 <- fixParTable(pt0, facinfix[i], "=~", fix[i], 1:ngroups, oldmarkerval)
                pt1 <- fixParTable(pt1, facinfix[i], "=~", fix[i], 1:ngroups, oldmarkerval)
                fixLoadingFac[[facinfix[i]]] <- fix[i]
              }
            } else {
              pt0 <- constrainParTable(pt0, facinfix[i], "=~", fix[i], 1:ngroups)
              pt1 <- constrainParTable(pt1, facinfix[i], "=~", fix[i], 1:ngroups)
            }
          }
        }
      }
      if (!is.null(free)) {
        facinfree <- findFactor(free, facList)
        for (i in seq_along(free)) {
          # Need to change marker variable if fixed
          oldmarker <- fixLoadingFac[[facinfree[i]]]
          if (length(oldmarker) > 0 && oldmarker == free[i]) {
            oldmarkerval <- pt1$ustart[pt1$lhs == facinfix[i] & pt1$op == "=~" & pt1$rhs == oldmarker & pt1$group == 1]
            candidatemarker <- setdiff(facList[[facinfree[i]]], free[i])[1]
            pt0 <- freeParTable(pt0, facinfree[i], "=~", free[i], 1:ngroups)
            pt1 <- freeParTable(pt1, facinfree[i], "=~", free[i], 1:ngroups)
            pt0 <- fixParTable(pt0, facinfix[i], "=~", candidatemarker, 1:ngroups, oldmarkerval)
            pt1 <- fixParTable(pt1, facinfix[i], "=~", candidatemarker, 1:ngroups, oldmarkerval)
            fixLoadingFac[[facinfix[i]]] <- candidatemarker
          } else {
            pt0 <- freeParTable(pt0, facinfree[i], "=~", free[i], 1:ngroups)
            pt1 <- freeParTable(pt1, facinfree[i], "=~", free[i], 1:ngroups)
          }
        }
      }
      namept1 <- paramNameFromPt(pt1)
      namept0 <- paramNameFromPt(pt0)
      fit0 <- refit(pt0, fit0)
      fit1 <- refit(pt1, fit1)
      beta <- lavaan::coef(fit1)
      beta0 <- lavaan::coef(fit0)
      waldMat <- matrix(0, ngroups - 1, length(beta))
      varfree <- setdiff(varfree, c(free, fix))
    }

    estimates <- matrix(NA, length(varfree), ngroups + 1)
    stdestimates <- matrix(NA, length(varfree), ngroups)
    colnames(estimates) <- c("poolest", paste0("load:", grouplab))
    colnames(stdestimates) <- paste0("std:", grouplab)
    esstd <- esz <- matrix(NA, length(varfree), ngroups - 1)
    colnames(esstd) <- paste0("diff_std:", complab)
    colnames(esz) <- paste0("q:", complab)
    fixCon <- freeCon <- matrix(NA, length(varfree), 4)
    waldCon <- matrix(NA, length(varfree), 3)
    colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
    colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
    colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")
    index <- which((pt1$rhs %in% varfree) & (pt1$op == "=~") & (pt1$group == 1))
    facinfix <- findFactor(fix, facList)
    varinfixvar <- unlist(facList[facinfix])
    varinfixvar <- setdiff(varinfixvar, setdiff(varinfixvar, varfree))
    indexfixvar <- which((pt1$rhs %in% varinfixvar) & (pt1$op == "=~") & (pt1$group == 1))
    varnonfixvar <- setdiff(varfree, varinfixvar)
    indexnonfixvar <- setdiff(index, indexfixvar)

    pos <- 1
    for (i in seq_along(indexfixvar)) {
      runnum <- indexfixvar[i]
      temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
      tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
      if (!is(tryresult, "try-error")) {
        compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
        if (!is(compresult, "try-error"))  fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
      }
      listFixCon <- c(listFixCon, tryresult)
      temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
      estimates[pos, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
      tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
      if (!is(tryresult0, "try-error")) {
        compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
        if (!is(compresult0, "try-error"))  freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
        loadVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
        estimates[pos, 2:ncol(estimates)] <- loadVal
        facVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~~", pt0$lhs[runnum], 1:ngroups)
        totalVal <- sapply(thetaImpliedTotalVar(tempfit0), function(x, v) x[v, v], v = pt0$rhs[runnum])
        names(facVal) <- names(totalVal) <- grouplab
        ifelse(poolvar, refFacVal <- poolVariance(facVal, neach), refFacVal <- facVal[refgroup])
        ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
        stdLoadVal <- loadVal * sqrt(refFacVal) / sqrt(refTotalVal)
        stdestimates[pos,] <- stdLoadVal
        stdLoadVal <- stdLoadVal[grouporder]
        esstd[pos,] <- stdLoadVal[2:ngroups] - stdLoadVal[1]
        if (any(abs(stdLoadVal) > 0.9999))
          warning(paste("Standardized Loadings of", pt0$rhs[runnum],
                        "in some groups are less than -1 or over 1. The",
                        " standardized loadings used in Fisher z",
                        " transformation are changed to -0.9999 or 0.9999."))
        stdLoadVal[stdLoadVal > 0.9999] <- 0.9999
        stdLoadVal[stdLoadVal < -0.9999] <- -0.9999
        zLoadVal <- atanh(stdLoadVal)
        esz[pos,] <- zLoadVal[2:ngroups] - zLoadVal[1]
      }
      listFreeCon <- c(listFreeCon, tryresult0)
      waldCon[pos,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
      pos <- pos + 1
    }

    facinvarfree <- findFactor(varnonfixvar, facList)
    for (i in seq_along(indexnonfixvar)) {
      runnum <- indexnonfixvar[i]
      # Need to change marker variable if fixed
      oldmarker <- fixLoadingFac[[facinvarfree[i]]]
      if (length(oldmarker) > 0 && oldmarker == varnonfixvar[i]) {
        candidatemarker <- setdiff(facList[[facinvarfree[i]]], varnonfixvar[i])[1]
        temp <- freeParTable(pt1, facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups)
        temp <- constrainParTable(temp, facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups)
        temp <- fixParTable(temp, facinvarfree[i], "=~", candidatemarker, 1:ngroups)
        newparent <- freeParTable(pt1, facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups)
        newparent <- fixParTable(newparent, facinvarfree[i], "=~", candidatemarker, 1:ngroups)
        newparentresult <- try(newparentfit <- refit(newparent, fit1), silent = TRUE)
        if (!is(newparentresult, "try-error")) {
          tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
          if (!is(tryresult, "try-error")) {
            compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, newparentfit, method = method), silent = TRUE)
            if (!is(compresult, "try-error")) fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(newparentfit, tempfit))
          }
          waldCon[pos,] <- waldConstraint(newparentfit, newparent, waldMat, cbind(facinvarfree[i], "=~", varnonfixvar[i], 1:ngroups))
        }
      } else {
        temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
        tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
        if (!is(tryresult, "try-error")) {
          compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
          if (!is(compresult, "try-error"))  fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
        }
        waldCon[pos,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
      }
      listFixCon <- c(listFixCon, tryresult)
      if (length(oldmarker) > 0 && oldmarker == varnonfixvar[i]) {
        temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
      } else {
        temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
      }
      estimates[pos, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
      tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
      if (!is(tryresult0, "try-error")) {
        compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
        if (!is(compresult0, "try-error")) freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
        loadVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
        estimates[pos, 2:ncol(estimates)] <- loadVal
        facVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~~", pt0$lhs[runnum], 1:ngroups)
        totalVal <- sapply(thetaImpliedTotalVar(tempfit0), function(x, v) x[v, v], v = pt0$rhs[runnum])
        names(facVal) <- names(totalVal) <- grouplab
        ifelse(poolvar, refFacVal <- poolVariance(facVal, neach), refFacVal <- facVal[refgroup])
        ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
        stdLoadVal <- loadVal * sqrt(refFacVal) / sqrt(refTotalVal)
        stdestimates[pos,] <- stdLoadVal
        stdLoadVal <- stdLoadVal[grouporder]
        esstd[pos,] <- stdLoadVal[2:ngroups] - stdLoadVal[1]
        if (any(abs(stdLoadVal) > 0.9999))
          warning(paste("Standardized Loadings of", pt0$rhs[runnum],
                        "in some groups are less than -1 or over 1. The",
                        " standardized loadings used in Fisher z",
                        " transformation are changed to -0.9999 or 0.9999."))
        stdLoadVal[stdLoadVal > 0.9999] <- 0.9999
        stdLoadVal[stdLoadVal < -0.9999] <- -0.9999
        zLoadVal <- atanh(stdLoadVal)
        esz[pos,] <- zLoadVal[2:ngroups] - zLoadVal[1]
      }
      listFreeCon <- c(listFreeCon, tryresult0)
      pos <- pos + 1
    }
    freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
    fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
    waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)

    rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[c(indexfixvar, indexnonfixvar)]
    estimates <- cbind(estimates, stdestimates, esstd, esz)
    result <- cbind(freeCon, fixCon, waldCon)
  } else if (numType == 2) {
    if (!is.null(free) | !is.null(fix)) {
      if (!is.null(fix)) {
        facinfix <- findFactor(fix, facList)
        dup <- duplicated(facinfix)
        for (i in seq_along(fix)) {
          numfixthres <- numThreshold[fix[i]]
          if (numfixthres > 1) {
            if (dup[i]) {
              for (s in 2:numfixthres) {
                pt0 <- constrainParTable(pt0, fix[i], "|", paste0("t", s), 1:ngroups)
                pt1 <- constrainParTable(pt1, fix[i], "|", paste0("t", s), 1:ngroups)
              }
            } else {
              oldmarker <- fixIntceptFac[[facinfix[i]]]
              numoldthres <- numThreshold[oldmarker]
              if (length(oldmarker) > 0) {
                if (oldmarker == fix[i]) {
                  for (s in 2:numfixthres) {
                    pt0 <- constrainParTable(pt0, fix[i], "|", paste0("t", s), 1:ngroups)
                    pt1 <- constrainParTable(pt1, fix[i], "|", paste0("t", s), 1:ngroups)
                  }
                } else {
                  for (r in 2:numoldthres) {
                    pt1 <- freeParTable(pt1, oldmarker, "|", paste0("t", r), 1:ngroups)
                  }
                  for (s in 2:numfixthres) {
                    pt0 <- constrainParTable(pt0, fix[i], "|", paste0("t", s), 1:ngroups)
                    pt1 <- constrainParTable(pt1, fix[i], "|", paste0("t", s), 1:ngroups)
                  }
                  fixIntceptFac[[facinfix[i]]] <- fix[i]
                }
              } else {
                for (s in 2:numfixthres) {
                  pt0 <- constrainParTable(pt0, fix[i], "|", paste0("t", s), 1:ngroups)
                  pt1 <- constrainParTable(pt1, fix[i], "|", paste0("t", s), 1:ngroups)
                }
              }
            }
          }
        }
      }
      if (!is.null(free)) {
        facinfree <- findFactor(free, facList)
        for (i in seq_along(free)) {
          numfreethres <- numThreshold[free[i]]
          # Need to change marker variable if fixed
          oldmarker <- fixIntceptFac[[facinfree[i]]]
          numoldthres <- numThreshold[oldmarker]
          if (length(oldmarker) > 0 && oldmarker == free[i]) {
            candidatemarker <- setdiff(facList[[facinfree[i]]], free[i])
            candidatemarker <- candidatemarker[numThreshold[candidatemarker] > 1][1]
            numcandidatethres <- numThreshold[candidatemarker]
            pt0 <- constrainParTable(pt0, candidatemarker, "|", "t2", 1:ngroups)
            pt1 <- constrainParTable(pt1, candidatemarker, "|", "t2", 1:ngroups)
            for (s in 2:numfixthres) {
              pt0 <- freeParTable(pt0, free[i], "|", paste0("t", s), 1:ngroups)
              pt1 <- freeParTable(pt1, free[i], "|", paste0("t", s), 1:ngroups)
            }
            fixIntceptFac[[facinfix[i]]] <- candidatemarker
          } else {
            for (s in 2:numfixthres) {
              pt0 <- freeParTable(pt0, free[i], "|", paste0("t", s), 1:ngroups)
              pt1 <- freeParTable(pt1, free[i], "|", paste0("t", s), 1:ngroups)
            }
          }
        }
      }
      namept1 <- paramNameFromPt(pt1)
      namept0 <- paramNameFromPt(pt0)
      fit0 <- refit(pt0, fit0)
      fit1 <- refit(pt1, fit1)
      beta <- lavaan::coef(fit1)
      beta0 <- lavaan::coef(fit0)
      waldMat <- matrix(0, ngroups - 1, length(beta))
      varfree <- setdiff(varfree, c(free, fix))
    }

    maxcolumns <- max(numThreshold[varfree]) - 1
    tname <- paste0("t", 2:(maxcolumns + 1))
    estimates <- matrix(NA, length(varfree), (ngroups * length(tname)) + length(tname))
    stdestimates <- matrix(NA, length(varfree), ngroups * length(tname))
    tnameandlab <- expand.grid(tname, grouplab)
    colnames(estimates) <- c(paste0("pool:", tname), paste0(tnameandlab[,1], ":", tnameandlab[,2]))
    colnames(stdestimates) <- paste0("std:", tnameandlab[,1], ":", tnameandlab[,2])
    esstd <- matrix(NA, length(varfree), (ngroups - 1)* length(tname))
    tnameandcomplab <- expand.grid(tname, complab)
    colnames(esstd) <- paste0("diff_std:", tnameandcomplab[,1], ":", tnameandcomplab[,2])
    fixCon <- freeCon <- matrix(NA, length(varfree), 4)
    waldCon <- matrix(NA, length(varfree), 3)
    colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
    colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
    colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")

    facinfix <- findFactor(fix, facList)
    varinfixvar <- unlist(facList[facinfix])
    varinfixvar <- setdiff(varinfixvar, setdiff(varinfixvar, varfree))
    varnonfixvar <- setdiff(varfree, varinfixvar)

    pos <- 1
    for (i in seq_along(varinfixvar)) {
      temp <- pt1
      for (s in 2:numThreshold[varinfixvar[i]]) {
        runnum <- which((pt1$lhs == varfree[i]) & (pt1$op == "|") & (pt1$rhs == paste0("t", s)) & (pt1$group == 1))
        temp <- constrainParTable(temp, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
      }
      tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
      if (!is(tryresult, "try-error")) {
        compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
        if (!is(compresult, "try-error"))  fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
      }
      listFixCon <- c(listFixCon, tryresult)
      temp0 <- pt0
      for (s in 2:numThreshold[varinfixvar[i]]) {
        runnum <- which((pt0$lhs == varfree[i]) & (pt0$op == "|") & (pt0$rhs == paste0("t", s)) & (pt0$group == 1))
        temp0 <- freeParTable(temp0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
        estimates[pos, s - 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
      }
      tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
      if (!is(tryresult0, "try-error")) {
        compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
        if (!is(compresult0, "try-error"))  freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
        for (s in 2:numThreshold[varinfixvar[i]]) {
          runnum <- which((pt0$lhs == varfree[i]) & (pt0$op == "|") & (pt0$rhs == paste0("t", s)) & (pt0$group == 1))
          thresVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
          estimates[pos, maxcolumns*(1:ngroups) + (s - 1)] <- thresVal
          totalVal <- sapply(thetaImpliedTotalVar(tempfit0), function(x, v) x[v, v], v = pt0$lhs[runnum])
          ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
          stdIntVal <- thresVal / sqrt(refTotalVal)
          stdestimates[pos, maxcolumns*(1:ngroups - 1) + (s - 1)] <- stdIntVal
          stdIntVal <- stdIntVal[grouporder]
          esstd[pos, maxcolumns*(1:length(complab) - 1) + (s - 1)] <- stdIntVal[2:ngroups] - stdIntVal[1]
        }
      }
      listFreeCon <- c(listFreeCon, tryresult0)
      args <- list(fit1, pt1, waldMat)
      for (s in 2:numThreshold[varinfixvar[i]]) {
        runnum <- which((pt1$lhs == varfree[i]) & (pt1$op == "|") & (pt1$rhs == paste0("t", s)) & (pt1$group == 1))
        args <- c(args, list(cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)))
      }
      waldCon[pos,] <- do.call(waldConstraint, args)
      pos <- pos + 1
    }

    facinvarfree <- findFactor(varnonfixvar, facList)
    for (i in seq_along(varnonfixvar)) {
      # Need to change marker variable if fixed
      oldmarker <- fixIntceptFac[[facinvarfree[i]]]
      if (length(oldmarker) > 0 && oldmarker == varfree[i]) {
        candidatemarker <- setdiff(facList[[facinvarfree[i]]], varnonfixvar[i])
        candidatemarker <- candidatemarker[numThreshold[candidatemarker] > 1][1]
        numcandidatethres <- numThreshold[candidatemarker]
        newparent <- constrainParTable(pt1, candidatemarker, "|", "t2", 1:ngroups)
        for (s in 2:numcandidatethres) {
          newparent <- freeParTable(newparent, varnonfixvar[i], "|", paste0("t", s), 1:ngroups)
        }
        temp <- newparent
        for (s in 2:numThreshold[varnonfixvar[i]]) {
          runnum <- which((newparent$lhs == varnonfixvar[i]) & (newparent$op == "|") & (newparent$rhs == paste0("t", s)) & (newparent$group == 1))
          temp <- constrainParTable(temp, newparent$lhs[runnum], newparent$op[runnum], newparent$rhs[runnum], 1:ngroups)
        }
        newparentresult <- try(newparentfit <- refit(newparent, fit1), silent = TRUE)
        if (!is(newparentresult, "try-error")) {
          tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
          if (!is(tryresult, "try-error")) {
            compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, newparentfit, method = method), silent = TRUE)
            if (!is(compresult, "try-error")) fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(newparentfit, tempfit))
          }
          args <- list(newparentfit, newparent, waldMat)
          for (s in 2:numThreshold[varnonfixvar[i]]) {
            runnum <- which((newparent$lhs == varnonfixvar[i]) & (newparent$op == "|") & (newparent$rhs == paste0("t", s)) & (newparent$group == 1))
            args <- c(args, list(cbind(newparent$lhs[runnum], newparent$op[runnum], newparent$rhs[runnum], 1:ngroups)))
          }
          waldCon[pos,] <- do.call(waldConstraint, args)
        }
      } else {
        temp <- pt1
        for (s in 2:numThreshold[varnonfixvar[i]]) {
          runnum <- which((pt1$lhs == varfree[i]) & (pt1$op == "|") & (pt1$rhs == paste0("t", s)) & (pt1$group == 1))
          temp <- constrainParTable(temp, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
        }
        tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
        if (!is(tryresult, "try-error")) {
          compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
          if (!is(compresult, "try-error"))  fixCon[pos,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
        }
        args <- list(fit1, pt1, waldMat)
        for (s in 2:numThreshold[varnonfixvar[i]]) {
          runnum <- which((pt1$lhs == varfree[i]) & (pt1$op == "|") & (pt1$rhs == paste0("t", s)) & (pt1$group == 1))
          args <- c(args, list(cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)))
        }
        waldCon[pos,] <- do.call(waldConstraint, args)
      }
      listFixCon <- c(listFixCon, tryresult)

      temp0 <- pt0
      for (s in 2:numThreshold[varnonfixvar[i]]) {
        runnum <- which((pt0$lhs == varfree[i]) & (pt0$op == "|") & (pt0$rhs == paste0("t", s)) & (pt0$group == 1))
        temp0 <- freeParTable(temp0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
        estimates[pos, s - 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
      }
      tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
      if (!is(tryresult0, "try-error")) {
        compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
        if (!is(compresult0, "try-error"))  freeCon[pos,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
        for (s in 2:numThreshold[varnonfixvar[i]]) {
          runnum <- which((pt0$lhs == varfree[i]) & (pt0$op == "|") & (pt0$rhs == paste0("t", s)) & (pt0$group == 1))
          thresVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
          estimates[pos, maxcolumns*(1:ngroups) + (s - 1)] <- thresVal
          totalVal <- sapply(thetaImpliedTotalVar(tempfit0), function(x, v) x[v, v], v = pt0$lhs[runnum])
          ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
          stdIntVal <- thresVal / sqrt(refTotalVal)
          stdestimates[pos, maxcolumns*(1:ngroups - 1) + (s - 1)] <- stdIntVal
          stdIntVal <- stdIntVal[grouporder]
          esstd[pos, maxcolumns*(1:length(complab) - 1) + (s - 1)] <- stdIntVal[2:ngroups] - stdIntVal[1]
        }
      }
      listFreeCon <- c(listFreeCon, tryresult0)
      pos <- pos + 1
    }
    freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
    fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
    waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)
    rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- paste0(c(varinfixvar, varnonfixvar), "|")
    estimates <- cbind(estimates, stdestimates, esstd)
    result <- cbind(freeCon, fixCon, waldCon)
  } else if (numType == 3) {
    if (!is.null(free) | !is.null(fix)) {
      if (!is.null(fix)) {
        for (i in seq_along(fix)) {
          pt0 <- constrainParTable(pt0, fix[i], "~~", fix[i], 1:ngroups)
          pt1 <- constrainParTable(pt1, fix[i], "~~", fix[i], 1:ngroups)
        }
      }
      if (!is.null(free)) {
        for (i in seq_along(free)) {
          pt0 <- freeParTable(pt0, free[i], "~~", free[i], 1:ngroups)
          pt1 <- freeParTable(pt1, free[i], "~~", free[i], 1:ngroups)
        }
      }
      namept1 <- paramNameFromPt(pt1)
      namept0 <- paramNameFromPt(pt0)
      fit0 <- refit(pt0, fit0)
      fit1 <- refit(pt1, fit1)
      beta <- lavaan::coef(fit1)
      beta0 <- lavaan::coef(fit0)
      waldMat <- matrix(0, ngroups - 1, length(beta))
      varfree <- setdiff(varfree, c(free, fix))
    }

    estimates <- matrix(NA, length(varfree), ngroups + 1)
    stdestimates <- matrix(NA, length(varfree), ngroups)
    colnames(estimates) <- c("poolest", paste0("errvar:", grouplab))
    colnames(stdestimates) <- paste0("std:", grouplab)
    esstd <- esz <- matrix(NA, length(varfree), ngroups - 1)
    colnames(esstd) <- paste0("diff_std:", complab)
    colnames(esz) <- paste0("h:", complab)
    fixCon <- freeCon <- matrix(NA, length(varfree), 4)
    waldCon <- matrix(NA, length(varfree), 3)
    colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
    colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
    colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")
    index <- which((pt1$lhs %in% varfree) & (pt1$op == "~~") & (pt1$lhs == pt1$rhs) & (pt1$group == 1))
    for (i in seq_along(index)) {
      runnum <- index[i]
      ustart <- getValue(pt1, beta, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1)
      temp <- fixParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 2:ngroups, ustart)
      tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
      if (!is(tryresult, "try-error")) {
        compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
        if (!is(compresult, "try-error"))  fixCon[i,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
      }
      listFixCon <- c(listFixCon, tryresult)
      temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
      estimates[i, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
      tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
      if (!is(tryresult0, "try-error")) {
        compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
        if (!is(compresult0, "try-error"))  freeCon[i,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
        errVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
        estimates[i, 2:ncol(estimates)] <- errVal
        totalVal <- sapply(thetaImpliedTotalVar(tempfit0), function(x, v) x[v, v], v = pt0$rhs[runnum])
        ifelse(poolvar, refTotalVal <- poolVariance(totalVal, neach), refTotalVal <- totalVal[refgroup])
        stdErrVal <- errVal / sqrt(refTotalVal)
        stdestimates[i,] <- stdErrVal
        stdErrVal <- stdErrVal[grouporder]
        esstd[i,] <- stdErrVal[2:ngroups] - stdErrVal[1]
        if (any(abs(stdErrVal) > 0.9999))
          warning(paste("The uniqueness of", pt0$rhs[runnum],
                        "in some groups are over 1. The uniqueness used in",
                        " arctan transformation are changed to 0.9999."))
        stdErrVal[stdErrVal > 0.9999] <- 0.9999
        zErrVal <- asin(sqrt(stdErrVal))
        esz[i,] <- zErrVal[2:ngroups] - zErrVal[1]
      }
      listFreeCon <- c(listFreeCon, tryresult0)
      waldCon[i,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
    }
    freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
    fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
    waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)
    rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[index]
    estimates <- cbind(estimates, stdestimates, esstd, esz)
    result <- cbind(freeCon, fixCon, waldCon)
  } else if (numType == 4) {
    varfree <- facnames
    if (!is.null(free) | !is.null(fix)) {
      if (!is.null(fix)) {
        for (i in seq_along(fix)) {
          pt0 <- constrainParTable(pt0, fix[i], "~1", "", 1:ngroups)
          pt1 <- constrainParTable(pt1, fix[i], "~1", "", 1:ngroups)
        }
      }
      if (!is.null(free)) {
        for (i in seq_along(free)) {
          pt0 <- freeParTable(pt0, free[i], "~1", "", 1:ngroups)
          pt1 <- freeParTable(pt1, free[i], "~1", "", 1:ngroups)
        }
      }
      namept1 <- paramNameFromPt(pt1)
      namept0 <- paramNameFromPt(pt0)
      fit0 <- refit(pt0, fit0)
      fit1 <- refit(pt1, fit1)
      beta <- lavaan::coef(fit1)
      beta0 <- lavaan::coef(fit0)
      waldMat <- matrix(0, ngroups - 1, length(beta))
      varfree <- setdiff(varfree, c(free, fix))
    }

    estimates <- matrix(NA, length(varfree), ngroups + 1)
    stdestimates <- matrix(NA, length(varfree), ngroups)
    colnames(estimates) <- c("poolest", paste0("mean:", grouplab))
    colnames(stdestimates) <- paste0("std:", grouplab)
    esstd <- matrix(NA, length(varfree), ngroups - 1)
    colnames(esstd) <- paste0("diff_std:", complab)
    fixCon <- freeCon <- matrix(NA, length(varfree), 4)
    waldCon <- matrix(NA, length(varfree), 3)
    colnames(fixCon) <- c("fix.chi", "fix.df", "fix.p", "fix.cfi")
    colnames(freeCon) <- c("free.chi", "free.df", "free.p", "free.cfi")
    colnames(waldCon) <- c("wald.chi", "wald.df", "wald.p")
    index <- which((pt1$lhs %in% varfree) & (pt1$op == "~1") & (pt1$group == 1))
    for (i in seq_along(index)) {
      runnum <- index[i]
      isfree <- pt1$free[runnum] != 0
      if (isfree) {
        temp <- constrainParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups)
      } else {
        temp <- fixParTable(pt1, pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 2:ngroups, ustart = pt1$ustart[runnum])
      }
      tryresult <- try(tempfit <- refit(temp, fit1), silent = TRUE)
      if (!is(tryresult, "try-error")) {
        compresult <- try(modelcomp <- lavaan::lavTestLRT(tempfit, fit1, method = method), silent = TRUE)
        if (!is(compresult, "try-error"))  fixCon[i,] <- c(unlist(modelcomp[2,5:7]), deltacfi(fit1, tempfit))
      }
      listFixCon <- c(listFixCon, tryresult)
      isfree0 <- pt0$free[runnum] != 0
      if (isfree0) {
        temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
      } else {
        temp0 <- freeParTable(pt0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 2:ngroups)
      }
      estimates[i, 1] <- getValue(pt0, beta0, pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1)
      tryresult0 <- try(tempfit0 <- refit(temp0, fit0), silent = TRUE)
      if (!is(tryresult0, "try-error")) {
        compresult0 <- try(modelcomp0 <- lavaan::lavTestLRT(tempfit0, fit0, method = method), silent = TRUE)
        if (!is(compresult0, "try-error"))  freeCon[i,] <- c(unlist(modelcomp0[2,5:7]), deltacfi(tempfit0, fit0))
        meanVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], pt0$op[runnum], pt0$rhs[runnum], 1:ngroups)
        estimates[i, 2:ncol(estimates)] <- meanVal
        facVal <- getValue(temp0, lavaan::coef(tempfit0), pt0$lhs[runnum], "~~", pt0$lhs[runnum], 1:ngroups)
        ifelse(poolvar, refFacVal <- poolVariance(facVal, neach), refFacVal <- facVal[refgroup])
        stdMeanVal <- meanVal / sqrt(refFacVal)
        stdestimates[i,] <- stdMeanVal
        stdMeanVal <- stdMeanVal[grouporder]
        esstd[i,] <- stdMeanVal[2:ngroups] - stdMeanVal[1]
      }
      listFreeCon <- c(listFreeCon, tryresult0)
      waldCon[i,] <- waldConstraint(fit1, pt1, waldMat, cbind(pt1$lhs[runnum], pt1$op[runnum], pt1$rhs[runnum], 1:ngroups))
    }
    freeCon[,3] <- stats::p.adjust(freeCon[,3], p.adjust)
    fixCon[,3] <- stats::p.adjust(fixCon[,3], p.adjust)
    waldCon[,3] <- stats::p.adjust(waldCon[,3], p.adjust)
    rownames(fixCon) <- names(listFixCon) <- rownames(freeCon) <- names(listFreeCon) <- rownames(waldCon) <- rownames(estimates) <- namept1[index]
    estimates <- cbind(estimates, stdestimates, esstd)
    result <- cbind(freeCon, fixCon, waldCon)
  }
  if (return.fit) {
    return(invisible(list(estimates = estimates, results = result, models = list(free = listFreeCon, fix = listFixCon, nested = fit0, parent = fit1))))
  } else {
    return(list(estimates = estimates, results = result))
  }
}


## ----------------
## Hidden Functions
## ----------------

findFactor <- function(var, facList) {
	tempfac <- lapply(facList, intersect, var)
	facinvar <- rep(names(tempfac), sapply(tempfac, length))
	facinvar[match(unlist(tempfac), var)]
}

## Terry moved here from wald.R so that wald() could be removed (redundant with lavaan::lavTestWald)
## FIXME: Update WaldConstraint to rely on lavaan::lavTestWald instead
#' @importFrom stats pchisq
waldContrast <- function(object, contrast) {
  beta <- lavaan::coef(object)
  acov <- lavaan::vcov(object)
  chisq <- t(contrast %*% beta) %*%  solve(contrast %*% as.matrix(acov) %*% t(contrast)) %*% (contrast %*% beta)
  df <- nrow(contrast)
  p <- pchisq(chisq, df, lower.tail=FALSE)
  c(chisq = chisq, df = df, p = p)
}

#' @importFrom lavaan parTable
waldConstraint <- function(fit, pt, mat, ...) {
	dotdotdot <- list(...)
	overallMat <- NULL
	for(i in seq_along(dotdotdot)) {
		target <- dotdotdot[[i]]
		tempMat <- mat
		element <- apply(target, 1, matchElement, parTable = pt)
		freeIndex  <- pt$free[element]
		tempMat[,freeIndex[1]] <- -1
		for(m in 2:length(freeIndex)) {
			tempMat[m - 1, freeIndex[m]] <- 1
		}
		overallMat <- rbind(overallMat, tempMat)
	}
	result <- rep(NA, 3)
	if(!any(apply(overallMat, 1, sum) != 0)) {
		try(result <- waldContrast(fit, overallMat), silent = TRUE)
	}
	return(result)
}

poolVariance <- function(var, n) {
	nm <- n - 1
	sum(var * nm) / sum(nm)
}

deltacfi <- function(parent, nested) lavaan::fitmeasures(nested)["cfi"] - lavaan::fitmeasures(parent)["cfi"]

## For categorical.
## FIXME: Why is this even necessary?
##        Did Sunthud not know implied Sigma is available?
#' @importFrom lavaan lavInspect
thetaImpliedTotalVar <- function(object) {
  # param <- lavInspect(object, "est")
  # ngroup <- lavInspect(object, "ngroups")
  # name <- names(param)
  # if(ngroup == 1) {
  # 	ly <- param[name == "lambda"]
  # } else {
  # 	ly <- lapply(param, "[[", "lambda")
  # }
  # ps <- lavInspect(object, "cov.lv")
  # if(ngroup == 1) ps <- list(ps)
  # if(ngroup == 1) {
  # 	te <- param[name == "theta"]
  # } else {
  # 	te <- lapply(param, "[[", "theta")
  # }
  # result <- list()
  # for(i in 1:ngroup) {
  # 	result[[i]] <- ly[[i]] %*% ps[[i]] %*% t(ly[[i]]) + te[[i]]
  # }
  # result
  if (lavInspect(object, "ngroups") == 1L) return(list(lavInspect(object, "cov.ov")))
  lavInspect(object, "cov.ov")
}




## MOVED FROM lonInvariance.R when it was removed from semTools 0.5-8 (9 Feb 2026)


# constrainParTable: Impose equality constraints in any set of elements in the parameter table
constrainParTable <- function(parTable, lhs, op, rhs, group) {
  parTable$start <- parTable$est <- parTable$se <- NULL
  target <- cbind(lhs, op, rhs, group)
  element <- apply(target, 1, matchElement, parTable=parTable)

  #     id lhs  op rhs user group free ustart exo label plabel  start
  for (i in 2:length(element)) {
    len <- length(parTable$id)
    newline <- list(lhs = parTable$plabel[element[1]], op = "==",
                    rhs = parTable$plabel[element[i]])
    if (!any(parTable$lhs == newline$lhs & parTable$op == newline$op &
             parTable$rhs == newline$rhs)) parTable <- patMerge(pt1 = parTable, pt2 = newline)
  }
  parTable
}

# matchElement: Find the number of row that have the specification in vec (lhs, op, rhs, group)
matchElement <- function(parTable, vec) {
  if (is.null(parTable$group)) {
    return(which((parTable$lhs == vec[1]) & (parTable$op == vec[2]) & (parTable$rhs == vec[3])))
  } else {
    return(which((parTable$lhs == vec[1]) & (parTable$op == vec[2]) & (parTable$rhs == vec[3]) & (parTable$group == vec[4])))
  }
}

getValue <- function(parTable, est, lhs, op, rhs, group) {
  target <- cbind(lhs, op, rhs, group)
  element <- apply(target, 1, matchElement, parTable = parTable)
  free <- parTable$free[element]
  out <- parTable$ustart[element]
  out[free != 0] <- est[free[free != 0]]
  out
}



