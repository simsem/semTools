### Sunthud Pornprasertmanit & Yves Rosseel
### Last updated: 10 January 2021
### Deprecated for years, removed 9 February 2026 (semTools 0.5-8)

##' Measurement Invariance Tests Within Person
##'
##' Testing measurement invariance across timepoints (longitudinal) or any
##' context involving the use of the same scale in one case (e.g., a dyad case
##' with husband and wife answering the same scale). The measurement invariance
##' uses a typical sequence of model comparison tests. This function currently
##' works with only one scale, and only with continuous indicators.
##'
##' If `strict = FALSE`, the following four models are tested in order:
##' \enumerate{
##' \item Model 1: configural invariance. The same factor structure is
##'   imposed on all units.
##' \item Model 2: weak invariance. The factor loadings are constrained to be
##'  equal across units.
##' \item Model 3: strong invariance. The factor loadings and intercepts are
##'  constrained to be equal across units.
##' \item Model 4: The factor loadings, intercepts and means are constrained to
##'  be equal across units.
##' }
##'
##' Each time a more restricted model is fitted, a \eqn{\Delta\chi^2} test is
##' reported, comparing the current model with the previous one, and comparing
##' the current model to the baseline model (Model 1). In addition, the
##' difference in CFA is also reported (\eqn{\Delta}CFI).
##'
##' If `strict = TRUE`, the following five models are tested in order:
##'
##' \enumerate{
##' \item Model 1: configural invariance. The same factor structure is imposed
##'  on all units.
##' \item Model 2: weak invariance. The factor loadings are constrained to be
##'  equal across units.
##' \item Model 3: strong invariance. The factor loadings and intercepts are
##'  constrained to be equal across units.
##' \item Model 4: strict invariance. The factor loadings, intercepts and
##'  residual variances are constrained to be equal across units.
##' \item Model 5: The factor loadings, intercepts, residual variances and
##'  means are constrained to be equal across units.
##' }
##'
##' Note that if the \eqn{\chi^2} test statistic is scaled (eg. a Satorra-Bentler
##' or Yuan-Bentler test statistic), a special version of the \eqn{\Delta\chi^2}
##' test is used as described in <http://www.statmodel.com/chidiff.shtml>
##'
##'
##' @param model lavaan syntax or parameter table
##' @param varList A list containing indicator names of factors used in the
##'   invariance testing, such as the list that the first element is the vector
##'   of indicator names in the first timepoint and the second element is the
##'   vector of indicator names in the second timepoint. The order of indicator
##'   names should be the same (but measured in different times or different
##'   units).
##' @param auto The order of autocorrelation on the measurement errors on the
##'   similar items across factor (e.g., Item 1 in Time 1 and Time 2). If 0 is
##'   specified, the autocorrelation will be not imposed. If 1 is specified,
##'   the autocorrelation will imposed for the adjacent factor listed in
##'   `varList`. The maximum number can be specified is the number of
##'   factors specified minus 1. If `"all"` is specified, the maximum
##'   number of order will be used.
##' @param constrainAuto If `TRUE`, the function will equate the
##'   auto-*covariance* to be equal within the same item across factors.
##'   For example, the covariance of item 1 in time 1 and time 2 is equal to
##'   the covariance of item 1 in time 2 and time 3.
##' @param fixed.x See [lavaan::lavOptions()].
##' @param std.lv See [lavaan::lavOptions()].
##' @param group See [lavaan::lavaan()].
##' @param group.equal See [lavaan::lavOptions()].
##' @param group.partial See [lavaan::lavOptions()].
##' @param strict If `TRUE`, the sequence requires strict invariance. See
##' @param warn See See [lavaan::lavOptions()].
##' @param debug See See [lavaan::lavOptions()]. details for more information.
##' @param quiet If `FALSE` (default), a summary is printed out containing
##'   an overview of the different models that are fitted, together with some
##'   model comparison tests. If `TRUE`, no summary is printed.
##' @param fit.measures Fit measures used to calculate the differences between
##'   nested models.
##' @param baseline.model custom baseline model passed to
##'  [lavaan::fitMeasures()]
##' @param method The method used to calculate likelihood ratio test. See
##'   [lavaan::lavTestLRT()] for available options
##' @param ... Additional arguments in the [lavaan::lavaan()]
##'   function. See also [lavaan::lavOptions()]
##'
##' @return Invisibly, all model fits in the sequence are returned as a list.
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##'  Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
##'
##'  Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references Vandenberg, R. J., and Lance, C. E. (2000). A review and
##'   synthesis of the measurement invariance literature: Suggestions,
##'   practices, and recommendations for organizational research.
##'   *Organizational Research Methods, 3*(1), 4--70.
##'   \doi{10.1177/109442810031002}
##'
##' @examples
##'
##' model <- ' f1t1 =~ y1t1 + y2t1 + y3t1
##'            f1t2 =~ y1t2 + y2t2 + y3t2
##' 			      f1t3 =~ y1t3 + y2t3 + y3t3 '
##'
##' ## Create list of variables
##' var1 <- c("y1t1", "y2t1", "y3t1")
##' var2 <- c("y1t2", "y2t2", "y3t2")
##' var3 <- c("y1t3", "y2t3", "y3t3")
##' constrainedVar <- list(var1, var2, var3)
##'
##' ## Invariance of the same factor across timepoints
##' longInvariance(model, auto = 1, constrainAuto = TRUE,
##'                varList = constrainedVar, data = exLong)
##'
##' ## Invariance of the same factor across timepoints and groups
##' longInvariance(model, auto = 1, constrainAuto = TRUE,
##'                varList = constrainedVar, data = exLong, group = "sex",
##' 	              group.equal = c("loadings", "intercepts"))
##'
##' @name longInvariance-deprecated
##' @usage
##' longInvariance(model, varList, auto = "all", constrainAuto = FALSE,
##'                fixed.x = TRUE, std.lv = FALSE, group = NULL,
##'                group.equal = "", group.partial = "", strict = FALSE,
##'                warn = TRUE, debug = FALSE, quiet = FALSE,
##'                fit.measures = "default", baseline.model = NULL,
##'                method = "satorra.bentler.2001", ...)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL


##' @rdname semTools-deprecated
##'
##' @export
longInvariance <- function(model, varList, auto = "all", constrainAuto = FALSE,
                           fixed.x = TRUE, std.lv = FALSE, group = NULL,
                           group.equal = "", group.partial = "", strict = FALSE,
                           warn = TRUE, debug = FALSE, quiet = FALSE,
                           fit.measures = "default", baseline.model = NULL,
                           method = "satorra.bentler.2001", ...) {

  .Deprecated(msg = c("The longInvariance function is deprecated, and ",
                      "it will cease to be included in future versions of ",
                      "semTools. See help('semTools-deprecated) for details."))

  List <- list(...)

	# Find the number of groups
	ngroups <- 1
	if (!is.null(group)) {
		if (!is.null(List$data)) {
			ngroups <- length(unique(List$data[,group]))
		} else if (!is.null(List$sample.cov)) {
			ngroups <- length(List$sample.cov)
		} else {
			stop("Cannot find the specifying variable name in the 'group' argument.")
		}
	}

	# Get the lavaan parameter table
    if (is.character(model)) {
      lavaanParTable <-
        lavaan::lavaanify(model           = model,
                          meanstructure   = TRUE,
                          int.ov.free     = TRUE,
                          int.lv.free     = FALSE,
                          orthogonal      = FALSE,
                          fixed.x         = fixed.x,
                          std.lv          = std.lv,

                          auto.fix.first  = ifelse(std.lv, FALSE, TRUE),
                          auto.fix.single = TRUE,
                          auto.var        = TRUE,
                          auto.cov.lv.x   = TRUE,
                          auto.cov.y      = TRUE,

                          ngroups         = ngroups,
                          group.equal     = group.equal,
                          group.partial   = group.partial,
                          debug           = debug,
                          warn            = warn,
                          as.data.frame.  = TRUE)
    } else if (is.list(model)) {
      if (!is.null(model$lhs) && !is.null(model$op)  &&
          !is.null(model$rhs) && !is.null(model$free)) {
        lavaanParTable <- model
      } else if (is.character(model[[1]])) {
        stop("lavaan ERROR: model is a list, but not a parameterTable?")
      }
    } else {
      cat("model type: ", class(model), "\n")
      stop("lavaan ERROR: model is not of type character or list")
    }

	# Error checking on the varList argument and get the factor name corresponding to each elements of the list
	facName <- lapply(varList,
	                  function(vec, pt) pt$lhs[(pt$op == "=~") & (pt$rhs %in% vec)],
	                  pt = lavaanParTable)
	if (any(sapply(facName, function(x) length(unique(x)) > 1)))
	  stop("The factor names of the same element of the 'varList' are not the same.")
	if (length(unique(sapply(facName, function(x) length(x)))) > 1)
	  stop("The numbers of variables in each element are not equal.")
	facName <- unlist(lapply(facName, unique))

	# Impose the autocorrelation in the parameter table
	if (auto != 0) {
		if (is.numeric(auto) && auto >= length(varList))
		  stop("The number of lag in auto-correlation is not possible in the current number of timepoints.")
		if (auto == "all") auto <- length(varList) - 1
		for (k in 1:ngroups) {
			for (i in 1:length(varList[[1]])) {
				name <- sapply(varList, function(x, element) x[element], element = i)
				for (j in 1:auto) {
					vec <- 1:(length(varList) - j)
					lavaanParTable <- freeParTable(lavaanParTable, name[vec], "~~", name[vec + j], k, ustart = NA)
					if (constrainAuto & (length(vec) > 1))
					  lavaanParTable <- constrainParTable(lavaanParTable, name[vec], "~~", name[vec + j], k)
				}
			}
		}
	}

	# Fit configural invariance
	fitConfigural <- try(lavaan::lavaan(lavaanParTable, ..., group = group,
	                                    group.equal = group.equal,
	                                    group.partial = group.partial,
	                                    warn = TRUE, debug = FALSE),
	                     silent = TRUE)

	# Create the parameter table for metric invariance
	ptMetric <- lavaanParTable
	if (std.lv) {
		for (k in 1:ngroups) {
			# Free variances of factor 2, 3, ...
			ptMetric <- freeParTable(ptMetric, facName[-1], "~~", facName[-1], k, ustart = NA)

			# Constrain factor loadings
			for (i in 1:length(varList[[1]])) {
				ptMetric <- constrainParTable(ptMetric, facName, "=~", sapply(varList, function(x, element) x[element], element = i), k)
			}
		}
		ptMetric$ustart[(ptMetric$op == "=~") & (ptMetric$rhs %in% sapply(varList, function(x, element) x[element], element = 1))] <- 1

	} else {
		for (k in 1:ngroups) {
			# Constrain factor loadings but keep marker variables
			for (i in 2:length(varList[[1]])) {
				ptMetric <- constrainParTable(ptMetric, facName, "=~", sapply(varList, function(x, element) x[element], element = i), k)
			}
		}
	}
	fitMetric <- try(lavaan::lavaan(ptMetric, ..., group = group,
	                                group.equal = group.equal,
	                                group.partial = group.partial,
	                                warn = TRUE, debug = FALSE),
	                 silent = TRUE)

	# Create the parameter table for scalar invariance
	ptScalar <- ptMetric
	for (k in 1:ngroups) {
		# Free means of factors 2, 3, ...
		ptScalar <- freeParTable(ptScalar, facName[-1], "~1", "", k, ustart = NA)

		# Constrain measurement intercepts
		for (i in 1:length(varList[[1]])) {
			ptScalar <- constrainParTable(ptScalar, sapply(varList, function(x, element) x[element], element = i), "~1", "", k)
		}
	}
	ptScalar$ustart[(ptMetric$op == "~1") & (ptMetric$rhs %in% facName)] <- 0
	fitScalar <- try(lavaan::lavaan(ptScalar, ..., group = group,
	                                group.equal = group.equal,
	                                group.partial = group.partial,
	                                warn = TRUE, debug = FALSE),
	                 silent = TRUE)

	ptMeans <- ptScalar

	# Create the parameter table for strict invariance if specified
	ptStrict <- ptScalar
	fitStrict <- NULL
	if (strict) {
		ptStrict <- ptScalar
		for (k in 1:ngroups) {
			# Constrain measurement error variances
			for (i in 1:length(varList[[1]])) {
				name <- sapply(varList, function(x, element) x[element], element = i)
				ptStrict <- constrainParTable(ptStrict, name, "~~", name, k)
			}
		}
		fitStrict <- try(lavaan::lavaan(ptStrict, ..., group = group,
		                                group.equal = group.equal,
		                                group.partial = group.partial,
		                                warn = TRUE, debug = FALSE),
		                 silent = TRUE)
		ptMeans <- ptStrict
	}

	# Create the parameter table for mean equality

	# Constrain factor means to be equal
	for (k in 1:ngroups) {
		ptMeans <- fixParTable(ptMeans, facName[-1], "~1", "", k, ustart = 0)
	}
	fitMeans <- try(lavaan::lavaan(ptMeans, ..., group = group,
	                               group.equal = group.equal,
	                               group.partial = group.partial,
	                               warn = TRUE, debug = FALSE),
	                silent = TRUE)

	FIT <- invisible(list(fit.configural = fitConfigural, fit.loadings = fitMetric,
	                      fit.intercepts = fitScalar, fit.residuals = fitStrict,
	                      fit.means = fitMeans))
	FIT <- FIT[!sapply(FIT, is.null)]

	if (!quiet) printInvarianceResult(FIT, fit.measures, baseline.model, method)

	# Modify these functions from measurementInvariance function
	# if(!quiet) {
		# cat("\n#################### Measurement invariance tests ####################\n")
		# cat("\nThe order of autocorrelation: ", auto, "\n")
		# cat("\n#################### Model 1: configural invariance:\n")
		# printFitLine(fitConfigural)

		# cat("\n#################### Model 2: weak invariance (equal loadings):\n")
		# printFitLine(fitMetric)

		# cat("\n[Model 1 versus model 2]\n")
		# difftest(fitConfigural, fitMetric)

		# cat("\n#################### Model 3: strong invariance (equal loadings + intercepts):\n")
		# printFitLine(fitScalar)
		# cat("\n[Model 1 versus model 3]\n")
		# difftest(fitConfigural, fitScalar)
		# cat("\n[Model 2 versus model 3]\n")
		# difftest(fitMetric, fitScalar)
		# if(strict) {
            # cat("\n#################### Model 4: strict invariance (equal loadings + intercepts + residuals):\n")
            # printFitLine(fitStrict)
            # cat("\n[Model 1 versus model 4]\n")
            # difftest(fitConfigural, fitStrict)
            # cat("\n[Model 2 versus model 4]\n")
            # difftest(fitMetric, fitStrict)
            # cat("\n[Model 3 versus model 4]\n")
            # difftest(fitScalar, fitStrict)

            # cat("\n#################### Model 5: equal loadings + intercepts + residuals + means:\n")
            # printFitLine(fitMeans, horizontal=TRUE)
            # cat("\n[Model 1 versus model 5]\n")
            # difftest(fitConfigural, fitMeans)
            # cat("\n[Model 2 versus model 5]\n")
            # difftest(fitMetric, fitMeans)
            # cat("\n[Model 3 versus model 5]\n")
            # difftest(fitScalar, fitMeans)
            # cat("\n[Model 4 versus model 5]\n")
            # difftest(fitStrict, fitMeans)
        # } else {
            # cat("\n#################### Model 4: equal loadings + intercepts + means:\n")
            # printFitLine(fitMeans)
            # cat("\n[Model 1 versus model 4]\n")
            # difftest(fitConfigural, fitMeans)
            # cat("\n[Model 2 versus model 4]\n")
            # difftest(fitMetric, fitMeans)
            # cat("\n[Model 3 versus model 4]\n")
            # difftest(fitScalar, fitMeans)
        # }
	# }
	# return(invisible(list(fit.configural = fitConfigural, fit.loadings = fitMetric, fit.intercepts = fitScalar, fit.residuals = fitStrict, fit.means = fitMeans)))
	invisible(FIT)
}




