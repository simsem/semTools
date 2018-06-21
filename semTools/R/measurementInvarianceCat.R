### Sunthud Pornprasertmanit, Yves Rosseel, & Terrence D. Jorgensen
### Last updated: 14 May 2018
### automate measurement invariance tests for categorical indicators


#' Measurement Invariance Tests for Categorical Items
#'
#' Testing measurement invariance across groups using a typical sequence of
#' model comparison tests.
#'
#' Theta parameterization is used to represent SEM for categorical items.  That
#' is, residual variances are modeled instead of the total variance of
#' underlying normal variate for each item.  Five models can be tested based on
#' different constraints across groups.
#' \enumerate{
#'  \item Model 1: configural invariance. The same factor structure is imposed
#'   on all groups.
#'  \item Model 2: weak invariance. The factor loadings are constrained to be
#'   equal across groups.
#'  \item Model 3: strong invariance. The factor loadings and thresholds are
#'   constrained to be equal across groups.
#'  \item Model 4: strict invariance. The factor loadings, thresholds and
#'   residual variances are constrained to be equal across groups.
#'   For categorical variables, all residual variances are fixed as 1.
#'  \item Model 5: The factor loadings, threshoulds, residual variances and
#'   means are constrained to be equal across groups.
#' }
#'
#' However, if all items have two items (dichotomous), scalar invariance and
#' weak invariance cannot be separated because thresholds need to be equal
#' across groups for scale identification. Users can specify \code{strict}
#' option to include the strict invariance model for the invariance testing.
#' See the further details of scale identification and different
#' parameterization in Millsap and Yun-Tein (2004).
#'
#' @importFrom lavaan lavInspect parTable
#'
#' @param ... The same arguments as for any lavaan model.  See
#' \code{\link{cfa}} for more information.
#' @param std.lv If \code{TRUE}, the fixed-factor method of scale
#' identification is used. If \code{FALSE}, the first variable for each factor
#' is used as marker variable.
#' @param strict If \code{TRUE}, the sequence requires `strict' invariance.
#' See details for more information.
#' @param quiet If \code{FALSE} (default), a summary is printed out containing
#' an overview of the different models that are fitted, together with some
#' model comparison tests. If \code{TRUE}, no summary is printed.
#' @param fit.measures Fit measures used to calculate the differences between
#' nested models.
#' @param baseline.model custom baseline model passed to
#'  \code{\link[lavaan]{fitMeasures}}
#' @param method The method used to calculate likelihood ratio test. See
#' \code{\link[lavaan]{lavTestLRT}} for available options
#' @return Invisibly, all model fits in the sequence are returned as a list.
#' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#'
#'  Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
#'
#'  Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
#' @seealso \code{\link{measurementInvariance}} for measurement invariance for
#' continuous variables; \code{\link{longInvariance}} For the measurement
#' invariance test within person with continuous variables;
#' \code{partialInvariance} for the automated function for finding partial
#' invariance models
#' @references Millsap, R. E., & Yun-Tein, J. (2004). Assessing factorial
#' invariance in ordered-categorical measures. \emph{Multivariate Behavioral
#' Research, 39}(3), 479--515. doi:10.1207/S15327906MBR3903_4
#' @examples
#'
#' \dontrun{
#' syntax <- ' f1 =~ u1 + u2 + u3 + u4'
#'
#' measurementInvarianceCat(model = syntax, data = datCat, group = "g",
#'                          parameterization = "theta", estimator = "wlsmv",
#'                          ordered = c("u1", "u2", "u3", "u4"))
#' }
#'
#' @export
measurementInvarianceCat <- function(..., std.lv = FALSE, strict = FALSE,
                                     quiet = FALSE, fit.measures = "default",
                                     baseline.model = NULL,
                                     method = "default") {
	List <- list(...)
	if (is.null(List$model)) stop('all lavaan() and lavOptions() arguments must ',
	                              'named, including the "model=" argument.')
	lavaancfa <- function(...) { lavaan::cfa(...) }
	lavaanlavaan <- function(...) { lavaan::lavaan(...) }
	if (!is.null(List$parameterization) && tolower(List$parameterization) != "theta")
	  warning("The parameterization is set to 'theta' by default.")
	List$parameterization <- "theta"

	# Find the number of groups
	if (is.null(List$group)) stop("Please specify the group variable")

	# Get the lavaan parameter table
	template <- do.call(lavaancfa, c(List, do.fit = FALSE))
	lavaanParTable <- parTable(template)

	# Find the number of groups
	ngroups <- max(lavaanParTable$group)

	# Check whether all variables are categorical
	sampstat <- lavInspect(template, "samp")[[1]]
	meanname <- names(sampstat$mean)
	thname <- names(sampstat$th)
	if (any(is.na(charmatch(meanname, thname))))
	  stop("Some variables in your model are not identified as categorical.")

	varList <- lavaanParTable$rhs[lavaanParTable$op == "=~"]
	facName <- lavaanParTable$lhs[(lavaanParTable$op == "=~") & (lavaanParTable$rhs %in% varList)]
	if (length(unique(sapply(facName, function(x) length(x)))) > 1)
	  stop("The numbers of variables in each element are not equal.")
	varList <- unique(varList)
	facName <- unique(facName)

	# Check whether the factor configuration is the same across gorups
	groupParTable <- split(lavaanParTable, lavaanParTable$group)
	group1pt <- groupParTable[[1]]
	groupParTable <- lapply(groupParTable, "[", c("lhs", "op", "rhs"))
	if (!multipleAllEqualList(lapply(groupParTable, function(x) sapply(x, "[", x$op == "=~"))))
	  stop("Factor configuration is not the same across groups")

	# Extract the number of thresholds
	numThreshold <- table(sapply(group1pt, "[", group1pt$op == "|")[,"lhs"])

	# Find the indicators of each factor
	group1facload <- sapply(group1pt, "[", group1pt$op == "=~")
	factorRep <- split(group1facload[,"rhs"], group1facload[,"lhs"])

	# Find marker variables
	marker <- rep(NA, length(factorRep))
	numThresholdMarker <- rep(NA, length(factorRep))
	for (i in seq_along(factorRep)) {
		temp <- sapply(group1pt, "[", group1pt$rhs %in% factorRep[[i]] & group1pt$op == "=~" & group1pt$lhs == names(factorRep)[i])
		marker[i] <- temp[!is.na(temp[,"ustart"]), "rhs"]
		numThresholdMarker[i] <- numThreshold[marker[i]]
	}

	numThresholdFactorRep <- lapply(factorRep, function(x) numThreshold[x])
	constraintSecondThreshold <- unlist(lapply(numThresholdFactorRep, function(x) names(which(x > 1)[1])))
	constraintSecondThreshold <- constraintSecondThreshold[!is.na(constraintSecondThreshold)]
			# Find the marker variable of each facto

	for (i in names(numThreshold)) {
		lavaanParTable <- constrainParTable(lavaanParTable, i, "|", "t1", 1:ngroups)
	}

	if (length(constraintSecondThreshold) > 0) {
		for (i in constraintSecondThreshold) {
			lavaanParTable <- constrainParTable(lavaanParTable, i, "|", "t2", 1:ngroups)
		}
	}

	# Group 1
	for (i in facName) {
		lavaanParTable <- fixParTable(lavaanParTable, i, "~1", "", 1, 0) # Fix factor means as 0
		if (std.lv) {
			lavaanParTable <- fixParTable(lavaanParTable, i, "~~", i, 1, 1)
		} else {
			lavaanParTable <- freeParTable(lavaanParTable, i, "~~", i, 1, NA) # Free factor variances
		}
		# Assuming that all factor covariances are freeParTable
	}
	for (i in varList) {
		lavaanParTable <- fixParTable(lavaanParTable, i, "~~", i, 1, 1)
	}

	# Other groups
	for (k in 2:ngroups) {
		for (i in facName) {
			lavaanParTable <- freeParTable(lavaanParTable, i, "~1", "", k, NA)
			if (std.lv) {
				lavaanParTable <- fixParTable(lavaanParTable, i, "~~", i, k, 1)
			} else {
				lavaanParTable <- freeParTable(lavaanParTable, i, "~~", i, k, NA)
			}
		}
		for (i in varList) {
			lavaanParTable <- freeParTable(lavaanParTable, i, "~~", i, k, NA)
		}
		# Fix the indicator variances of marker variables with two categories as 1
		for (i in seq_along(marker)) {
			if (numThresholdMarker[i] == 1)  lavaanParTable <- fixParTable(lavaanParTable, marker[i], "~~", marker[i], k, 1)
		}
	}

	if (std.lv) {
		for (i in seq_along(factorRep)) {
			lavaanParTable <- freeParTable(lavaanParTable, names(factorRep)[i], "=~", marker[i], 1:ngroups, NA)
		}
	}
	# Fit configural invariance
	ListConfigural <- List
	ListConfigural$model <- lavaanParTable
	fitConfigural <- try(do.call(lavaanlavaan, ListConfigural), silent = TRUE)

	# Create the parameter table for metric invariance
	ptMetric <- lavaanParTable

	for (i in seq_along(factorRep)) {
		varwithin <- factorRep[[i]]
		if (!std.lv) {
			varwithin <- setdiff(varwithin, marker[i])
		}
		for (j in seq_along(varwithin)) {
			ptMetric <- constrainParTable(ptMetric, names(factorRep)[i], "=~", varwithin[j], 1:ngroups)
		}
	}
	if (std.lv) {
		for (k in 2:ngroups) {
			for (i in facName) {
				ptMetric <- freeParTable(ptMetric, i, "~~", i, k, NA)
			}
		}
	}

	ListMetric <- List
	ListMetric$model <- ptMetric
	fitMetric <- try(do.call(lavaanlavaan, ListMetric), silent = TRUE)

	ptMeans <- ptStrict <- ptMetric

	nonMarker <- setdiff(names(numThreshold), marker)
	nonDichoMarker <- numThreshold[which(numThreshold[nonMarker] > 1)]
	scalar <- length(nonDichoMarker) > 0
	if (scalar) {
		ptScalar <- ptMetric
		for (i in seq_along(numThreshold)) {
			thresholdName <- paste0("t", 1:numThreshold[i])
			for(j in seq_along(thresholdName)) {
				ptScalar <- constrainParTable(ptScalar, names(numThreshold)[i], "|", thresholdName[j], 1:ngroups)
			}
		}
		ListScalar <- List
		ListScalar$model <- ptScalar
		fitScalar <- try(do.call(lavaanlavaan, ListScalar), silent = TRUE)
		ptMeans <- ptStrict <- ptScalar
	} else fitScalar <- NULL


	fitStrict <- NULL
	# Create the parameter table for strict invariance if specified
	if (strict) {
	  if (scalar) ptStrict <- ptScalar
		for (k in 2:ngroups) {
			# Constrain measurement error variances
			for (i in varList) {
				ptStrict <- fixParTable(ptStrict, i, "~~", i, k, 1)
			}
		}
		ListStrict <- List
		ListStrict$model <- ptStrict
		fitStrict <- try(do.call(lavaanlavaan, ListStrict), silent = TRUE)
		ptMeans <- ptStrict
	}

	# Create the parameter table for mean equality

	# Constrain factor means to be equal
	for (k in 2:ngroups) {
		ptMeans <- fixParTable(ptMeans, facName, "~1", "", k, ustart = 0)
	}
	ListMeans <- List
	ListMeans$model <- ptMeans
	fitMeans <- try(do.call(lavaanlavaan, ListMeans), silent = TRUE)

	FIT <- invisible(list(fit.configural = fitConfigural, fit.loadings = fitMetric,
	                      fit.thresholds = fitScalar, fit.residuals = fitStrict,
	                      fit.means = fitMeans))
	FIT <- FIT[!sapply(FIT, is.null)]

	if (!quiet) {
        printInvarianceResult(FIT, fit.measures, baseline.model, method)
    }

    invisible(FIT)
}



## ----------------
## Hidden Functions
## ----------------

multipleAllEqual <- function(...) {
  obj <- list(...)
  multipleAllEqualList(obj)
}

multipleAllEqualList <- function(obj) {
  for (i in 2:length(obj)) {
    for (j in 1:(i - 1)) {
      temp <- isTRUE(all.equal(obj[[i]], obj[[j]]))
      if (!temp)
        return(FALSE)
    }
  }
  return(TRUE)
}

multipleAnyEqual <- function(...) {
  obj <- list(...)
  multipleAnyEqualList(obj)
}

multipleAnyEqualList <- function(obj) {
  for (i in 2:length(obj)) {
    for (j in 1:(i - 1)) {
      temp <- isTRUE(all.equal(obj[[i]], obj[[j]]))
      if (temp)
        return(TRUE)
    }
  }
  return(FALSE)
}


