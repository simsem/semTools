### Sunthud Pornprasertmanit, Yves Rosseel, and Terrence D. Jorgensen
### Last updated: 1 September 2018


##' Measurement Invariance Tests
##'
##' Testing measurement invariance across groups using a typical sequence of
##' model comparison tests.
##'
##' If \code{strict = FALSE}, the following four models are tested in order:
##' \enumerate{
##'  \item Model 1: configural invariance. The same factor structure
##' is imposed on all groups.
##'  \item Model 2: weak invariance. The factor loadings are constrained to
##'   be equal across groups.
##'  \item Model 3: strong invariance. The factor loadings and intercepts
##'   are constrained to be equal across groups.
##'  \item Model 4: The factor loadings, intercepts and means are constrained
##'   to be equal across groups.
##' }
##'
##' Each time a more restricted model is fitted, a \eqn{\Delta\chi^2} test is
##' reported, comparing the current model with the previous one, and comparing
##' the current model to the baseline model (Model 1). In addition, the
##' difference in CFI is also reported (\eqn{\Delta}CFI).
##'
##' If \code{strict = TRUE}, the following five models are tested in order:
##' \enumerate{
##'  \item Model 1: configural invariance. The same factor structure
##'    is imposed on all groups.
##'  \item Model 2: weak invariance. The factor loadings are constrained to be
##'    equal across groups.
##'  \item Model 3: strong invariance. The factor loadings and intercepts are
##'    constrained to be equal across groups.
##'  \item Model 4: strict invariance. The factor loadings, intercepts and
##'    residual variances are constrained to be equal across groups.
##'  \item Model 5: The factor loadings, intercepts, residual variances and means
##'    are constrained to be equal across groups.
##' }
##'
##' Note that if the \eqn{\chi^2} test statistic is scaled (e.g., a Satorra-Bentler
##' or Yuan-Bentler test statistic), a special version of the \eqn{\Delta\chi^2}
##' test is used as described in \url{http://www.statmodel.com/chidiff.shtml}
##'
##' @importFrom lavaan parTable
##'
##' @param ... The same arguments as for any lavaan model.  See
##'   \code{\link{cfa}} for more information.
##' @param std.lv If \code{TRUE}, the fixed-factor method of scale
##'   identification is used. If \code{FALSE}, the first variable for each factor
##'   is used as marker variable.
##' @param strict If \code{TRUE}, the sequence requires `strict' invariance.
##'   See details for more information.
##' @param quiet If \code{FALSE} (default), a summary is printed out containing
##'   an overview of the different models that are fitted, together with some
##'   model comparison tests. If \code{TRUE}, no summary is printed.
##' @param fit.measures Fit measures used to calculate the differences between
##'   nested models.
##' @param baseline.model custom baseline model passed to
##'   \code{\link[lavaan]{fitMeasures}}
##' @param method The method used to calculate likelihood ratio test. See
##'   \code{\link[lavaan]{lavTestLRT}} for available options
##'
##' @return Invisibly, all model fits in the sequence are returned as a list.
##'
##' @author Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
##'
##' Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
##'
##' @references
##'   Vandenberg, R. J., and Lance, C. E. (2000). A review and synthesis of the
##'   measurement invariance literature: Suggestions, practices, and
##'   recommendations for organizational research. \emph{Organizational
##'   Research Methods, 3,} 4--70.
##'
##' @examples
##'
##' HW.model <- ' visual =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed =~ x7 + x8 + x9 '
##'
##' measurementInvariance(model = HW.model, data = HolzingerSwineford1939,
##'                       group = "school", fit.measures = c("cfi","aic"))
##'
##' @name measurementInvariance-deprecated
##' @usage
##' measurementInvariance(..., std.lv = FALSE, strict = FALSE, quiet = FALSE,
##'                       fit.measures = "default", baseline.model = NULL,
##'                       method = "satorra.bentler.2001")
##' @seealso \code{\link{semTools-deprecated}}
##' @keywords internal
NULL


##' @rdname semTools-deprecated
##' @section Previous measurement-invariance functions:
##' The \code{measurementInvariance}, \code{measurementInvarianceCat}, and
##' \code{longInvariance} functions will no longer be supported. Instead, use
##' the \code{\link{measEq.syntax}} function, which is much more flexible and
##' supports a wider range of data (e.g., any mixture of \code{numeric} and
##' \code{ordered} indicators, any combination of multiple groups and repeated
##' measures, models fit to multiple imputations with \code{\link{runMI}}).
##'
##' @export
measurementInvariance <- function(..., std.lv = FALSE, strict = FALSE,
                                  quiet = FALSE,  fit.measures = "default",
                                  baseline.model = NULL,
                                  method = "satorra.bentler.2001") {

    .Deprecated(msg = c("The measurementInvariance function is deprecated, and",
                        " it will cease to be included in future versions of ",
                        "semTools. See help('semTools-deprecated) for details."))

	lavaancfa <- function(...) { lavaan::cfa(...) }
  ## check for a group.equal argument in ...
  dotdotdot <- list(...)
  if (is.null(dotdotdot$model)) stop('all lavaan() and lavOptions() arguments must',
                                     ' named, including the "model=" argument.')
  if (!is.null(dotdotdot$group.equal))
      stop("lavaan ERROR: group.equal argument should not be used")
  ## and a model
  if (names(dotdotdot)[1] == "") names(dotdotdot)[1] <- "model"

  res <- list()
  ## base-line model: configural invariance

	configural <- dotdotdot
	configural$group.equal <- ""
	template <- try(do.call(lavaancfa, configural), silent = TRUE)
	if (class(template) == "try-error") stop('Configural model did not converge.')
	pttemplate <- parTable(template)
	varnames <- unique(pttemplate$rhs[pttemplate$op == "=~"])
	facnames <- unique(pttemplate$lhs[(pttemplate$op == "=~") & (pttemplate$rhs %in% varnames)])
	ngroups <- max(pttemplate$group)
	if (ngroups <= 1) stop("Well, the number of groups is 1. Measurement",
	                       " invariance across 'groups' cannot be done.")

	if (std.lv) {
		for (i in facnames) {
			pttemplate <- fixParTable(pttemplate, i, "~~", i, 1:ngroups, 1)
		}
		fixloadings <- which(pttemplate$op == "=~" & pttemplate$free == 0)
		for (i in fixloadings) {
			pttemplate <- freeParTable(pttemplate, pttemplate$lhs[i], "=~",
			                           pttemplate$rhs[i], pttemplate$group[i])
		}
		dotdotdot$model <- pttemplate
		res$fit.configural <- try(do.call(lavaancfa, dotdotdot), silent = TRUE)
	} else {
		res$fit.configural <- template
	}

  ## fix loadings across groups
	if (std.lv) {
		findloadings <- which(pttemplate$op == "=~" & pttemplate$free != 0 & pttemplate$group == 1)
		for (i in findloadings) {
			pttemplate <- constrainParTable(pttemplate, pttemplate$lhs[i],
			                                "=~", pttemplate$rhs[i], 1:ngroups)
		}
		for (i in facnames) {
			pttemplate <- freeParTable(pttemplate, i, "~~", i, 2:ngroups)
		}
		dotdotdot$model <- pttemplate
		res$fit.loadings <- try(do.call(lavaancfa, dotdotdot), silent = TRUE)
	} else {
		loadings <- dotdotdot
		loadings$group.equal <- c("loadings")
		res$fit.loadings <- try(do.call(lavaancfa, loadings), silent = TRUE)
	}

  ## fix loadings + intercepts across groups
	if (std.lv) {
		findintcepts <- which(pttemplate$op == "~1" & pttemplate$lhs %in% varnames &
		                        pttemplate$free != 0 & pttemplate$group == 1)
		for (i in findintcepts) {
			pttemplate <- constrainParTable(pttemplate,
			                                pttemplate$lhs[i], "~1", "", 1:ngroups)
		}
		for (i in facnames) {
			pttemplate <- freeParTable(pttemplate, i, "~1", "", 2:ngroups)
		}
		dotdotdot$model <- pttemplate
		res$fit.intercepts <- try(do.call(lavaancfa, dotdotdot), silent = TRUE)
	} else {
		intercepts <- dotdotdot
		intercepts$group.equal <- c("loadings", "intercepts")
		res$fit.intercepts <- try(do.call(lavaancfa, intercepts), silent = TRUE)
	}

  if (strict) {
		if (std.lv) {
			findresiduals <- which(pttemplate$op == "~~" &
			                         pttemplate$lhs %in% varnames &
			                         pttemplate$rhs == pttemplate$lhs &
			                         pttemplate$free != 0 & pttemplate$group == 1)
			for (i in findresiduals) {
				pttemplate <- constrainParTable(pttemplate, pttemplate$lhs[i], "~~",
				                                pttemplate$rhs[i], 1:ngroups)
			}
			dotdotdot$model <- pttemplate
			res$fit.residuals <- try(do.call(lavaancfa, dotdotdot), silent = TRUE)
			for (i in facnames) {
				pttemplate <- fixParTable(pttemplate, i, "~1", "", 1:ngroups, 0)
			}
			dotdotdot$model <- pttemplate
			res$fit.means <- try(do.call(lavaancfa, dotdotdot), silent = TRUE)
		} else {
			# fix loadings + intercepts + residuals
			residuals <- dotdotdot
			residuals$group.equal <- c("loadings", "intercepts", "residuals")
			res$fit.residuals <- try(do.call(lavaancfa, residuals), silent = TRUE)

			# fix loadings + residuals + intercepts + means
			means <- dotdotdot
			means$group.equal <- c("loadings", "intercepts", "residuals", "means")
			res$fit.means <- try(do.call(lavaancfa, means), silent = TRUE)
		}
  } else {
		if (std.lv) {
			for (i in facnames) {
				pttemplate <- fixParTable(pttemplate, i, "~1", "", 1:ngroups, 0)
			}
		  dotdotdot$model <- pttemplate
		  res$fit.means <- try(do.call(lavaancfa, dotdotdot), silent = TRUE)
		} else {
			# fix loadings + intercepts + means
			means <- dotdotdot
			means$group.equal <- c("loadings", "intercepts", "means")
			res$fit.means <- try(do.call(lavaancfa, means), silent = TRUE)
		}
  }

	if (!quiet) printInvarianceResult(res, fit.measures, baseline.model, method)
  invisible(res)
}



## ----------------
## Hidden Functions
## ----------------

##' @importFrom lavaan lavInspect
##' @importMethodsFrom lavaan fitMeasures
printInvarianceResult <- function(FIT, fit.measures, baseline.model, method) {
  ## check whether models converged
  NAMES <- names(FIT)
  nonconv <- which(sapply(FIT, class) == "try-error")
  if (length(nonconv)) {
    message('The following model(s) did not converge: \n', paste(NAMES[nonconv], sep = "\n"))
    FIT <- FIT[-nonconv]
    NAMES <- NAMES[-nonconv]
  }
	names(FIT) <- NULL
	## compare models
	lavaanLavTestLRT <- function(...) lavaan::lavTestLRT(...)
	TABLE <- do.call(lavaanLavTestLRT, c(FIT, list(model.names = NAMES,
	                                               method = method)))

	if (length(fit.measures) == 1L && fit.measures == "default") {
		## scaled test statistic?
		if (length(lavInspect(FIT[[1]], "test")) > 1L) {
		  if (lavInspect(FIT[[1]], "test")[[2]]$test %in% c("satorra.bentler", "yuan.bentler")) {
		    fit.measures <- c("cfi.robust", "rmsea.robust")
		  } else fit.measures <- c("cfi.scaled", "rmsea.scaled")
		} else {
			fit.measures <- c("cfi", "rmsea")
		}
	}

	## add some fit measures
	if (length(fit.measures)) {

		FM <- lapply(FIT, fitMeasures,
		             fit.measures = fit.measures, baseline.model = baseline.model)
		FM.table1 <- sapply(fit.measures, function(x) sapply(FM, "[[", x))
		if (length(FM) == 1L) {
			FM.table1 <- rbind( rep(as.numeric(NA), length(fit.measures)), FM.table1)
		}
		if (length(FM) > 1L) {
			FM.table2 <- rbind(as.numeric(NA),
							   abs(apply(FM.table1, 2, diff)))
			colnames(FM.table2) <- paste(colnames(FM.table2), ".delta", sep = "")
			FM.TABLE <- as.data.frame(cbind(FM.table1, FM.table2))
		} else {
			FM.TABLE <- as.data.frame(FM.table1)
		}
		rownames(FM.TABLE) <- rownames(TABLE)
		class(FM.TABLE) <- c("lavaan.data.frame", "data.frame")
	}
	cat("\n")
	cat("Measurement invariance models:\n\n")
	cat(paste(paste("Model", seq_along(FIT), ":", NAMES), collapse = "\n"))
	cat("\n\n")

	print(TABLE)
	if (length(fit.measures)) {
		cat("\n\n")
		cat("Fit measures:\n\n")
		print(FM.TABLE)
		cat("\n")
		return(list(anova = TABLE, fitMeasures = FM.TABLE))
	}
	TABLE
}


