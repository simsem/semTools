### Terrence D. Jorgensen
### Last updated: 25 June 2018
### permutation randomization test for measurement equivalence and DIF


## -----------------
## Class and Methods
## -----------------

##' Class for the Results of Permutation Randomization Tests of Measurement
##' Equivalence and DIF
##'
##' This class contains the results of tests of Measurement Equivalence and
##' Differential Item Functioning (DIF).
##'
##'
##' @name permuteMeasEq-class
##' @aliases permuteMeasEq-class show,permuteMeasEq-method
##' summary,permuteMeasEq-method hist,permuteMeasEq-method
##' @docType class
##'
##' @slot PT A \code{data.frame} returned by a call to
##'   \code{\link[lavaan]{parTable}} on the constrained model
##' @slot modelType A character indicating the specified \code{modelType} in the
##'   call to \code{permuteMeasEq}
##' @slot ANOVA A \code{numeric} vector indicating the results of the observed
##'   (\eqn{\Delta})\eqn{\chi^2} test, based on the central \eqn{\chi^2}
##'   distribution
##' @slot AFI.obs A vector of observed (changes in) user-selected fit measures
##' @slot AFI.dist The permutation distribution(s) of user-selected fit measures.
##'   A \code{data.frame} with \code{n.Permutations} rows and one column for each
##'   \code{AFI.obs}.
##' @slot AFI.pval A vector of \emph{p} values (one for each element in slot
##'   \code{AFI.obs}) calculated using slot \code{AFI.dist}, indicating the
##'   probability of observing a change at least as extreme as \code{AFI.obs}
##'   if the null hypothesis were true
##' @slot MI.obs A \code{data.frame} of observed Lagrange Multipliers
##'   (modification indices) associated with the equality constraints or fixed
##'   parameters specified in the \code{param} argument. This is a subset of the
##'   output returned by a call to \code{\link[lavaan]{lavTestScore}} on the
##'   constrained model.
##' @slot MI.dist The permutation distribution of the maximum modification index
##'   (among those seen in slot \code{MI.obs$X2}) at each permutation of group
##'   assignment or of \code{covariates}
##' @slot extra.obs If \code{permuteMeasEq} was called with an \code{extra}
##'   function, the output when applied to the original data is concatenated
##'   into this vector
##' @slot extra.dist A \code{data.frame}, each column of which contains the
##'   permutation distribution of the corresponding statistic in slot
##'   \code{extra.obs}
##' @slot n.Permutations An \code{integer} indicating the number of permutations
##'   requested by the user
##' @slot n.Converged An \code{integer} indicating the number of permuation
##'   iterations which yielded a converged solution
##' @slot n.nonConverged An \code{integer} vector of length
##'   \code{n.Permutations} indicating how many times group assignment was
##'   randomly permuted (at each iteration) before converging on a solution
##' @slot n.Sparse Only relevant with \code{ordered} indicators when
##'   \code{modelType == "mgcfa"}. An \code{integer} vector of length
##'   \code{n.Permutations} indicating how many times group assignment was
##'   randomly permuted (at each iteration) before obtaining a sample with all
##'   categories observed in all groups.
##' @slot oldSeed An \code{integer} vector storing the value of
##'   \code{.Random.seed} before running \code{permuteMeasEq}. Only relevant
##'   when using a parallel/multicore option and the original
##'   \code{RNGkind() != "L'Ecuyer-CMRG"}. This enables users to restore their
##'   previous \code{.Random.seed} state, if desired, by running:
##'   \code{.Random.seed[-1] <- permutedResults@oldSeed[-1]}
##' @section Objects from the Class: Objects can be created via the
##'   \code{\link[semTools]{permuteMeasEq}} function.
##'
##' @return
##' \itemize{
##' \item The \code{show} method prints a summary of the multiparameter
##'   omnibus test results, using the user-specified AFIs. The parametric
##'  (\eqn{\Delta})\eqn{\chi^2} test is also displayed.
##' \item The \code{summary} method prints the same information from the
##'   \code{show} method, but when \code{extra = FALSE} (the default) it also
##'   provides a table summarizing any requested follow-up tests of DIF using
##'   modification indices in slot \code{MI.obs}. The user can also specify an
##'   \code{alpha} level for flagging modification indices as significant, as
##'   well as \code{nd} (the number of digits displayed). For each modification
##'   index, the \emph{p} value is displayed using a central \eqn{\chi^2}
##'   distribution with the \emph{df} shown in that column. Additionally, a
##'   \emph{p} value is displayed using the permutation distribution of the
##'   maximum index, which controls the familywise Type I error rate in a manner
##'   similar to Tukey's studentized range test. If any indices are flagged as
##'   significant using the \code{tukey.p.value}, then a message is displayed for
##'   each flagged index. The invisibly returned \code{data.frame} is the
##'   displayed table of modification indices, unless
##'   \code{\link[semTools]{permuteMeasEq}} was called with \code{param = NULL},
##'   in which case the invisibly returned object is \code{object}. If
##'   \code{extra = TRUE}, the permutation-based \emph{p} values for each
##'   statistic returned by the \code{extra} function are displayed and returned
##'   in a \code{data.frame} instead of the modification indices requested in the
##'   \code{param} argument.
##' \item The \code{hist} method returns a list of \code{length == 2},
##'    containing the arguments for the call to \code{hist} and the arguments
##'    to the call for \code{legend}, respectively. This list may facilitate
##'    creating a customized histogram of \code{AFI.dist}, \code{MI.dist}, or
##'    \code{extra.dist}
##' }
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' @seealso \code{\link[semTools]{permuteMeasEq}}
##'
##' @examples
##'
##' # See the example from the permuteMeasEq function
##'
setClass("permuteMeasEq", slots = c(PT = "data.frame",
                                    modelType = "character",
                                    ANOVA = "vector",
                                    AFI.obs = "vector",
                                    AFI.dist = "data.frame",
                                    AFI.pval = "vector",
                                    MI.obs = "data.frame",
                                    MI.dist = "vector",
                                    extra.obs = "vector",
                                    extra.dist = "data.frame",
                                    n.Permutations = "integer",
                                    n.Converged = "integer",
                                    n.nonConverged = "vector",
                                    n.Sparse = "vector",
                                    oldSeed = "integer"))


##' @rdname permuteMeasEq-class
##' @aliases show,permuteMeasEq-method
##' @export
setMethod("show", "permuteMeasEq", function(object) {
  ## print warning if there are nonConverged permutations
  if (object@n.Permutations != object@n.Converged) {
    warning(paste("Only", object@n.Converged, "out of",
                  object@n.Permutations, "models converged within",
                  max(object@n.nonConverged), "attempts per permutation.\n\n"))
  }
  ## print ANOVA
  cat("Omnibus p value based on parametric chi-squared difference test:\n\n")
  print(round(object@ANOVA, digits = 3))
  ## print permutation results
  cat("\n\nOmnibus p values based on nonparametric permutation method: \n\n")
  AFI <- data.frame(AFI.Difference = object@AFI.obs, p.value = object@AFI.pval)
  class(AFI) <- c("lavaan.data.frame","data.frame")
  print(AFI, nd = 3)
  invisible(object)
})

##' @rdname permuteMeasEq-class
##' @aliases summary,permuteMeasEq-method
##' @export
setMethod("summary", "permuteMeasEq", function(object, alpha = .05, nd = 3,
                                               extra = FALSE) {
  ## print warning if there are nonConverged permutations
  if (object@n.Permutations != object@n.Converged) {
    warning(paste("Only", object@n.Converged, "out of",
                  object@n.Permutations, "models converged within",
                  max(object@n.nonConverged), "attempts per permutation.\n\n"))
  }
  ## print ANOVA
  cat("Omnibus p value based on parametric chi-squared difference test:\n\n")
  print(round(object@ANOVA, digits = nd))
  ## print permutation results
  cat("\n\nOmnibus p values based on nonparametric permutation method: \n\n")
  AFI <- data.frame(AFI.Difference = object@AFI.obs, p.value = object@AFI.pval)
  class(AFI) <- c("lavaan.data.frame","data.frame")
  print(AFI, nd = nd)

  ## print extras or DIF test results, if any were requested
  if (extra && length(object@extra.obs)) {
    cat("\n\nUnadjusted p values of extra statistics,\n",
        "based on permutation distribution of each statistic: \n\n")
    MI <- data.frame(Statistic = object@extra.obs)
    class(MI) <- c("lavaan.data.frame","data.frame")
    MI$p.value <- sapply(names(object@extra.dist), function(nn) {
      mean(abs(object@extra.dist[,nn]) >= abs(object@extra.obs[nn]), na.rm = TRUE)
    })
    MI$flag <- ifelse(MI$p.value < alpha, "*   ", "")
    print(MI, nd = nd)
  } else if (length(object@MI.dist)) {
    cat("\n\n Modification indices for equality constrained parameter estimates,\n",
        "with unadjusted 'p.value' based on chi-squared distribution and\n",
        "adjusted 'tukey.p.value' based on permutation distribution of the\n",
        "maximum modification index per iteration: \n\n")
    MI <- do.call(paste("summ", object@modelType, sep = "."),
                  args = list(object = object, alpha = alpha))
    print(MI, nd = nd)

    ## print messages about potential DIF
    if (all(MI$tukey.p.value > alpha)) {
      cat("\n\n No equality constraints were flagged as significant.\n\n")
      return(invisible(MI))
    }
    if (object@modelType == "mgcfa") {
      cat("\n\nThe following equality constraints were flagged as significant:\n\n")
      for (i in which(MI$tukey.p.value < alpha)) {
        cat("Parameter '", MI$parameter[i], "' may differ between Groups '",
            MI$group.lhs[i], "' and '", MI$group.rhs[i], "'.\n", sep = "")
      }
      cat("\nUse lavTestScore(..., epc = TRUE) on your constrained model to",
          "display expected parameter changes for these equality constraints\n\n")
    }

  } else return(invisible(object))

  invisible(MI)
})

summ.mgcfa <- function(object, alpha) {
  MI <- object@MI.obs
  class(MI) <- c("lavaan.data.frame","data.frame")
  PT <- object@PT
  eqPar <- rbind(PT[PT$plabel %in% MI$lhs, ], PT[PT$plabel %in% MI$rhs, ])
  MI$flag <- ""
  MI$parameter <- ""
  MI$group.lhs <- ""
  MI$group.rhs <- ""
  for (i in 1:nrow(MI)) {
    par1 <- eqPar$par[ eqPar$plabel == MI$lhs[i] ]
    par2 <- eqPar$par[ eqPar$plabel == MI$rhs[i] ]
    MI$parameter[i] <- par1
    MI$group.lhs[i] <- eqPar$group.label[ eqPar$plabel == MI$lhs[i] ]
    MI$group.rhs[i] <- eqPar$group.label[ eqPar$plabel == MI$rhs[i] ]
    if (par1 != par2) {
      myMessage <- paste0("Constraint '", MI$lhs[i], "==", MI$rhs[i],
                          "' refers to different parameters: \n'",
                          MI$lhs[i], "' is '", par1, "' in group '",
                          MI$group.lhs[i], "'\n'",
                          MI$rhs[i], "' is '", par2, "' in group '",
                          MI$group.rhs[i], "'\n")
      warning(myMessage)
    }
    if (MI$tukey.p.value[i] < alpha) MI$flag[i] <- "*  -->"
  }
  MI
}

summ.mimic <- function(object, alpha) {
  MI <- object@MI.obs
  class(MI) <- c("lavaan.data.frame","data.frame")
  MI$flag <- ifelse(MI$tukey.p.value < alpha, "*   ", "")
  MI
}


##' @rdname permuteMeasEq-class
##' @aliases hist,permuteMeasEq-method
##' @importFrom stats qchisq dchisq quantile
##' @param object,x object of class \code{permuteMeasEq}
##' @param ... Additional arguments to pass to \code{\link[graphics]{hist}}
##' @param AFI \code{character} indicating the fit measure whose permutation
##'  distribution should be plotted
##' @param alpha alpha level used to draw confidence limits in \code{hist} and
##'   flag significant statistics in \code{summary} output
##' @param nd number of digits to display
##' @param extra \code{logical} indicating whether the \code{summary} output
##'   should return permutation-based \emph{p} values for each statistic returned
##'   by the \code{extra} function.  If \code{FALSE} (default), \code{summary}
##'   will return permutation-based \emph{p} values for each modification index.
##' @param printLegend \code{logical}. If \code{TRUE} (default), a legend will
##'  be printed with the histogram
##' @param legendArgs \code{list} of arguments passed to the
##'  \code{\link[graphics]{legend}} function.  The default argument is a list
##'  placing the legend at the top-left of the figure.
##' @export
setMethod("hist", "permuteMeasEq", function(x, ..., AFI, alpha = .05, nd = 3,
                                            printLegend = TRUE,
                                            legendArgs = list(x = "topleft")) {
  histArgs <- list(...)
  histArgs$x <- x@AFI.dist[[AFI]]
  if (is.null(histArgs$col)) histArgs$col <- "grey69"
  histArgs$freq <- !grepl("chi", AFI)
  histArgs$ylab <- if (histArgs$freq) "Frequency" else "Probability Density"

  if (printLegend) {
    if (is.null(legendArgs$box.lty)) legendArgs$box.lty <- 0
    if (nd < length(strsplit(as.character(1 / alpha), "")[[1]]) - 1) {
      warning(paste0("The number of digits argument (nd = ", nd ,
                     ") is too low to display your p value at the",
                     " same precision as your requested alpha level (alpha = ",
                     alpha, ")"))
    }
    if (x@AFI.pval[[AFI]] < (1 / 10^nd)) {
      pVal <- paste(c("< .", rep(0, nd - 1),"1"), collapse = "")
    } else {
      pVal <- paste("=", round(x@AFI.pval[[AFI]], nd))
    }
  }

  delta <- length(x@MI.dist) > 0L && x@modelType == "mgcfa"
  if (grepl("chi", AFI)) {   ####################################### Chi-squared
    ChiSq <- x@AFI.obs[AFI]
    DF <- x@ANOVA[2]
    histArgs$xlim <- range(c(ChiSq, x@AFI.dist[[AFI]], qchisq(c(.01, .99), DF)))
    xVals <- seq(histArgs$xlim[1], histArgs$xlim[2], by = .1)
    theoDist <- dchisq(xVals, df = DF)
    TheoCrit <- round(qchisq(p = alpha, df = DF, lower.tail = FALSE), 2)
    Crit <- quantile(histArgs$x, probs = 1 - alpha)
    if (ChiSq > histArgs$xlim[2]) histArgs$xlim[2] <- ChiSq
    if (delta) {
      histArgs$main <- expression(Permutation~Distribution~of~Delta*chi^2)
      histArgs$xlab <- expression(Delta*chi^2)
      if (printLegend) {
        legendArgs$legend <- c(bquote(Theoretical~Delta*chi[Delta*.(paste("df =", DF))]^2 ~ Distribution),
                               bquote(Critical~chi[alpha~.(paste(" =", alpha))]^2 == .(round(TheoCrit, nd))),
                               bquote(.(paste("Permuted Critical Value =", round(Crit, nd)))),
                               bquote(Observed~Delta*chi^2 == .(round(ChiSq, nd))),
                               expression(paste("")),
                               bquote(Permuted~italic(p)~.(pVal)))
      }
    } else {
      histArgs$main <- expression(Permutation~Distribution~of~chi^2)
      histArgs$xlab <- expression(chi^2)
      if (printLegend) {
        legendArgs$legend <- c(bquote(Theoretical~chi[.(paste("df =", DF))]^2 ~ Distribution),
                               bquote(Critical~chi[alpha~.(paste(" =", alpha))]^2 == .(round(TheoCrit, nd))),
                               bquote(.(paste("Permuted Critical Value =", round(Crit, nd)))),
                               bquote(Observed~chi^2 == .(round(ChiSq, nd))),
                               expression(paste("")),
                               bquote(Permuted~italic(p)~.(pVal)))
      }
    }
    H <- do.call(hist, c(histArgs["x"], plot = FALSE))
    histArgs$ylim <- c(0, max(H$density, theoDist))
    if (printLegend) {
      legendArgs <- c(legendArgs, list(lty = c(2, 2, 1, 1, 0, 0),
                                       lwd = c(2, 2, 2, 3, 0, 0),
                                       col = c("black","black","black","red","","")))
    }
  } else {        ################################################### other AFIs
    badness <- grepl(pattern = "fmin|aic|bic|rmr|rmsea|cn|sic|hqc",
                     x = AFI, ignore.case = TRUE)
    if (badness) {
      Crit <- quantile(histArgs$x, probs = 1 - alpha)
    } else {
      Crit <- quantile(histArgs$x, probs = alpha)
    }
    histArgs$xlim <- range(histArgs$x, x@AFI.obs[AFI])
    if (delta) {
      histArgs$main <- bquote(~Permutation~Distribution~of~Delta*.(toupper(AFI)))
      histArgs$xlab <- bquote(~Delta*.(toupper(AFI)))
      if (printLegend) {
        legendArgs$legend <- c(bquote(Critical~Delta*.(toupper(AFI))[alpha~.(paste(" =", alpha))] == .(round(Crit, nd))),
                               bquote(Observed~Delta*.(toupper(AFI)) == .(round(x@AFI.obs[AFI], nd))),
                               expression(paste("")),
                               bquote(Permuted~italic(p)~.(pVal)))

      }
    } else {
      histArgs$main <- paste("Permutation Distribution of", toupper(AFI))
      histArgs$xlab <- toupper(AFI)
      if (printLegend) {
        legendArgs$legend <- c(bquote(Critical~.(toupper(AFI))[alpha~.(paste(" =", alpha))] == .(round(Crit, nd))),
                               bquote(Observed~.(toupper(AFI)) == .(round(x@AFI.obs[AFI], nd))),
                               expression(paste("")),
                               bquote(Permuted~italic(p)~.(pVal)))

      }
    }
    if (printLegend) {
      legendArgs <- c(legendArgs, list(lty = c(1, 1, 0, 0),
                                       lwd = c(2, 3, 0, 0),
                                       col = c("black","red","","")))
    }
  }
  ## print histogram (and optionally, print legend)
  suppressWarnings({
    do.call(hist, histArgs)
    if (grepl("chi", AFI)) {
      lines(x = xVals, y = theoDist, lwd = 2, lty = 2)
      abline(v = TheoCrit, col = "black", lwd = 2, lty = 2)
    }
    abline(v = Crit, col = "black", lwd = 2)
    abline(v = x@AFI.obs[AFI], col = "red", lwd = 3)
    if (printLegend) do.call(legend, legendArgs)
  })
  ## return arguments to create histogram (and optionally, legend)
  invisible(list(hist = histArgs, legend = legendArgs))
})



## --------------------
## Constructor Function
## --------------------

##' Permutation Randomization Tests of Measurement Equivalence and Differential
##' Item Functioning (DIF)
##'
##' The function \code{permuteMeasEq} provides tests of hypotheses involving
##' measurement equivalence, in one of two frameworks: multigroup CFA or MIMIC
##' models.
##'
##'
##' The function \code{permuteMeasEq} provides tests of hypotheses involving
##' measurement equivalence, in one of two frameworks:
##' \enumerate{
##'   \item{1} For multiple-group CFA models, provide a pair of nested lavaan objects,
##'   the less constrained of which (\code{uncon}) freely estimates a set of
##'   measurement parameters (e.g., factor loadings, intercepts, or thresholds;
##'   specified in \code{param}) in all groups, and the more constrained of which
##'   (\code{con}) constrains those measurement parameters to equality across
##'   groups. Group assignment is repeatedly permuted and the models are fit to
##'   each permutation, in order to produce an empirical distribution under the
##'   null hypothesis of no group differences, both for (a) changes in
##'   user-specified fit measures (see \code{AFIs} and \code{moreAFIs}) and for
##'   (b) the maximum modification index among the user-specified equality
##'   constraints. Configural invariance can also be tested by providing that
##'   fitted lavaan object to \code{con} and leaving \code{uncon = NULL}, in which
##'   case \code{param} must be \code{NULL} as well.
##'
##'   \item{2} In MIMIC models, one or a set of continuous and/or discrete
##'   \code{covariates} can be permuted, and a constrained model is fit to each
##'   permutation in order to provide a distribution of any fit measures (namely,
##'   the maximum modification index among fixed parameters in \code{param}) under
##'   the null hypothesis of measurement equivalence across levels of those
##'   covariates.
##' }
##'
##' In either framework, modification indices for equality constraints or fixed
##' parameters specified in \code{param} are calculated from the constrained
##' model (\code{con}) using the function \code{\link[lavaan]{lavTestScore}}.
##'
##' For multiple-group CFA models, the multiparameter omnibus null hypothesis of
##' measurement equivalence/invariance is that there are no group differences in
##' any measurement parameters (of a particular type). This can be tested using
##' the \code{anova} method on nested \code{lavaan} objects, as seen in the
##' output of \code{\link[semTools]{measurementInvariance}}, or by inspecting
##' the change in alternative fit indices (AFIs) such as the CFI. The
##' permutation randomization method employed by \code{permuteMeasEq} generates
##' an empirical distribution of any \code{AFIs} under the null hypothesis, so
##' the user is not restricted to using fixed cutoffs proposed by Cheung &
##' Rensvold (2002), Chen (2007), or Meade, Johnson, & Braddy (2008).
##'
##' If the multiparameter omnibus null hypothesis is rejected, partial
##' invariance can still be established by freeing invalid equality constraints,
##' as long as equality constraints are valid for at least two indicators per
##' factor. Modification indices can be calculated from the constrained model
##' (\code{con}), but multiple testing leads to inflation of Type I error rates.
##' The permutation randomization method employed by \code{permuteMeasEq}
##' creates a distribution of the maximum modification index if the null
##' hypothesis is true, which allows the user to control the familywise Type I
##' error rate in a manner similar to Tukey's \emph{q} (studentized range)
##' distribution for the Honestly Significant Difference (HSD) post hoc test.
##'
##' For MIMIC models, DIF can be tested by comparing modification indices of
##' regression paths to the permutation distribution of the maximum modification
##' index, which controls the familywise Type I error rate. The MIMIC approach
##' could also be applied with multiple-group models, but the grouping variable
##' would not be permuted; rather, the covariates would be permuted separately
##' within each group to preserve between-group differences. So whether
##' parameters are constrained or unconstrained across groups, the MIMIC
##' approach is only for testing null hypotheses about the effects of
##' \code{covariates} on indicators, controlling for common factors.
##'
##' In either framework, \code{\link[lavaan]{lavaan}}'s \code{group.label}
##' argument is used to preserve the order of groups seen in \code{con} when
##' permuting the data.
##'
##'
##' @importFrom lavaan lavInspect parTable
##'
##' @param nPermute An integer indicating the number of random permutations used
##' to form empirical distributions under the null hypothesis.
##' @param modelType A character string indicating type of model employed:
##' multiple-group CFA (\code{"mgcfa"}) or MIMIC (\code{"mimic"}).
##' @param con The constrained \code{lavaan} object, in which the parameters
##' specified in \code{param} are constrained to equality across all groups when
##' \code{modelType = "mgcfa"}, or which regression paths are fixed to zero when
##' \code{modelType = "mimic"}. In the case of testing \emph{configural}
##' invariance when \code{modelType = "mgcfa"}, \code{con} is the configural
##' model (implicitly, the unconstrained model is the saturated model, so use
##' the defaults \code{uncon = NULL} and \code{param = NULL}). When
##' \code{modelType = "mimic"}, \code{con} is the MIMIC model in which the
##' covariate predicts the latent construct(s) but no indicators (unless they
##' have already been identified as DIF items).
##' @param uncon Optional.  The unconstrained \code{lavaan} object, in which the
##' parameters specified in \code{param} are freely estimated in all groups.
##' When \code{modelType = "mgcfa"}, only in the case of testing
##' \emph{configural} invariance should \code{uncon = NULL}. When
##' \code{modelType = "mimic"}, any non-\code{NULL uncon} is silently set to
##' \code{NULL}.
##' @param null Optional.  A \code{lavaan} object, in which an alternative null
##' model is fit (besides the default independence model specified by
##' \code{lavaan}) for the calculation of incremental fit indices. See Widamin &
##' Thompson (2003) for details. If \code{NULL}, \code{lavaan}'s default
##' independence model is used.
##' @param param An optional character vector or list of character vectors
##' indicating which parameters the user would test for DIF following a
##' rejection of the omnibus null hypothesis tested using
##' (\code{more})\code{AFIs}. Note that \code{param} does not guarantee certain
##' parameters \emph{are} constrained in \code{con}; that is for the user to
##' specify when fitting the model. If users have any "anchor items" that they
##' would never intend to free across groups (or levels of a covariate), these
##' should be excluded from \code{param}; exceptions to a type of parameter can
##' be specified in \code{freeParam}. When \code{modelType = "mgcfa"},
##' \code{param} indicates which parameters of interest are constrained across
##' groups in \code{con} and are unconstrained in \code{uncon}. Parameter names
##' must match those returned by \code{names(coef(con))}, but omitting any
##' group-specific suffixes (e.g., \code{"f1~1"} rather than \code{"f1~1.g2"})
##' or user-specified labels (that is, the parameter names must follow the rules
##' of lavaan's \code{\link[lavaan]{model.syntax}}). Alternatively (or
##' additionally), to test all constraints of a certain type (or multiple types)
##' of parameter in \code{con}, \code{param} may take any combination of the
##' following values: \code{"loadings"}, \code{"intercepts"},
##' \code{"thresholds"}, \code{"residuals"}, \code{"residual.covariances"},
##' \code{"means"}, \code{"lv.variances"}, and/or \code{"lv.covariances"}. When
##' \code{modelType = "mimic"}, \code{param} must be a vector of individual
##' parameters or a list of character strings to be passed one-at-a-time to
##' \code{\link[lavaan]{lavTestScore}}\code{(object = con, add = param[i])},
##' indicating which (sets of) regression paths fixed to zero in \code{con} that
##' the user would consider freeing (i.e., exclude anchor items). If
##' \code{modelType = "mimic"} and \code{param} is a list of character strings,
##' the multivariate test statistic will be saved for each list element instead
##' of 1-\emph{df} modification indices for each individual parameter, and
##' \code{names(param)} will name the rows of the \code{MI.obs} slot (see
##' \linkS4class{permuteMeasEq}). Set \code{param = NULL} (default) to avoid
##' collecting modification indices for any follow-up tests.
##' @param freeParam An optional character vector, silently ignored when
##' \code{modelType = "mimic"}. If \code{param} includes a type of parameter
##' (e.g., \code{"loadings"}), \code{freeParam} indicates exceptions (i.e.,
##' anchor items) that the user would \emph{not} intend to free across groups
##' and should therefore be ignored when calculating \emph{p} values adjusted
##' for the number of follow-up tests. Parameter types that are already
##' unconstrained across groups in the fitted \code{con} model (i.e., a
##' \emph{partial} invariance model) will automatically be ignored, so they do
##' not need to be specified in \code{freeParam}. Parameter names must match
##' those returned by \code{names(coef(con))}, but omitting any group-specific
##' suffixes (e.g., \code{"f1~1"} rather than \code{"f1~1.g2"}) or
##' user-specified labels (that is, the parameter names must follow the rules of
##' lavaan \code{\link[lavaan]{model.syntax}}).
##' @param covariates An optional character vector, only applicable when
##' \code{modelType = "mimic"}. The observed data are partitioned into columns
##' indicated by \code{covariates}, and the rows are permuted simultaneously for
##' the entire set before being merged with the remaining data.  Thus, the
##' covariance structure is preserved among the covariates, which is necessary
##' when (e.g.) multiple dummy codes are used to represent a discrete covariate
##' or when covariates interact. If \code{covariates = NULL} when
##' \code{modelType = "mimic"}, the value of \code{covariates} is inferred by
##' searching \code{param} for predictors (i.e., variables appearing after the
##' "\code{~}" operator).
##' @param AFIs A character vector indicating which alternative fit indices (or
##' chi-squared itself) are to be used to test the multiparameter omnibus null
##' hypothesis that the constraints specified in \code{con} hold in the
##' population. Any fit measures returned by \code{\link[lavaan]{fitMeasures}}
##' may be specified (including constants like \code{"df"}, which would be
##' nonsensical). If both \code{AFIs} and \code{moreAFIs} are \code{NULL}, only
##' \code{"chisq"} will be returned.
##' @param moreAFIs Optional. A character vector indicating which (if any)
##' alternative fit indices returned by \code{\link[semTools]{moreFitIndices}}
##' are to be used to test the multiparameter omnibus null hypothesis that the
##' constraints specified in \code{con} hold in the population.
##' @param maxSparse Only applicable when \code{modelType = "mgcfa"} and at
##' least one indicator is \code{ordered}. An integer indicating the maximum
##' number of consecutive times that randomly permuted group assignment can
##' yield a sample in which at least one category (of an \code{ordered}
##' indicator) is unobserved in at least one group, such that the same set of
##' parameters cannot be estimated in each group. If such a sample occurs, group
##' assignment is randomly permuted again, repeatedly until a sample is obtained
##' with all categories observed in all groups. If \code{maxSparse} is exceeded,
##' \code{NA} will be returned for that iteration of the permutation
##' distribution.
##' @param maxNonconv An integer indicating the maximum number of consecutive
##' times that a random permutation can yield a sample for which the model does
##' not converge on a solution. If such a sample occurs, permutation is
##' attempted repeatedly until a sample is obtained for which the model does
##' converge. If \code{maxNonconv} is exceeded, \code{NA} will be returned for
##' that iteration of the permutation distribution, and a warning will be
##' printed when using \code{show} or \code{summary}.
##' @param showProgress Logical. Indicating whether to display a progress bar
##' while permuting. Silently set to \code{FALSE} when using parallel options.
##' @param warn Sets the handling of warning messages when fitting model(s) to
##' permuted data sets. See \code{\link[base]{options}}.
##' @param datafun An optional function that can be applied to the data
##' (extracted from \code{con}) after each permutation, but before fitting the
##' model(s) to each permutation. The \code{datafun} function must have an
##' argument named \code{data} that accepts a \code{data.frame}, and it must
##' return a \code{data.frame} containing the same column names. The column
##' order may differ, the values of those columns may differ (so be careful!),
##' and any additional columns will be ignored when fitting the model, but an
##' error will result if any column names required by the model syntax do not
##' appear in the transformed data set. Although available for any
##' \code{modelType}, \code{datafun} may be useful when using the MIMIC method
##' to test for nonuniform DIF (metric/weak invariance) by using product
##' indicators for a latent factor representing the interaction between a factor
##' and one of the \code{covariates}, in which case the product indicators would
##' need to be recalculated after each permutation of the \code{covariates}. To
##' access other R objects used within \code{permuteMeasEq}, the arguments to
##' \code{datafun} may also contain any subset of the following: \code{"con"},
##' \code{"uncon"}, \code{"null"}, \code{"param"}, \code{"freeParam"},
##' \code{"covariates"}, \code{"AFIs"}, \code{"moreAFIs"}, \code{"maxSparse"},
##' \code{"maxNonconv"}, and/or \code{"iseed"}. The values for those arguments
##' will be the same as the values supplied to \code{permuteMeasEq}.
##' @param extra An optional function that can be applied to any (or all) of the
##' fitted lavaan objects (\code{con}, \code{uncon}, and/or \code{null}). This
##' function will also be applied after fitting the model(s) to each permuted
##' data set. To access the R objects used within \code{permuteMeasEq}, the
##' arguments to \code{extra} must be any subset of the following: \code{"con"},
##' \code{"uncon"}, \code{"null"}, \code{"param"}, \code{"freeParam"},
##' \code{"covariates"}, \code{"AFIs"}, \code{"moreAFIs"}, \code{"maxSparse"},
##' \code{"maxNonconv"}, and/or \code{"iseed"}. The values for those arguments
##' will be the same as the values supplied to \code{permuteMeasEq}. The
##' \code{extra} function must return a named \code{numeric} vector or a named
##' \code{list} of scalars (i.e., a \code{list} of \code{numeric} vectors of
##' \code{length == 1}). Any unnamed elements (e.g., \code{""} or \code{NULL})
##' of the returned object will result in an error.
##' @param parallelType The type of parallel operation to be used (if any). The
##' default is \code{"none"}. Forking is not possible on Windows, so if
##' \code{"multicore"} is requested on a Windows machine, the request will be
##' changed to \code{"snow"} with a message.
##' @param ncpus Integer: number of processes to be used in parallel operation.
##' If \code{NULL} (the default) and \code{parallelType %in%
##' c("multicore","snow")}, the default is one less than the maximum number of
##' processors detected by \code{\link[parallel]{detectCores}}. This default is
##' also silently set if the user specifies more than the number of processors
##' detected.
##' @param cl An optional \pkg{parallel} or \pkg{snow} cluster for use when
##' \code{parallelType = "snow"}.  If \code{NULL}, a \code{"PSOCK"} cluster on
##' the local machine is created for the duration of the \code{permuteMeasEq}
##' call. If a valid \code{\link[parallel]{makeCluster}} object is supplied,
##' \code{parallelType} is silently set to \code{"snow"}, and \code{ncpus} is
##' silently set to \code{length(cl)}.
##' @param iseed Integer: Only used to set the states of the RNG when using
##' parallel options, in which case \code{\link[base]{RNGkind}} is set to
##' \code{"L'Ecuyer-CMRG"} with a message. See
##' \code{\link[parallel]{clusterSetRNGStream}} and Section 6 of
##' \code{vignette("parallel", "parallel")} for more details. If user supplies
##' an invalid value, \code{iseed} is silently set to the default (12345). To
##' set the state of the RNG when not using parallel options, call
##' \code{\link[base]{set.seed}} before calling \code{permuteMeasEq}.
##'
##' @return The \linkS4class{permuteMeasEq} object representing the results of
##' testing measurement equivalence (the multiparameter omnibus test) and DIF
##' (modification indices), as well as diagnostics and any \code{extra} output.
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##' \email{TJorgensen314@@gmail.com})
##'
##' @seealso \code{\link[stats]{TukeyHSD}}, \code{\link[lavaan]{lavTestScore}},
##' \code{\link[semTools]{measurementInvariance}},
##' \code{\link[semTools]{measurementInvarianceCat}}
##'
##' @references
##'
##' \bold{Papers about permutation tests of measurement equivalence:}
##'
##' Jorgensen, T. D., Kite, B. A., Chen, P.-Y., & Short, S. D. (in press).
##' Permutation randomization methods for testing measurement equivalence and
##' detecting differential item functioning in multiple-group confirmatory
##' factor analysis. \emph{Psychological Methods}. doi:10.1037/met0000152
##'
##' Kite, B. A., Jorgensen, T. D., & Chen, P.-Y. (in press). Random permutation
##' testing applied to measurement invariance testing with ordered-categorical
##' indicators. \emph{Structural Equation Modeling}.
##' doi:10.1080/10705511.2017.1421467
##'
##' Jorgensen, T. D. (2017). Applying permutation tests and multivariate
##' modification indices to configurally invariant models that need
##' respecification. \emph{Frontiers in Psychology, 8}(1455).
##' doi:10.3389/fpsyg.2017.01455
##'
##' \bold{Additional reading:}
##'
##' Chen, F. F. (2007). Sensitivity of goodness of fit indexes to
##' lack of measurement invariance.  \emph{Structural Equation Modeling, 14}(3),
##' 464--504. doi:10.1080/10705510701301834
##'
##' Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit indexes
##' for testing measurement invariance. \emph{Structural Equation Modeling,
##' 9}(2), 233--255. doi:10.1207/S15328007SEM0902_5
##'
##' Meade, A. W., Johnson, E. C., & Braddy, P. W. (2008). Power and sensitivity
##' of alternative fit indices in tests of measurement invariance. \emph{Journal
##' of Applied Psychology, 93}(3), 568--592. doi:10.1037/0021-9010.93.3.568
##'
##' Widamin, K. F., & Thompson, J. S. (2003). On specifying the null model for
##' incremental fit indices in structural equation modeling. \emph{Psychological
##' Methods, 8}(1), 16--37. doi:10.1037/1082-989X.8.1.16
##' @examples
##'
##' \dontrun{
##'
##' ########################
##' ## Multiple-Group CFA ##
##' ########################
##'
##' ## create 3-group data in lavaan example(cfa) data
##' HS <- lavaan::HolzingerSwineford1939
##' HS$ageGroup <- ifelse(HS$ageyr < 13, "preteen",
##'                       ifelse(HS$ageyr > 13, "teen", "thirteen"))
##'
##' ## specify and fit an appropriate null model for incremental fit indices
##' mod.null <- c(paste0("x", 1:9, " ~ c(T", 1:9, ", T", 1:9, ", T", 1:9, ")*1"),
##'               paste0("x", 1:9, " ~~ c(L", 1:9, ", L", 1:9, ", L", 1:9, ")*x", 1:9))
##' fit.null <- cfa(mod.null, data = HS, group = "ageGroup")
##'
##' ## fit target model with varying levels of measurement equivalence
##' mod.config <- '
##' visual  =~ x1 + x2 + x3
##' textual =~ x4 + x5 + x6
##' speed   =~ x7 + x8 + x9
##' '
##' miout <- measurementInvariance(mod.config, data = HS, std.lv = TRUE,
##'                                group = "ageGroup")
##'
##' (fit.config <- miout[["fit.configural"]])
##' (fit.metric <- miout[["fit.loadings"]])
##' (fit.scalar <- miout[["fit.intercepts"]])
##'
##'
##' ####################### Permutation Method
##'
##' ## fit indices of interest for multiparameter omnibus test
##' myAFIs <- c("chisq","cfi","rmsea","mfi","aic")
##' moreAFIs <- c("gammaHat","adjGammaHat")
##'
##' ## Use only 20 permutations for a demo.  In practice,
##' ## use > 1000 to reduce sampling variability of estimated p values
##'
##' ## test configural invariance
##' set.seed(12345)
##' out.config <- permuteMeasEq(nPermute = 20, con = fit.config)
##' out.config
##'
##' ## test metric equivalence
##' set.seed(12345) # same permutations
##' out.metric <- permuteMeasEq(nPermute = 20, uncon = fit.config, con = fit.metric,
##'                             param = "loadings", AFIs = myAFIs,
##'                             moreAFIs = moreAFIs, null = fit.null)
##' summary(out.metric, nd = 4)
##'
##' ## test scalar equivalence
##' set.seed(12345) # same permutations
##' out.scalar <- permuteMeasEq(nPermute = 20, uncon = fit.metric, con = fit.scalar,
##'                             param = "intercepts", AFIs = myAFIs,
##'                             moreAFIs = moreAFIs, null = fit.null)
##' summary(out.scalar)
##'
##' ## Not much to see without significant DIF.
##' ## Try using an absurdly high alpha level for illustration.
##' outsum <- summary(out.scalar, alpha = .50)
##'
##' ## notice that the returned object is the table of DIF tests
##' outsum
##'
##' ## visualize permutation distribution
##' hist(out.config, AFI = "chisq")
##' hist(out.metric, AFI = "chisq", nd = 2, alpha = .01,
##'      legendArgs = list(x = "topright"))
##' hist(out.scalar, AFI = "cfi", printLegend = FALSE)
##'
##'
##' ####################### Extra Output
##'
##' ## function to calculate expected change of Group-2 and -3 latent means if
##' ## each intercept constraint were released
##' extra <- function(con) {
##'   output <- list()
##'   output["x1.vis2"] <- lavTestScore(con, release = 19:20, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[70]
##'   output["x1.vis3"] <- lavTestScore(con, release = 19:20, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[106]
##'   output["x2.vis2"] <- lavTestScore(con, release = 21:22, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[70]
##'   output["x2.vis3"] <- lavTestScore(con, release = 21:22, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[106]
##'   output["x3.vis2"] <- lavTestScore(con, release = 23:24, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[70]
##'   output["x3.vis3"] <- lavTestScore(con, release = 23:24, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[106]
##'   output["x4.txt2"] <- lavTestScore(con, release = 25:26, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[71]
##'   output["x4.txt3"] <- lavTestScore(con, release = 25:26, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[107]
##'   output["x5.txt2"] <- lavTestScore(con, release = 27:28, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[71]
##'   output["x5.txt3"] <- lavTestScore(con, release = 27:28, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[107]
##'   output["x6.txt2"] <- lavTestScore(con, release = 29:30, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[71]
##'   output["x6.txt3"] <- lavTestScore(con, release = 29:30, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[107]
##'   output["x7.spd2"] <- lavTestScore(con, release = 31:32, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[72]
##'   output["x7.spd3"] <- lavTestScore(con, release = 31:32, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[108]
##'   output["x8.spd2"] <- lavTestScore(con, release = 33:34, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[72]
##'   output["x8.spd3"] <- lavTestScore(con, release = 33:34, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[108]
##'   output["x9.spd2"] <- lavTestScore(con, release = 35:36, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[72]
##'   output["x9.spd3"] <- lavTestScore(con, release = 35:36, univariate = FALSE,
##'                                     epc = TRUE, warn = FALSE)$epc$epc[108]
##'   output
##' }
##'
##' ## observed EPC
##' extra(fit.scalar)
##'
##' ## permutation results, including extra output
##' set.seed(12345) # same permutations
##' out.scalar <- permuteMeasEq(nPermute = 20, uncon = fit.metric, con = fit.scalar,
##'                             param = "intercepts", AFIs = myAFIs,
##'                             moreAFIs = moreAFIs, null = fit.null, extra = extra)
##' ## summarize extra output
##' summary(out.scalar, extra = TRUE)
##'
##'
##' ###########
##' ## MIMIC ##
##' ###########
##'
##' ## Specify Restricted Factor Analysis (RFA) model, equivalent to MIMIC, but
##' ## the factor covaries with the covariate instead of being regressed on it.
##' ## The covariate defines a single-indicator construct, and the
##' ## double-mean-centered products of the indicators define a latent
##' ## interaction between the factor and the covariate.
##' mod.mimic <- '
##' visual  =~ x1 + x2 + x3
##' age =~ ageyr
##' age.by.vis =~ x1.ageyr + x2.ageyr + x3.ageyr
##'
##' x1 ~~ x1.ageyr
##' x2 ~~ x2.ageyr
##' x3 ~~ x3.ageyr
##' '
##'
##' HS.orth <- indProd(var1 = paste0("x", 1:3), var2 = "ageyr", match = FALSE,
##'                    data = HS[ , c("ageyr", paste0("x", 1:3))] )
##' fit.mimic <- cfa(mod.mimic, data = HS.orth, meanstructure = TRUE)
##' summary(fit.mimic, stand = TRUE)
##'
##' ## Whereas MIMIC models specify direct effects of the covariate on an indicator,
##' ## DIF can be tested in RFA models by specifying free loadings of an indicator
##' ## on the covariate's construct (uniform DIF, scalar invariance) and the
##' ## interaction construct (nonuniform DIF, metric invariance).
##' param <- as.list(paste0("age + age.by.vis =~ x", 1:3))
##' names(param) <- paste0("x", 1:3)
##' # param <- as.list(paste0("x", 1:3, " ~ age + age.by.vis")) # equivalent
##'
##' ## test both parameters simultaneously for each indicator
##' do.call(rbind, lapply(param, function(x) lavTestScore(fit.mimic, add = x)$test))
##' ## or test each parameter individually
##' lavTestScore(fit.mimic, add = as.character(param))
##'
##'
##' ####################### Permutation Method
##'
##' ## function to recalculate interaction terms after permuting the covariate
##' datafun <- function(data) {
##'   d <- data[, !names(data) %in% paste0("x", 1:3, ".ageyr")]
##'   indProd(var1 = paste0("x", 1:3), var2 = "ageyr", match = FALSE, data = d)
##' }
##'
##' set.seed(12345)
##' perm.mimic <- permuteMeasEq(nPermute = 20, modelType = "mimic",
##'                             con = fit.mimic, param = param,
##'                             covariates = "ageyr", datafun = datafun)
##' summary(perm.mimic)
##'
##' }
##'
##' @export
permuteMeasEq <- function(nPermute, modelType = c("mgcfa","mimic"),
                          con, uncon = NULL, null = NULL,
                          param = NULL, freeParam = NULL, covariates = NULL,
                          AFIs = NULL, moreAFIs = NULL,
                          maxSparse = 10, maxNonconv = 10, showProgress = TRUE,
                          warn = -1, datafun, extra,
                          parallelType = c("none","multicore","snow"),
                          ncpus = NULL, cl = NULL, iseed = 12345) {

  ## save arguments from call
  availableArgs <- as.list(formals(permuteMeasEq))
  argNames <- names(availableArgs)
  if (missing(datafun)) argNames <- setdiff(argNames, "datafun")
  if (missing(extra)) argNames <- setdiff(argNames, "extra")
  for (aa in argNames) {
    if (!is.null(eval(as.name(aa))))
      suppressWarnings(availableArgs[[aa]] <- eval(as.name(aa)))
  }
  ## check and return them
  fullCall <- do.call(checkPermArgs, availableArgs)
  ## assign them to workspace (also adds old_RNG & oldSeed to workspace)
  for (aa in names(fullCall)) assign(aa, fullCall[[aa]])

  ###################### SAVE OBSERVED RESULTS ##########################
  AFI.obs <- do.call(getAFIs, fullCall)
  ## save modification indices if !is.null(param)
  if (is.null(param)) {
    MI.obs <- data.frame(NULL)
  } else MI.obs <- do.call(getMIs, fullCall)

  ## anything extra?
  if (!missing(extra)) {
    extraArgs <- formals(extra)
    neededArgs <- intersect(names(extraArgs), names(fullCall))
    extraArgs <- do.call(c, lapply(neededArgs, function(nn) fullCall[nn]))
    extraOut <- do.call(extra, extraArgs)
    ## check that extra() returns a named list of scalars
    if (!is.list(extraOut)) extraOut <- as.list(extraOut)
    wrongFormat <- paste('Function "extra" must return a numeric vector or a',
                         'list of scalars, with each element named.')
    if (!all(sapply(extraOut, is.numeric))) stop(wrongFormat)
    if (!all(sapply(extraOut, length) == 1L)) stop(wrongFormat)
    if (is.null(names(extraOut)) | any(names(extraOut) == "")) stop(wrongFormat)
    extra.obs <- do.call(c, extraOut)
  } else extra.obs <- numeric(length = 0L)

  ######################### PREP DATA ##############################
  argList <- fullCall[c("con","uncon","null","param","freeParam","covariates",
                        "AFIs","moreAFIs","maxSparse","maxNonconv","warn","iseed")]
  argList$G <- lavInspect(con, "group")
    ## check for categorical variables
    # catVars <- lavaan::lavNames(con, type = "ov.ord")
    # numVars <- lavaan::lavNames(con, type = "ov.num")
    # latentVars <- lavaan::lavNames(con, type = "lv.regular")
  ## assemble data to which the models were fit
  if (length(argList$G)) {
    dataList <- mapply(FUN = function(x, g, n) {
      y <- data.frame(as.data.frame(x), g, stringsAsFactors = FALSE)
      names(y) <- c(n, argList$G)
      y
    }, SIMPLIFY = FALSE,
    x = lavInspect(con, "data"), g = lavInspect(con, "group.label"),
    n = lavaan::lavNames(con, type = "ov",
                         group = seq_along(lavInspect(con, "group.label"))))
    argList$d <- do.call(rbind, dataList)
  } else {
    argList$d <- as.data.frame(lavInspect(con, "data"))
    names(argList$d) <- lavaan::lavNames(con, type = "ov")
  }
  ## check that covariates are actual variables
  if (modelType == "mimic") {
    if (length(covariates) && !all(covariates %in% names(argList$d)))
      stop('These specified covariates are not columns in the data.frame:\n',
           paste(setdiff(covariates, names(argList$d)), collapse = ", "))
  }
  ## anything extra?
  if (!missing(extra)) argList$extra <- extra
  if (!missing(datafun)) argList$datafun <- datafun

  ###################### PERMUTED RESULTS ###########################
  ## permute and return distributions of (delta)AFIs, largest MI, and extras
  if (showProgress) {
    mypb <- utils::txtProgressBar(min = 1, max = nPermute, initial = 1,
                                  char = "=", width = 50, style = 3, file = "")
    permuDist <- list()
    for (j in 1:nPermute) {
      permuDist[[j]] <- do.call(paste("permuteOnce", modelType, sep = "."),
                                args = c(argList, i = j))
      utils::setTxtProgressBar(mypb, j)
    }
    close(mypb)
  } else if (parallelType == "multicore") {
    if (length(iseed)) set.seed(iseed)
    argList$FUN <- paste("permuteOnce", modelType, sep = ".")
    argList$X <- 1:nPermute
    argList$mc.cores <- ncpus
    argList$mc.set.seed <- TRUE
	pmcl <- function(...) { parallel::mclapply(...) }
    permuDist <- do.call(pmcl, args = argList)
    ## restore old RNG type
    if (fullCall$old_RNG[1] != "L'Ecuyer-CMRG") RNGkind(fullCall$old_RNG[1])
  } else if (parallelType == "snow") {
    stopTheCluster <- FALSE
    if (is.null(cl)) {
      stopTheCluster <- TRUE
      cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
    }
    parallel::clusterSetRNGStream(cl, iseed = iseed)
    argList$cl <- cl
    argList$X <- 1:nPermute
    argList$fun <- paste("permuteOnce", modelType, sep = ".")
    parallel::clusterExport(cl, varlist = c(argList$fun, "getAFIs","getMIs")) #FIXME: need update?
	tempppl <- function(...) { parallel::parLapply(...) }
    permuDist <- do.call(tempppl, args = argList)
    if (stopTheCluster) parallel::stopCluster(cl)
    ## restore old RNG type
    if (fullCall$old_RNG[1] != "L'Ecuyer-CMRG") RNGkind(fullCall$old_RNG[1])
  } else {
    argList$X <- 1:nPermute
    argList$FUN <- paste("permuteOnce", modelType, sep = ".")
    permuDist <- do.call(lapply, args = argList)
  }

  ## extract AFI distribution
  if (length(AFI.obs) > 1) {
    AFI.dist <- as.data.frame(t(sapply(permuDist, function(x) x$AFI)))
  }
  if (length(AFI.obs) == 1L) {
    AFI.dist <- data.frame(sapply(permuDist, function(x) x$AFI))
    colnames(AFI.dist) <- names(AFI.obs)
  }
  ## identify badness-of-fit measures
  badness <- grepl(pattern = "fmin|chi|aic|bic|rmr|rmsea|cn|sic|hqc",
                   x = names(AFI.obs), ignore.case = TRUE)
  ## calculate all one-directional p-values
  AFI.pval <- mapply(FUN = function(x, y, b) {
      if (b) return(mean(x >= y, na.rm = TRUE))
      mean(x <= y, na.rm = TRUE)
    }, x = unclass(AFI.dist), y = AFI.obs, b = badness)

  ## extract distribution of maximum modification indices
  MI.dist <- as.numeric(unlist(lapply(permuDist, function(x) x$MI)))
  ## calculate Tukey-adjusted p values for modification indices
  if (!is.null(param)) {
    MI.obs$tukey.p.value <- sapply(MI.obs$X2,
                                   function(i) mean(i <= MI.dist, na.rm = TRUE))
    MI.obs <- as.data.frame(unclass(MI.obs))
    rownames(MI.obs) <- names(param)
  }

  ## anything extra?
  if (!missing(extra)) {
    extra.dist <- do.call(rbind, lapply(permuDist, function(x) x$extra))
  } else extra.dist <- data.frame(NULL)

  ## save parameter table for show/summary methods
  PT <- as.data.frame(parTable(con))
  PT$par <- paste0(PT$lhs, PT$op, PT$rhs)
  if (length(lavInspect(con, "group")))
    PT$group.label[PT$group > 0] <- lavInspect(con, "group.label")[PT$group[PT$group > 0] ]

  ## return observed results, permutation p values, and ANOVA results
  if (is.null(uncon)) {
    delta <- lavaan::anova(con)
  } else {
    delta <- lavaan::anova(uncon, con)
  }
  ANOVA <- sapply(delta[,c("Chisq diff","Df diff","Pr(>Chisq)")], function(x) x[2])
  out <- new("permuteMeasEq", PT = PT, modelType = modelType, ANOVA = ANOVA,
             AFI.obs = AFI.obs, AFI.dist = AFI.dist, AFI.pval = AFI.pval,
             MI.obs = MI.obs, MI.dist = MI.dist,
             extra.obs = extra.obs, extra.dist = extra.dist,
             n.Permutations = nPermute, n.Converged = sum(!is.na(AFI.dist[,1])),
             n.nonConverged = sapply(permuDist, function(x) x$n.nonConverged),
             n.Sparse = sapply(permuDist, function(x) x$n.Sparse),
             oldSeed = fullCall$oldSeed)
  out
}



## ----------------
## Hidden Functions
## ----------------


## function to check validity of arguments to permuteMeasEq()
#' @importFrom lavaan lavInspect parTable
checkPermArgs <- function(nPermute, modelType, con, uncon, null,
                          param, freeParam, covariates, AFIs, moreAFIs,
                          maxSparse, maxNonconv, showProgress, warn,
                          datafun, extra, parallelType, ncpus, cl, iseed) {
  fixedCall <- as.list(match.call())[-1]

  fixedCall$nPermute <- as.integer(nPermute[1])
  fixedCall$modelType <- modelType[1]
  if (!fixedCall$modelType %in% c("mgcfa","mimic","long"))
    stop('modelType must be one of c("mgcfa","mimic","long")')
  if (fixedCall$modelType == "long") stop('modelType "long" is not yet available.')
  if (fixedCall$modelType == "mgcfa" && lavInspect(con, "ngroups") == 1L)
    stop('modelType = "mgcfa" applies only to multigroup models.')
  if (fixedCall$modelType == "mimic") {
    uncon <- NULL
    fixedCall$uncon <- NULL
    fixedCall <- c(fixedCall, list(uncon = NULL))
  }
  ## strip white space
  if (is.list(param)) {
    fixedCall$param <- lapply(param, function(cc) gsub("[[:space:]]+", "", cc))
  } else if (!is.null(param)) fixedCall$param <- gsub("[[:space:]]+", "", param)
  if (!is.null(freeParam)) fixedCall$freeParam <- gsub("[[:space:]]+", "", freeParam)
  if (fixedCall$modelType == "mimic") {
    # PT <- lavaan::lavaanify(fixedCall$param)
    # checkCovs <- unique(PT$rhs[PT$op == "~"])
    # if (is.null(covariates)) covariates <- checkCovs
    # if (length(setdiff(covariates, checkCovs)))
    #   warning('Argument "covariates" includes predictors not in argument "param"')
    ##### ordVars <- lavaan::lavNames(con, type = "ov.ord")
    fixedCall$covariates <- as.character(covariates)
  }
  fixedCall$maxSparse <- as.integer(maxSparse[1])
  fixedCall$maxNonconv <- as.integer(maxNonconv[1])
  fixedCall$showProgress <- as.logical(showProgress[1])
  fixedCall$warn <- as.integer(warn[1])
  fixedCall$oldSeed <- as.integer(NULL)
  parallelType <- as.character(parallelType[1])
  if (!parallelType %in% c("none","multicore","snow")) parallelType <- "none"
  if (!is.null(cl)) {
    if (!is(cl, "cluster")) stop("Invalid cluster object.  Check class(cl)")
    parallelType <- "snow"
    ncpus <- length(cl)
  }
  if (parallelType == "multicore" && .Platform$OS.type == "windows") {
    parallelType <- "snow"
    message("'multicore' option unavailable on Windows. Using 'snow' instead.")
  }
  ## parallel settings, adapted from boot::boot()
  if (parallelType != "none") {
    if (is.null(ncpus) || ncpus > parallel::detectCores()) {
      ncpus <- parallel::detectCores() - 1
    }
    if (ncpus <= 1L) {
      parallelType <- "none"
    } else {
      fixedCall$showProgress <- FALSE
      fixedCall$old_RNG <- RNGkind()
      fixedCall$oldSeed <- .Random.seed
      if (fixedCall$old_RNG[1] != "L'Ecuyer-CMRG") {
        RNGkind("L'Ecuyer-CMRG")
        message("Your RNGkind() was changed from ", fixedCall$old_RNG[1],
                " to L'Ecuyer-CMRG, which is required for reproducibility ",
                " in parallel jobs.  Your RNGkind() has been returned to ",
                fixedCall$old_RNG[1], " but the seed has not been set. ",
                " The state of your previous RNG is saved in the slot ",
                " named 'oldSeed', if you want to restore it using ",
                " the syntax:\n",
                ".Random.seed[-1] <- permuteMeasEqObject@oldSeed[-1]")
      }
      fixedCall$iseed <- as.integer(iseed[1])
      if (is.na(fixedCall$iseed)) fixedCall$iseed <- 12345
    }
  }
  fixedCall$parallelType <- parallelType
  if (is.null(ncpus)) {
    fixedCall$ncpus <- NULL
    fixedCall <- c(fixedCall, list(ncpus = NULL))
  } else fixedCall$ncpus <- ncpus

  ## check that "param" is NULL if uncon is NULL, and check for lavaan class
  notLavaan <- "Non-NULL 'con', 'uncon', or 'null' must be fitted lavaan object."
  if (is.null(uncon)) {
    if (!is.null(fixedCall$param) && fixedCall$modelType == "mgcfa") {
      message(c(" When 'uncon = NULL', only configural invariance is tested.",
                "\n So the 'param' argument was changed to NULL."))
      fixedCall$param <- NULL
      fixedCall <- c(fixedCall, list(param = NULL))
    }
    if (class(con) != "lavaan") stop(notLavaan)
  } else {
    if (class(con) != "lavaan") stop(notLavaan)
    if (class(uncon) != "lavaan") stop(notLavaan)
  }
  if (!is.null(null)) {
    if (class(null) != "lavaan") stop(notLavaan)
  }

  ############ FIXME: check that lavInspect(con, "options")$conditional.x = FALSE (find defaults for continuous/ordered indicators)
  if (!is.null(fixedCall$param)) {
    ## Temporarily warn about testing thresholds without necessary constraints.   FIXME: check for binary indicators
    if ("thresholds" %in% fixedCall$param | any(grepl("\\|", fixedCall$param))) {
      warning(c("This function is not yet optimized for testing thresholds.\n",
                "Necessary identification contraints might not be specified."))
    }
    ## collect parameter types for "mgcfa"
    if (fixedCall$modelType != "mimic") {
      ## save all estimates from constrained model
      PT <- parTable(con)[ , c("lhs","op","rhs","group","plabel")]
      ## extract parameters of interest
      paramTypes <- c("loadings","intercepts","thresholds","residuals","means",
                      "residual.covariances","lv.variances","lv.covariances")
      params <- PT[paste0(PT$lhs, PT$op, PT$rhs) %in% setdiff(fixedCall$param,
                                                              paramTypes), ]
      ## add parameters by type, if any are specified
      types <- intersect(fixedCall$param, paramTypes)
      ov.names <- lavaan::lavNames(con, "ov")
      isOV <- PT$lhs %in% ov.names
      lv.names <- con@pta$vnames$lv[[1]]
      isLV <- PT$lhs %in% lv.names & PT$rhs %in% lv.names
      if ("loadings" %in% types) params <- rbind(params, PT[PT$op == "=~", ])
      if ("intercepts" %in% types) {
        params <- rbind(params, PT[isOV & PT$op == "~1", ])
      }
      if ("thresholds" %in% types) params <- rbind(params, PT[PT$op == "|", ])
      if ("residuals" %in% types) {
        params <- rbind(params, PT[isOV & PT$lhs == PT$rhs & PT$op == "~~", ])
      }
      if ("residual.covariances" %in% types) {
        params <- rbind(params, PT[isOV & PT$lhs != PT$rhs & PT$op == "~~", ])
      }
      if ("means" %in% types) {
        params <- rbind(params, PT[PT$lhs %in% lv.names & PT$op == "~1", ])
      }
      if ("lv.variances" %in% types) {
        params <- rbind(params, PT[isLV & PT$lhs == PT$rhs & PT$op == "~~", ])
      }
      if ("lv.covariances" %in% types) {
        params <- rbind(params, PT[isLV & PT$lhs != PT$rhs & PT$op == "~~", ])
      }
      ## remove parameters specified by "freeParam" argument
      params <- params[!paste0(params$lhs, params$op, params$rhs) %in% fixedCall$freeParam, ]
      fixedCall$param <- paste0(params$lhs, params$op, params$rhs)
    }
  }


  if (is.null(AFIs) & is.null(moreAFIs)) {
    message("No AFIs were selected, so only chi-squared will be permuted.\n")
    fixedCall$AFIs <- "chisq"
    AFIs <- "chisq"
  }
  if ("ecvi" %in% AFIs & lavInspect(con, "ngroups") > 1L)
    stop("ECVI is not available for multigroup models.")

  ## check estimators
  leastSq <- grepl("LS", lavInspect(con, "options")$estimator)
  if (!is.null(uncon)) {
    if (uncon@Options$estimator != lavInspect(con, "options")$estimator)
      stop("Models must be fit using same estimator.")
  }
  if (!is.null(null)) {
    if (lavInspect(null, "options")$estimator != lavInspect(con, "options")$estimator)
      stop("Models must be fit using same estimator.")
  }

  ## check extra functions, if any
  restrictedArgs <- c("con","uncon","null","param","freeParam","covariates",
                      "AFIs","moreAFIs","maxSparse","maxNonconv","iseed")
  if (!missing(datafun)) {
    if (!is.function(datafun)) stop('Argument "datafun" must be a function.')
    extraArgs <- formals(datafun)
    if (!all(names(extraArgs) %in% c(restrictedArgs, "data")))
      stop('The user-supplied function "datafun" can only have any among the ',
           'following arguments:\n', paste(restrictedArgs, collapse = ", "))
  }
  if (!missing(extra)) {
    if (!is.function(extra)) stop('Argument "extra" must be a function.')
    extraArgs <- formals(extra)
    if (!all(names(extraArgs) %in% restrictedArgs))
      stop('The user-supplied function "extra" can only have any among the ',
           'following arguments:\n', paste(restrictedArgs, collapse = ", "))
  }

  ## return evaluated list of other arguments
  lapply(fixedCall, eval)
}


## function to extract fit measures
#' @importFrom lavaan lavInspect
getAFIs <- function(...) {
  dots <- list(...)

  AFI1 <- list()
  AFI0 <- list()
  leastSq <- grepl("LS", lavInspect(dots$con, "options")$estimator)
  ## check validity of user-specified AFIs, save output
  if (!is.null(dots$AFIs)) {
    IC <- grep("ic|logl", dots$AFIs, value = TRUE)
    if (leastSq & length(IC)) {
      stop(paste("Argument 'AFIs' includes invalid options:",
                 paste(IC, collapse = ", "),
                 "Information criteria unavailable for least-squares estimators.",
                 sep = "\n"))
    }
    if (!is.null(dots$uncon))
      AFI1[[1]] <- lavaan::fitMeasures(dots$uncon, fit.measures = dots$AFIs,
                                       baseline.model = dots$null)
    AFI0[[1]] <- lavaan::fitMeasures(dots$con, fit.measures = dots$AFIs,
                                     baseline.model = dots$null)
  }
  ## check validity of user-specified moreAFIs
  if (!is.null(dots$moreAFIs)) {
    IC <- grep("ic|hqc", dots$moreAFIs, value = TRUE)
    if (leastSq & length(IC)) {
      stop(paste("Argument 'moreAFIs' includes invalid options:",
                 paste(IC, collapse = ", "),
                 "Information criteria unavailable for least-squares estimators.",
                 sep = "\n"))
    }
    if (!is.null(dots$uncon))
      AFI1[[2]] <- moreFitIndices(dots$uncon, fit.measures = dots$moreAFIs)
    AFI0[[2]] <- moreFitIndices(dots$con, fit.measures = dots$moreAFIs)
  }

  ## save observed AFIs or delta-AFIs
  if (is.null(dots$uncon)) {
    AFI.obs <- unlist(AFI0)
  } else {
    AFI.obs <- unlist(AFI0) - unlist(AFI1)
  }
  AFI.obs
}

## Function to extract modification indices for equality constraints
#' @importFrom lavaan parTable
getMIs <- function(...) {
  dots <- list(...)

  if (dots$modelType == "mgcfa") {
    ## save all estimates from constrained model
    PT <- parTable(dots$con)[ , c("lhs","op","rhs","group","plabel")]
    ## extract parameters of interest
    params <- PT[paste0(PT$lhs, PT$op, PT$rhs) %in% dots$param, ]
    ## return modification indices for specified constraints (param)
    MIs <- lavaan::lavTestScore(dots$con)$uni
    MI.obs <- MIs[MIs$lhs %in% params$plabel, ]
  } else if (dots$modelType == "mimic") {
    if (is.list(dots$param)) {
      MI <- lapply(dots$param, function(x) lavaan::lavTestScore(dots$con, add = x)$test)
      MI.obs <- do.call(rbind, MI)
    } else MI.obs <- lavaan::lavTestScore(dots$con, add = dots$param)$uni
  } else if (dots$modelType == "long") {
    ## coming soon
  }

  MI.obs
}

## Functions to find delta-AFIs & maximum modification index in one permutation
#' @importFrom lavaan lavInspect
permuteOnce.mgcfa <- function(i, d, G, con, uncon, null, param, freeParam,
                              covariates, AFIs, moreAFIs, maxSparse, maxNonconv,
                              iseed, warn, extra = NULL, datafun = NULL) {
  old_warn <- options()$warn
  options(warn = warn)
  ## save arguments from call
  argNames <- names(formals(permuteOnce.mgcfa))
  availableArgs <- lapply(argNames, function(x) eval(as.name(x)))
  names(availableArgs) <- argNames

  nSparse <- 0L
  nTries <- 1L
  while ( (nSparse <= maxSparse) & (nTries <= maxNonconv) ) {
    ## permute grouping variable
    d[ , G] <- sample(d[ , G])
    ## transform data?
    if (!is.null(datafun)) {
      extraArgs <- formals(datafun)
      neededArgs <- intersect(names(extraArgs), names(availableArgs))
      extraArgs <- do.call(c, lapply(neededArgs, function(nn) availableArgs[nn]))
      extraArgs$data <- d
      originalNames <- colnames(d)
      d <- do.call(datafun, extraArgs)
      ## coerce extraOut to data.frame
      if (!is.data.frame(d)) stop('Argument "datafun" did not return a data.frame')
      if (!all(originalNames %in% colnames(d)))
        stop('The data.frame returned by argument "datafun" did not contain ',
             'column names required by the model:\n',
             paste(setdiff(originalNames, colnames(d)), collapse = ", "))
    }

    ## for ordered indicators, check that groups have same observed categories
    ordVars <- lavaan::lavNames(con, type = "ov.ord")
    if (length(ordVars) > 0) {
      try(onewayTables <- lavaan::lavTables(d, dimension = 1L,
                                            categorical = ordVars, group = G),
          silent = TRUE)
      if (exists("onewayTables")) {
        if (any(onewayTables$obs.prop == 1)) {
          nSparse <- nSparse + 1L
          next
        }
      } else {
        ## no "onewayTables" probably indicates empty categories in 1+ groups
        nSparse <- nSparse + 1L
        next
      }
    }
    ## fit null model, if it exists
    if (!is.null(null)) {
      out.null <- lavaan::update(null, data = d, group.label = lavInspect(con, "group.label"))
    }

    ## fit constrained model, check for convergence
    try(out0 <- lavaan::update(con, data = d, group.label = lavInspect(con, "group.label")))
    if (!exists("out0")) {
      nTries <- nTries + 1L
      next
    }
    if (!lavInspect(out0, "converged")) {
      nTries <- nTries + 1L
      next
    }

    ## fit unconstrained model (unless NULL), check for convergence
    if (!is.null(uncon)) {
      try(out1 <- lavaan::update(uncon, data = d, group.label = lavInspect(con, "group.label")))
      if (!exists("out1")) {
        nTries <- nTries + 1L
        next
      }
      if (!lavInspect(out1, "converged")) {
        nTries <- nTries + 1L
        next
      }

    }
    ## If you get this far, everything converged, so break WHILE loop
    break
  }
  ## if WHILE loop ended before getting results, return NA
  if ( (nSparse == maxSparse) | (nTries == maxNonconv) ) {
    allAFIs <- c(AFIs, moreAFIs)
    AFI <- rep(NA, sum(!is.na(allAFIs)))
    names(AFI) <- allAFIs[!is.na(allAFIs)]
    MI <- if (is.null(param)) NULL else NA
    extra.obs <- NA
    nTries <- nTries + 1L
  } else {
    availableArgs$con <- out0
    if (exists("out1")) availableArgs$uncon <- out1
    if (exists("out.null")) availableArgs$null <- out.null
    AFI <- do.call(getAFIs, availableArgs)
    ## save max(MI) if !is.null(param)
    if (is.null(param)) {
      MI <- NULL
    } else {
      MI <- max(do.call(getMIs, c(availableArgs, modelType = "mgcfa"))$X2)
    }
    ## anything extra?
    if (!is.null(extra)) {
      extraArgs <- formals(extra)
      neededArgs <- intersect(names(extraArgs), names(availableArgs))
      extraArgs <- do.call(c, lapply(neededArgs, function(nn) availableArgs[nn]))
      extraOut <- do.call(extra, extraArgs)
      ## coerce extraOut to data.frame
      if (!is.list(extraOut)) extraOut <- as.list(extraOut)
      extra.obs <- data.frame(extraOut)
    } else extra.obs <- data.frame(NULL)
  }
  options(warn = old_warn)
  list(AFI = AFI, MI = MI, extra = extra.obs,
       n.nonConverged = nTries - 1L, n.Sparse = nSparse)
}

#' @importFrom lavaan lavInspect
permuteOnce.mimic <- function(i, d, G, con, uncon, null, param, freeParam,
                              covariates, AFIs, moreAFIs, maxSparse, maxNonconv,
                              iseed, warn, extra = NULL, datafun = NULL) {
  old_warn <- options()$warn
  options(warn = warn)
  ## save arguments from call
  argNames <- names(formals(permuteOnce.mimic))
  availableArgs <- lapply(argNames, function(x) eval(as.name(x)))
  names(availableArgs) <- argNames

  nTries <- 1L
  while (nTries <= maxNonconv) {
    ## permute covariate(s) within each group
    if (length(G)) {
      for (gg in lavInspect(con, "group.label")) {
        dG <- d[ d[[G]] == gg, ]
        N <- nrow(dG)
        newd <- dG[sample(1:N, N), covariates, drop = FALSE]
        for (COV in covariates) d[d[[G]] == gg, COV] <- newd[ , COV]
      }
    } else {
      N <- nrow(d)
      newd <- d[sample(1:N, N), covariates, drop = FALSE]
      for (COV in covariates) d[ , COV] <- newd[ , COV]
    }
    ## transform data?
    if (!is.null(datafun)) {
      extraArgs <- formals(datafun)
      neededArgs <- intersect(names(extraArgs), names(availableArgs))
      extraArgs <- do.call(c, lapply(neededArgs, function(nn) availableArgs[nn]))
      extraArgs$data <- d
      originalNames <- colnames(d)
      d <- do.call(datafun, extraArgs)
      ## coerce extraOut to data.frame
      if (!is.data.frame(d)) stop('Argument "datafun" did not return a data.frame')
      if (!all(originalNames %in% colnames(d)))
        stop('The data.frame returned by argument "datafun" did not contain ',
             'column names required by the model:\n',
             paste(setdiff(originalNames, colnames(d)), collapse = ", "))
    }


    ## fit null model, if it exists
    if (!is.null(null)) {
      out.null <- lavaan::update(null, data = d, group.label = lavInspect(con, "group.label"))
    }

    ## fit constrained model
    try(out0 <- lavaan::update(con, data = d, group.label = lavInspect(con, "group.label")))
    ## check for convergence
    if (!exists("out0")) {
      nTries <- nTries + 1L
      next
    }
    if (!lavInspect(out0, "converged")) {
      nTries <- nTries + 1L
      next
    }
    ## If you get this far, everything converged, so break WHILE loop
    break
  }
  ## if WHILE loop ended before getting results, return NA
  if (nTries == maxNonconv) {
    allAFIs <- c(AFIs, moreAFIs)
    AFI <- rep(NA, sum(!is.na(allAFIs)))
    names(AFI) <- allAFIs[!is.na(allAFIs)]
    MI <- if (is.null(param)) NULL else NA
    extra.obs <- NA
    nTries <- nTries + 1L
  } else {
    availableArgs$con <- out0
    if (exists("out.null")) availableArgs$null <- out.null
    AFI <- do.call(getAFIs, availableArgs)
    if (is.null(param)) {
      MI <- NULL
    } else {
      MI <- max(do.call(getMIs, c(availableArgs, modelType = "mimic"))$X2)
    }
    ## anything extra?
    if (!is.null(extra)) {
      extraArgs <- formals(extra)
      neededArgs <- intersect(names(extraArgs), names(availableArgs))
      extraArgs <- do.call(c, lapply(neededArgs, function(nn) availableArgs[nn]))
      extraOut <- do.call(extra, extraArgs)
      ## coerce extraOut to data.frame
      if (!is.list(extraOut)) extraOut <- as.list(extraOut)
      extra.obs <- data.frame(extraOut)
    } else extra.obs <- data.frame(NULL)
  }
  options(warn = old_warn)
  list(AFI = AFI, MI = MI, extra = extra.obs,
       n.nonConverged = nTries - 1L, n.Sparse = integer(length = 0))
}



