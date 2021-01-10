### Terrence D. Jorgensen
### Last updated: 10 January 2021
### Class and Methods for lavaan.mi object, returned by runMI()


##' Class for a lavaan Model Fitted to Multiple Imputations
##'
##' This class extends the \code{\linkS4class{lavaanList}} class, created by
##' fitting a lavaan model to a list of data sets. In this case, the list of
##' data sets are multiple imputations of missing data.
##'
##'
##' @name lavaan.mi-class
##' @importClassesFrom lavaan lavaanList
##' @aliases lavaan.mi-class show,lavaan.mi-method summary,lavaan.mi-method
##'   fitMeasures,lavaan.mi-method fitmeasures,lavaan.mi-method
##'   anova,lavaan.mi-method nobs,lavaan.mi-method coef,lavaan.mi-method
##'   vcov,lavaan.mi-method fitted,lavaan.mi-method fitted.values,lavaan.mi-method
##'   residuals,lavaan.mi-method resid,lavaan.mi-method
##' @docType class
##'
##' @slot coefList \code{list} of estimated coefficients in matrix format (one
##'   per imputation) as output by \code{\link[lavaan]{lavInspect}(fit, "est")}
##' @slot phiList \code{list} of model-implied latent-variable covariance
##'   matrices (one per imputation) as output by
##'   \code{\link[lavaan]{lavInspect}(fit, "cov.lv")}
##' @slot miList \code{list} of modification indices output by
##'   \code{\link[lavaan]{modindices}}
##' @slot seed \code{integer} seed set before running imputations
##' @slot lavListCall call to \code{\link[lavaan]{lavaanList}} used to fit the
##'   model to the list of imputed data sets in \code{@@DataList}, stored as a
##'   \code{list} of arguments
##' @slot imputeCall call to imputation function (if used), stored as a
##'   \code{list} of arguments
##' @slot convergence \code{list} of \code{logical} vectors indicating whether,
##'   for each imputed data set, (1) the model converged on a solution, (2)
##'   \emph{SE}s could be calculated, (3) the (residual) covariance matrix of
##'   latent variables (\eqn{\Psi}) is non-positive-definite, and (4) the
##'   residual covariance matrix of observed variables (\eqn{\Theta}) is
##'   non-positive-definite.
##' @slot lavaanList_slots All remaining slots are from
##'   \code{\linkS4class{lavaanList}}, but \code{\link{runMI}} only populates a
##'   subset of the \code{list} slots, two of them with custom information:
##' @slot DataList The \code{list} of imputed data sets
##' @slot SampleStatsList List of output from
##'   \code{\link[lavaan]{lavInspect}(fit, "sampstat")} applied to each fitted
##'   model
##' @slot ParTableList See \code{\linkS4class{lavaanList}}
##' @slot vcovList See \code{\linkS4class{lavaanList}}
##' @slot testList See \code{\linkS4class{lavaanList}}
##' @slot h1List See \code{\linkS4class{lavaanList}}. An additional element is
##'   added to the \code{list}: \code{$PT} is the "saturated" model's parameter
##'   table, returned by \code{\link[lavaan]{lav_partable_unrestricted}}.
##' @slot baselineList See \code{\linkS4class{lavaanList}}
##'
##' @param object An object of class \code{lavaan.mi}
##' @param se,ci,level,standardized,rsquare,header,output See
##'        \code{\link[lavaan]{parameterEstimates}}. \code{output}
##'        can also be passed to \code{\link[lavaan]{fitMeasures}}.
##' @param fmi \code{logical} indicating whether to include the Fraction Missing
##'        Information (FMI) for parameter estimates in the \code{summary}
##'        output (see \bold{Value} section).
##' @param asymptotic \code{logical}. If \code{FALSE} (typically a default, but
##'        see \bold{Value} section for details using various methods), pooled
##'        tests (of fit or pooled estimates) will be \emph{F} or \emph{t}
##'        statistics with associated degrees of freedom (\emph{df}). If
##'        \code{TRUE}, the (denominator) \emph{df} are assumed to be
##'        sufficiently large for a \emph{t} statistic to follow a normal
##'        distribution, so it is printed as a \emph{z} statisic; likewise,
##'        \emph{F} times its numerator \emph{df} is printed, assumed to follow
##'        a \eqn{\chi^2} distribution.
##' @param scale.W \code{logical}. If \code{TRUE} (default), the \code{vcov}
##'        method will calculate the pooled covariance matrix by scaling the
##'        within-imputation component by the ARIV (see Enders, 2010, p. 235,
##'        for definition and formula). Otherwise, the pooled matrix is
##'        calculated as the weighted sum of the within-imputation and
##'        between-imputation components (see Enders, 2010, ch. 8, for details).
##'        This in turn affects how the \code{summary} method calcualtes its
##'        pooled standard errors, as well as the Wald test
##'        (\code{\link{lavTestWald.mi}}).
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
##'        imputations warrant further investigation. Specific imputation
##'        numbers can also be included in this argument, in case users want to
##'        apply their own custom omission criteria (or simulations can use
##'        different numbers of imputations without redundantly refitting the
##'        model).
##' @param labels \code{logical} indicating whether the \code{coef} output
##'        should include parameter labels. Default is \code{TRUE}.
##' @param total \code{logical} (default: \code{TRUE}) indicating whether the
##'        \code{nobs} method should return the total sample size or (if
##'        \code{FALSE}) a vector of group sample sizes.
##' @param type The meaning of this argument varies depending on which method it
##'        it used for. Find detailed descriptions in the \bold{Value} section
##'        under \code{coef}, \code{vcov}, and \code{residuals}.
##' @param fit.measures,baseline.model See \code{\link[lavaan]{fitMeasures}}.
##'        \code{summary(object, fit.measures = TRUE)} will print (but not
##'        return) a table of fit measures to the console.
##' @param ... Additional arguments passed to \code{\link{lavTestLRT.mi}}, or
##'        subsequently to \code{\link[lavaan]{lavTestLRT}}.
##'
##' @return
##'
##' \item{coef}{\code{signature(object = "lavaan.mi", type = "free",
##'   labels = TRUE, omit.imps = c("no.conv","no.se"))}:
##'   See \code{\linkS4class{lavaan}}. Returns the pooled point estimates (i.e.,
##'   averaged across imputed data sets; see Rubin, 1987).}
##'
##' \item{vcov}{\code{signature(object = "lavaan.mi", scale.W = TRUE,
##'   omit.imps = c("no.conv","no.se"),
##'   type = c("pooled","between","within","ariv"))}:  By default, returns the
##'   pooled covariance matrix of parameter estimates (\code{type = "pooled"}),
##'   the within-imputations covariance matrix (\code{type = "within"}), the
##'   between-imputations covariance matrix (\code{type = "between"}), or the
##'   average relative increase in variance (\code{type = "ariv"}) due to
##'   missing data.}
##'
##' \item{fitted.values}{\code{signature(object = "lavaan.mi",
##'   omit.imps = c("no.conv","no.se"))}: See \code{\linkS4class{lavaan}}.
##'   Returns model-implied moments, evaluated at the pooled point estimates.}
##' \item{fitted}{alias for \code{fitted.values}}
##'
##' \item{residuals}{\code{signature(object = "lavaan.mi",
##'   type = c("raw","cor"), omit.imps = c("no.conv","no.se"))}:
##'   See \code{\linkS4class{lavaan}}. By default (\code{type = "raw"}), returns
##'   the difference between the model-implied moments from \code{fitted.values}
##'   and the pooled observed moments (i.e., averaged across imputed data sets).
##'   Standardized residuals are also available, using Bollen's
##'   (\code{type = "cor"} or \code{"cor.bollen"}) or Bentler's
##'   (\code{type = "cor.bentler"}) formulas.}
##' \item{resid}{alias for \code{residuals}}
##'
##' \item{nobs}{\code{signature(object = "lavaan.mi", total = TRUE)}: either
##'   the total (default) sample size or a vector of group sample sizes
##'   (\code{total = FALSE}).}
##'
##' \item{anova}{\code{signature(object = "lavaan.mi", ...)}:
##'   Returns a test of model fit for a single model (\code{object}) or test(s)
##'   of the difference(s) in fit between nested models passed via \code{...}.
##'   See \code{\link{lavTestLRT.mi}} and \code{\link{compareFit}} for details.}
##'
##' \item{fitMeasures}{\code{signature(object = "lavaan.mi",
##'   fit.measures = "all", baseline.model = NULL, output = "vector",
##'   omit.imps = c("no.conv","no.se"), ...)}: See lavaan's
##'   \code{\link[lavaan]{fitMeasures}} for details. Pass additional arguments
##'   to \code{\link{lavTestLRT.mi}} via \code{...}.}
##' \item{fitmeasures}{alias for \code{fitMeasures}.}
##'
##' \item{show}{\code{signature(object = "lavaan.mi")}: returns a message about
##'  convergence rates and estimation problems (if applicable) across imputed
##'  data sets.}
##'
##' \item{summary}{\code{signature(object = "lavaan.mi", se = TRUE, ci = FALSE,
##'  level = .95, standardized = FALSE, rsquare = FALSE, fmi = FALSE,
##'  scale.W = !asymptotic, omit.imps = c("no.conv","no.se"), asymptotic = FALSE,
##'   header = TRUE, output = "text", fit.measures = FALSE, ...)}:
##'  see \code{\link[lavaan]{parameterEstimates}} for details.
##'  By default, \code{summary} returns pooled point and \emph{SE}
##'  estimates, along with \emph{t} test statistics and their associated
##'  \emph{df} and \emph{p} values. If \code{ci = TRUE}, confidence intervales
##'  are returned with the specified confidence \code{level} (default 95\% CI).
##'  If \code{asymptotic = TRUE}, \emph{z} instead of \emph{t} tests are
##'  returned. \code{standardized} solution(s) can also be requested by name
##'  (\code{"std.lv"} or \code{"std.all"}) or both are returned with \code{TRUE}.
##'  \emph{R}-squared for endogenous variables can be requested, as well as the
##'  Fraction Missing Information (FMI) for parameter estimates. By default, the
##'  output will appear like \code{lavaan}'s \code{summary} output, but if
##'  \code{output == "data.frame"}, the returned \code{data.frame} will resemble
##'  the \code{parameterEstimates} output. The \code{scale.W} argument is
##'  passed to \code{vcov} (see description above).
##'  Setting \code{fit.measures=TRUE} will additionally print fit measures to
##'  the console, but they will not be returned; additional arguments may be
##'  passed via \code{...} to \code{\link[lavaan]{fitMeasures}} and
##'  subsequently to \code{\link{lavTestLRT.mi}}.}
##'
##' @section Objects from the Class: See the \code{\link{runMI}} function for
##'   details. Wrapper functions include \code{\link{lavaan.mi}},
##'   \code{\link{cfa.mi}}, \code{\link{sem.mi}}, and \code{\link{growth.mi}}.
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Asparouhov, T., & Muthen, B. (2010). \emph{Chi-square statistics
##'   with multiple imputation}. Technical Report. Retrieved from
##'   \url{http://www.statmodel.com/}
##'
##'   Enders, C. K. (2010). \emph{Applied missing data analysis}. New York, NY:
##'   Guilford.
##'
##'   Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
##'   Significance levels from repeated \emph{p}-values with multiply-imputed
##'   data. \emph{Statistica Sinica, 1}(1), 65--92. Retrieved from
##'   \url{https://www.jstor.org/stable/24303994}
##'
##'   Meng, X.-L., & Rubin, D. B. (1992). Performing likelihood ratio tests with
##'   multiply-imputed data sets. \emph{Biometrika, 79}(1), 103--111.
##'   \doi{10.2307/2337151}
##'
##'   Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}.
##'   New York, NY: Wiley.
##'
##' @examples
##'
##' ## See ?runMI help page
##'
setClass("lavaan.mi", contains = "lavaanList",
         slots = c(coefList = "list",     # coefficients in matrix format
                   phiList = "list",      # list of model-implied latent covariance matrices
                   miList = "list",       # modification indices
                   seed = "integer",      # seed set before running imputations
                   lavListCall = "list",  # store actual call to lavaanList
                   imputeCall = "list",   # store call from imputation, if used
                   convergence = "list")) # also check SEs and Heywood cases



##' @name lavaan.mi-class
##' @aliases show,lavaan.mi-method
##' @export
setMethod("show", "lavaan.mi", function(object) {
  nData <- object@meta$ndat

  useImps <- sapply(object@convergence, "[[", i = "converged")
  nConverged <- sum(useImps)

  SE <- sapply(object@convergence, "[[", "SE")
  SE[is.na(SE)] <- FALSE

  Heywood.ov <- sapply(object@convergence, "[[", "Heywood.ov")
  Heywood.ov[is.na(Heywood.ov)] <- FALSE

  Heywood.lv <- sapply(object@convergence, "[[", "Heywood.lv")
  Heywood.lv[is.na(Heywood.lv)] <- FALSE

  cat('lavaan.mi object based on ', nData, ' imputed data sets. \n',
      'See class?lavaan.mi help page for available methods. \n\n',
      'Convergence information:\n', 'The model converged on ',
      nConverged, ' imputed data sets \n\n', sep = "")

  if (!all(SE)) cat('Standard errors could not be computed for data set(s)',
                    paste(which(!SE), collapse = ", "), '\nTry fitting the',
                    'model to the individual data set(s) to diagnose',
                    'problems. If they cannot be fixed, try inspecting the',
                    'imputations. It may be necessary to reimpute the data',
                    'with some restrictions imposed. \n\n')

  if (any(Heywood.ov | Heywood.lv))
    cat('Heywood cases detected for data set(s)',
        paste(which(Heywood.ov | Heywood.lv), collapse = ", "),
        '\nThese are not necessarily a cause for concern, unless a pooled',
        'estimate is also a Heywood case. \n\n')

  object
})


##' @importFrom stats pt qt pnorm qnorm
##' @importFrom lavaan lavListInspect parTable lavNames
##' @importFrom methods getMethod
summary.lavaan.mi <- function(object, se = TRUE, ci = FALSE, level = .95,
                              standardized = FALSE, rsquare = FALSE,
                              fmi = FALSE, scale.W = !asymptotic,
                              omit.imps = c("no.conv","no.se"),
                              asymptotic = FALSE, header = TRUE,
                              output = "text", fit.measures = FALSE, ...) {
  useImps <- rep(TRUE, length(object@DataList))
  if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
  if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
  if ("no.npd" %in% omit.imps) {
    Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
    Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
    useImps <- useImps & !(Heywood.lv | Heywood.ov)
  }
  ## custom removal by imputation number
  rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
  if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
  ## whatever is left
  m <- sum(useImps)
  if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
  useImps <- which(useImps)

  lavoptions <- lavListInspect(object, "options")

  ## extract parameter table with attributes for printing
  PT <- parTable(object)
  myCols <- c("lhs","op","rhs","exo")
  if (lavListInspect(object, "ngroups") > 1L) myCols <- c(myCols,"block","group")
  if (lavListInspect(object, "nlevels") > 1L) myCols <- c(myCols,"block","level")
  PE <- PT[ , unique(myCols)]
  free <- PT$free > 0L | PT$op == ":="
  STDs <- !(PT$op %in% c("==","<",">")) # which rows can be standardized

  PE$est <- getMethod("coef","lavaan.mi")(object, type = "all",
                                          omit.imps = omit.imps)

  if (lavoptions$se == "none") {
    warning('pooled variances and tests unavailable when se="none" is requested')
    se <- FALSE
  }
  if (!se) fmi <- FALSE
  messPool <- paste0("Rubin's (1987) rules were used to pool point",
                     if (se) " and SE",
                     " estimates across ", m, " imputed data sets",
                     if (se & !asymptotic) ", and to calculate degrees of",
                     if (se & !asymptotic) " freedom for each parameter's t",
                     if (se & !asymptotic) " test and CI.",
                     "\n")
  if (se) {
    VCOV <- getMethod("vcov","lavaan.mi")(object, scale.W = scale.W,
                                          omit.imps = omit.imps)
    PE$se <- lavaan::lav_model_vcov_se(object@Model, VCOV = VCOV,
                                       lavpartable = object@ParTable)
    W <- rowMeans(sapply(object@ParTableList[useImps], "[[", i = "se")^2)
    B <- apply(sapply(object@ParTableList[useImps], "[[", i = "est"), 1, var)
    Bm <- B + B/m
    Tot <- W + Bm
    if (asymptotic) {
      PE$z[free] <- PE$est[free] / PE$se[free]
      PE$pvalue <- pnorm(-abs(PE$z))*2
      crit <- qnorm(1 - (1 - level) / 2)
    } else {
      PE$t[free] <- PE$est[free] / PE$se[free]
      ## calculate df for t test
      ## can't do finite-sample correction because Wald z tests have no df
      ## (see Enders, 2010, p. 231, eq. 8.13 & 8.14)
      PE$df[free] <- (m - 1) * (1 + W[free] / Bm[free])^2
      ## if DF are obscenely large, set them to infinity for pretty printing
      PE$df <- ifelse(PE$df > 9999, Inf, PE$df)
      PE$pvalue <- pt(-abs(PE$t), df = PE$df)*2
      crit <- qt(1 - (1 - level) / 2, df = PE$df)
    }
    if (ci) {
      PE$ci.lower <- PE$est - crit * PE$se
      PE$ci.upper <- PE$est + crit * PE$se
      PE$ci.lower[!free] <- PE$ci.upper[!free] <- PE$est[!free]
    }
  }

  if (is.logical(standardized)) {
    if (standardized) {
      standardized <- c("std.lv","std.all")
      if (length(lavNames(object, "ov.x")) && lavoptions$fixed.x) {
        standardized <- c(standardized, "std.nox")
      }
    } else standardized <- NULL
  } else standardized <- tolower(as.character(standardized))

  if (length(standardized) || rsquare) {
    ## pooled estimates for standardizedSolution()
    est <- getMethod("coef", "lavaan.mi")(object, omit.imps = omit.imps)
    ## updates @Model@GLIST for standardizedSolution(..., GLIST=)
    object@Model <- lavaan::lav_model_set_parameters(object@Model, x = est)
  }

  if ("std.lv" %in% standardized) {
    PE$std.lv[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                    type = "std.lv",
                                                    GLIST = object@Model@GLIST,
                                                    est = PE$est)$est.std
  }
  if ("std.all" %in% standardized) {
    PE$std.all[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                     type = "std.all",
                                                     GLIST = object@Model@GLIST,
                                                     est = PE$est)$est.std
  }
  if ("std.nox" %in% standardized) {
    PE$std.nox[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                     type = "std.nox",
                                                     GLIST = object@Model@GLIST,
                                                     est = PE$est)$est.std
  }

  if (fmi) {
    PE$fmi[free] <- Bm[free] / Tot[free]
    PE$riv[free] <- Bm[free] / W[free] # (Enders, 2010, p. 226, eq. 8.10)
    # == PE$riv[free] <- PE$fmi1[free] / (1 - PE$fmi1[free])
    messRIV <- paste("The RIV will exceed 1 whenever between-imputation",
                     "variance exceeds within-imputation variance",
                     "(when FMI(1) > 50%).\n\n")
  }
  ## fancy or not?
  if (output == "text") {
    PE$label <- PT$label
    #FIXME: no longer needed?  PE$exo <- 0L
    class(PE) <- c("lavaan.parameterEstimates","lavaan.data.frame","data.frame")
    attr(PE, "information") <- lavoptions$information[1]
    attr(PE, "information.meat") <- lavoptions$information.meat
    attr(PE, "se") <- lavoptions$se
    attr(PE, "group.label") <- lavListInspect(object, "group.label")
    attr(PE, "level.label") <- c("within", lavListInspect(object, "cluster"))
    attr(PE, "bootstrap") <- lavoptions$bootstrap
    attr(PE, "bootstrap.successful") <- 0L #FIXME: assumes none. Implement Wei & Fan's mixing method?
    attr(PE, "missing") <- lavoptions$missing
    attr(PE, "observed.information") <- lavoptions$observed.information[1]
    attr(PE, "h1.information") <- lavoptions$h1.information[1]
    attr(PE, "h1.information.meat") <- lavoptions$h1.information.meat
    attr(PE, "header") <- header
    # FIXME: lavaan may add more!!
    if (fmi) cat("\n", messRIV, sep = "")
  } else {
    PE$exo <- NULL
    class(PE) <- c("lavaan.data.frame","data.frame")
  }
  ## requested R-squared?
  endoNames <- c(lavNames(object, "ov.nox"), lavNames(object, "lv.nox"))
  if (rsquare & length(endoNames)) {
    isEndo <- sapply(PE$lhs, function(x) x %in% endoNames)
    rsqPE <- PE[PE$lhs == PE$rhs & PE$op == "~~" & isEndo, ]
    rsqPE$op <- "r2"
    for (i in which(!sapply(colnames(PE),
                            function(x) x %in% c("lhs","op","rhs","block",
                                                 "level","group","est","exo")))) {
      rsqPE[ , i] <- NA
    }
    STD <- lavaan::standardizedSolution(object, se = FALSE, type = "std.all",
                                        GLIST = object@Model@GLIST, est = PE$est)
    isEndoSTD <- sapply(STD$lhs, function(x) x %in% endoNames)
    std.all <- STD$est.std[STD$lhs == STD$rhs & STD$op == "~~" & isEndoSTD]
    rsqPE$est <- ifelse(std.all < 0, NA, 1 - std.all) # negative variances
    if (output == "text") rsqPE$label <- ""
    PE <- rbind(PE, rsqPE)
  }

  if (output == "data.frame") PE <- PE[!(PE$op %in% c("==","<",">")), ]
  rownames(PE) <- NULL

  if (output == "text") {
    getMethod("show", "lavaan.mi")(object)
    cat(messPool)
  }
  if (fit.measures) {
    indices <- c("chisq","df","pvalue","cfi","tli","rmsea","srmr")
    FITS <- suppressWarnings(fitMeasures(object, fit.measures = indices,
                                         output = "text", ...))
    try(print(FITS, add.h0 = TRUE), silent = TRUE)
  }

  PE
}
##' @name lavaan.mi-class
##' @aliases summary,lavaan.mi-method
##' @export
setMethod("summary", "lavaan.mi", summary.lavaan.mi)


##' @name lavaan.mi-class
##' @aliases nobs,lavaan.mi-method
##' @importFrom lavaan lavListInspect
##' @export
setMethod("nobs", "lavaan.mi", function(object, total = TRUE) {
  if (total) return(lavListInspect(object, "ntotal"))
  #FIXME: cluster N for multilevel?
  N <- lavListInspect(object, "norig")
  if (length(N) > 1L) names(N) <- lavListInspect(object, "group.label")
  N
})



##' @importFrom lavaan parTable
coef.lavaan.mi <- function(object, type = "free", labels = TRUE,
                           omit.imps = c("no.conv","no.se")) {
  useImps <- rep(TRUE, length(object@DataList))
  if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
  if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
  if ("no.npd" %in% omit.imps) {
    Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
    Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
    useImps <- useImps & !(Heywood.lv | Heywood.ov)
  }
  ## custom removal by imputation number
  rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
  if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
  ## whatever is left
  m <- sum(useImps)
  if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
  useImps <- which(useImps)

  PT <- parTable(object)
  if (type == "user" || type == "all") {
    type <- "user"
    idx <- 1:length(PT$lhs)
  } else if (type == "free") {
    ## FIXME: duplicated leftover from old way of handling EQ constraints?
    idx <- which(PT$free > 0L & !duplicated(PT$free))
  }
  ## extract coefficients for converged models
  coefList <- lapply(object@ParTableList[useImps], "[[", i = "est")
  out <- colMeans(do.call(rbind, coefList))[idx]
  ## attach names, set class
  if (labels) names(out) <- lavaan::lav_partable_labels(PT, type = type)
  class(out) <- c("lavaan.vector","numeric")
  out
}
##' @name lavaan.mi-class
##' @aliases coef,lavaan.mi-method
##' @export
setMethod("coef", "lavaan.mi", coef.lavaan.mi)



##' @importFrom stats cov
##' @importFrom lavaan lavListInspect parTable
vcov.lavaan.mi <- function(object, type = c("pooled","between","within","ariv"),
                           scale.W = TRUE, omit.imps = c("no.conv","no.se")) {
  useImps <- rep(TRUE, length(object@DataList))
  if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
  if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
  if ("no.npd" %in% omit.imps) {
    Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
    Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
    useImps <- useImps & !(Heywood.lv | Heywood.ov)
  }
  ## custom removal by imputation number
  rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
  if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
  ## whatever is left
  m <- sum(useImps)
  if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
  useImps <- which(useImps)

  if (lavListInspect(object, "options")$se == "none") {
    warning('requested se="none", so only between-imputation (co)variance can',
            ' be computed')
    type <- "between"
  }
  type <- tolower(type[1])
  if (!(type %in% c("pooled","between","within","ariv")))
    stop("'", type, "' is not a valid option for 'type'")

  PT <- parTable(object)
  ncon <- sum(PT$op == "==")
  npar <- max(PT$free) - ncon

  coefList <- lapply(object@ParTableList[useImps], "[[", i = "est")
  B <- cov(do.call(rbind, coefList)[ , PT$free > 0L & !duplicated(PT$free)])
  class(B) <- c("lavaan.matrix.symmetric","matrix")
  rownames(B) <- colnames(B) <- lavaan::lav_partable_labels(PT, type = "free")
  if (type == "between") return(B)

  W <- Reduce("+", lapply(object@vcovList[useImps], function(x) x$vcov)) / m
  class(W) <- c("lavaan.matrix.symmetric","matrix")
  dimnames(W) <- dimnames(B)
  if (type == "within") return(W)

  ## check whether equality constraints prevent inversion of W
  if (scale.W || type == "ariv") {
    inv.W <- if (ncon == 0) try(solve(W), silent = TRUE) else MASS::ginv(W)
    if (inherits(inv.W, "try-error")) {
      if (ncon == 0) {
        warning("Could not invert within-imputation covariance matrix. ",
                "Generalized inverse used instead.\nIt may be ",
                "safer to set `scale.W = FALSE' (and `asymptotic = TRUE').")
      }
      inv.W <- MASS::ginv(W)
    }
    ## relative increase in variance due to missing data
    r <- (1 + 1/m)/npar * sum(diag(B %*% inv.W)) # Enders (2010, p. 235) eqs. 8.20-21
    if (type == "ariv") return(r)
    Total <- (1 + r) * W # FIXME: asked Yves for a hack, says it can't be inverted back to infoMat
  } else {
    ## less reliable, but constraints prevent inversion of W
    Total <- W + B + (1/m)*B ## Enders (2010, p. 235) eq. 8.19
  }
  ## return pooled variance
  Total
}
##' @name lavaan.mi-class
##' @aliases vcov,lavaan.mi-method
##' @export
setMethod("vcov", "lavaan.mi", vcov.lavaan.mi)


##' @importFrom lavaan lavListInspect lavTestLRT
anova.lavaan.mi <- function(object, ...) {
  ## save model names
  objname <- deparse(substitute(object))
  dotnames <- as.character(sapply(substitute(list(...))[-1], deparse))

  ## check class
  if (!inherits(object, "lavaan.mi")) stop("object is not class 'lavaan.mi'")
  ## check for additional arguments
  dots <- list(...)
  if (length(dots)) {
    ## separate lavaan.mi objects from other lavTestLRT.mi() arguments
    idx.mi <- which(sapply(dots, inherits, what = "lavaan.mi"))
    if (length(idx.mi)) {
      mods <- dots[idx.mi]
      dots <- dots[-idx.mi]
      ## save names for mods, so compareFit() doesn't break
      modnames <- dotnames[idx.mi]
      nonames <- which(names(mods) == "")
      names(mods)[nonames] <- modnames[nonames]
    } else {
      mods <- NULL
      modnames <- NULL
    }
    LRT.names <- intersect(names(dots),
                           union(names(formals(lavTestLRT)),
                                 names(formals(lavTestLRT.mi))))
    dots <- if (length(LRT.names)) dots[LRT.names] else NULL
    if (!is.null(dots$h1)) {
      #FIXME: this shouldn't be necessary: mods <- c(mods, list(h1 = dots$h1))
      dots$h1 <- NULL
    }
  } else mods <- NULL

  ## run compareFit if length(idx.mi) > 1L
  if (length(mods) == 0L) {
    argList <- c(list(object = object), dots)
    results <- do.call(lavTestLRT.mi, argList)
  } else if (length(mods) == 1L) {
    argList <- c(list(object = object, h1 = mods[[1]]), dots)
    results <- do.call(lavTestLRT.mi, argList)
  } else if (length(mods) > 1L) {
    modList <- c(list(object), mods)
    names(modList) <- c(objname, modnames)
    argList <- c(modList, list(argsLRT = dots, indices = FALSE))
    results <- do.call(compareFit, argList)@nested
    class(results) <- c("lavaan.data.frame","data.frame")
    attr(results, "header") <- "Nested Model Comparisons:"
  }

  results
}
##' @name lavaan.mi-class
##' @aliases anova,lavaan.mi-method
##' @export
setMethod("anova", "lavaan.mi", anova.lavaan.mi)


##' @importFrom lavaan lavListInspect lavNames
##' @importFrom methods getMethod
## utility function called within fitMeasures.mi()
getSRMR <- function(object, type = "cor.bentler", level = "within",
                    include.mean = TRUE, omit.imps = c("no.conv","no.se")) {
  conditional.x <- lavListInspect(object, "options")$conditional.x
  include.mean <- include.mean && lavListInspect(object, "meanstructure")
  include.diag <- type %in% c("cor.bentler","raw")
  mplus <- type == "mplus"
  if (mplus) type <- "cor.bollen"

  ## how many blocks to loop over
  nG <- lavListInspect(object, "ngroups")
  nlevels <- lavListInspect(object, "nlevels")
  ## save relevant sample sizes
  if (nlevels > 1L && level != "within") {
    n.per.group <- lavListInspect(object, "nclusters") #FIXME: only works for 2 levels
    N <- sum(n.per.group)
  } else {
    n.per.group <- lavListInspect(object, "nobs")
    N <- lavListInspect(object, "ntotal")
  }

  ## grab residuals
  R <- getMethod("resid", "lavaan.mi")(object, type = type,
                                       omit.imps = omit.imps)
  if (mplus) Rd <- getMethod("resid", "lavaan.mi")(object, type = "cor.bentler",
                                                   omit.imps = omit.imps)
  ## restructure, if necessary
  if (nG == 1L) {
    loopBlocks <- 1L

    ## extract relevant level
    if (nlevels > 1L) {
      R <- R[[level]]
      if (mplus) Rd <- Rd[[level]]
    }
    ## to loop over blocks
    R <- list(R)
    if (mplus) Rd <- list(Rd)


  ## multiple groups AND multilevel
  } else if (nlevels > 1L) {
    loopBlocks <- 2*(1:nG)
    if (level == "within") loopBlocks <- loopBlocks - 1L
    R <- R[loopBlocks]
    if (mplus) Rd <- Rd[loopBlocks]

  } else loopBlocks <- 1:nG # no restructure necessary for multigroup 1-level models


  ## store vector of squared residuals
  RR <- vector("list", nG)
  for (b in loopBlocks) {
    index <- if (conditional.x) "res.cov" else "cov"

    RR[[b]] <- R[[b]][[index]][lower.tri(R[[b]][[index]], diag = FALSE)]^2
    ## only capture means/variances of numeric modeled variables (not conditional.x)
    vv <- intersect(lavNames(object, type = "ov.num", block = b),
                    lavNames(object, type = "ov.model", block = b))
    if (include.diag)  RR[[b]] <- c(RR[[b]], diag(R[[b]][[index]])[vv]^2)
    if (mplus)  RR[[b]] <- c(RR[[b]], diag(Rd[[b]][[index]])[vv]^2)

    if (include.mean) {
      index <- if (conditional.x) "res.int" else "mean"
      RR[[b]] <- c(RR[[b]], R[[b]][[index]][vv]^2)
    }
  }

  ## take weighted average of group means
  as.numeric( (n.per.group %*% sqrt(sapply(RR, mean))) / N )
}
##' @importFrom lavaan lavNames lavListInspect
##' @importFrom stats pchisq uniroot
fitMeasures.mi <- function(object, fit.measures = "all", baseline.model = NULL,
                           output = "vector", omit.imps = c("no.conv","no.se"),
                           ...) {

  useImps <- rep(TRUE, length(object@DataList))
  if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
  if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
  if ("no.npd" %in% omit.imps) {
    Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
    Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
    useImps <- useImps & !(Heywood.lv | Heywood.ov)
  }
  ## custom removal by imputation number
  rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
  if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
  ## whatever is left
  m <- sum(useImps)
  if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
  useImps <- which(useImps)

  lavoptions <- lavListInspect(object, "options")

  fit.measures <- tolower(fit.measures)
  if (length(fit.measures) == 0L) fit.measures <- "all"
  ## narrow down fit indices
  incremental <- c("cfi","tli","nnfi","rfi","nfi","pnfi","ifi","rni")
  if ("all" %in% fit.measures) {
    indices <- c("chisq","df","pvalue","scaling", incremental,
                 "rmsea","rmr","mfi","gammahat")
  } else {
    indices <- grep(pattern = paste(c("chisq","df","pvalue","scaling",
                                      incremental, "mfi","rmsea",
                                      "gammahat","rmr"), collapse = "|"),
                    x = fit.measures, ignore.case = TRUE, value = TRUE)
  }

  ## CHI-SQUARED-BASED FIT INDICES
  notest <- length(lavoptions$test) == 1L && lavoptions$test == "none"
  if (notest || any(!grepl(pattern = "rmr", x = indices))) {

    ## check for additional arguments
    dots <- list(...)
    if (length(dots)) {
      LRT.names <- intersect(names(dots),
                             union(names(formals(lavTestLRT)),
                                   names(formals(lavTestLRT.mi))))
      dots <- if (length(LRT.names)) dots[LRT.names] else list(asymptotic = TRUE)
    } else dots <- list(asymptotic = TRUE)

    ## check test options (adapted from lavTestLRT.mi, limits duplicate warnings)
    test <- dots$test
    if (is.null(test)) {
      test <- "d3" # default
    } else test <- tolower(test[1])
    if (tolower(test) %in% c("mr","meng.rubin","likelihood","lrt","mplus","d3")) test <- "D3"
    if (tolower(test) %in% c("lmrr","li.et.al","pooled.wald","d2")) test <- "D2"
    if (test == "D3" && !lavoptions$estimator %in% c("ML","PML","FML")) {
      message('"D3" only available using maximum likelihood estimation. ',
              'Changed test to "D2".')
      test <- "D2"
    }

    ## check for robust
    test.names <- lavoptions$test
    # lavaan 0.6-5: for now, we only acknowledge the first non-standard @test
    if (length(test.names) > 1L) {
      ## remove standard and any bootstrapped tests
      rm.idx <- which(test.names %in% c("standard","bootstrap","bollen.stine"))
      if (length(rm.idx) > 0L) {
        test.names <- test.names[-rm.idx]
      }
      ## only acknowledge the first scaled test statistic
      if (length(test.names) > 1L) {
        test.names <- test.names[1]
      }
    }

    robust <- any(test.names %in% c("satorra.bentler","yuan.bentler",
                                    "yuan.bentler.mplus","scaled.shifted",
                                    "mean.var.adjusted","satterthwaite"))
    if (robust) {
      ## assign pool.robust option to object
      if (is.null(dots$pool.robust)) {
        pool.robust <- formals(lavTestLRT.mi)$pool.robust # default value
      } else {
        pool.robust <- dots$pool.robust # user-specified value
      }
    } else dots$pool.robust <- pool.robust <- FALSE

    scaleshift <- any(test.names == "scaled.shifted")
    if (scaleshift) {
      if (test == "D3") {
        message("If test = 'scaled.shifted' (estimator = 'WLSMV' or 'MLMV'), ",
                "model evaluation is only available by (re)setting .",
                "test = 'D2'.\nControl more options by passing arguments to ",
                "lavTestLRT() via the '...' argument.\n")
        test <- 'D2'
      }
    }

    if (pool.robust && test == "D3") {
      message('pool.robust = TRUE is only applicable when test = "D2". ',
              'Changed test to "D2".')
      test <- "D2"
    }

    dots$test <- test


    ## pooled test statistic(s)
    argList <- c(list(object = object), dots)
    argList$asymptotic <- TRUE # in case it wasn't set in list(...)
    argList$omit.imps <- omit.imps
    out <- do.call(lavTestLRT.mi, argList)
    ## check for scaled test statistic (if not, set robust=FALSE)
    if (robust && is.na(out["chisq.scaled"])) robust <- FALSE

    ## fit baseline model if necessary
    if (any(indices %in% incremental)) {
      if (inherits(baseline.model, "lavaan.mi")) {
        baseFit <- baseline.model
      } else if (inherits(object@external$baseline.model, "lavaan.mi")) {
        baseFit <- object@external$baseline.model

        ## MUST fit PTb for "D3" likelihoods, but for "D2" use @baselineList
      } else if (test == "D2") {
        ## length(baseImps) == m, not just length(useImps)
        baseImps <- object@meta$baseline.ok
        if (!all(baseImps[useImps])) warning('The default independence model ',
                                             'did not converge for data set(s): ',
                                             which(!baseImps))
        ## only use imputations that converged for both
        baseImps <- intersect(useImps, which(baseImps))

        w <- sapply(object@baselineList[baseImps],
                    function(x) x$test$standard[["stat"]])
        if (is.list(w)) {
          #TODO: figure out why this happens!
          w <- unlist(w)
          DF <- mean(unlist(sapply(object@baselineList[baseImps],
                                   function(x) x$test$standard[["df"]])))
        } else {
          DF <- mean(sapply(object@baselineList[baseImps],
                            function(x) x$test$standard[["df"]]))
        }
        baseOut <- calculate.D2(w, DF, asymptotic = TRUE)
        if (robust) {
          if (pool.robust) {
            w.r <- sapply(object@baselineList[baseImps],
                          function(x) x$test[[ test.names[1] ]][["stat"]])
            if (is.list(w.r)) {
              w.r <- unlist(w.r)
              DF.r <- mean(unlist(sapply(object@baselineList[baseImps],
                                         function(x) x$test[[ test.names[1] ]][["df"]])))
            } else {
              DF.r <- mean(sapply(object@baselineList[baseImps],
                                  function(x) x$test[[ test.names[1] ]][["df"]]))
            }
            base.robust <- calculate.D2(w.r, DF.r, asymptotic = TRUE)
            names(base.robust) <- paste0(names(base.robust), ".scaled")
            baseOut <- c(baseOut, base.robust)
          } else {
            baseOut <- robustify(ChiSq = baseOut, object = object,
                                 baseline = TRUE, useImps = baseImps)
          }
        }
        baseFit <- NULL # for later checking, to avoid unnecessary calls

      } else {
        PTb <- object@baselineList[[ useImps[1] ]]$partable
        PTb[c("est","se")] <- NULL
        # FIXME: shouldn't need this line, but lav_partable_merge() fails when
        #        lavaan:::lav_object_extended() returns a NULL slot instead of "plabel"
        PTb$plabel <- paste0(".p", PTb$id, ".")
        group <- lavListInspect(object, "group")
        if (length(group) == 0L) group <- NULL
        cluster <- lavListInspect(object, "cluster")
        if (length(cluster) == 0L) cluster <- NULL
        baseFit <- runMI(model = PTb, data = object@DataList[useImps],
                         group = group, cluster = cluster,
                         test = lavoptions$test, estimator = lavoptions$estimator,
                         fixed.x = lavoptions$fixed.x, se = "none", # to save time
                         conditional.x = lavoptions$conditional.x,
                         ordered = lavListInspect(object, "ordered"),
                         parameterization = lavoptions$parameterization)
      }

      if (!is.null(baseFit)) {
        ## length(baseImps) is only as long as length(useImps), not the original m
        baseImps <- sapply(baseFit@convergence, "[[", i = "converged")
        if (!all(baseImps)) warning('baseline.model did not converge for data set(s): ',
                                    useImps[!baseImps])
        argList <- c(list(object = baseFit), dots)
        argList$asymptotic <- TRUE # in case it wasn't set in list(...)
        argList$omit.imps <- setdiff(omit.imps, "no.se") # se="none" in baseFit
        baseOut <- do.call(lavTestLRT.mi, argList)
      }
      # else { already used "D2" with @baselineList info to make baseOut }

    }


    X2 <- out[["chisq"]]
    DF <- out[["df"]]
    if (robust) {
      X2.sc <- out[["chisq.scaled"]]
      DF.sc <- out[["df.scaled"]] ## for mean.var.adjusted, mean DF across imputations
      if (!pool.robust) ch <- out[["chisq.scaling.factor"]] ## mean c_hat across imputations
      if (X2 < .Machine$double.eps && DF == 0) ch <- 0
      ## for RMSEA
      if ("rmsea" %in% indices) {
        d <- mean(sapply(object@testList[useImps],
                         function(x) sum(x[[ test.names[1] ]][["trace.UGamma"]])))
        if (is.na(d) || d == 0) d <- NA # FIXME: only relevant when mean.var.adjusted?
      }
    }

    ## for CFI, TLI, etc.
    if (any(indices %in% incremental)) {
      bX2 <- baseOut[["chisq"]]
      bDF <- baseOut[["df"]]
      out <- c(out, baseline.chisq = bX2, baseline.df = bDF,
               baseline.pvalue = baseOut[["pvalue"]])
      if (robust) {
        out["baseline.chisq.scaled"] <- bX2.sc <- baseOut[["chisq.scaled"]]
        out["baseline.df.scaled"]    <- bDF.sc <- baseOut[["df.scaled"]]
        out["baseline.pvalue.scaled"] <- baseOut[["pvalue.scaled"]]
        if (!pool.robust) {
          cb <- baseOut[["chisq.scaling.factor"]]
          out["baseline.chisq.scaling.factor"] <- cb
          if (scaleshift) out["baseline.chisq.shift.parameters"] <- baseOut[["chisq.shift.parameters"]]
        }
      }
    }

    if ("cfi" %in% indices) {
      t1 <- max(X2 - DF, 0)
      t2 <- max(X2 - DF, bX2 - bDF, 0)
      out["cfi"] <- if(t1 == 0 && t2 == 0) 1 else 1 - t1/t2
      if (robust) {
        ## scaled
        t1 <- max(X2.sc - DF.sc, 0)
        t2 <- max(X2.sc - DF.sc, bX2.sc - bDF.sc, 0)
        if (is.na(t1) || is.na(t2)) {
          out["cfi.scaled"] <- NA
        } else if (t1 == 0 && t2 == 0) {
          out["cfi.scaled"] <- 1
        } else out["cfi.scaled"] <- 1 - t1/t2
        ## Brosseau-Liard & Savalei MBR 2014, equation 15
        if (!pool.robust & test.names[1] %in%
            c("satorra.bentler","yuan.bentler","yuan.bentler.mplus")) {
          t1 <- max(X2 - ch*DF, 0)
          t2 <- max(X2 - ch*DF, bX2 - cb*bDF, 0)
          if (is.na(t1) || is.na(t2)) {
            out["cfi.robust"] <- NA
          } else if (t1 == 0 && t2 == 0) {
            out["cfi.robust"] <- 1
          } else out["cfi.robust"] <- 1 - t1/t2
        }
      }
    }
    if ("rni" %in% indices) {
      t1 <- X2 - DF
      t2 <- bX2 - bDF
      out["rni"] <- if (t2 == 0) NA else 1 - t1/t2
      if (robust) {
        ## scaled
        t1 <- X2.sc - DF.sc
        t2 <- bX2.sc - bDF.sc
        if (is.na(t1) || is.na(t2)) {
          out["rni.scaled"] <- NA
        } else if (t2 == 0) {
          out["rni.scaled"] <- NA
        } else out["rni.scaled"] <- 1 - t1/t2
        ## Brosseau-Liard & Savalei MBR 2014, equation 15
        if (!pool.robust & test.names[1] %in%
            c("satorra.bentler","yuan.bentler","yuan.bentler.mplus")) {
          t1 <- X2 - ch*DF
          t2 <- bX2 - cb*bDF
          if (is.na(t1) || is.na(t2)) {
            out["rni.robust"] <- NA
          } else if (t1 == 0 && t2 == 0) {
            out["rni.robust"] <- NA
          } else out["rni.robust"] <- 1 - t1/t2
        }
      }
    }
    if (any(indices %in% c("tli","nnfi"))) {
      t1 <- (X2 - DF)*bDF
      t2 <- (bX2 - bDF)*DF
      out["tli"] <- out["nnfi"] <- if (DF > 0) 1 - t1/t2 else 1
      if (robust) {
        ## scaled
        t1 <- (X2.sc - DF.sc)*bDF.sc
        t2 <- (bX2.sc - bDF.sc)*DF.sc
        if (is.na(t1) || is.na(t2)) {
          out["tli.scaled"] <- out["nnfi.scaled"] <- NA
        } else if (DF > 0 && t2 != 0) {
          out["tli.scaled"] <- out["nnfi.scaled"] <- 1 - t1/t2
        } else {
          out["tli.scaled"] <- out["nnfi.scaled"] <- 1
        }
        ## Brosseau-Liard & Savalei MBR 2014, equation 15
        if (!pool.robust & test.names[1] %in%
            c("satorra.bentler","yuan.bentler","yuan.bentler.mplus")) {
          t1 <- (X2 - ch*DF)*bDF
          t2 <- (bX2 - cb*bDF)*DF
          if (is.na(t1) || is.na(t2)) {
            out["tli.robust"] <- out["nnfi.robust"] <- NA
          } else if (t1 == 0 && t2 == 0) {
            out["tli.robust"] <- out["nnfi.robust"] <- 1 - t1/t2
          } else out["tli.robust"] <- out["nnfi.robust"] <- 1
        }
      }
    }
    if ("rfi" %in% indices) {
      if (DF > 0) {
        t2 <- bX2 / bDF
        t1 <- t2 - X2/DF
        out["rfi"] <- if (t1 < 0 || t2 < 0) 1 else t1/t2
      } else out["rfi"] <- 1
      if (robust) {
        if (DF > 0) {
          t2 <- bX2.sc / bDF.sc
          t1 <- t2 - X2.sc/DF.sc
          out["rfi.scaled"] <- if (t1 < 0 || t2 < 0) 1 else t1/t2
        } else out["rfi.scaled"] <- 1
      }
    }
    if ("nfi" %in% indices) {
      if (DF > 0) {
        t1 <- bX2 - X2
        t2 <- bX2
        out["nfi"] <- t1 / t2
      } else out["nfi"] <- 1
      if (robust) out["nfi.scaled"] <- (bX2.sc - X2.sc) / bX2.sc
    }
    if ("pnfi" %in% indices) {
      t1 <- bX2 - X2
      t2 <- bX2
      out["pnfi"] <- (DF / bDF) * t1/t2
      if (robust) {
        t1 <- bX2.sc - X2.sc
        t2 <- bX2.sc
        out["pnfi.scaled"] <- (DF / bDF) * t1/t2
      }
    }
    if ("ifi" %in% indices) {
      t1 <- bX2 - X2
      t2 <- bX2 - DF
      out["ifi"] <- if (t2 < 0) 1 else t1/t2
      if (robust) {
        t1 <- bX2.sc - X2.sc
        t2 <- bX2.sc - DF.sc
        if (is.na(t2)) {
          out["ifi.scaled"] <- NA
        } else if (t2 < 0) {
          out["ifi.scaled"] <- 1
        } else out["ifi.scaled"] <- t1/t2
      }
    }

    N <- lavListInspect(object, "ntotal")
    Ns <- lavListInspect(object, "nobs") # N per group
    nG <- lavListInspect(object, "ngroups")
    nVars <- length(lavNames(object))
    if (!(lavoptions$likelihood == "normal" |
          lavoptions$estimator %in% c("ML","PML","FML"))) {
      N <- N - nG
      Ns <- Ns - 1
    }

    if ("mfi" %in% indices) {
      out["mfi"] <- exp(-0.5 * (X2 - DF) / N)
    }

    if ("rmsea" %in% indices) {
      N.RMSEA <- max(N, X2*4) # FIXME: good strategy??

      if (is.na(X2) || is.na(DF)) {
        out["rmsea"] <- as.numeric(NA)
      } else if (DF > 0) {
        getLambda <- function(lambda, chi, df, p) pchisq(chi, df, ncp=lambda) - p

        out["rmsea"] <- sqrt( max(0, (X2/N)/DF - 1/N) ) * sqrt(nG)
        ## lower confidence limit
        if (getLambda(0, X2, DF, .95) < 0.0) out["rmsea.ci.lower"] <- 0 else {
          lambda.l <- try(uniroot(f = getLambda, chi = X2, df = DF, p = .95,
                                  lower = 0, upper = X2)$root, silent = TRUE)
          if (inherits(lambda.l, "try-error")) lambda.l <- NA
          out["rmsea.ci.lower"] <- sqrt( lambda.l/(N*DF) ) * sqrt(nG)
        }
        ## upper confidence limit
        if (getLambda(N.RMSEA, X2, DF, .05) > 0 || getLambda(0, X2, DF, .05) < 0) {
          out["rmsea.ci.upper"] <- 0
        } else {
          lambda.u <- try(uniroot(f = getLambda, chi = X2, df = DF, p = .05,
                                  lower = 0, upper = N.RMSEA)$root, silent = TRUE)
          if (inherits(lambda.u, "try-error")) lambda.u <- NA
          out["rmsea.ci.upper"] <- sqrt( lambda.u/(N*DF) ) * sqrt(nG)
        }
        ## p value
        out["rmsea.pvalue"] <- pchisq(X2, DF, ncp = N*DF*0.05^2/nG,
                                      lower.tail = FALSE)

        ## Scaled versions (naive and robust)
        if (robust & !scaleshift) {
          ## naive
          out["rmsea.scaled"] <- sqrt( max(0, (X2/N)/d - 1/N) ) * sqrt(nG)
          ## lower confidence limit
          if (DF < 1 || d < 1 || getLambda(0, X2, d, .95) < 0.0) {
            out["rmsea.ci.lower.scaled"] <- 0
          } else {
            lambda.l <- try(uniroot(f = getLambda, chi = X2, df = d, p = .95,
                                    lower = 0, upper = X2)$root, silent = TRUE)
            if (inherits(lambda.l, "try-error")) lambda.l <- NA
            out["rmsea.ci.lower.scaled"] <- sqrt( lambda.l/(N*DF) ) * sqrt(nG)
          }
          ## upper confidence limit
          if (DF < 1|| d < 1 || getLambda(0, X2, d, .95) < 0.0 || getLambda(N.RMSEA, X2, d, .05) > 0.0) {
            out["rmsea.ci.upper.scaled"] <- 0
          } else {
            lambda.u <- try(uniroot(f = getLambda, chi = X2, df = d, p = .05,
                                    lower = 0, upper = N.RMSEA)$root, silent = TRUE)
            if (inherits(lambda.u, "try-error")) lambda.u <- NA
            out["rmsea.ci.upper.scaled"] <- sqrt( lambda.u/(N*DF) ) * sqrt(nG)
          }
          ## p value
          out["rmsea.pvalue.scaled"] <- pchisq(X2, d, ncp = N*d*0.05^2/nG,
                                               lower.tail = FALSE)

          if (!pool.robust & test.names[1] %in%
              c("satorra.bentler","yuan.bentler","yuan.bentler.mplus")) {
            ## robust
            out["rmsea.robust"] <- sqrt( max(0, (X2/N)/DF - ch/N ) ) * sqrt(nG)
            ## lower confidence limit
            if (DF.sc < 1 | getLambda(0, X2.sc, DF.sc, .95) < 0.0) {
              out["rmsea.ci.lower.robust"] <- 0
            } else {
              lambda.l <- try(uniroot(f = getLambda, chi = X2.sc, df = DF.sc, p = .95,
                                      lower = 0, upper = X2.sc)$root, silent = TRUE)
              if (inherits(lambda.l, "try-error")) lambda.l <- NA
              out["rmsea.ci.lower.robust"] <- sqrt( (ch*lambda.l)/(N*DF.sc) ) * sqrt(nG)
            }
            ## upper confidence limit
            if (DF.sc < 1 | getLambda(N.RMSEA, X2.sc, DF.sc, .05) > 0.0) {
              out["rmsea.ci.upper.robust"] <- 0
            } else {
              lambda.u <- try(uniroot(f = getLambda, chi = X2.sc, df = DF.sc, p = .05,
                                      lower = 0, upper = N.RMSEA)$root, silent = TRUE)
              if (inherits(lambda.u, "try-error")) lambda.u <- NA
              out["rmsea.ci.upper.robust"] <- sqrt( (ch*lambda.u)/(N*DF.sc) ) * sqrt(nG)
            }
            ## p value
            ########## To be discovered?
          }
        } else if (robust & scaleshift) {
          ## naive only
          out["rmsea.scaled"] <- sqrt( max(0, (X2.sc/N)/DF - 1/N) ) * sqrt(nG)
          ## lower confidence limit
          if (DF < 1 | getLambda(0, X2.sc, DF, .95) < 0.0) {
            out["rmsea.ci.lower.scaled"] <- 0
          } else {
            lambda.l <- try(uniroot(f = getLambda, chi = X2.sc, df = DF, p = .95,
                                    lower = 0, upper = X2.sc)$root, silent = TRUE)
            if (inherits(lambda.l, "try-error")) lambda.l <- NA
            out["rmsea.ci.lower.scaled"] <- sqrt( lambda.l/(N*DF) ) * sqrt(nG)
          }
          ## upper confidence limit
          if (DF < 1 | getLambda(N.RMSEA, X2.sc, DF, .05) > 0.0) {
            out["rmsea.ci.upper.scaled"] <- 0
          } else {
            lambda.u <- try(uniroot(f = getLambda, chi = X2.sc, df = DF, p = .05,
                                    lower = 0, upper = N.RMSEA)$root, silent = TRUE)
            if (inherits(lambda.u, "try-error")) lambda.u <- NA
            out["rmsea.ci.upper.scaled"] <- sqrt( lambda.u/(N*DF) ) * sqrt(nG)
          }
          ## p value
          out["rmsea.pvalue.scaled"] <- pchisq(X2.sc, DF, ncp = N*DF*0.05^2/nG,
                                               lower.tail = FALSE)
        }
      }
    }

    if ("gammahat" %in% indices) {
      out["gammaHat"] <- nVars / (nVars + 2*((X2 - DF) / N))
      out["adjGammaHat"] <- 1 - (((nG * nVars * (nVars + 1)) / 2) / DF) * (1 - out["gammaHat"])
      if (robust) {
        out["gammaHat.scaled"] <- nVars / (nVars + 2*((X2.sc - DF.sc) / N))
        out["adjGammaHat.scaled"] <- 1 - (((nG * nVars * (nVars + 1)) / 2) / DF.sc) * (1 - out["gammaHat.scaled"])
      }
    }

    ## END CHI-SQUARED-BASED FIT INDICES
  } else out <- numeric(0)


  ## RESIDUALS-BASED FIT INDICES

  if (any(grepl(pattern = "rmr", x = indices))) {
    if (lavListInspect(object, "nlevels") > 1L) {
      out["srmr"] <- NA # to preserve the order in lavaan output
      out["srmr_within"] <- getSRMR(object, type = "cor", include.mean = FALSE,
                                    level = "within", omit.imps = omit.imps)
      out["srmr_between"] <- getSRMR(object, type = "cor", include.mean = FALSE,
                                     level = lavListInspect(object, "cluster"),
                                     omit.imps = omit.imps)
      out["srmr"] <- out["srmr_within"] + out["srmr_between"]
    } else {
      out["rmr"] <- getSRMR(object, type = "raw", include.mean = TRUE,
                            omit.imps = omit.imps)
      out["rmr_nomean"] <- getSRMR(object, type = "raw", include.mean = FALSE,
                                   omit.imps = omit.imps)
      out["srmr_bentler"] <- out["srmr"] <- getSRMR(object, type = "cor.bentler",
                                                    include.mean = TRUE,
                                                    omit.imps = omit.imps)
      out["srmr_bentler_nomean"] <- getSRMR(object, type = "cor.bentler",
                                            include.mean = FALSE,
                                            omit.imps = omit.imps)
      out["crmr"] <- getSRMR(object, type = "cor.bollen", include.mean = TRUE,
                             omit.imps = omit.imps)
      out["crmr_nomean"] <- getSRMR(object, type = "cor.bollen",
                                    include.mean = FALSE, omit.imps = omit.imps)
      out["srmr_mplus"] <- getSRMR(object, type = "mplus", include.mean = TRUE,
                                   omit.imps = omit.imps)
      out["srmr_mplus_nomean"] <- getSRMR(object, type = "mplus",
                                          include.mean = FALSE,
                                          omit.imps = omit.imps)
    }
    ## END RESIDUALS-BASED FIT INDICES
  }


  ## return requested measures (loosely matched)
  if ("all" %in% fit.measures) {
    fits <- out
  } else {
    fits <- out[grepl(pattern = paste(fit.measures, collapse = "|"),
                      x = names(out), ignore.case = TRUE)]
    fits <- fits[which(!is.na(names(fits)))]
  }
  class(fits) <- c("lavaan.vector","numeric")
  if (output == "text") class(fits) <- c("lavaan.fitMeasures", class(fits))
  fits
}
##' @name lavaan.mi-class
##' @aliases fitMeasures,lavaan.mi-method
##' @importFrom lavaan fitMeasures
##' @export
setMethod("fitMeasures", "lavaan.mi", fitMeasures.mi)
## lowercase 'm'
##' @name lavaan.mi-class
##' @aliases fitmeasures,lavaan.mi-method
##' @importFrom lavaan fitmeasures
##' @export
setMethod("fitmeasures", "lavaan.mi", fitMeasures.mi)


##' @importFrom lavaan lavListInspect lavNames
##' @importFrom methods getMethod
fitted.lavaan.mi <- function(object, omit.imps = c("no.conv","no.se")) {
  useImps <- rep(TRUE, length(object@DataList))
  if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
  if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
  if ("no.npd" %in% omit.imps) {
    Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
    Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
    useImps <- useImps & !(Heywood.lv | Heywood.ov)
  }
  ## custom removal by imputation number
  rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
  if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
  ## whatever is left
  m <- sum(useImps)
  if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
  useImps <- which(useImps)

  ## how many blocks to loop over
  nG <- lavListInspect(object, "ngroups")
  nlevels <- lavListInspect(object, "nlevels")
  nBlocks <- nG * nlevels #FIXME: always?
  group.label <- if (nG > 1L) lavListInspect(object, "group.label") else NULL
  clus.label <- if (nlevels > 1L) c("within", lavListInspect(object, "cluster")) else NULL
  if (nBlocks > 1L) {
    block.label <- paste(rep(group.label, each = nlevels), clus.label,
                         sep = if (nG > 1L && nlevels > 1L) "_" else "")
  }

  est <- getMethod("coef", "lavaan.mi")(object, omit.imps = omit.imps)
  setpar <- lavaan::lav_model_set_parameters(object@Model, x = est)
  impMats <- lavaan::lav_model_implied(setpar)
  if (lavListInspect(object, "categorical")) {
    th.idx <- lavListInspect(object, "th.idx") # to select $(res.)th
    if (nBlocks == 1L) th.idx <- list(th.idx)  # to loop over
    #FIXME when multilevel accepts categorical
  }

  #TODO: adapt to multilevel, multigroup, or both
  ## loop over (blocks and) moments
  Implied <- vector("list", nBlocks)
  for (b in 1:nBlocks) {
    for (nm in names(impMats)) {

      ## skip any empty objects
      if (is.null(impMats[[nm]][[b]])) next

      Implied[[b]][[nm]] <- impMats[[nm]][[b]]

      ## assign names and classes
      if (nm %in% c("cov","res.cov")) {
        NAMES <- lavNames(object, type = "ov.model", block = b)
        dimnames(Implied[[b]][[nm]]) <- list(NAMES, NAMES)
        class(Implied[[b]][[nm]]) <- c("lavaan.matrix.symmetric","matrix")

      } else if (nm %in% c("mean","res.int")) {
        Implied[[b]][[nm]] <- as.numeric(Implied[[b]][[nm]]) # remove matrix
        names(Implied[[b]][[nm]]) <- lavNames(object, type = "ov.model", block = b)
        class(Implied[[b]][[nm]]) <- c("lavaan.vector","numeric")

      } else if (nm %in% c("th","res.th")) {
        #FIXME: When lavaan allows multilevel categorical, thresholds only
        ##      apply once (not to each level, like for all groups).
        ##      Will lavaan return a vector of zeros for all but "within"?
        ##      If not, it will not exist for each block, so count over groups.
        Implied[[b]][[nm]] <- as.numeric(Implied[[b]][[nm]])[ th.idx[[b]] ] # remove matrix & numeric -means
        names(Implied[[b]][[nm]]) <- lavNames(object, type = "th",
                                              block = b) #FIXME?
        class(Implied[[b]][[nm]]) <- c("lavaan.vector","numeric")

      } else if (nm == "group.w") {
        ## Only for (D)WLS estimation, but when is it relevant?
        ## For now, assign no names/class


      ## The remaining only exist when conditional.x
      } else if (nm %in% c("slopes","res.slopes")) {
        dimnames(Implied[[b]][[nm]]) <- list(lavNames(object, type = "ov.nox", block = b),
                                             lavNames(object, type = "ov.x", block = b))
        class(Implied[[b]][[nm]]) <- c("lavaan.matrix","matrix")

      } else if (nm == "cov.x") {
        NAMES <- lavNames(object, type = "ov.x", block = b)
        dimnames(Implied[[b]][[nm]]) <- list(NAMES, NAMES)
        class(Implied[[b]][[nm]]) <- c("lavaan.matrix.symmetric","matrix")

      } else if (nm == "mean.x") {
        Implied[[b]][[nm]] <- as.numeric(Implied[[b]][[nm]]) # remove matrix
        names(Implied[[b]][[nm]]) <- lavNames(object, type = "ov.x", block = b)
        class(Implied[[b]][[nm]]) <- c("lavaan.vector","numeric")
      }

    ## end loops
    }
  }

  ## drop list for 1 block, or add labels for multiple
  if (nBlocks == 1L)  {
    Implied <- Implied[[1]]
  } else names(Implied) <- block.label

  Implied
}
##' @name lavaan.mi-class
##' @aliases fitted,lavaan.mi-method
##' @export
setMethod("fitted", "lavaan.mi", fitted.lavaan.mi)
##' @name lavaan.mi-class
##' @aliases fitted.values,lavaan.mi-method
##' @export
setMethod("fitted.values", "lavaan.mi", fitted.lavaan.mi)



##' @importFrom lavaan lavListInspect
##' @importFrom methods getMethod
##' @importFrom stats cov2cor
resid.lavaan.mi <- function(object, type = c("raw","cor"),
                            omit.imps = c("no.conv","no.se")) {
  ## @SampleStatsList is (for each imputation) output from:
  ##    getSampStats <- function(obj) lavInspect(obj, "sampstat")
  useImps <- rep(TRUE, length(object@DataList))
  if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
  if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
  if ("no.npd" %in% omit.imps) {
    Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
    Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
    useImps <- useImps & !(Heywood.lv | Heywood.ov)
  }
  ## custom removal by imputation number
  rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
  if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
  ## whatever is left
  m <- sum(useImps)
  if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
  useImps <- which(useImps)

  ## check type options
  type <- tolower(type[1])
  if (type %in% c("raw","rmr")) {
    type = "raw"
  } else if (type %in% c("cor","cor.bollen","crmr")) {
    type <- "cor.bollen"
  } else if (type %in% c("cor.bentler","cor.eqs","srmr")) {
    type <- "cor.bentler"
  } else stop('type="', type, '" not supported for lavaan.mi objects')

  ## how many blocks to loop over
  nG <- lavListInspect(object, "ngroups")
  nlevels <- lavListInspect(object, "nlevels")
  nBlocks <- nG * nlevels #FIXME: always?
  group.label <- if (nG > 1L) lavListInspect(object, "group.label") else NULL
  clus.label <- if (nlevels > 1L) c("within", lavListInspect(object, "cluster")) else NULL
  if (nBlocks > 1L) {
      block.label <- paste(rep(group.label, each = nlevels), clus.label,
                           sep = if (nG > 1L && nlevels > 1L) "_" else "")
  }

  if (lavListInspect(object, "categorical")) {
    th.idx <- lavListInspect(object, "th.idx") # to select $(res.)th
    if (nBlocks == 1L) th.idx <- list(th.idx)  # to loop over
    #FIXME when multilevel accepts categorical
  }

  ## H0-model-implied moments, already pooled
  ## (moments-list nested in block-list)
  Implied <- getMethod("fitted", "lavaan.mi")(object, omit.imps = omit.imps)
  if (nBlocks == 1L) Implied <- list(Implied) # store single block in a block-list

  ## template to store observed moments & residuals
  RES <- OBS <- vector("list", nBlocks)

  ## loop over (blocks and) moments
  for (b in 1:nBlocks) {
    for (nm in names(Implied[[b]])) {

      ## skip if Implied element is not part of the saturated list
      if (is.null(object@h1List[[ useImps[1] ]]$implied[[nm]][[b]])) next

      ## H1 (saturated model) implied moments
      ## (block-list nested in moments-list)
      momentList <- lapply(object@h1List[useImps],
                           function(x) x$implied[[nm]][[b]])
      OBS[[b]][[nm]] <- Reduce("+", momentList) / m
      #TODO: unnecessary calculation if standardized and nm %in% c("th","slopes")

      ## remove numeric -means from thresholds
      if (nm %in% c("th","res.th")) OBS[[b]][[nm]] <- as.numeric(OBS[[b]][[nm]])[ th.idx[[b]] ]

      ## calculate residuals
      if (type == "raw") {
        RES[[b]][[nm]] <- OBS[[b]][[nm]] - Implied[[b]][[nm]]
        class(RES[[b]][[nm]]) <- class(Implied[[b]][[nm]])


        ## correlation residuals
      } else if (type == "cor.bollen") {

        if (nm %in% c("cov","res.cov")) {
          RES[[b]][[nm]] <- cov2cor(OBS[[b]][[nm]]) - cov2cor(Implied[[b]][[nm]])
          class(RES[[b]][[nm]]) <- c("lavaan.matrix.symmetric","matrix")

          ## mean structure
        } else if (nm == "mean") {
          std.obs.M <- OBS[[b]][[nm]] / sqrt(diag(OBS[[b]]$cov))
          std.mod.M <- Implied[[b]][[nm]] / sqrt(diag(Implied[[b]]$cov))
          RES[[b]][[nm]] <- std.obs.M - std.mod.M
          class(RES[[b]][[nm]]) <- c("lavaan.vector","numeric")
        } else if (nm == "res.int") {
          std.obs.M <- OBS[[b]][[nm]] / sqrt(diag(OBS[[b]]$res.cov))
          std.mod.M <- Implied[[b]][[nm]] / sqrt(diag(Implied[[b]]$res.cov))
          RES[[b]][[nm]] <- std.obs.M - std.mod.M
          class(RES[[b]][[nm]]) <- c("lavaan.vector","numeric")

          ## thresholds, slopes, cov.x, mean.x
        } else {
          #FIXME: lavaan currently (0.6-4.1399) returns nothing
          next
        }


        ## standardized (by observed SDs) residuals
      } else if (type == "cor.bentler") {

        if (nm %in% c("cov","mean")) {
          SDs <- diag(sqrt(diag(OBS[[b]]$cov)))
          dimnames(SDs) <- dimnames(OBS[[b]][[nm]])
        } else if (nm %in% c("res.cov","res.int")) {
          SDs <- diag(sqrt(diag(OBS[[b]]$res.cov)))
          dimnames(SDs) <- dimnames(OBS[[b]][[nm]])
        } else {
          #FIXME: lavaan currently (0.6-4.1399) returns nothing for "th" or "slopes"
          next
        }


        if (nm %in% c("cov","res.cov")) {
          RES[[b]][[nm]] <- solve(SDs) %*% (OBS[[b]][[nm]] - Implied[[b]][[nm]]) %*% solve(SDs)
          class(RES[[b]][[nm]]) <- c("lavaan.matrix.symmetric","matrix")
        } else if (nm %in% c("mean","res.int")) {
          RES[[b]][[nm]] <- (OBS[[b]][[nm]] - Implied[[b]][[nm]]) / diag(SDs)
          class(RES[[b]][[nm]]) <- c("lavaan.vector","numeric")
        }

      }

      ## copy names from fitted() results
      if (is.null(dim(RES[[b]][[nm]]))) {
        names(RES[[b]][[nm]]) <- names(Implied[[b]][[nm]])
      } else dimnames(RES[[b]][[nm]]) <- dimnames(Implied[[b]][[nm]])

    ## end loop over moments
    }

    ## add type to beginning of each block's list
    RES[[b]] <- c(list(type = type), RES[[b]])

    #TODO: Rename (res.)cov to (res.)cor?  lavResiduals() does not

  }

  ## drop list for 1 block
  if (nBlocks == 1L) {
    RES <- RES[[1]]
  } else names(RES) <- block.label #FIXME: will lavaan do this in the future?

  RES
}
##' @name lavaan.mi-class
##' @aliases residuals,lavaan.mi-method
##' @export
setMethod("residuals", "lavaan.mi", resid.lavaan.mi)
##' @name lavaan.mi-class
##' @aliases resid,lavaan.mi-method
##' @export
setMethod("resid", "lavaan.mi", resid.lavaan.mi)



