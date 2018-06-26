### Terrence D. Jorgensen
### Last updated: 26 June 2018
### Class and Methods for lavaan.mi object, returned by runMI()


#' Class for a lavaan Model Fitted to Multiple Imputations
#'
#' This class extends the \code{\linkS4class{lavaanList}} class, created by
#' fitting a lavaan model to a list of data sets. In this case, the list of
#' data sets are multiple imputations of missing data.
#'
#'
#' @name lavaan.mi-class
#' @importClassesFrom lavaan lavaanList
#' @aliases lavaan.mi-class show,lavaan.mi-method summary,lavaan.mi-method
#' anova,lavaan.mi-method nobs,lavaan.mi-method coef,lavaan.mi-method
#' vcov,lavaan.mi-method fitted,lavaan.mi-method fitted.values,lavaan.mi-method
#' residuals,lavaan.mi-method resid,lavaan.mi-method
#' @docType class
#'
#' @slot coefList \code{list} of estimated coefficients in matrix format (one
#'  per imputation) as output by \code{\link[lavaan]{lavInspect}(fit, "est")}
#' @slot GLIST pooled \code{list} of coefficients in GLIST format
#' @slot miList \code{list} of modification indices output by
#'  \code{\link[lavaan]{modindices}}
#' @slot seed \code{integer} seed set before running imputations
#' @slot lavListCall call to \code{\link[lavaan]{lavaanList}} used to fit the
#'  model to the list of imputed data sets in \code{@@DataList}, stored as a
#'  \code{list} of arguments
#' @slot imputeCall call to imputation function (if used), stored as a
#'  \code{list} of arguments
#' @slot convergence \code{list} of \code{logical} vectors indicating whether,
#'  for each imputed data set, (1) the model converged on a solution, (2)
#'  \emph{SE}s could be calculated, (3) the (residual) covariance matrix of
#'  latent variables (\eqn{\Psi}) is non-positive-definite, and (4) the residual
#'  covariance matrix of observed variables (\eqn{\Theta}) is
#'  non-positive-definite.
#' @slot lavaanList_slots All remaining slots are from
#'  \code{\linkS4class{lavaanList}}, but \code{\link{runMI}} only populates a
#'  subset of the \code{list} slots, two of them with custom information:
#' @slot DataList The \code{list} of imputed data sets
#' @slot SampleStatsList List of output from
#'  \code{\link[lavaan]{lavInspect}(fit, "sampstat")} applied to each fitted
#'  model
#' @slot ParTableList See \code{\linkS4class{lavaanList}}
#' @slot vcovList See \code{\linkS4class{lavaanList}}
#' @slot testList See \code{\linkS4class{lavaanList}}
#'
#' @param object An object of class \code{lavaan.mi}
#' @param se,ci,level,standardized,rsquare,header,add.attributes See
#'        \code{\link[lavaan]{parameterEstimates}}.
#' @param fmi \code{logical} indicating whether to include the Fraction Missing
#'        Information (FMI) for parameter estimates in the \code{summary} output
#'        (see \bold{Value} section).
#' @param asymptotic \code{logical}. If \code{FALSE} (typically a default, but
#'       see \bold{Value} section for details using various methods), pooled
#'       tests (of fit or pooled estimates) will be \emph{F} or \emph{t}
#'       statistics with associated degrees of freedom (\emph{df}). If
#'       \code{TRUE}, the (denominator) \emph{df} are assumed to be sufficiently
#'       large for a \emph{t} statistic to follow a normal distribution, so it
#'       is printed as a \emph{z} statisic; likewise, \emph{F} times its
#'       numerator \emph{df} is printed, assumed to follow a \eqn{\chi^2}
#'       distribution.
#' @param scale.W \code{logical}. If \code{TRUE} (default), the \code{vcov}
#'       method will calculate the pooled covariance matrix by scaling the
#'       within-imputation component by the ARIV (see Enders, 2010, p. 235,
#'       for definition and formula). Otherwise, the pooled matrix is calculated
#'       as the weighted sum of the within-imputation and between-imputation
#'       components (see Enders, 2010, ch. 8, for details). This in turn affects
#'       how the \code{summary} method calcualtes its pooled standard errors, as
#'       well as the Wald test (\code{anova(..., test = "D1")}).
#' @param labels \code{logical} indicating whether the \code{coef} output should
#'        include parameter labels. Default is \code{TRUE}.
#' @param total \code{logical} (default: \code{TRUE}) indicating whether the
#'        \code{nobs} method should return the total sample size or (if
#'        \code{FALSE}) a vector of group sample sizes.
#' @param type The meaning of this argument varies depending on which method it
#'        it used for. Find detailed descriptions in the \bold{Value} section
#'        under \code{coef}, \code{vcov}, \code{residuals}, and \code{anova}.
#' @param h1 An object of class \code{lavaan.mi} in which \code{object} is
#'        nested, so that their difference in fit can be tested using
#'        \code{anova} (see \bold{Value} section for details).
#' @param test \code{character} indicating the method used to pool model-fit or
#'        model-comparison test statistics:
#'   \itemize{
#'      \item{\code{"D3": }}{The default test (\code{"D3"}, or any of
#'            \code{"mr", "Meng.Rubin", "likelihood", "LRT"}) is a pooled
#'            likeliehood-ratio test (see Enders, 2010, ch. 8).
#'            \code{test = "mplus"} implies \code{"D3"} and \code{asymptotic =
#'            TRUE} (see Asparouhov & Muthen, 2010). When using a non-likelihood
#'            estimator (e.g., DWLS for categorical outcomes), \code{"D3"} is
#'            unavailable, so the default is changed to \code{"D2"}.}
#'      \item{\code{"D2": }}{Returns a pooled test statistic, as described by
#'            Li, Meng, Raghunathan, & Rubin (1991) and Enders (2010, chapter 8).
#'            Aliases include \code{"lmrr", "Li.et.al", "pooled.wald"}).}
#'      \item{\code{"D1": }}{Returns a Wald test calculated for constraints on
#'            the pooled point estimates, using the pooled covariance matrix of
#'            parameter estimates; see \code{\link[lavaan]{lavTestWald}} for
#'            details. \code{h1} is ignored when \code{test = "D1"}, and
#'            \code{constraints} is ignored when \code{test != "D1"}. The
#'            \code{scale.W} argument is passed to the \code{vcov} method (see
#'            \bold{Value} section for details).}
#'    }
#' @param pool.robust \code{logical}. Ignored unless \code{test = "D2"} and a
#'      robust test was requested. If \code{pool.robust = TRUE}, the robust test
#'      statistic is pooled, whereas \code{pool.robust = FALSE} will pool
#'      the naive test statistic (or difference statistic) and apply the average
#'      scale/shift parameter to it (unavailable for mean- and variance-adjusted
#'      difference statistics, so \code{pool.robust} will be set \code{TRUE}).
#'      If \code{test = "D2"} and \code{pool.robust = TRUE}, further options
#'      can be passed to \code{\link[lavaan]{lavTestLRT}} (see below).
#' @param indices \code{logical}, or \code{character} vector naming fit indices
#'      to be printed with test of model fit.  Ignored \code{if (!is.null(h1))}.
#'      See description of \code{anova} in \bold{Value} section for details.
#' @param constraints See \code{\link[lavaan]{lavTestWald}}.
#' @param method,A.method,H1,scaled.shifted See \code{\link[lavaan]{lavTestLRT}}.
#' @param fit.measures,baseline.model See \code{\link[lavaan]{fitMeasures}}.
#'
#' @return
#' \item{coef}{\code{signature(object = "lavaan.mi", type = "free", labels = TRUE)}:
#'  See \code{\linkS4class{lavaan}}. Returns the pooled point estimates (i.e.,
#'  averaged across imputed data sets; see Rubin, 1987).}
#'
#' \item{vcov}{\code{signature(object = "lavaan.mi", scale.W = TRUE,
#'  type = c("pooled","between","within","ariv"))}:  By default, returns the
#'  pooled covariance matrix of parameter estimates (\code{type = "pooled"}),
#'  the within-imputations covariance matrix (\code{type = "within"}), the
#'  between-imputations covariance matrix (\code{type = "between"}), or the
#'  average relative increase in variance (\code{type = "ariv"}) due to missing
#'  data.}
#'
#' \item{fitted.values}{\code{signature(object = "lavaan.mi")}: See
#'  \code{\linkS4class{lavaan}}. Returns model-implied moments, evaluated at the
#'  pooled point estimates.}
#' \item{fitted}{\code{signature(object = "lavaan.mi")}:
#'   alias for \code{fitted.values}}
#'
#' \item{residuals}{\code{signature(object = "lavaan.mi", type = c("raw","cor"))}:
#'  See \code{\linkS4class{lavaan}}. By default (\code{type = "raw"}), returns
#'  the difference between the model-implied moments from \code{fitted.values}
#'  and the pooled observed moments (i.e., averaged across imputed data sets).
#'  Standardized residuals are also available, using Bollen's
#'  (\code{type = "cor"} or \code{"cor.bollen"}) or Bentler's
#'  (\code{type = "cor.bentler"}) formulas.}
#' \item{resid}{\code{signature(object = "lavaan.mi", type = c("raw","cor"))}:
#'  alias for \code{residuals}}
#'
#' \item{nobs}{\code{signature(object = "lavaan.mi", total = TRUE)}: either
#'  the total (default) sample size or a vector of group sample sizes
#'  (\code{total = FALSE}).}
#'
#' \item{anova}{\code{signature(object = "lavaan.mi", h1 = NULL,
#'   test = c("D3","D2","D1"), pool.robust = FALSE, scale.W = TRUE,
#'   asymptotic = FALSE, constraints = NULL, indices = FALSE, baseline.model = NULL,
#'   method = "default", A.method = "delta", H1 = TRUE, type = "Chisq")}:
#'   Returns a test of model fit if \code{h1} is \code{NULL}, or a test
#'   of the difference in fit between nested models if \code{h1} is another
#'   \code{lavaan.mi} object, assuming \code{object} is nested in \code{h1}. If
#'   \code{asymptotic}, the returned test statistic will follow a \eqn{\chi^2}
#'   distribution in sufficiently large samples; otherwise, it will follow an
#'   \emph{F} distribution. If a robust test statistic is detected in the
#'   \code{object} results (it is assumed the same was requested in \code{h1},
#'   if provided), then \code{asymptotic} will be set to \code{TRUE} and the
#'   pooled test statistic will be scaled using the average scaling factor (and
#'   average shift parameter or \emph{df}, if applicable) across imputations
#'   (unless \code{pool.robust = FALSE} and \code{test = "D2"}; see below).
#'
#'   When \code{indices = TRUE} and \code{is.null(h1)}, popular indices of
#'   approximate fit (CFI, TLI/NNFI, RMSEA with CI, and SRMR) will be returned
#'   for \code{object}; see \code{\link[lavaan]{fitMeasures}} for more details.
#'   Specific indices can be requested with a \code{character} vector (any of
#'   \code{"mfi", "rmsea", "gammaHat", "rmr", "srmr", "cfi", "tli", "nnfi",
#'   "rfi", "nfi", "pnfi", "ifi", "rni"}), or all available indices will be
#'   returned if \code{indices = "all"}. Users can specify a custom
#'   \code{baseline.model}, also fit using \code{runMI}, to calculate
#'   incremental fit indices (e.g., CFI, TLI). If \code{baseline.model = NULL},
#'   the default independence model will be used.}
#'
#' \item{fitMeasures}{\code{signature(object = "lavaan.mi",
#'   fit.measures = "all", baseline.model = NULL)}: arguments are consistent
#'   with lavaan's \code{\link[lavaan]{fitMeasures}}. This merely calls the
#'   \code{anova} method described above, with \code{indices = fit.measures}
#'   and \code{baseline.model = baseline.model}, and default values for the
#'   remaining arguments. The user has more control (e.g., over pooling methods)
#'   using \code{anova} directly.}
#' \item{fitmeasures}{alias for \code{fitMeasures}.}
#'
#' \item{show}{\code{signature(object = "lavaan.mi")}: returns a message about
#'  convergence rates and estimation problems (if applicable) across imputed
#'  data sets.}
#'
#' \item{summary}{\code{signature(object = "lavaan.mi", se = TRUE, ci = FALSE,
#'  level = .95, standardized = FALSE, rsquare = FALSE, fmi = FALSE,
#'  scale.W = FALSE, asymptotic = FALSE, add.attributes = TRUE)}: see
#'  \code{\link[lavaan]{parameterEstimates}} for details.
#'  By default, \code{summary} returns pooled point and \emph{SE}
#'  estimates, along with \emph{t} test statistics and their associated
#'  \emph{df} and \emph{p} values. If \code{ci = TRUE}, confidence intervales
#'  are returned with the specified confidence \code{level} (default 95\% CI).
#'  If \code{asymptotic = TRUE}, \emph{z} instead of \emph{t} tests are
#'  returned. \code{standardized} solution(s) can also be requested by name
#'  (\code{"std.lv"} or \code{"std.all"}) or both are returned with \code{TRUE}.
#'  \emph{R}-squared for endogenous variables can be requested, as well as the
#'  Fraction Missing Information (FMI) for parameter estimates. By default, the
#'  output will appear like \code{lavaan}'s \code{summary} output, but if
#'  \code{add.attributes = FALSE}, the returned \code{data.frame} will resemble
#'  the \code{parameterEstimates} output. The \code{scale.W} argument is
#'  passed to \code{vcov} (see description above).}
#'
#' @section Objects from the Class: See the \code{\link{runMI}} function for
#' details. Wrapper functions include \code{\link{lavaan.mi}},
#' \code{\link{cfa.mi}}, \code{\link{sem.mi}}, and \code{\link{growth.mi}}.
#' @author Terrence D. Jorgensen (University of Amsterdam;
#' \email{TJorgensen314@@gmail.com})
#' @references Asparouhov, T., & Muthen, B. (2010). \emph{Chi-square statistics
#' with multiple imputation}. Technical Report. Retrieved from
#' \url{www.statmodel.com}
#'
#' Enders, C. K. (2010). \emph{Applied missing data analysis}. New York, NY:
#' Guilford.
#'
#' Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
#' Significance levels from repeated \emph{p}-values with multiply-imputed data.
#' \emph{Statistica Sinica, 1}(1), 65--92. Retrieved from
#' \url{http://www.jstor.org/stable/24303994}
#'
#' Meng, X.-L., & Rubin, D. B. (1992). Performing likelihood ratio tests with
#' multiply-imputed data sets. \emph{Biometrika, 79}(1), 103--111. Retrieved
#' from \url{http://www.jstor.org/stable/2337151}
#'
#' Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}.
#' New York, NY: Wiley.
#' @examples
#'
#' ## See ?runMI help page
#'
setClass("lavaan.mi", contains = "lavaanList",
         slots = c(coefList = "list",     # coefficients in matrix format
                   GLIST = "list",        # list of pooled coefs in GLIST format
                   miList = "list",       # modification indices
                   seed = "integer",      # seed set before running imputations
                   lavListCall = "list",  # store actual call to lavaanList
                   imputeCall = "list",   # store call from imputation, if used
                   convergence = "list")) # also check SEs and Heywood cases



#' @name lavaan.mi-class
#' @aliases show,lavaan.mi-method
#' @export
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


#' @importFrom stats pt qt pnorm qnorm
#' @importFrom lavaan lavListInspect parTable
summary.lavaan.mi <- function(object, se = TRUE, ci = FALSE, level = .95,
                              standardized = FALSE, rsquare = FALSE,
                              fmi = FALSE, header = TRUE, scale.W = TRUE,
                              asymptotic = FALSE, add.attributes = TRUE) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  m <- sum(useImps)
  ## extract parameter table with attributes for printing
  PT <- parTable(object)
  myCols <- c("lhs","op","rhs")
  if (lavListInspect(object, "ngroups") > 1L) myCols <- c(myCols,"block","group")
  PE <- PT[ , myCols]
  free <- PT$free > 0L | PT$op == ":="
  STDs <- !(PT$op %in% c("==","<",">")) # which rows can be standardized

  PE$est <- rowMeans(sapply(object@ParTableList[useImps], "[[", i = "est"))

  if (lavListInspect(object, "options")$se == "none") {
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
    VCOV <- getMethod("vcov","lavaan.mi")(object, scale.W = scale.W)
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
      ## can't do finite-sample correction because Wald z tests have no df (see Enders, 2010, p. 231, eq. 8.13 & 8.14)
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
      PE$std.lv[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                      type = "std.lv",
                                                      GLIST = object@GLIST,
                                                      est = PE$est)$est.std
      PE$std.all[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                       type = "std.all",
                                                       GLIST = object@GLIST,
                                                       est = PE$est)$est.std
    }
  } else if (tolower(as.character(standardized)[1]) == "std.lv") {
    PE$std.lv[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                    type = "std.lv",
                                                    GLIST = object@GLIST,
                                                    est = PE$est)$est.std
  } else if (tolower(as.character(standardized)[1]) == "std.all") {
    PE$std.all[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                     type = "std.all",
                                                     GLIST = object@GLIST,
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
  if (add.attributes) {
    PE$label <- PT$label
    PE$exo <- 0L # because PT$exo must be when !fixed.x
    class(PE) <- c("lavaan.parameterEstimates","lavaan.data.frame","data.frame")
    lavops <- lavListInspect(object, "options")
    attr(PE, "information") <- lavops$information
    attr(PE, "se") <- lavops$se
    attr(PE, "group.label") <- lavListInspect(object, "group.label")
    attr(PE, "level.label") <- object@Data@level.label #FIXME: lavListInspect?
    attr(PE, "bootstrap") <- lavops$bootstrap
    attr(PE, "bootstrap.successful") <- 0L #FIXME: assumes none. Implement Wei & Fan's mixing method?
    attr(PE, "missing") <- lavops$missing
    attr(PE, "observed.information") <- lavops$observed.information
    attr(PE, "h1.information") <- lavops$h1.information
    attr(PE, "header") <- header
    # FIXME: lavaan may add more!!
    if (fmi) cat("\n", messRIV, sep = "")
  } else {
    class(PE) <- c("lavaan.data.frame","data.frame")
  }
  ## requested R-squared?
  endoNames <- c(lavaan::lavNames(object, "ov.nox"),
                 lavaan::lavNames(object, "lv.nox"))
  if (rsquare & length(endoNames)) {
    isEndo <- sapply(PE$lhs, function(x) x %in% endoNames)
    rsqPE <- PE[PE$lhs == PE$rhs & PE$op == "~~" & isEndo, ]
    rsqPE$op <- "r2"
    for (i in which(!sapply(colnames(PE),
                            function(x) x %in% c("lhs","op","rhs","block","group","est","exo")))) {
      rsqPE[ , i] <- NA
    }
    STD <- lavaan::standardizedSolution(object, se = FALSE, type = "std.all",
                                        GLIST = object@GLIST, est = PE$est)
    isEndoSTD <- sapply(STD$lhs, function(x) x %in% endoNames)
    std.all <- STD$est.std[STD$lhs == STD$rhs & STD$op == "~~" & isEndoSTD]
    rsqPE$est <- ifelse(std.all < 0, NA, 1 - std.all) # negative variances
    if (add.attributes) rsqPE$label <- ""
    PE <- rbind(PE, rsqPE)
  }

  if (!add.attributes) PE <- PE[!(PE$op %in% c("==","<",">")), ]
  rownames(PE) <- NULL
  if (add.attributes) {
    getMethod("show", "lavaan.mi")(object)
    cat(messPool)
  }
  ## FIXME: ask Yves to make this accessible somehow, or hack it?
  # if (fit.measures) lavaan:::print.fit.measures(fitMeasures(object))
  PE
}
#' @name lavaan.mi-class
#' @aliases summary,lavaan.mi-method
#' @export
setMethod("summary", "lavaan.mi", summary.lavaan.mi)


#' @name lavaan.mi-class
#' @aliases nobs,lavaan.mi-method
#' @importFrom lavaan lavListInspect
#' @export
setMethod("nobs", "lavaan.mi", function(object, total = TRUE) {
  if (total) return(lavListInspect(object, "ntotal"))
  N <- lavListInspect(object, "norig")
  if (length(N) > 1L) names(N) <- lavListInspect(object, "group.label")
  N
})



#' @importFrom lavaan parTable
coef.lavaan.mi <- function(object, type = "free", labels = TRUE) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
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
#' @name lavaan.mi-class
#' @aliases coef,lavaan.mi-method
#' @export
setMethod("coef", "lavaan.mi", coef.lavaan.mi)



#' @importFrom stats cov
#' @importFrom lavaan lavListInspect parTable
vcov.lavaan.mi <- function(object, type = c("pooled","between","within","ariv"),
                           scale.W = TRUE) {
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
  useImps <- sapply(object@convergence, "[[", i = "converged")
  m <- sum(useImps)

  useSE <- sapply(object@convergence, "[[", i = "SE")
  useSE[is.na(useSE)] <- FALSE

  coefList <- lapply(object@ParTableList[useImps], "[[", i = "est")
  B <- cov(do.call(rbind, coefList)[ , PT$free > 0L & !duplicated(PT$free)])
  class(B) <- c("lavaan.matrix.symmetric","matrix")
  rownames(B) <- colnames(B) <- lavaan::lav_partable_labels(PT, type = "free")
  if (type == "between") return(B)

  if (sum(useSE) == 0L) stop('Standard errors could not be computed in any ',
                             'imputations, so it is not possible to calculate ',
                             'the within-imputation portion of sampling variance.')
  W <- Reduce("+", lapply(object@vcovList[useSE], function(x) x$vcov)) / sum(useSE)
  class(W) <- c("lavaan.matrix.symmetric","matrix")
  dimnames(W) <- dimnames(B)
  if (type == "within") return(W)

  if (!all(useImps == useSE))
    warning('Between-imputation covariance matrix based on estimated parameters',
            ' from ', m, ' converged solutions, but the mean within-imputation',
            ' covariance matrix based on ', sum(useSE), ' solutions for which',
            ' standard errors could be calculated.  Pooled total covariance',
            ' matrix is therefore based on different imputed data sets.')

  ## check whether equality constraints prevent inversion of W
  if (scale.W || type == "ariv") {
    inv.W <- if (ncon == 0) try(solve(W), silent = TRUE) else MASS::ginv(W)
    if (inherits(inv.W, "try-error")) {
      if (ncon == 0) {
        warning("Could not invert within-imputation covariance matrix. ",
                "Generalized inverse used instead.\n",
                "It may be safer to set `scale.W = FALSE'.")
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
#' @name lavaan.mi-class
#' @aliases vcov,lavaan.mi-method
#' @export
setMethod("vcov", "lavaan.mi", vcov.lavaan.mi)


#' @importFrom stats pf pchisq
#' @importFrom lavaan parTable
D1 <- function(object, constraints = NULL, scale.W = FALSE,
               asymptotic = FALSE, verbose = FALSE) {
  ## "borrowed" lavTestWald()
  nImps <- sum(sapply(object@convergence, "[[", i = "converged"))
  if (nImps == 1L) stop("model did not converge on any imputations")
  if (is.null(constraints) || nchar(constraints) == 0L) stop("constraints are empty")

  # remove == constraints from parTable, save as list
  PT <- parTable(object)
  partable <- as.list(PT[PT$op != "==", ])

  # parse constraints
  FLAT <- lavaan::lavParseModelString( constraints )
  CON <- attr(FLAT, "constraints")
  LIST <- list()
  if (length(CON) > 0L) {
    lhs <- unlist(lapply(CON, "[[", i = "lhs"))
    op <- unlist(lapply(CON, "[[", i = "op"))
    rhs <- unlist(lapply(CON, "[[", i = "rhs"))
    LIST$lhs <- c(LIST$lhs, lhs) # FIXME: why concatenate with NULL?
    LIST$op  <- c(LIST$op,  op)
    LIST$rhs <- c(LIST$rhs, rhs)
  } else stop("no equality constraints found in constraints argument")

  # theta = free parameters only (equality-constrained allowed)
  theta <- getMethod("coef", "lavaan.mi")(object) #object@optim$x

  # build constraint function
  ceq.function <- lavaan::lav_partable_constraints_ceq(partable = partable,
                                                       con = LIST, debug = FALSE)
  # compute jacobian restrictions
  JAC <- try(lavaan::lav_func_jacobian_complex(func = ceq.function, x = theta),
             silent = TRUE)
  if (inherits(JAC, "try-error")) { # eg. pnorm()
    JAC <- lavaan::lav_func_jacobian_simple(func = ceq.function, x = theta)
  }
  if (verbose) {cat("Restriction matrix (jacobian):\n"); print(JAC); cat("\n")}

  # linear restriction
  theta.r <- ceq.function( theta )
  if (verbose) {cat("Restricted theta values:\n"); print(theta.r); cat("\n")}

  # get VCOV
  VCOV <- getMethod("vcov","lavaan.mi")(object, scale.W = scale.W)

  # restricted vcov
  info.r  <- JAC %*% VCOV %*% t(JAC)

  # Wald test statistic
  test.stat <- as.numeric(t(theta.r) %*% solve( info.r ) %*% theta.r)

  # number of constraints (k in Enders (2010, p. 235) eqs. 8.23-25)
  DF <- nrow(JAC)

  if (asymptotic) {
    out <- c("chisq" = test.stat, df = DF,
             pvalue = pchisq(test.stat, df = DF, lower.tail = FALSE))
  } else {
    W <- getMethod("vcov", "lavaan.mi")(object, type = "within")
    B <- getMethod("vcov", "lavaan.mi")(object, type = "between")
    #FIXME: only valid for linear constraints?
    ## restricted B & W components of VCOV
    W.r  <- JAC %*% W %*% t(JAC)
    B.r  <- JAC %*% B %*% t(JAC)
    ## relative increase in variance due to missing data
    W.inv <- MASS::ginv(W.r)
    ariv <- (1 + 1/nImps) * sum(diag(B.r %*% W.inv)) / DF
    ## calculate denominator DF for F statistic
    a <- DF*(nImps - 1)
    if (a > 4) {
      v2 <- 4 + (a - 4) * (1 + (1 - 2/a)*(1 / ariv))^2 # Enders (eq. 8.24)
    } else {
      v2 <- a*(1 + 1/DF) * (1 + 1/ariv)^2 / 2 # Enders (eq. 8.25)
    }
    out <- c("F" = test.stat / DF, df1 = DF, df2 = v2,
             pvalue = pf(test.stat / DF, df1 = DF, df2 = v2, lower.tail = FALSE))
  }

  class(out) <- c("lavaan.vector","numeric")
  out
}
#' @importFrom stats var pf pchisq
calculate.D2 <- function(w, DF, asymptotic = FALSE) {
  if (!length(w)) return(NA)
  nImps <- sum(!is.na(w))
  if (nImps == 0) return(NA)
  w_bar <- mean(w, na.rm = TRUE)
  ariv <- (1 + 1/nImps) * var(sqrt(w), na.rm = TRUE)
  test.stat <- (w_bar/DF - ((nImps + 1) * ariv / (nImps - 1))) / (1 + ariv)
  if (test.stat < 0) test.stat <- 0
  if (asymptotic) {
    out <- c("chisq" = test.stat * DF, df = DF,
             pvalue = pchisq(test.stat * DF, df = DF, lower.tail = FALSE))
  } else {
    v3 <- DF^(-3 / nImps) * (nImps - 1) * (1 + (1 / ariv))^2
    out <- c("F" = test.stat, df1 = DF, df2 = v3,
             pvalue = pf(test.stat, df1 = DF, df2 = v3, lower.tail = FALSE))
  }
  out
}
#' @importFrom lavaan lavListInspect parTable
D2 <- function(object, h1 = NULL, asymptotic = FALSE, pool.robust = FALSE,
               method = "default", A.method = "delta", H1 = TRUE,
               scaled.shifted = TRUE, type = "Chisq") {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  lavoptions <- lavListInspect(object, "options")

  if (pool.robust & !is.null(h1)) {
    PT1 <- parTable(h1)
    op1 <- lavListInspect(h1, "options")
    oldCall <- object@lavListCall #re-run lavaanList() and save DIFFTEST
    if (!is.null(oldCall$parallel)) {
      if (oldCall$parallel == "snow") {
        oldCall$parallel <- "no"
        oldCall$ncpus <- 1L
        if (lavoptions$warn) warning("Unable to pass lavaan::lavTestLRT() ",
                          "arguments when parallel = 'snow'.\n",
                          "Switching to parallel = 'no'.",
                          " Unless using Windows, parallel='multicore' works.")
      }
    }

    ## call lavaanList() again to run lavTestLRT() on each imputation
    oldCall$FUN <- function(obj) {
      fit1 <- try(lavaan::lavaan(PT1, slotOptions = op1, slotData = obj@Data),
                  silent = TRUE)
      if (inherits(fit1, "try-error")) return("fit failed")
      out <- try(lavaan::lavTestLRT(obj, fit1, H1 = H1, method = method,
                                    A.method = A.method, type = type,
                                    scaled.shifted = scaled.shifted),
                 silent = TRUE)
      if (inherits(out, "try-error")) return("lavTestLRT() failed")
      c(chisq = out[2, "Chisq diff"], df = out[2, "Df diff"])
    }
    FIT <- eval(as.call(oldCall))
    ## check if there are any results
    noFit <- sapply(FIT@funList, function(x) x[1] == "fit failed")
    noLRT <- sapply(FIT@funList, function(x) x[1] == "lavTestLRT() failed")
    if (all(noFit | noLRT)) stop("No success using lavTestScore() on any imputations.")

    chiList <- sapply(FIT@funList[useImps & !(noFit | noLRT)], "[[", i = "chisq")
    dfList <- sapply(FIT@funList[useImps & !(noFit | noLRT)], "[[", i = "df")
    out <- calculate.D2(chiList, DF = mean(dfList), asymptotic)
    names(out) <- paste0(names(out), ".scaled")
    class(out) <- c("lavaan.vector","numeric")
    return(out)
  }
  ## else, return model fit OR naive difference test to be robustified


  test <- if (pool.robust) 2L else 1L
  ## pool Wald tests
  if (is.null(h1)) {
    DF <- mean(sapply(object@testList[useImps], function(x) x[[test]][["df"]]))
    w <- sapply(object@testList[useImps], function(x) x[[test]][["stat"]])
  } else {
    ## this will not get run if !pool.robust because logic catches that first
    DF0 <- mean(sapply(object@testList[useImps], function(x) x[[1]][["df"]]))
    DF1 <- mean(sapply(h1@testList[useImps], function(x) x[[1]][["df"]]))
    DF <- DF0 - DF1
    w0 <- sapply(object@testList[useImps], function(x) x[[1]][["stat"]])
    w1 <- sapply(h1@testList[useImps], function(x) x[[1]][["stat"]])
    w <- w0 - w1
  }
  out <- calculate.D2(w, DF, asymptotic)
  ## add .scaled suffix
  if (pool.robust) names(out) <- paste0(names(out), ".scaled")
  ## for 1 model, add extra info (redundant if pool.robust)
  if (is.null(h1) & !pool.robust) {
    PT <- parTable(object)
    out <- c(out, npar = max(PT$free) - sum(PT$op == "=="),
             ntotal = lavListInspect(object, "ntotal"))
  }

  class(out) <- c("lavaan.vector","numeric")
  out
}
#' @importFrom lavaan parTable lavaan lavListInspect
#' @importFrom methods getMethod
getLLs <- function(object, saturated = FALSE) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  ## FIXME: lavaanList does not return info when fixed because no convergence!
  dataList <- object@DataList[useImps]
  lavoptions <- lavListInspect(object, "options")
  group <- lavListInspect(object, "group")
  if (length(group) == 0L) group <- NULL
  if (saturated) {
    fit <- lavaan(parTable(object), data = dataList[[ which(useImps)[1] ]],
                  slotOptions = lavoptions, group = group)
    ## use saturated parameter table as new model
    PT <- lavaan::lav_partable_unrestricted(fit)
    ## fit saturated parameter table to each imputation, return estimates
    satParams <- lapply(object@DataList[useImps], function(d) {
      parTable(lavaan(model = PT, data = d,
                      slotOptions = lavoptions, group = group))$est
    })
    ## set all parameters fixed
    PT$free <- 0L
    PT$user <- 1L
    ## fix them to pooled estimates
    PT$ustart <- colMeans(do.call(rbind, satParams))
    PT$start <- NULL
    PT$est <- NULL
    PT$se <- NULL
  } else {
    ## save parameter table as new model
    PT <- parTable(object)
    ## set all parameters fixed
    PT$free <- 0L
    PT$user <- 1L
    ## fix them to pooled estimates
    fixedValues <- getMethod("coef","lavaan.mi")(object, type = "user")
    PT$ustart <- fixedValues
    PT$start <- NULL
    PT$est <- NULL
    PT$se <- NULL
    ## omit (in)equality constraints and user-defined parameters
    params <- !(PT$op %in% c("==","<",">",":="))
    PT <- PT[params, ]
  }
  ## return log-likelihoods
  sapply(object@DataList[useImps], function(d) {
    lavaan::logLik(lavaan(PT, data = d, slotOptions = lavoptions, group = group))
  })
}
#' @importFrom stats pf pchisq
#' @importFrom lavaan lavListInspect parTable
D3 <- function(object, h1 = NULL, asymptotic = FALSE) {
  N <- lavListInspect(object, "ntotal")
  useImps <- sapply(object@convergence, "[[", i = "converged")
  nImps <- sum(useImps)
  # m <- length(object@testList)
  if (is.null(h1)) {
    DF <- object@testList[[ which(useImps)[1] ]][[1]][["df"]]
  } else {
    DF1 <- h1@testList[[ which(useImps)[1] ]][[1]][["df"]]
    DF0 <- object@testList[[ which(useImps)[1] ]][[1]][["df"]]
    DF <- DF0 - DF1
  }

  ## calculate m log-likelihoods under pooled H0 estimates
  LL0 <- getLLs(object)
  ## calculate m log-likelihoods under pooled H1 estimates
  LL1 <- if (is.null(h1)) getLLs(object, saturated = TRUE) else getLLs(h1)
  #FIXME: check whether LL1 or LL0 returned errors?  add try()?

  ## calculate average of m LRTs
  LRT_con <- mean(-2*(LL0 - LL1)) # getLLs() already applies [useImps]
  ## average chisq across imputations
  if (is.null(h1)) {
    LRT_bar <- mean(sapply(object@testList[useImps], function(x) x[[1]]$stat))
  } else {
    LRT_bar <- mean(sapply(object@testList[useImps], function(x) x[[1]]$stat) -
                      sapply(h1@testList[useImps], function(x) x[[1]]$stat))
  }
  ## calculate average relative increase in variance
  a <- DF*(nImps - 1)
  ariv <- ((nImps + 1) / a) * (LRT_bar - LRT_con)
  test.stat <- LRT_con / (DF*(1 + ariv))
  if (is.na(test.stat)) stop('D3 test statistic could not be calculated. ',
                             'Try the D2 pooling method.') #FIXME: check whether model-implied Sigma is NPD
  if (test.stat < 0) {
    message('Negative test statistic set to zero \n')
    test.stat <- 0
  }
  if (asymptotic) {
    out <- c("chisq" = test.stat * DF, df = DF,
             pvalue = pchisq(test.stat * DF, df = DF, lower.tail = FALSE))
  } else {
    ## F statistic
    if (a > 4) {
      v4 <- 4 + (a - 4) * (1 + (1 - (2 / a))*(1 / ariv))^2 # Enders (eq. 8.34)
    } else {
      v4 <- a*(1 + 1/DF)*(1 + 1/ariv)^2 / 2 # Enders (eq. 8.35)
      # v4 <- (DF + 1)*(m - 1)*(1 + (1 / ariv))^2 / 2 # Grund et al. (eq. 9)
    }
    out <- c("F" = test.stat, df1 = DF, df2 = v4,
             pvalue = pf(test.stat, df1 = DF, df2 = v4, lower.tail = FALSE))
  }
  ## add log-likelihood and AIC/BIC for target model
  if (is.null(h1)) {
    PT <- parTable(object)
    npar <- max(PT$free) - sum(PT$op == "==")
    out <- c(out, npar = npar, ntotal = lavListInspect(object, "ntotal"),
             logl = mean(LL0), unrestricted.logl = mean(LL1),
             aic = -2*mean(LL0) + 2*npar, bic = -2*mean(LL0) + npar*log(N),
             bic2 = -2*mean(LL0) + npar*log((N + 2) / 24))
    ## NOTE: Mplus reports the average of m likelihoods evaluated at the
    ##       m point estimates, not evaluated at the pooled point estimates.
    ##       Mplus also uses those to calcluate AIC and BIC.
  }

  class(out) <- c("lavaan.vector","numeric")
  out
}
#' @importFrom stats pchisq
#' @importFrom lavaan lavListInspect
robustify <- function(ChiSq, object, h1 = NULL) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  scaleshift <- lavListInspect(object, "options")$test == "scaled.shifted"

  d0 <- mean(sapply(object@testList[useImps], function(x) x[[2]][["df"]]))
  c0 <- mean(sapply(object@testList[useImps],
                    function(x) x[[2]][["scaling.factor"]]))
  if (!is.null(h1)) {
    d1 <- mean(sapply(h1@testList[useImps], function(x) x[[2]][["df"]]))
    c1 <- mean(sapply(h1@testList[useImps],
                      function(x) x[[2]][["scaling.factor"]]))
    delta_c <- (d0*c0 - d1*c1) / (d0 - d1)
    ChiSq["chisq.scaled"] <- ChiSq[["chisq"]] / delta_c
    ChiSq["df.scaled"] <- d0 - d1
    ChiSq["pvalue.scaled"] <- pchisq(ChiSq[["chisq.scaled"]],
                                     df = ChiSq[["df.scaled"]],
                                     lower.tail = FALSE)
    ChiSq["chisq.scaling.factor"] <- delta_c
  } else {
    ChiSq["chisq.scaled"] <- ChiSq[["chisq"]] / c0
    ChiSq["df.scaled"] <- d0
    if (scaleshift) {
      ## add average shift parameter (or average of sums, if nG > 1)
      shift <- mean(sapply(object@testList[useImps],
                           function(x) sum(x[[2]][["shift.parameter"]]) ))
      ChiSq["chisq.scaled"] <- ChiSq[["chisq.scaled"]] + shift
      ChiSq["pvalue.scaled"] <- pchisq(ChiSq[["chisq.scaled"]],
                                       df = ChiSq[["df.scaled"]],
                                       lower.tail = FALSE)
      ChiSq["chisq.scaling.factor"] <- c0
      ChiSq["chisq.shift.parameters"] <- shift
    } else {
      ChiSq["pvalue.scaled"] <- pchisq(ChiSq[["chisq.scaled"]],
                                       df = ChiSq[["df.scaled"]],
                                       lower.tail = FALSE)
      ChiSq["chisq.scaling.factor"] <- c0
    }
  }
  ChiSq
}
#' @importFrom stats pchisq uniroot
#' @importFrom lavaan lavListInspect
anova.lavaan.mi <- function(object, h1 = NULL, test = c("D3","D2","D1"),
                            pool.robust = FALSE, scale.W = FALSE,
                            asymptotic = FALSE, constraints = NULL,
                            indices = FALSE, baseline.model = NULL,
                            method = "default", A.method = "delta",
                            scaled.shifted = TRUE, H1 = TRUE, type = "Chisq") {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  nImps <- sum(useImps)
  ## check class
  if (!inherits(object, "lavaan.mi")) stop("object is not class 'lavaan.mi'")
  if (!is.null(h1) & !inherits(object, "lavaan.mi")) stop("h1 is not class 'lavaan.mi'")
  test <- as.character(test[1])
  ## check test options, backward compatibility?
  if (tolower(test) == "mplus") {
    test <- "D3"
    asymptotic <- TRUE
  }
  if (tolower(test) %in% c("mr","meng.rubin","likelihood","lrt","d3")) test <- "D3"
  if (tolower(test) %in% c("lmrr","li.et.al","pooled.wald","d2")) test <- "D2"
  if (toupper(test) == "D3" & !lavListInspect(object, "options")$estimator %in% c("ML","PML","FML")) {
    message('"D3" only available using maximum likelihood estimation. ',
            'Changed test to "D2".')
    test <- "D2"
  }

  ## Everything else obsolete if test = "D1"
  if (toupper(test) == "D1") {
    out <- D1(object, constraints = constraints, scale.W = scale.W,
              asymptotic = asymptotic)
    message('D1 (Wald test) calculated using pooled "',
            lavListInspect(object, "options")$se,
            '" asymptotic covariance matrix of model parameters')
    return(out)
  }

  ## check for robust
  robust <- lavListInspect(object, "options")$test != "standard"
  if (robust & !pool.robust) {
    if (!asymptotic)
      message('Robust correction can only be applied to pooled chi-squared',
              ' statistic, not F statistic. "asymptotic" was switched to TRUE.')
    asymptotic <- TRUE
  }
  scaleshift <- lavListInspect(object, "options")$test == "scaled.shifted"
  if (scaleshift & !is.null(h1)) {
    if (test == "D3" | !pool.robust)
      message("If test = 'scaled.shifted' (estimator = 'WLSMV' or 'MLMV'), ",
              "model comparison is only available by (re)setting test = 'D2' ",
              "and pool.robust = TRUE.\n",
              "Control more options by passing arguments to lavTestLRT().\n")
    pool.robust <- TRUE
    test <- 'D2'
  }


  ## check request for fit indices
  if (is.null(h1)) {
    incremental <- c("cfi","tli","nnfi","rfi","nfi","pnfi","ifi","rni")
  } else {
    indices <- FALSE
    incremental <- c("")
  }
  if (is.logical(indices)) {
    moreFit <- is.null(h1) & indices
    if (moreFit) indices <- c("cfi","tli","rmsea","srmr")
  } else if (is.character(indices)) {
    indices <- tolower(indices)
    moreFit <- is.null(h1) & any(indices %in% c(incremental, "all","mfi","rmr",
                                                "srmr","rmsea","gammaHat"))
    if (moreFit & any(indices == "all")) {
      indices <- c(incremental, "mfi","rmsea","gammaHat","rmr","srmr")
    }
  } else indices <- moreFit <- FALSE
  ## fit baseline model if necessary
  if (moreFit & any(indices %in% incremental)) {
    if (is.null(baseline.model)) {
      PTb <- lavaan::lav_partable_independence(lavdata = object@Data,
                                               lavoptions = lavListInspect(object, "options"))
      # FIXME: shouldn't need this line, but lav_partable_merge() fails when
      #        lavaan:::lav_object_extended() returns a NULL slot instead of "plabel"
      PTb$plabel <- paste0(".p", PTb$id, ".")
      group <- lavListInspect(object, "group")
      if (length(group) == 0L) group <- NULL
      baseFit <- runMI(model = PTb, data = object@DataList[useImps],
                       group = group, se = "none", # to save time
                       test = lavListInspect(object, "options")$test,
                       estimator = lavListInspect(object, "options")$estimator,
                       ordered = lavListInspect(object, "ordered"),
                       parameterization = lavListInspect(object,
                                                         "parameterization"))
    } else if (!inherits(baseline.model, "lavaan.mi")) {
      stop('User-supplied baseline.model must be "lavaan.mi" class fit',
           ' to the same imputed data')
    } else baseFit <- baseline.model
    baseImps <- sapply(baseFit@convergence, "[[", i = "converged")
    if (!all(baseImps)) warning('baseline.model did not converge for data set(s): ',
                                which(useImps)[!baseImps])
  }

  ## check DF
  DF0 <- object@testList[[ which(useImps)[1] ]][[1]][["df"]]
  if (!is.null(h1)) {
    if (!inherits(h1, "lavaan.mi")) stop("h1 is not class 'lavaan.mi'")
    DF1 <- h1@testList[[ which(useImps)[1] ]][[1]][["df"]]
    if (DF0 == DF1) stop("models have equal degrees of freedom")
    if (DF0 < DF1) {
      H0 <- h1
      h1 <- object
      object <- H0
      H0 <- DF1
      DF1 <- DF0
      DF0 <- H0
    }
    DF <- DF0 - DF1
  } else DF <- DF0
  if (DF == 0) indices <- moreFit <- FALSE # arbitrary perfect fit, no indices
  if (moreFit) asymptotic <- TRUE

  ## calculate pooled test
  if (test == "D3") {
    if (pool.robust & moreFit) stop('pool.robust = TRUE only applicable ',
                                    'when test = "D2".')
    out <- D3(object, h1 = h1, asymptotic = asymptotic)
    if (any(indices %in% incremental)) baseOut <- D3(baseFit, asymptotic = TRUE)
  } else if (test == "D2") {
    out <- D2(object, h1 = h1, asymptotic = asymptotic, pool.robust = FALSE)
    if (any(indices %in% incremental)) baseOut <- D2(baseFit, asymptotic = TRUE,
                                                     pool.robust = FALSE)
    if (robust & pool.robust) {
      out <- c(out,
               D2(object, h1 = h1, asymptotic = asymptotic, pool.robust = TRUE,
                  method = method, A.method = A.method,
                  scaled.shifted = scaled.shifted, H1 = H1, type = type))
      if (any(indices %in% incremental)) {
        baseOut <- c(baseOut, D2(baseFit, asymptotic = TRUE, pool.robust = TRUE,
                                 method = method, A.method = A.method, H1 = H1,
                                 scaled.shifted = scaled.shifted, type = type))
      }
    }
  } else stop("'", test, "' is an invalid option for the 'test' argument.")
  ## If test statistic is negative, return without any indices or robustness
  if (asymptotic & (moreFit | robust)) {
    if (out[["chisq"]] == 0) {
      message('Negative test statistic set to zero, so fit will appear to be ',
              'arbitrarily perfect.  Robust corrections and additional fit ',
              'indices are not returned because they are uninformative.\n')
      class(out) <- c("lavaan.vector","numeric")
      return(out)
    }
  }

  ## If robust statistics were not pooled above, robustify naive statistics
  if (robust & !pool.robust) {
    out <- robustify(ChiSq = out, object, h1)
    if (scaleshift) {
      extraWarn <- ' and shift parameter'
    } else if (lavListInspect(object, "options")$test == "mean.var.adjusted") {
      extraWarn <- ' and degrees of freedom'
    } else extraWarn <- ''
    message('Robust corrections are made by pooling the naive chi-squared ',
            'statistic across ', nImps, ' imputations for which the model ',
            'converged, then applying the average (across imputations) scaling',
            ' factor', extraWarn, ' to that pooled value. \n',
            'To instead pool the robust test statistics, set test = "D2" and ',
            'pool.robust = TRUE. \n')
  }

  ## add fit indices for single model
  if (moreFit) {
    X2 <- out[["chisq"]]
    # if (pool.robust) message('All fit indices are calculated using the pooled',
    #                          ' robust test statistic. \n')
    if (robust) {
      X2.sc <- out[["chisq.scaled"]]
      DF.sc <- out[["df.scaled"]] ## for mean.var.adjusted, mean DF across imputations
      if (!pool.robust) ch <- out[["chisq.scaling.factor"]] ## mean c_hat across imputations
      if (X2 < .Machine$double.eps && DF == 0) ch <- 0
      ## for RMSEA
      if ("rmsea" %in% indices) {
        d <- mean(sapply(object@testList[useImps],
                         function(x) sum(x[[2]][["trace.UGamma"]])))
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
        if (!pool.robust) baseOut <- robustify(ChiSq = baseOut, object = baseFit)
        out["baseline.chisq.scaled"] <- bX2.sc <- baseOut[["chisq.scaled"]]
        out["baseline.df.scaled"]    <- bDF.sc <- baseOut[["df.scaled"]]
        out["baseline.pvalue.scaled"] <- baseOut[["pvalue.scaled"]]
        if (!pool.robust) {
          cb <- baseOut[["chisq.scaling.factor"]]
          out["baseline.chisq.scaling.factor"] <- cb
          if (scaleshift) {
            out["baseline.chisq.shift.parameters"] <- baseOut[["chisq.shift.parameters"]]
          }
        }
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
      if (!pool.robust & lavListInspect(object, "options")$test %in%
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
      if (!pool.robust & lavListInspect(object, "options")$test %in%
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
      if (!pool.robust & lavListInspect(object, "options")$test %in%
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
  Ns <- lavListInspect(object, "nobs")
  nG <- lavListInspect(object, "ngroups")
  nVars <- length(lavaan::lavNames(object))
  if (!(lavListInspect(object, "options")$likelihood == "normal" |
        lavListInspect(object, "options")$estimator %in% c("ML","PML","FML"))) {
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
        if (DF.sc < 1 | getLambda(0, X2, DF.sc, .95) < 0.0) {
          out["rmsea.ci.lower.scaled"] <- 0
        } else {
          lambda.l <- try(uniroot(f = getLambda, chi = X2, df = DF.sc, p = .95,
                                  lower = 0, upper = X2)$root, silent = TRUE)
          if (inherits(lambda.l, "try-error")) lambda.l <- NA
          out["rmsea.ci.lower.scaled"] <- sqrt( lambda.l/(N*DF) ) * sqrt(nG)
        }
        ## upper confidence limit
        if (DF.sc < 1 | getLambda(N.RMSEA, X2, DF.sc, .05) > 0.0) {
          out["rmsea.ci.upper.scaled"] <- 0
        } else {
          lambda.u <- try(uniroot(f = getLambda, chi = X2, df = DF.sc, p = .05,
                                  lower = 0, upper = N.RMSEA)$root, silent = TRUE)
          if (inherits(lambda.u, "try-error")) lambda.u <- NA
          out["rmsea.ci.upper.scaled"] <- sqrt( lambda.u/(N*DF) ) * sqrt(nG)
        }
        ## p value
        out["rmsea.pvalue.scaled"] <- pchisq(X2, DF.sc, ncp = N*DF.sc*0.05^2/nG,
                                             lower.tail = FALSE)

        if (!pool.robust & lavListInspect(object, "options")$test %in%
            c("satorra.bentler","yuan.bentler","yuan.bentler.mplus")) {
          ## robust
          out["rmsea.robust"] <- sqrt( max(0, (X2/N)/DF - ch/N ) ) * sqrt(nG)
          ## lower confidence limit
          if (DF.sc < 1 | getLambda(0, X2.sc, DF.sc, .95) < 0.0) {
            out["rmsea.ci.lower.robust"] <- 0
          } else {
            lambda.l <- try(uniroot(f = getLambda, chi = X2.sc, df = DF.sc, p = .95,
                                    lower = 0, upper = X2)$root, silent = TRUE)
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
      } else if (scaleshift) {
        ## naive only
        out["rmsea.scaled"] <- sqrt( max(0, (X2.sc/N)/DF - 1/N) ) * sqrt(nG)
        ## lower confidence limit
        if (DF.sc < 1 | getLambda(0, X2.sc, DF.sc, .95) < 0.0) {
          out["rmsea.ci.lower.scaled"] <- 0
        } else {
          lambda.l <- try(uniroot(f = getLambda, chi = X2.sc, df = DF.sc, p = .95,
                                  lower = 0, upper = X2.sc)$root, silent = TRUE)
          if (inherits(lambda.l, "try-error")) lambda.l <- NA
          out["rmsea.ci.lower.scaled"] <- sqrt( lambda.l/(N*DF.sc) ) * sqrt(nG)
        }
        ## upper confidence limit
        if (DF.sc < 1 | getLambda(N.RMSEA, X2.sc, DF.sc, .05) > 0.0) {
          out["rmsea.ci.upper.scaled"] <- 0
        } else {
          lambda.u <- try(uniroot(f = getLambda, chi = X2.sc, df = DF.sc, p = .05,
                                  lower = 0, upper = N.RMSEA)$root, silent = TRUE)
          if (inherits(lambda.u, "try-error")) lambda.u <- NA
          out["rmsea.ci.upper.scaled"] <- sqrt( lambda.u/(N*DF.sc) ) * sqrt(nG)
        }
        ## p value
        out["rmsea.pvalue.scaled"] <- pchisq(X2.sc, DF.sc, ncp = N*DF.sc*0.05^2/nG,
                                             lower.tail = FALSE)
      }
    }
  }

  if ("gammaHat" %in% indices) {
    out["gammaHat"] <- nVars / (nVars + 2*((X2 - DF) / N))
    out["adjGammaHat"] <- 1 - (((nG * nVars * (nVars + 1)) / 2) / DF) * (1 - out["gammaHat"])
    if (robust) {
      out["gammaHat.scaled"] <- nVars / (nVars + 2*((X2.sc - DF.sc) / N))
      out["adjGammaHat.scaled"] <- 1 - (((nG * nVars * (nVars + 1)) / 2) / DF.sc) * (1 - out["gammaHat.scaled"])
    }
  }

  getSRMR <- function(object, type) {
    vv <- lavaan::lavNames(object, type = "ov.num")
    R <- getMethod("resid", "lavaan.mi")(object, type = type)
    index <- if (type == "raw") "cov" else "cor"
    if (nG > 1L) {
      RR <- list()
      for (g in 1:nG) {
        RR[[g]] <- c(R[[g]][[index]][lower.tri(R[[g]][[index]], diag = FALSE)]^2,
                     diag(R[[g]][[index]])[vv]^2)
      }
    } else RR <- c(R[[index]][lower.tri(R[[index]], diag = FALSE)]^2,
                   diag(R[[index]])[vv]^2)

    if (lavListInspect(object, "meanstructure")) {
      if (nG > 1L) {
        for (g in 1:nG) RR[[g]] <- c(RR[[g]], R[[g]]$mean[vv]^2)
      } else RR <- c(RR, R$mean[vv]^2)
    }

    SS <- if (nG > 1L) sqrt(sapply(RR, mean)) else sqrt(mean(RR))
    as.numeric( (lavListInspect(object, "nobs") %*% SS) / lavListInspect(object, "ntotal") )
  }
  if("rmr" %in% indices) out["rmr"] <- getSRMR(object, type = "raw")
  if("srmr" %in% indices) {
    out["srmr_bollen"] <- getSRMR(object, type = "cor.bollen")
    out["srmr_bentler"] <- getSRMR(object, type = "cor.bentler")
  }

  class(out) <- c("lavaan.vector","numeric")
  out # FIXME: in future, accept more than 2 models, arrange sequentially by DF
}
#' @name lavaan.mi-class
#' @aliases anova,lavaan.mi-method
#' @export
setMethod("anova", "lavaan.mi", anova.lavaan.mi)


#' @name lavaan.mi-class
#' @aliases fitMeasures,lavaan.mi-method
#' @importFrom lavaan fitMeasures
#' @export
setMethod("fitMeasures", "lavaan.mi", function(object, fit.measures = "all",
                                               baseline.model = NULL) {
  if (!is.character(fit.measures)) stop("'fit.measures' must be a character ",
                                        "string specifying name(s) of desired ",
                                        "fit indices.")
  message('anova() provides more control over options for pooling chi-squared',
          ' before calculating fit indices from multiple imputations. ',
          'See the class?lavaan.mi help page for details.\n\n')
  fits <- anova.lavaan.mi(object, indices = "all", baseline.model = baseline.model)
  if ("all" %in% fit.measures) return(fits)
  out <- fits[grepl(paste(fit.measures, collapse = "|"),
                    names(fits), ignore.case = TRUE)]
  out <- out[which(!is.na(names(out)))]
  class(out) <- c("lavaan.vector","numeric")
  out
})
# lowercase 'm'
#' @name lavaan.mi-class
#' @aliases fitmeasures,lavaan.mi-method
#' @importFrom lavaan fitmeasures
#' @export
setMethod("fitmeasures", "lavaan.mi", function(object, fit.measures = "all",
                                               baseline.model = NULL) {
  if (!is.character(fit.measures)) stop("'fit.measures' must be a character ",
                                        "string specifying name(s) of desired ",
                                        "fit indices.")
  message('anova() provides more control over options for pooling chi-squared',
          ' before calculating fit indices from multiple imputations. ',
          'See the class?lavaan.mi help page for details.\n\n')
  fits <- anova.lavaan.mi(object, indices = "all", baseline.model = baseline.model)
  if ("all" %in% fit.measures) return(fits)
  out <- fits[grepl(paste(fit.measures, collapse = "|"),
                    names(fits), ignore.case = TRUE)]
  out <- out[which(!is.na(names(out)))]
  class(out) <- c("lavaan.vector","numeric")
  out
})


## function to pool each group's list of sample stats
sampstat.lavaan.mi <- function(lst, means = FALSE, categ = FALSE, m = m) {
  ## average sample stats across imputations
  out <- list(cov = Reduce("+", lapply(lst, "[[", i = "cov")) / m)
  if (means) out$mean <- Reduce("+", lapply(lst, "[[", i = "mean")) / m
  if (categ) out$th <- Reduce("+", lapply(lst, "[[", i = "th")) / m
  out
}
#' @importFrom lavaan lavListInspect
fitted.lavaan.mi <- function(object) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  m <- sum(useImps)
  meanstructure <- lavListInspect(object, "meanstructure")
  categ <- lavListInspect(object, "categorical")
  nG <- lavListInspect(object, "ngroups")
  ov.names <- lavaan::lavNames(object)

  est <- getMethod("coef", "lavaan.mi")(object)
  imp <- lavaan::lav_model_implied(lavaan::lav_model_set_parameters(object@Model,
                                                                    x = est))
  out <- list()
  if (nG > 1L) {
    group.label <- lavListInspect(object, "group.label")
    for (i in seq_along(imp)) names(imp[[i]]) <- group.label
    for (g in group.label) {
      out[[g]]$cov <- imp$cov[[g]]
      dimnames(out[[g]]$cov) <- list(ov.names, ov.names)
      class(out[[g]]$cov) <- c("lavaan.matrix.symmetric","matrix")
      if (meanstructure) {
        out[[g]]$mean <- as.numeric(imp$mean[[g]])
        names(out[[g]]$mean) <- ov.names
        class(out[[g]]$mean) <- c("lavaan.vector","numeric")
      } else {
        out[[g]]$mean <- sampstat.lavaan.mi(lapply(object@SampleStatsList[useImps], "[[", g),
                                            means = TRUE, categ = categ, m = m)$mean
      }
      if (categ) {
        out[[g]]$th <- imp$th[[g]]
        names(out[[g]]$th) <- lavaan::lavNames(object, "th")
        class(out[[g]]$th) <- c("lavaan.vector","numeric")
      }
    }
  } else {
    out$cov <- imp$cov[[1]]
    dimnames(out$cov) <- list(ov.names, ov.names)
    class(out$cov) <- c("lavaan.matrix.symmetric","matrix")
    if (meanstructure) {
      out$mean <- as.numeric(imp$mean[[1]])
      names(out$mean) <- ov.names
      class(out$mean) <- c("lavaan.vector","numeric")
    } else {
      out$mean <- sampstat.lavaan.mi(object@SampleStatsList[useImps],
                                     means = TRUE, categ = categ, m = m)$mean
    }
    if (categ) {
      out$th <- imp$th[[1]]
      names(out$th) <- lavaan::lavNames(object, "th")
      class(out$th) <- c("lavaan.vector","numeric")
    }
  }
  out
}
#' @name lavaan.mi-class
#' @aliases fitted,lavaan.mi-method
#' @export
setMethod("fitted", "lavaan.mi", fitted.lavaan.mi)
#' @name lavaan.mi-class
#' @aliases fitted.values,lavaan.mi-method
#' @export
setMethod("fitted.values", "lavaan.mi", fitted.lavaan.mi)



## function to calculate residuals for one group
#' @importFrom stats cov2cor
gp.resid.lavaan.mi <- function(Observed, N, Implied, type,
                               means = FALSE, categ = FALSE, m) {
  obsMats <- sampstat.lavaan.mi(Observed, means = means, categ = categ, m = m)
  ## average sample stats across imputations
  S_mean <- if (is.null(N)) obsMats$cov else (obsMats$cov * ((N - 1L) / N))
  if (means) M_mean <- obsMats$mean
  if (categ) Th_mean <- obsMats$th

  if (type == "raw") {
    out <- list(cov = S_mean - Implied$cov)
    if (means) out$mean <- M_mean - Implied$mean else {
      out$mean <- rep(0, nrow(out$cov))
      names(out$mean) <- rownames(out$cov)
    }
    if (categ) out$th <- Th_mean - Implied$th
    return(out)
  } else if (type == "cor.bollen") {
    out <- list(cor = cov2cor(S_mean) - cov2cor(Implied$cov))
    if (!means) {
      out$mean <- rep(0, nrow(out$cor))
      names(out$mean) <- rownames(out$cor)
    } else {
      std.obs.M <- M_mean / sqrt(diag(S_mean))
      std.mod.M <- Implied$mean / sqrt(diag(Implied$cov))
      out$mean <- std.obs.M - std.mod.M
    }
  } else if (type == "cor.bentler") {
    SDs <- diag(sqrt(diag(S_mean)))
    dimnames(SDs) <- dimnames(S_mean)
    out <- list(cor = solve(SDs) %*% (S_mean - Implied$cov) %*% solve(SDs))
    class(out$cor) <- c("lavaan.matrix.symmetric","matrix")
    if (!means) {
      out$mean <- rep(0, nrow(out$cor))
      names(out$mean) <- rownames(out$cor)
    } else out$mean <- (M_mean - Implied$mean) / diag(SDs)
  } else stop("argument 'type' must be 'raw', 'cor', 'cor.bollen', ",
              "or 'cor.bentler'.")
  if (categ) out$th <- Th_mean - Implied$th
  out
}
#' @importFrom lavaan lavListInspect
resid.lavaan.mi <- function(object, type = c("raw","cor")) {
  ## @SampleStatsList is (for each imputation) output from:
  ##    getSampStats <- function(obj) lavInspect(obj, "sampstat")
  useImps <- sapply(object@convergence, "[[", i = "converged")
  m <- sum(useImps)
  rescale <- lavListInspect(object, "options")$sample.cov.rescale
  meanstructure <- lavListInspect(object, "meanstructure")
  categ <- lavListInspect(object, "categorical")
  type <- tolower(type[1])
  ## check for type = "cor" ("cor.bollen") or "cor.bentler"
  if (type == "cor") type <- "cor.bollen"
  ## model-implied moments, already pooled
  Implied <- getMethod("fitted", "lavaan.mi")(object)
  ## Calculate residuals
  nG <- lavListInspect(object, "ngroups")
  N <- lavListInspect(object, "nobs")
  if (nG > 1L) {
    group.label <- names(Implied)
    if (is.null(group.label)) group.label <- 1:length(Implied) else names(N) <- group.label
    out <- list()
    for (g in group.label) {
      out[[g]] <- gp.resid.lavaan.mi(Observed = lapply(object@SampleStatsList[useImps], "[[", g),
                                     N = if (rescale) N[g] else NULL,
                                     Implied = Implied[[g]], type = type,
                                     means = meanstructure, m = m, categ = categ)
    }
  } else {
    out <- gp.resid.lavaan.mi(Observed = object@SampleStatsList[useImps],
                              N = if (rescale) N else NULL,
                              Implied = Implied, type = type,
                              means = meanstructure, m = m, categ = categ)
  }
  out
}
#' @name lavaan.mi-class
#' @aliases residuals,lavaan.mi-method
#' @export
setMethod("residuals", "lavaan.mi", resid.lavaan.mi)
#' @name lavaan.mi-class
#' @aliases resid,lavaan.mi-method
#' @export
setMethod("resid", "lavaan.mi", resid.lavaan.mi)



