### Terrence D. Jorgensen
### Last updated: 15 September 2018
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
##' @slot GLIST pooled \code{list} of coefficients in GLIST format
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
##'
##' @param object An object of class \code{lavaan.mi}
##' @param se,ci,level,standardized,rsquare,header,add.attributes See
##'        \code{\link[lavaan]{parameterEstimates}}.
##' @param fmi \code{logical} indicating whether to include the Fraction Missing
##'        Information (FMI) for parameter estimates in the \code{summary}
##'        output (see \bold{Value} section).
##' @param asymptotic \code{logical}. If \code{FALSE} (typically a default, but
##'       see \bold{Value} section for details using various methods), pooled
##'       tests (of fit or pooled estimates) will be \emph{F} or \emph{t}
##'       statistics with associated degrees of freedom (\emph{df}). If
##'       \code{TRUE}, the (denominator) \emph{df} are assumed to be sufficiently
##'       large for a \emph{t} statistic to follow a normal distribution, so it
##'       is printed as a \emph{z} statisic; likewise, \emph{F} times its
##'       numerator \emph{df} is printed, assumed to follow a \eqn{\chi^2}
##'       distribution.
##' @param scale.W \code{logical}. If \code{TRUE} (default), the \code{vcov}
##'       method will calculate the pooled covariance matrix by scaling the
##'       within-imputation component by the ARIV (see Enders, 2010, p. 235,
##'       for definition and formula). Otherwise, the pooled matrix is
##'       calculated as the weighted sum of the within-imputation and
##'       between-imputation components (see Enders, 2010, ch. 8, for details).
##'       This in turn affects how the \code{summary} method calcualtes its
##'       pooled standard errors, as well as the Wald test
##'       (\code{\link{lavTestWald.mi}}).
##' @param labels \code{logical} indicating whether the \code{coef} output
##'        should include parameter labels. Default is \code{TRUE}.
##' @param total \code{logical} (default: \code{TRUE}) indicating whether the
##'        \code{nobs} method should return the total sample size or (if
##'        \code{FALSE}) a vector of group sample sizes.
##' @param type The meaning of this argument varies depending on which method it
##'        it used for. Find detailed descriptions in the \bold{Value} section
##'        under \code{coef}, \code{vcov}, and \code{residuals}.
##' @param fit.measures,baseline.model See \code{\link[lavaan]{fitMeasures}}.
##' @param ... Additional arguments passed to \code{\link{lavTestLRT.mi}}, or
##'   subsequently to \code{\link[lavaan]{lavTestLRT}}.
##'
##' @return
##'
##' \item{coef}{\code{signature(object = "lavaan.mi", type = "free", labels = TRUE)}:
##'   See \code{\linkS4class{lavaan}}. Returns the pooled point estimates (i.e.,
##'   averaged across imputed data sets; see Rubin, 1987).}
##'
##' \item{vcov}{\code{signature(object = "lavaan.mi", scale.W = TRUE,
##'   type = c("pooled","between","within","ariv"))}:  By default, returns the
##'   pooled covariance matrix of parameter estimates (\code{type = "pooled"}),
##'   the within-imputations covariance matrix (\code{type = "within"}), the
##'   between-imputations covariance matrix (\code{type = "between"}), or the
##'   average relative increase in variance (\code{type = "ariv"}) due to
##'   missing data.}
##'
##' \item{fitted.values}{\code{signature(object = "lavaan.mi")}: See
##'   \code{\linkS4class{lavaan}}. Returns model-implied moments, evaluated at
##'   the pooled point estimates.}
##' \item{fitted}{\code{signature(object = "lavaan.mi")}:
##'   alias for \code{fitted.values}}
##'
##' \item{residuals}{\code{signature(object = "lavaan.mi", type = c("raw","cor"))}:
##'   See \code{\linkS4class{lavaan}}. By default (\code{type = "raw"}), returns
##'   the difference between the model-implied moments from \code{fitted.values}
##'   and the pooled observed moments (i.e., averaged across imputed data sets).
##'   Standardized residuals are also available, using Bollen's
##'   (\code{type = "cor"} or \code{"cor.bollen"}) or Bentler's
##'   (\code{type = "cor.bentler"}) formulas.}
##' \item{resid}{\code{signature(object = "lavaan.mi", type = c("raw","cor"))}:
##'   alias for \code{residuals}}
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
##'   fit.measures = "all", baseline.model = NULL)}: See lavaan's
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
##'  scale.W = FALSE, asymptotic = FALSE, add.attributes = TRUE)}: see
##'  \code{\link[lavaan]{parameterEstimates}} for details.
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
##'  \code{add.attributes = FALSE}, the returned \code{data.frame} will resemble
##'  the \code{parameterEstimates} output. The \code{scale.W} argument is
##'  passed to \code{vcov} (see description above).}
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
##'   \url{www.statmodel.com}
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
##'   multiply-imputed data sets. \emph{Biometrika, 79}(1), 103--111. Retrieved
##'   from \url{https://www.jstor.org/stable/2337151}
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
                   GLIST = "list",        # list of pooled coefs in GLIST format
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
summary.lavaan.mi <- function(object, se = TRUE, ci = FALSE, level = .95,
                              standardized = FALSE, rsquare = FALSE,
                              fmi = FALSE, header = TRUE, scale.W = TRUE,
                              asymptotic = FALSE, add.attributes = TRUE) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  m <- sum(useImps)
  ## extract parameter table with attributes for printing
  PT <- parTable(object)
  myCols <- c("lhs","op","rhs","exo")
  if (lavListInspect(object, "ngroups") > 1L) myCols <- c(myCols,"block","group")
  PE <- PT[ , myCols]
  free <- PT$free > 0L | PT$op == ":="
  STDs <- !(PT$op %in% c("==","<",">")) # which rows can be standardized

  # PE$est <- rowMeans(sapply(object@ParTableList[useImps], "[[", i = "est"))
  PE$est <- getMethod("coef","lavaan.mi")(object, type = "all")

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
    #FIXME: no longer needed?  PE$exo <- 0L
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
  N <- lavListInspect(object, "norig")
  if (length(N) > 1L) names(N) <- lavListInspect(object, "group.label")
  N
})



##' @importFrom lavaan parTable
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
##' @name lavaan.mi-class
##' @aliases coef,lavaan.mi-method
##' @export
setMethod("coef", "lavaan.mi", coef.lavaan.mi)



##' @importFrom stats cov
##' @importFrom lavaan lavListInspect parTable
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
##' @name lavaan.mi-class
##' @aliases vcov,lavaan.mi-method
##' @export
setMethod("vcov", "lavaan.mi", vcov.lavaan.mi)


##' @importFrom lavaan lavListInspect lavTestLRT
anova.lavaan.mi <- function(object, ...) {

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
    } else mods <- NULL
    LRT.names <- intersect(names(dots),
                           union(names(formals(lavTestLRT)),
                                 names(formals(lavTestLRT.mi))))
    dots <- if (length(LRT.names)) dots[LRT.names] else NULL
    if (!is.null(dots$h1)) {
      mods <- c(mods, list(dots$h1))
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
    argList <- c(list(object), mods, list(argsLRT = dots, indices = FALSE))
    out <- do.call(compareFit, argList)
    results <- getMethod("summary", "FitDiff")(out)$test.statistics
  }

  results
}
##' @name lavaan.mi-class
##' @aliases anova,lavaan.mi-method
##' @export
setMethod("anova", "lavaan.mi", anova.lavaan.mi)

##' @importFrom lavaan lavNames
##' @importFrom stats pchisq uniroot
fitMeasures.mi <- function(object, fit.measures = "all", #FIXME: lavaan's generic needs "..."
                           baseline.model = NULL) {

  useImps <- sapply(object@convergence, "[[", i = "converged")
  lavoptions <- lavListInspect(object, "options")
  robust <- lavoptions$test != "standard" #TODO: check for bootstrap test
  scaleshift <- lavoptions$test == "scaled.shifted"

  if (!is.character(fit.measures)) stop("'fit.measures' must be a character ",
                                        "string specifying name(s) of desired ",
                                        "fit indices.")
  if (length(fit.measures) == 0L) fit.measures <- "all"
  ## narrow down fit indices
  incremental <- c("cfi","tli","nnfi","rfi","nfi","pnfi","ifi","rni")
  if ("all" %in% tolower(fit.measures)) {
    indices <- c(incremental, "mfi","rmsea","gammaHat","rmr")
  } else {
    indices <- grep(pattern = paste(c(incremental, "mfi","rmsea",
                                      "gammaHat","rmr"), collapse = "|"),
                    x = fit.measures, ignore.case = TRUE, value = TRUE)
  }

  ## check for additional arguments
  dots <- NULL #FIXME: list(...) once Yves accepts pull request
  if (length(dots)) {
    LRT.names <- intersect(names(dots),
                           union(names(formals(lavTestLRT)),
                                 names(formals(lavTestLRT.mi))))
    dots <- if (length(LRT.names)) dots[LRT.names] else list(asymptotic = TRUE)
  } else dots <- list(asymptotic = TRUE)
  if (robust) {
    if (is.null(dots$pool.robust)) {
      pool.robust <- formals(lavTestLRT.mi)$pool.robust # default value
    } else {
      pool.robust <- dots$pool.robust # user-specified value
    }
  }

  ## pooled test statistic(s)
  argList <- c(list(object = object), dots) #FIXME: make sure asymptotic = TRUE
  out <- do.call(lavTestLRT.mi, argList)

  ## fit baseline model if necessary
  if (any(indices %in% incremental)) {
    if (inherits(baseline.model, "lavaan.mi")) {
      baseFit <- baseline.model
    } else if (inherits(object@external$baseline.model, "lavaan.mi")) {
      baseFit <- object@external$baseline.model
    } else {
      PTb <- lavaan::lav_partable_independence(lavdata = object@Data,
                                               lavoptions = lavoptions)
      # FIXME: shouldn't need this line, but lav_partable_merge() fails when
      #        lavaan:::lav_object_extended() returns a NULL slot instead of "plabel"
      PTb$plabel <- paste0(".p", PTb$id, ".")
      group <- lavListInspect(object, "group")
      if (length(group) == 0L) group <- NULL
      baseFit <- runMI(model = PTb, data = object@DataList[useImps],
                       group = group, se = "none", # to save time
                       test = lavoptions$test, estimator = lavoptions$estimator,
                       ordered = lavListInspect(object, "ordered"),
                       parameterization = lavoptions$parameterization)
    }

    baseImps <- sapply(baseFit@convergence, "[[", i = "converged")
    if (!all(baseImps)) warning('baseline.model did not converge for data set(s): ',
                                which(useImps)[!baseImps])
  }

  ## pooled test statistic(s) for baseline model
  if (any(indices %in% incremental)) {
    argList <- c(list(object = baseFit), dots)
    baseOut <- do.call(lavTestLRT.mi, argList)
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
      if (!pool.robust & lavoptions$test %in%
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
      if (!pool.robust & lavoptions$test %in%
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
      if (!pool.robust & lavoptions$test %in%
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
  nlevels <- object@Data@nlevels #FIXME: lavListInspect(object, "nlevels")
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

        if (!pool.robust & lavoptions$test %in%
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

  getSRMR <- function(object, type, level = "within") {
    meanstructure <- lavListInspect(object, "meanstructure")
    N <- lavListInspect(object, "ntotal")
    nG <- lavListInspect(object, "ngroups")
    nlevels <- object@Data@nlevels #FIXME: lavListInspect(object, "nlevels")
    #TODO: ov.names(.x) should account for conditional.x (res.cov, res.int, etc.)

    R <- getMethod("resid", "lavaan.mi")(object, type = type)
    index <- if (type == "raw") "cov" else "cor"
    include.diag <- type != "cor.bollen"

    if (nG > 1L) {
      vv.g <- object@Data@ov.names #FIXME: assumes never nG > 1 && nlevels > 1
      RR <- list()
      for (g in 1:nG) {
        vv <- vv.g[[g]]
        RR[[g]] <- R[[g]][[index]][lower.tri(R[[g]][[index]], diag = FALSE)]^2
        if (include.diag)  RR[[g]] <- c(RR[[g]], diag(R[[g]][[index]])[vv]^2)
        if (meanstructure) RR[[g]] <- c(RR[[g]], R[[g]]$mean[vv]^2)
      }
      n.per.group <- lavListInspect(object, "nobs")

    } else if (nlevels > 1L) { #FIXME: needs to allow multiple levels & groups
      vv.l <- object@Data@ov.names.l[[1]]
      names(vv.l) <- c("within", lavListInspect(object, "cluster")) #FIXME: only works for 2 levels
      vv <- vv.l[[level]]

      RR <- R[[level]][[index]][lower.tri(R[[level]][[index]], diag = FALSE)]^2
      if (include.diag)  RR <- c(RR, diag(R[[level]][[index]])[vv]^2)
      if (meanstructure) RR <- c(RR, R[[level]]$mean[vv]^2)
      nclusters <- object@Data@Lp[[1]]$nclusters #TODO: add this to lavInspect
      names(nclusters) <- c("within", lavListInspect(object, "cluster")) #FIXME: only works for 2 levels
      n.per.group <- nclusters[[level]]

    } else {
      vv <- lavNames(object, type = "ov.num") #FIXME: not "ov" because always ignore correlation==1?
      RR <- R[[index]][lower.tri(R[[index]], diag = FALSE)]^2
      RR <- c(RR, diag(R[[index]])[vv]^2)
      if (meanstructure) RR <- c(RR, R$mean[vv]^2)
      n.per.group <- 1L
    }

    SS <- if (nG > 1L) sqrt(sapply(RR, mean)) else sqrt(mean(RR))
    as.numeric( (n.per.group %*% SS) / N )
  }

  if (any(c("rmr","srmr","crmr") %in% indices)) {
    if (nlevels > 1L) {
      out["srmr"] <- NA # to preserve the order in lavaan output
      out["srmr_within"] <- getSRMR(object, type = "cor", level = "within")
      out["srmr_between"] <- getSRMR(object, type = "cor",
                                     level = lavListInspect(object, "cluster"))
      out["srmr"] <- out["srmr_within"] + out["srmr_between"]
    } else {
      out["rmr"] <- getSRMR(object, type = "raw")
      out["crmr"] <- getSRMR(object, type = "cor.bollen")
      out["srmr"] <- getSRMR(object, type = "cor.bentler")
    }
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


## function to pool each group's list of sample stats
sampstat.lavaan.mi <- function(lst, means = FALSE, categ = FALSE, m = m) {
  ## average sample stats across imputations
  out <- list(cov = Reduce("+", lapply(lst, "[[", i = "cov")) / m)
  if (means) out$mean <- Reduce("+", lapply(lst, "[[", i = "mean")) / m
  if (categ) out$th <- Reduce("+", lapply(lst, "[[", i = "th")) / m
  #TODO: add others for conditional.x (e.g., slopes) ONLY if necessary
  out
}
##' @importFrom lavaan lavListInspect lavNames
fitted.lavaan.mi <- function(object) {
  useImps <- sapply(object@convergence, "[[", i = "converged")
  m <- sum(useImps)
  meanstructure <- lavListInspect(object, "meanstructure")
  categ <- lavListInspect(object, "categorical")
  nG <- lavListInspect(object, "ngroups")
  nlevels <- object@Data@nlevels #FIXME: lavListInspect(object, "nlevels")
  #TODO: account for fixed.x and conditional.x (res.cov, res.int, etc.)
  if (nG > 1L) {
    ov.names <- object@Data@ov.names #FIXME: assumes never nG > 1 && nlevels > 1
  } else if (nlevels > 1L) {
    ov.names <- object@Data@ov.names.l[[1]] # for first group (implies levels within groups?)
  } else ov.names <- lavNames(object)

  est <- getMethod("coef", "lavaan.mi")(object)
  imp <- lavaan::lav_model_implied(lavaan::lav_model_set_parameters(object@Model,
                                                                    x = est))

  #TODO: adapt to multilevel, multigroup, or both

  out <- list()
  if (nG > 1L || nlevels > 1L) {
    #FIXME: assumes never nG > 1 && nlevels > 1
    if (nG > 1L) {
      group.label <- lavListInspect(object, "group.label")
    }
    if (nlevels > 1L) {
      group.label <- c("within", lavListInspect(object, "cluster")) #FIXME: only works for 2 levels
    }
    names(ov.names) <- group.label
    for (i in seq_along(imp)) names(imp[[i]]) <- group.label
    for (g in group.label) {
      out[[g]]$cov <- imp$cov[[g]]
      dimnames(out[[g]]$cov) <- list(ov.names[[g]], ov.names[[g]])
      class(out[[g]]$cov) <- c("lavaan.matrix.symmetric","matrix")
      if (meanstructure) {
        out[[g]]$mean <- as.numeric(imp$mean[[g]])
        names(out[[g]]$mean) <- ov.names[[g]]
        class(out[[g]]$mean) <- c("lavaan.vector","numeric")
      }
      #TODO: omit "else" to match new lavaan::fitted output, which excludes $mean if !meanstructure
      # else {
      #   out[[g]]$mean <- sampstat.lavaan.mi(lapply(object@SampleStatsList[useImps], "[[", g),
      #                                       means = TRUE, categ = categ, m = m)$mean
      # }
      if (categ) {
        out[[g]]$th <- imp$th[[g]]
        names(out[[g]]$th) <- lavNames(object, "th")
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
    }
    #TODO: omit "else" to match new lavaan::fitted output, which excludes $mean if !meanstructure
    # else {
    #   out$mean <- sampstat.lavaan.mi(object@SampleStatsList[useImps],
    #                                  means = TRUE, categ = categ, m = m)$mean
    # }
    if (categ) {
      out$th <- imp$th[[1]]
      names(out$th) <- lavNames(object, "th")
      class(out$th) <- c("lavaan.vector","numeric")
    }
  }
  out
}
##' @name lavaan.mi-class
##' @aliases fitted,lavaan.mi-method
##' @export
setMethod("fitted", "lavaan.mi", fitted.lavaan.mi)
##' @name lavaan.mi-class
##' @aliases fitted.values,lavaan.mi-method
##' @export
setMethod("fitted.values", "lavaan.mi", fitted.lavaan.mi)



## function to calculate residuals for one group
##' @importFrom stats cov2cor
gp.resid.lavaan.mi <- function(Observed, N, Implied, type,
                               means = FALSE, categ = FALSE, m) {
  obsMats <- sampstat.lavaan.mi(Observed, means = means, categ = categ, m = m)
  ## average sample stats across imputations
  S_mean <- if (is.null(N)) obsMats$cov else (obsMats$cov * ((N - 1L) / N))
  if (means) M_mean <- obsMats$mean
  if (categ) Th_mean <- obsMats$th

  if (type == "raw") {
    out <- list(cov = S_mean - Implied$cov)
    if (means) out$mean <- M_mean - Implied$mean
    #TODO: omit "else" to match new lavaan::fitted output, which excludes $mean if !meanstructure
    # else {
    #   out$mean <- rep(0, nrow(out$cov))
    #   names(out$mean) <- rownames(out$cov)
    # }
    if (categ) out$th <- Th_mean - Implied$th
    return(out)
  } else if (type == "cor.bollen") {
    out <- list(cor = cov2cor(S_mean) - cov2cor(Implied$cov))
    if (means) {
      std.obs.M <- M_mean / sqrt(diag(S_mean))
      std.mod.M <- Implied$mean / sqrt(diag(Implied$cov))
      out$mean <- std.obs.M - std.mod.M
    }
  } else if (type == "cor.bentler") {
    SDs <- diag(sqrt(diag(S_mean)))
    dimnames(SDs) <- dimnames(S_mean)
    out <- list(cor = solve(SDs) %*% (S_mean - Implied$cov) %*% solve(SDs))
    class(out$cor) <- c("lavaan.matrix.symmetric","matrix")
    if (means) out$mean <- (M_mean - Implied$mean) / diag(SDs)
  } else stop("argument 'type' must be 'raw', 'cor', 'cor.bollen', ",
              "or 'cor.bentler'.")
  if (categ) out$th <- Th_mean - Implied$th
  out
}
##' @importFrom lavaan lavListInspect
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
  nlevels <- object@Data@nlevels #FIXME: lavListInspect(object, "nlevels")
  if (nG > 1L || nlevels > 1L) {
    if (nG > 1L) {
      group.label <- lavListInspect(object, "group.label")
      if (rescale) {
        N <- lavListInspect(object, "nobs")
        names(N) <- group.label
      } else N <- NULL
    }
    if (nlevels > 1L) {
      group.label <- c("within", lavListInspect(object, "cluster")) #FIXME: only works for 2 levels
      N <- NULL #FIXME: likelihood="wishart" does not change chisq in 0.6-3.1297
    }
    out <- list()
    for (g in group.label) {
      out[[g]] <- gp.resid.lavaan.mi(Observed = lapply(object@SampleStatsList[useImps], "[[", g),
                                     N = N, Implied = Implied[[g]], type = type,
                                     means = meanstructure, m = m, categ = categ)
    }
  } else {
    out <- gp.resid.lavaan.mi(Observed = object@SampleStatsList[useImps],
                              N = N, Implied = Implied, type = type,
                              means = meanstructure, m = m, categ = categ)
  }
  out
}
##' @name lavaan.mi-class
##' @aliases residuals,lavaan.mi-method
##' @export
setMethod("residuals", "lavaan.mi", resid.lavaan.mi)
##' @name lavaan.mi-class
##' @aliases resid,lavaan.mi-method
##' @export
setMethod("resid", "lavaan.mi", resid.lavaan.mi)



