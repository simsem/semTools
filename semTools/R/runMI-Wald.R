### Terrence D. Jorgensen & Yves Rosseel
### Last updated: 10 January 2021
### Pooled Wald test for multiple imputations
### Borrowed source code from lavaan/R/lav_test_Wald.R


##' Wald Test for Multiple Imputations
##'
##' Wald test for testing a linear hypothesis about the parameters of lavaan
##' models fitted to multiple imputed data sets. Statistics for constraining
##' one or more free parameters in a model can be calculated from the pooled
##' point estimates and asymptotic covariance matrix of model parameters
##' using Rubin's (1987) rules, or by pooling the Wald  test statistics
##' across imputed data sets (Li, Meng, Raghunathan, & Rubin, 1991).
##'
##' The constraints are specified using the \code{"=="} operator.
##' Both the left-hand side and the right-hand side of the equality can contain
##' a linear combination of model parameters, or a constant (like zero).
##' The model parameters must be specified by their user-specified labels from
##' the \code{link[lavaan]{model.syntax}}. Names of defined parameters
##' (using the ":=" operator) can be included too.
##'
##' @aliases lavTestWald.mi
##' @importFrom lavaan parTable lavListInspect
##' @importFrom stats pchisq pf
##' @importFrom methods getMethod
##'
##' @param object An object of class \code{\linkS4class{lavaan.mi}}.
##' @param constraints A \code{character} string (typically between single
##'   quotes) containing one or more equality constraints.
##'   See examples for more details
##' @param test \code{character} indicating which pooling method to use.
##'   \code{"D1"} or \code{"Rubin"} (default) indicates Rubin's (1987) rules
##'   will be applied to the point estimates and the asymptotic covariance
##'   matrix of model parameters, and those pooled values will be used to
##'   calculate the Wald test in the usual manner. \code{"D2"}, \code{"LMRR"},
##'   or \code{"Li.et.al"} indicate that the complete-data Wald test statistic
##'   should be calculated using each imputed data set, which will then be
##'   pooled across imputations, as described in Li, Meng, Raghunathan, & Rubin
##'   (1991) and Enders (2010, chapter 8).
##' @param asymptotic \code{logical}. If \code{FALSE} (default), the pooled test
##'   will be returned as an \emph{F}-distributed statistic with numerator
##'   (\code{df1}) and denominator (\code{df2}) degrees of freedom.
##'   If \code{TRUE}, the pooled \emph{F} statistic will be multiplied by its
##'   \code{df1} on the assumption that its \code{df2} is sufficiently large
##'   enough that the statistic will be asymptotically \eqn{\chi^2} distributed
##'   with \code{df1}.
##' @param scale.W \code{logical}. If \code{FALSE}, the pooled
##'   asymptotic covariance matrix of model parameters is calculated as the
##'   weighted sum of the within-imputation and between-imputation components.
##'   Otherwise, the pooled asymptotic covariance matrix of model parameters is
##'   calculated by scaling the within-imputation component by the
##'   average relative increase in variance (ARIV; see Enders, 2010, p. 235),
##'   which is \emph{only} consistent when requesting the \emph{F} test (i.e.,
##'   \code{asymptotic = FALSE}.  Ignored (irrelevant) if \code{test = "D2"}.
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'   imputations from pooled results.  Can include any of
##'   \code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
##'   default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (\code{"no.npd"}) would exclude any imputations which
##'   yielded a nonpositive definite covariance matrix for observed or
##'   latent variables, which would include any "improper solutions" such
##'   as Heywood cases. Specific imputation numbers can also be included in this
##'   argument, in case users want to  apply their own custom omission criteria
##'   (or simulations can use different numbers of imputations without
##'   redundantly refitting the model).
##' @param verbose \code{logical}. If \code{TRUE}, print the restriction
##'   matrix and the estimated restricted values.
##' @param warn \code{logical}. If \code{TRUE}, print warnings if they occur.
##'
##' @return
##'   A vector containing the Wald test statistic (either an \code{F} or
##'   \eqn{\chi^2} statistic, depending on the \code{asymptotic} argument),
##'   the degrees of freedom (numerator and denominator, if
##'   \code{asymptotic = FALSE}), and a \emph{p} value. If
##'   \code{asymptotic = FALSE}, the relative invrease in variance (RIV, or
##'   average for multiparameter tests: ARIV) used to calculate the denominator
##'   \emph{df} is also returned as a missing-data diagnostic, along with the
##'   fraction missing information (FMI = ARIV / (1 + ARIV)).
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##'   Adapted from \pkg{lavaan} source code, written by
##'   Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
##'
##' @references
##'   Enders, C. K. (2010). \emph{Applied missing data analysis}.
##'   New York, NY: Guilford.
##'
##'   Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
##'   Significance levels from repeated \emph{p}-values with multiply-imputed
##'   data. \emph{Statistica Sinica, 1}(1), 65--92. Retrieved from
##'   \url{https://www.jstor.org/stable/24303994}
##'
##'   Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}.
##'   New York, NY: Wiley.
##'
##' @seealso \code{\link[lavaan]{lavTestWald}}
##'
##' @examples
##'  \dontrun{
##' ## impose missing data for example
##' HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
##'                                       "ageyr","agemo","school")]
##' set.seed(12345)
##' HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
##' age <- HSMiss$ageyr + HSMiss$agemo/12
##' HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
##'
##' ## impute missing data
##' library(Amelia)
##' set.seed(12345)
##' HS.amelia <- amelia(HSMiss, m = 20, noms = "school", p2s = FALSE)
##' imps <- HS.amelia$imputations
##'
##' ## specify CFA model from lavaan's ?cfa help page
##' HS.model <- '
##'   visual  =~ x1 + b1*x2 + x3
##'   textual =~ x4 + b2*x5 + x6
##'   speed   =~ x7 + b3*x8 + x9
##' '
##'
##' fit <- cfa.mi(HS.model, data = imps)
##'
##' ## Testing whether a single parameter equals zero yields the 'chi-square'
##' ## version of the Wald z statistic from the summary() output, or the
##' ## 'F' version of the t statistic from the summary() output, depending
##' ## whether asymptotic = TRUE or FALSE
##' lavTestWald.mi(fit, constraints = "b1 == 0")      # default D1 statistic
##' lavTestWald.mi(fit, constraints = "b1 == 0", test = "D2") # D2 statistic
##'
##' ## The real advantage is simultaneously testing several equality
##' ## constraints, or testing more complex constraints:
##' con <- '
##'    2*b1 == b3
##'    b2 - b3 == 0
##' '
##' lavTestWald.mi(fit, constraints = con) # default F statistic
##' lavTestWald.mi(fit, constraints = con, asymptotic = TRUE) # chi-squared
##'
##' }
##'
##' @export
lavTestWald.mi <- function(object, constraints = NULL, test = c("D1","D2"),
                           asymptotic = FALSE, scale.W = !asymptotic,
                           omit.imps = c("no.conv","no.se"),
                           verbose = FALSE, warn = TRUE) {
  stopifnot(inherits(object, "lavaan.mi"))

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

  test <- tolower(test[1])
  if (test %in% c("d2", "lmrr", "li.et.al")) test <- "D2"
  if (test %in% c("d1", "rubin")) test <- "D1"
  if (!test %in% c("D1","D2")) stop('Invalid choice of "test" argument.')

  message('\nWald test calculated using se = "',
          lavListInspect(object, "options")$se, '"\n')

  if (test == "D2") {

    oldCall <- object@lavListCall

    if (!is.null(oldCall$parallel)) {
      if (oldCall$parallel == "snow") {
        oldCall$parallel <- "no"
        oldCall$ncpus <- 1L
        if (warn) warning("Unable to pass lavaan::lavTestWald() arguments ",
                          "when parallel='snow'. Switching to parallel='no'.",
                          " Unless using Windows, parallel='multicore' works.")
      }
    }

    ## call lavaanList() again to run lavTestWald() on each imputation
    oldCall$FUN <- function(obj) {
      out <- try(lavaan::lavTestWald(object = obj, constraints = constraints,
                                     verbose = FALSE), silent = TRUE)
      if (inherits(out, "try-error")) return(NULL)
      do.call(c, out[1:2])
    }
    FIT <- eval(as.call(oldCall))
    ## check if there are any results
    noStats <- sapply(FIT@funList, is.null)
    if (all(noStats)) stop("No success using lavTestWald() on any imputations.")

    ## template to fill in pooled values

    ## at a minimum, pool the total score test
    chiList <- sapply(FIT@funList[ intersect(useImps, which(!noStats)) ],
                      "[[", i = 1)
    DF <- FIT@funList[[ intersect(useImps, which(!noStats))[1] ]][[2]]
    out <- calculate.D2(chiList, DF = DF, asymptotic)
    class(out) <- c("lavaan.vector","numeric")
    return(out)
  } # else test == "D1", making 'scale.W=' relevant


  ## "borrowed" lavTestWald()
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
  theta <- getMethod("coef", "lavaan.mi")(object, omit.imps = omit.imps) #object@optim$x

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
  VCOV <- getMethod("vcov","lavaan.mi")(object, scale.W = scale.W,
                                        omit.imps = omit.imps)
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
    W <- getMethod("vcov", "lavaan.mi")(object, type = "within",
                                        omit.imps = omit.imps)
    B <- getMethod("vcov", "lavaan.mi")(object, type = "between",
                                        omit.imps = omit.imps)
    #FIXME: only valid for linear constraints?
    ## restricted B & W components of VCOV
    W.r  <- JAC %*% W %*% t(JAC)
    B.r  <- JAC %*% B %*% t(JAC)
    ## relative increase in variance due to missing data
    W.inv <- MASS::ginv(W.r)
    ariv <- (1 + 1/m) * sum(diag(B.r %*% W.inv)) / DF
    ## calculate denominator DF for F statistic
    a <- DF*(m - 1)
    if (a > 4) {
      v2 <- 4 + (a - 4) * (1 + (1 - 2/a)*(1 / ariv))^2 # Enders (eq. 8.24)
    } else {
      v2 <- a*(1 + 1/DF) * (1 + 1/ariv)^2 / 2 # Enders (eq. 8.25)
    }
    out <- c("F" = test.stat / DF, df1 = DF, df2 = v2,
             pvalue = pf(test.stat / DF, df1 = DF, df2 = v2, lower.tail = FALSE),
             ariv = ariv, fmi = ariv / (1 + ariv))
  }

  class(out) <- c("lavaan.vector","numeric")
  out
}



