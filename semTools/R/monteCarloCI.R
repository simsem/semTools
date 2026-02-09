### Terrence D. Jorgensen
### Last updated: 9 February 2026

## from https://www.da.ugent.be/cvs/pages/en/Presentations/Presentation%20Yves%20Rosseel.pdf
# dd <- read.table("https://www.statmodel.com/examples/shortform/4cat%20m.dat",
#                  col.names = c("intention", "intervention", "ciguse", "w"))
# myData <- do.call(rbind, lapply(1:nrow(dd), function(RR) {
#   data.frame(rep(1, dd$w[RR]) %*% as.matrix(dd[RR, 1:3]))
# }))
# model <- '
#     ciguse ~ c*intervention + b*intention
#     intention ~ a*intervention
# # label threshold for ciguse
#     ciguse | b0*t1
# # biased SEs
#     naive.indirect := a*b
#     naive.direct := c
# # correct
#     probit11 := (-b0+c+b*a)/sqrt(b^2+1)
#     probit10 := (-b0+c )/sqrt(b^2+1)
#     probit00 := (-b0 )/sqrt(b^2+1)
#     indirect := pnorm(probit11) - pnorm(probit10)
#     direct := pnorm(probit10) - pnorm(probit00)
#     OR.indirect := (pnorm(probit11)/(1-pnorm(probit11))) / (pnorm(probit10)/(1-pnorm(probit10)))
#     OR.direct := (pnorm(probit10)/(1-pnorm(probit10))) / (pnorm(probit00)/(1-pnorm(probit00)))
# '
# fit <- sem(model, data = myData, ordered = c("ciguse","intention"))
# summary(fit, ci = TRUE)



##' Monte Carlo Confidence Intervals to Test Functions of Parameter Estimates
##'
##' Robust confidence intervals for functions of parameter estimates,
##' based on empirical sampling distributions of estimated model parameters.
##'
##' This function implements the Monte Carlo method of obtaining an empirical
##' sampling distribution of estimated model parameters, as described by
##' MacKinnon et al. (2004) for testing indirect effects in mediation models.
##' This is essentially a parametric bootstrap method, which (re)samples
##' parameters (rather than raw data) from a multivariate-normal distribution
##' with mean vector equal to estimates in `coef()` and covariance matrix
##' equal to the asymptotic covariance matrix `vcov()` of estimated parameters.
##'
##' The easiest way to use the function is to fit a SEM to data with
##' [lavaan::lavaan()], using the `:=` operator in the
##' [lavaan::model.syntax()] to specify user-defined parameters.
##' All information is then available in the resulting
##' [lavaan::lavaan-class] object.  Alternatively (especially when using
##' external SEM software to fit the model), the expression(s) can be explicitly
##' passed to the function, along with the vector of estimated model parameters
##' and their associated asymptotic sampling covariance matrix (ACOV).
##' For further information on the Monte Carlo method, see MacKinnon et al.
##' (2004) and Preacher & Selig (2012).
##'
##' The asymptotic covariance matrix can be obtained easily from many popular
##' SEM software packages.
##' \itemize{
##'  \item{LISREL: Including the EC option on the OU line will print the ACM
##'    to a seperate file. The file contains the lower triangular elements of
##'    the ACM in free format and scientific notation.}
##'  \item{M*plus*: Include the command TECH3; in the OUTPUT section.
##'    The ACM will be printed in the output.}
##'  \item{`lavaan`: Use the [vcov()] method on the fitted [lavaan::lavaan-class]
##'    object to return the ACM.}
##' }
##'
##'
##' @importFrom stats quantile
##' @importFrom methods getMethod
##' @importFrom lavaan parTable lavInspect
##'
##' @param object A object of class [lavaan::lavaan-class] in which
##'   functions of parameters have already been defined using the `:=`
##'   operator in `lavaan`'s [lavaan::model.syntax()]. When
##'   `NULL`, users must specify `expr`, `coefs`, and `ACM`.
##' @param expr Optional `character` vector specifying functions of model
##'   parameters (e.g., an indirect effect). Ideally, the vector should have
##'   names, which is necessary if any user-defined parameters refer to other
##'   user-defined parameters defined earlier in the vector (order matters!).
##'   All parameters appearing in the vector must be provided in `coefs`,
##'   or defined (as functions of `coefs`) earlier in `expr`. If
##'   `length(expr) > 1L`, `nRep` samples will be drawn
##'   simultaneously from a single multivariate distribution; thus,
##'   `ACM` must include all parameters in `coefs`.
##' @param coefs `numeric` vector of parameter estimates used in
##'   `expr`. Ignored when `object` is used.
##' @param ACM Symmetric `matrix` representing the asymptotic sampling
##'   covariance matrix (ACOV) of the parameter estimates in `coefs`.
##'   Ignored when `object` is used. Information on how to obtain the ACOV
##'   in popular SEM software is described in **Details**.
##' @param nRep `integer`. The number of samples to draw, to obtain an
##'   empirical sampling distribution of model parameters. Many thousand are
##'   recommended to minimize Monte Carlo error of the estimated CIs.
##' @param standardized `logical` indicating whether to obtain CIs for the
##'   fully standardized (`"std.all"`) estimates, using their asymptotic
##'   sampling covariance matrix.
##' @param fast `logical` indicating whether to use a fast algorithm that
##'   assumes all functions of parameters (in `object` or `expr`) use
##'   standard operations. Set to `FALSE` if using (e.g.) [c()]
##'   to concatenate parameters in the definition, which would have unintended
##'   consequences when vectorizing functions in `expr` across sampled
##'   parameters.
##' @param level `numeric` confidence level, between 0--1
##' @param na.rm `logical` passed to [stats::quantile()]
##' @param append.samples `logical` indicating whether to return the
##'   simulated empirical sampling distribution of parameters (in `coefs`)
##'   and functions (in `expr`) in a `list` with the results. This
##'   could be useful to calculate more precise highest-density intervals (see
##'   examples).
##' @param plot `logical` indicating whether to plot the empirical sampling
##'   distribution of each function in `expr`
##' @param ask whether to prompt user before printing each plot
##' @param \dots arguments passed to [graphics::hist()] when
##'   `plot = TRUE`.
##'
##' @return A `lavaan.data.frame` (to use lavaan's `print` method)
##'   with point estimates and confidence limits of each requested function of
##'   parameters in `expr` is returned. If `append.samples = TRUE`,
##'   output will be a `list` with the same `$Results` along with a
##'   second `data.frame` with the `$Samples` (in rows) of each
##'   parameter (in columns), and an additional column for each requested
##'   function of those parameters.
##'
##' @author Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##' MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence limits
##' for the indirect effect: Distribution of the product and resampling methods.
##' *Multivariate Behavioral Research, 39*(1) 99--128.
##' \doi{10.1207/s15327906mbr3901_4}
##'
##' Preacher, K. J., & Selig, J. P. (2010, July). Monte Carlo method
##' for assessing multilevel mediation: An interactive tool for creating
##' confidence intervals for indirect effects in 1-1-1 multilevel models.
##' Computer software available from <https://quantpsy.org/>.
##'
##' Preacher, K. J., & Selig, J. P. (2012). Advantages of Monte Carlo confidence
##' intervals for indirect effects. *Communication Methods and Measures, 6*(2),
##' 77--98. \doi{10.1080/19312458.2012.679848}
##'
##' Selig, J. P., & Preacher, K. J. (2008, June). Monte Carlo method for
##' assessing mediation: An interactive tool for creating confidence intervals
##' for indirect effects. Computer software available from
##' <https://quantpsy.org/>.
##'
##' @aliases monteCarloCI monteCarloMed
##'
##' @examples
##'
##' ## From the mediation tutorial:
##' ## https://lavaan.ugent.be/tutorial/mediation.html
##'
##' set.seed(1234)
##' X <- rnorm(100)
##' M <- 0.5*X + rnorm(100)
##' Y <- 0.7*M + rnorm(100)
##' dat <- data.frame(X = X, Y = Y, M = M)
##'
##' mod <- ' # direct effect
##'   Y ~ c*X
##'   # mediator
##'   M ~ a*X
##'   Y ~ b*M
##'   # indirect effect (a*b)
##'   ind := a*b
##'   # total effect
##'   total := ind + c
##' '
##' fit <- sem(mod, data = dat)
##' summary(fit, ci = TRUE) # print delta-method CIs
##'
##' ## Automatically extract information from lavaan object
##' set.seed(1234)
##' monteCarloCI(fit) # CIs more robust than delta method in smaller samples
##'
##' ## delta method for standardized solution
##' standardizedSolution(fit)
##' ## compare to Monte Carlo CIs:
##' set.seed(1234)
##' monteCarloCI(fit, standardized = TRUE)
##'
##' \donttest{
##' ## save samples to calculate more precise intervals:
##' set.seed(1234)
##' foo <- monteCarloCI(fit, append.samples = TRUE)
##' # library(HDInterval) # not a dependency; must be installed
##' # hdi(foo$Samples)
##' }
##' ## Parameters can also be obtained from an external analysis
##' myParams <- c("a","b","c")
##' (coefs <- coef(fit)[myParams]) # names must match those in the "expression"
##' ## Asymptotic covariance matrix from an external analysis
##' (AsyCovMat <- vcov(fit)[myParams, myParams])
##' ## Compute CI, include a plot
##' set.seed(1234)
##' monteCarloCI(expr = c(ind = 'a*b', total = 'ind + c',
##'                       ## other arbitrary functions are also possible
##'                       meaningless = 'sqrt(a)^b / log(abs(c))'),
##'              coefs = coefs, ACM = AsyCovMat,
##'              plot = TRUE, ask = TRUE) # print a plot for each
##'
##' @export
monteCarloCI <- function(object = NULL, expr, coefs, ACM, nRep = 2e4,
                         standardized = FALSE, fast = TRUE, level = 0.95,
                         na.rm = TRUE, append.samples = FALSE, plot = FALSE,
                         ask = getOption("device.ask.default"), ...) {

  if (inherits(object, c("lavaan","lavaan.mi"))) {
    ## extract user-defined parameters from parTable (order of user-defined
    PT <- parTable(object) # parameters must be correct for model to be fitted)
    ## create expression vector
    expr <- PT$rhs[PT$op == ":="]
    names(expr) <- PT$lhs[PT$op == ":="]
    if (length(expr) == 0L) stop('No user-defined parameters found.')
  }
  ## provide names if there are none
  if (is.null(names(expr))) names(expr) <- expr

  ## Get names and the number of unique variables in the expression
  funcVars <- unique(do.call("c", lapply(paste("~", expr), function(x) {
    all.vars(stats::as.formula(x))
  })))

  ## isolate names of model parameters (not user-defined), which get sampled
  if (inherits(object, "lavaan")) {
    if (standardized) {
      STD <- lavaan::standardizedSolution(object)
      coefRows <- !(STD$op %in% c(":=","==","<",">","<=",">="))
      coefs <- STD$est.std[coefRows]
      names(coefs) <- lavaan::lav_partable_labels(STD[coefRows, ])
    } else coefs <- lavaan::coef(object)

  } else if (inherits(object, "lavaan.mi")) {
    if (!"package:lavaan.mi" %in% search()) attachNamespace("lavaan.mi")

    if (standardized) {
      ## only available after this lavaan version:
      if ( utils::packageDescription("lavaan", fields = "Version") < "0.6-19" ||
          (utils::packageDescription("lavaan", fields = "Version") > "0.6-19" &&
           utils::packageDescription("lavaan", fields = "Version") < "0.6-19.2148") ) {
        stop("standardized=TRUE for a lavaan.mi-class object requires lavaan ",
             "version >= 0.6-19 from CRAN, or development version >= ",
             "0.6-19.2148 from GitHub")
      }
      STD <- lavaan.mi::standardizedSolution.mi(object)
      coefRows <- !(STD$op %in% c(":=","==","<",">","<=",">="))
      coefs <- STD$est.std[coefRows]
      names(coefs) <- lavaan::lav_partable_labels(STD[coefRows, ])
    } else coefs <- getMethod(f = "coef", signature = "lavaan.mi",
                              where = getNamespace("lavaan.mi"))(object)
  }
  sampVars <- intersect(names(coefs), funcVars)

  ## If a lavaan(.mi) object is provided, extract coefs and ACM
  if (inherits(object, "lavaan")) {
    coefs <- coefs[sampVars]
    if (standardized) {
      ACM <- lavInspect(object, "vcov.std.all")[sampVars, sampVars]
    } else {
      ACM <- lavaan::vcov(object)[sampVars, sampVars]
    }

  } else if (inherits(object, "lavaan.mi")) {
    coefs <- coefs[sampVars]
    if (standardized) {
      ACM <- lavaan.mi::standardizedSolution.mi(object, return.vcov = TRUE,
                                                type = "std.all")[sampVars, sampVars]
    } else {
      ACM <- getMethod(f = "vcov", signature = "lavaan.mi",
                       where = getNamespace("lavaan.mi"))(object)[sampVars, sampVars]
    }
  }

  ## Apply the expression(s) to POINT ESTIMATES
  estList <- as.list(coefs)
  for (i in seq_along(expr)) {
    estList[names(expr[i])] <- eval(parse(text = expr[i]), envir = estList)
  }
  EST <- data.frame(est = do.call("c", estList[names(expr)]))
  ## old, buggy code (see issue #142)
  # estList <- within(as.list(coefs), expr = {
  #   for (i in seq_along(expr)) assign(names(expr[i]), eval(parse(text = expr[i])))
  # })[names(expr)]
  # EST <- data.frame(est = do.call("c", estList))
  rownames(EST) <- names(expr)
  if (standardized && inherits(object, "lavaan")) colnames(EST) <- "est.std"

  ## Matrix of sampled values
  # dat <-
  samples <- data.frame(mnormt::rmnorm(n = nRep, mean = coefs, varcov = ACM))
  colnames(samples) <- names(coefs)
  ## Apply the expression(s) to VECTORS of ESTIMATES
  if (fast) {
    for (i in seq_along(expr)) {
      samples[names(expr[i])] <- eval(parse(text = expr[i]), envir = samples)
    }
    ## old, buggy code (see issue #142)
    # samples <- within(dat, expr = {
    #   for (i in seq_along(expr)) assign(names(expr[i]), eval(parse(text = expr[i])))
    # })[c(sampVars, names(expr))]
  } else {
    ## SLOWER: only necessary if expr creates objects using (e.g.) c(), which
    ##         would concatenate parameters ACROSS samples as well as WITHIN
    datList <- lapply(1:nRep, function(Rep) {
      samples[Rep,] # dat[Rep,]
    })
    samples <- do.call(rbind, lapply(datList, function(Rep) {
      for (i in seq_along(expr)) {
        Rep[names(expr[i])] <- eval(parse(text = expr[i]), envir = Rep)
      }
      ## old, buggy code (see issue #142)
      #   within(Rep, expr = {
      #     for (i in seq_along(expr)) assign(names(expr[i]), eval(parse(text = expr[i])))
      #   })
      Rep
    })) # [c(sampVars, names(expr))]
  }

  ## Get the CI(s)
  halfAlpha <- (1-level)/2
  Limits <- lapply(samples[names(expr)], quantile, na.rm = na.rm,
                   probs = c(halfAlpha, 1 - halfAlpha))
  CIs <- data.frame(do.call("rbind", Limits))
  rownames(CIs) <- names(expr)
  colnames(CIs) <- c("ci.lower","ci.upper")

  ## Switch for outputting a plot
  if (plot) {
    if (length(expr) > 1L && ask) {
      opar <- grDevices::devAskNewPage()
      grDevices::devAskNewPage(ask = TRUE)
    }
    for (i in seq_along(expr)) {
      histArgs <- list(...)
      histArgs$x <- samples[[ names(expr)[i] ]]
      if (is.null(histArgs$breaks)) histArgs$breaks <- "FD"
      if (is.null(histArgs$xlab)) histArgs$xlab <- paste0(level*100, '% Confidence Interval')
      if (is.null(histArgs$main)) histArgs$main <- paste('Distribution of', names(expr)[i])
      do.call("hist", histArgs)
      abline(v = EST[i,1], lwd = 3)
      abline(v = CIs[i,1:2], lwd = 2, lty = "dashed")
    }
    if (length(expr) > 1L && ask) grDevices::devAskNewPage(ask = opar)
  }

  ## Always return point and interval estimates
  out <- cbind(EST, CIs)
  class(out) <- c("lavaan.data.frame","data.frame")
  ## also simulated values? (e.g., to calculate highest-density intervals)
  if (append.samples) return(list(Results = out, Samples = samples))
  out
}

#FIXME: Remove after a few version updates
monteCarloMed <- function(...) {
  .Defunct("monteCarloCI",
           msg = "monteCarloMed() has been replaced by monteCarloCI()")
}
