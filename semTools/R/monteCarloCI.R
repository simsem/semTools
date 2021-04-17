### Terrence D. Jorgensen
### Last updated: 16 April 2021

## from http://www.da.ugent.be/cvs/pages/en/Presentations/Presentation%20Yves%20Rosseel.pdf
# dd <- read.table("http://www.statmodel.com/examples/shortform/4cat%20m.dat",
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
##' sampling distriution of estimated model parameters, as described by
##' MacKinnon et al. (2004) for testing indirect effects in mediation models.
##' The easiest way to use the function is to fit a SEM to data with
##' \code{\link[lavaan]{lavaan}}, using the \code{:=} operator in the
##' \code{\link[lavaan]{model.syntax}} to specify user-defined parameters.
##' All information is then available in the resulting
##' \code{\linkS4class{lavaan}} object.  Alternatively (especially when using
##' external SEM software to fit the model), the expression(s) can be explicitly
##' passed to the function, along with the vector of estimated model parameters
##' and their associated asymptotic sampling covariance matrix (ACOV).
##' For further information on the Monte Carlo method, see MacKinnon et al.
##' (2004) and Preacher & Selig (2012).
##'
##' The asymptotic covariance matrix can be obtained easily from many popular
##' SEM software packages.
##' \itemize{
##'  \item LISREL: Including the EC option on the OU line will print the ACM
##'    to a seperate file. The file contains the lower triangular elements of
##'    the ACM in free format and scientific notation
##'  \item Mplus Include the command TECH3; in the OUTPUT section. The ACM will
##'    be printed in the output.
##'  \item \code{lavaan}: Use the \code{vcov} method on the fitted
##'    \code{\linkS4class{lavaan}} object to return the ACM.
##' }
##'
##'
##' @importFrom stats quantile
##' @importFrom lavaan parTable
##'
##' @param object A object of class \code{\linkS4class{lavaan}} in which
##'   functions of parameters have already been defined using the \code{:=}
##'   operator in \code{lavaan}'s \code{\link[lavaan]{model.syntax}}. When
##'   \code{NULL}, users must specify \code{expr}, \code{coefs}, and \code{ACM}.
##' @param expr Optional \code{character} vector specifying functions of model
##'   parameters (e.g., an indirect effect). Ideally, the vector should have
##'   names, which is necessary if any user-defined parameters refer to other
##'   user-defined parameters defined earlier in the vector (order matters!).
##'   All parameters appearing in the vector must be provided in \code{coefs},
##'   or defined (as functions of \code{coefs}) earlier in \code{expr}. If
##'   \code{length(expr) > 1L}, \code{nRep} samples will be drawn
##'   simultaneously from a single multivariate distribution; thus,
##'   \code{ACM} must include all parameters in \code{coefs}.
##' @param coefs \code{numeric} vector of parameter estimates used in
##'   \code{expr}. Ignored when \code{object} is used.
##' @param ACM Symmetric \code{matrix} representing the asymptotic sampling
##'   covariance matrix (ACOV) of the parameter estimates in \code{coefs}.
##'   Ignored when \code{object} is used. Information on how to obtain the ACOV
##'   in popular SEM software is described in \strong{Details}.
##' @param nRep \code{integer}. The number of samples to draw, to obtain an
##'   empirical sampling distribution of model parameters. Many thousand are
##'   recommended to minimize Monte Carlo error of the estimated CIs.
##' @param fast \code{logical} indicating whether to use a fast algorithm that
##'   assumes all functions of parameters (in \code{object} or \code{expr}) use
##'   standard operations. Set to \code{FALSE} if using (e.g.) \code{\link{c}()}
##'   to concatenate parameters in the definition, which would have unintended
##'   consequences when vectorizing functions in \code{expr} across sampled
##'   parameters.
##' @param level \code{numeric} confidence level, between 0--1
##' @param na.rm \code{logical} passed to \code{\link[stats]{quantile}}
##' @param return.samples \code{logical} indicating whether to return the
##'   simulated empirical sampling distribution of parameters (in \code{coefs})
##'   and functions (in \code{expr})
##' @param plot \code{logical} indicating whether to plot the empirical sampling
##'   distribution of each function in \code{expr}
##' @param ask whether to prompt user before printing each plot
##' @param \dots arguments passed to \code{\link[graphics]{hist}} when
##'   \code{plot = TRUE}.
##'
##' @return A \code{lavaan.data.frame} (to use lavaan's \code{print} method).
##'   By default, a \code{data.frame} with point estimates and confidence limits
##'   of each requested function of parameters in \code{expr} is returned.
##'   If \code{return.samples = TRUE}, output will be a \code{data.frame} with
##'   the samples (in rows) of each parameter (in columns), and an additional
##'   column for each requested function of those parameters.
##'
##' @author Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##' MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence limits
##' for the indirect effect: Distribution of the product and resampling methods.
##' \emph{Multivariate Behavioral Research, 39}(1) 99--128.
##' \doi{10.1207/s15327906mbr3901_4}
##'
##' Preacher, K. J., & Selig, J. P. (2010, July). Monte Carlo method
##' for assessing multilevel mediation: An interactive tool for creating
##' confidence intervals for indirect effects in 1-1-1 multilevel models
##' [Computer software]. Available from \url{http://quantpsy.org/}.
##'
##' Preacher, K. J., & Selig, J. P. (2012). Advantages of Monte Carlo confidence
##' intervals for indirect effects. \emph{Communication Methods and Measures,
##' 6}(2), 77--98. \doi{10.1080/19312458.2012.679848}
##'
##' Selig, J. P., & Preacher, K. J. (2008, June). Monte Carlo method for
##' assessing mediation: An interactive tool for creating confidence intervals
##' for indirect effects [Computer software]. Available from
##' \url{http://quantpsy.org/}.
##'
##' @examples
##'
##' ## From the mediation tutorial:
##' ## http://lavaan.ugent.be/tutorial/mediation.html
##'
##' set.seed(1234)
##' X <- rnorm(100)
##' M <- 0.5*X + rnorm(100)
##' Y <- 0.7*M + rnorm(100)
##' dat <- data.frame(X = X, Y = Y, M = M)
##' mod <- ' # direct effect
##' Y ~ c*X
##' # mediator
##' M ~ a*X
##' Y ~ b*M
##' # indirect effect (a*b)
##' ind := a*b
##' # total effect
##' total := ind + c
##' '
##' fit <- sem(mod, data = dat)
##' summary(fit, ci = TRUE) # print delta-method CIs
##'
##' ## Automatically extract information from lavaan object
##' set.seed(1234)
##' monteCarloCI(fit) # CIs more robust than delta method in smaller samples
##'
##'
##' ## Parameter can also be obtained from an external analysis
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
monteCarloCI <- function(object = NULL, expr, coefs, ACM, nRep = 2e5, fast = TRUE,
                         level = 0.95, na.rm = TRUE, return.samples = FALSE,
                         plot = FALSE, ask = getOption("device.ask.default"),
                         ...) {

  if (class(object) == "lavaan") {
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
  if (class(object) == "lavaan") coefs <- lavaan::coef(object)
  sampVars <- intersect(names(coefs), funcVars)

  ## If a lavaan object is provided, extract coefs and ACM
  if (class(object) == "lavaan") {
    coefs <- coefs[sampVars]
    ACM <- lavaan::vcov(object)[sampVars, sampVars]
  }

  ## Apply the expression(s) to POINT ESTIMATES
  estList <- within(as.list(coefs), expr = {
    for (i in seq_along(expr)) assign(names(expr[i]), eval(parse(text = expr[i])))
  })[names(expr)]
  EST <- data.frame(est = do.call("c", estList))
  rownames(EST) <- names(expr)

  ## Matrix of sampled values
  dat <- data.frame(MASS::mvrnorm(n = nRep, mu = coefs, Sigma = ACM))
  ## Apply the expression(s) to VECTORS of ESTIMATES
  if (fast) {
    samples <- within(dat, expr = {
      for (i in seq_along(expr)) assign(names(expr[i]), eval(parse(text = expr[i])))
    })[c(sampVars, names(expr))]
  } else {
    ## SLOWER: only necessary if expr creates objects using (e.g.) c(), which
    ##         would concatenate parameters ACROSS samples as well as WITHIN
    datList <- lapply(1:nRep, function(Rep) dat[Rep,])
    samples <- do.call(rbind, lapply(datList, function(Rep) {
      within(Rep, expr = {
        for (i in seq_along(expr)) assign(names(expr[i]), eval(parse(text = expr[i])))
      })
    }))[c(sampVars, names(expr))]
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

  ## Return simulated values OR point and interval estimates
  out <- if (return.samples) samples else cbind(EST, CIs)
  class(out) <- c("lavaan.data.frame","data.frame")
  out
}
