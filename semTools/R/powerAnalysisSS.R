### Alexander M. Schoemann & Terrence D. Jorgensen
### Last updated: 15 September 2018
### Function to apply Satorra & Saris method for chi-squared power analysis


## Steps:
##   1. Specify model (use lavaan syntax based on simulateData)
##   2. get model implied covariance matrix
##   3. Fit model with parameter constrained to 0 (or take a model specification for multiparameter tests?)
##   4. Use chi square from step 3 as non-centrality parameter to get power.
## Alternatively, combine steps 1 and 2 by providng population moments directly


##' Power for model parameters
##'
##' Apply Satorra & Saris (1985) method for chi-squared power analysis.
##'
##' Specify all non-zero parameters in a population model, either by using
##' lavaan syntax (\code{popModel}) or by submitting a population covariance
##' matrix (\code{Sigma}) and optional mean vector (\code{mu}) implied by the
##' population model. Then specify an analysis model that constrains at least
##' one nonzero parameter to an incorrect value. Note the number in the
##' \code{nparam} argument.
##'
##'
##' @importFrom stats qchisq pchisq
##'
##' @param powerModel lavaan \code{\link[lavaan]{model.syntax}} for the model to
##'   be analyzed. This syntax should constrain at least one nonzero parameter
##'   to 0 (or another number).
##' @param n \code{integer}. Sample size used in power calculation, or a vector
##'   of sample sizes if analyzing a multigroup model. If
##'   \code{length(n) < length(Sigma)} when \code{Sigma} is a list, \code{n} will
##'   be recycled.
##' @param nparam \code{integer}. Number of invalid constraints in \code{powerModel}.
##' @param popModel lavaan \code{\link[lavaan]{model.syntax}} specifying the
##'   data-generating model. This syntax should specify values for all nonzero
##'   paramters in the model. If \code{length(n) > 1}, the same population
##'   values will be used for each group. Different population values per group
##'   can only be specified by utilizing \code{Sigma} (and \code{mu}).
##' @param mu numeric or list. For a single-group model, a vector of population
##'   means. For a multigroup model, a list of vectors (one per group). If
##'   \code{mu} and \code{popModel} are missing, mean structure will be excluded
##'   from the analysis.
##' @param Sigma matrix or list. For a single-group model, a population covariance
##'   matrix. For a multigroup model, a list of matrices (one per group). If
##'   missing, popModel will be used to generate a model-implied Sigma.
##' @param fun character. Name of lavaan function used to fit \code{powerModel}
##'   (i.e., \code{"cfa"}, \code{"sem"}, \code{"growth"}, or \code{"lavaan"}).
##' @param alpha Type I error rate used to set a criterion for rejecting H0.
##' @param ... additional arguments to pass to \code{\link[lavaan]{lavaan}}.
##'
##' @author
##' Alexander M. Schoemann (East Carolina University; \email{schoemanna@@ecu.edu})
##'
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##'  Satorra, A., & Saris, W. E. (1985). Power of the likelihood ratio
##'  test in covariance structure analysis. \emph{Psychometrika, 50}, 83--90.
##'  doi:10.1007/BF02294150
##'
##' @examples
##' ## Specify population values. Note every paramter has a fixed value.
##' modelP <- '
##'   f1 =~ .7*V1 + .7*V2 + .7*V3 + .7*V4
##'   f2 =~ .7*V5 + .7*V6 + .7*V7 + .7*V8
##'   f1 ~~ .3*f2
##'   f1 ~~ 1*f1
##'   f2 ~~ 1*f2
##'   V1 ~~ .51*V1
##'   V2 ~~ .51*V2
##'   V3 ~~ .51*V3
##'   V4 ~~ .51*V4
##'   V5 ~~ .51*V5
##'   V6 ~~ .51*V6
##'   V7 ~~ .51*V7
##'   V8 ~~ .51*V8
##' '
##' ## Specify analysis model. Note parameter of interest f1~~f2 is fixed to 0.
##' modelA <- '
##'   f1 =~ V1 + V2 + V3 + V4
##'   f2 =~ V5 + V6 + V7 + V8
##'   f1 ~~ 0*f2
##' '
##' ## Calculate power
##' SSpower(powerModel = modelA, popModel = modelP, n = 150, nparam = 1,
##'         std.lv = TRUE)
##'
##' ## Get power for a range of sample sizes
##'
##' Ns <- seq(100, 500, 40)
##' Power <- rep(NA, length(Ns))
##' for(i in 1:length(Ns)) {
##'   Power[i] <- SSpower(powerModel = modelA, popModel = modelP,
##'                       n = Ns[i], nparam = 1, std.lv = TRUE)
##' }
##' plot(x = Ns, y = Power, type = "l", xlab = "Sample Size")
##'
##' ## Specify second population to calculate power for multigroup model
##'
##' popMoments1 <- fitted(cfa(modelP))
##' modelP2 <- '
##'   f1 =~ .7*V1 + .7*V2 + .7*V3 + .7*V4
##'   f2 =~ .7*V5 + .7*V6 + .7*V7 + .7*V8
##'   f1 ~~ .5*f2     ## higher correlation in Group 2
##'   f1 ~~ 1*f1
##'   f2 ~~ 1*f2
##'   V1 ~~ .51*V1
##'   V2 ~~ .51*V2
##'   V3 ~~ .51*V3
##'   V4 ~~ .51*V4
##'   V5 ~~ .51*V5
##'   V6 ~~ .51*V6
##'   V7 ~~ .51*V7
##'   V8 ~~ .51*V8
##' '
##' popMoments2 <- fitted(cfa(modelP2))
##' modelA2 <- '
##'   f1 =~ V1 + V2 + V3 + V4
##'   f2 =~ V5 + V6 + V7 + V8
##'   f1 ~~ c(0, 0)*f2
##' '
##' mu <- list(popMoments1$mean, popMoments2$mean) # ignored if NULL
##' Sigma <- list(popMoments1$cov, popMoments2$cov)
##' SSpower(powerModel = modelA2, mu = mu, Sigma = Sigma,
##'         n = c(60, 65), nparam = 2)
##'
##' @export
SSpower <- function(powerModel, n, nparam, popModel, mu, Sigma,
                    fun = "cfa", alpha = .05, ...) {
  if (missing(Sigma)) {

    ## Two item list, first item is covariance matrix, second item is mean vector
    popMoments <- lavaan::fitted(do.call(fun, list(model = popModel)))
    ## without data, can't apply fitted() to multigroup model syntax, so
    ## save the same fitted moments for each group
    if (length(n) > 1L) {
      Sigma <- list()
      mu <- if (!is.null(popMoments[[1]]$mean)) list() else NULL
      for (g in 1:length(n)) {
        Sigma[[g]] <- popMoments$cov
        if (!is.null(popMoments[[g]]$mean)) mu[[g]] <- popMoments$mean
      }
    } else {
      ## single group
      Sigma <- popMoments$cov
      mu <- popMoments$mean
    }

  } else {
    ## otherwise, user-supplied moments

    if (is.list(Sigma)) {
      nG <- length(Sigma)
      if (length(n) < nG) n <- rep(n, length.out = nG)
      if (length(n) > nG) n <- n[1:nG]
      no.mu <- missing(mu)
      if (!no.mu) null.mu <- any(sapply(mu, is.null))
      if (no.mu || null.mu) {
        mu <- NULL
      }

    } else if (is.matrix(Sigma)) {
      n <- n[[1]]

      if (missing(mu)) {
        mu <- NULL
      } else if (!is.numeric(mu) || !!is.null(mu)) {
        stop('mu must be a numeric vector, or a list (one vector per group)')
      }

    } else stop('Sigma must be a covariance matrix, or a list (one matrix per group)')

  }

  ## Fit (probably constrained) analysis model
  dots <- list(...)
  funArgs <- list(model = powerModel, sample.cov = Sigma,
                  sample.mean = mu, sample.nobs = n)
  useArgs <- c(funArgs, dots[setdiff(names(dots), names(funArgs))])
  fit <- do.call(fun, useArgs)

  ## get NCP from chi square
  ncp <- lavaan::fitmeasures(fit)["chisq"]
  ## critical value under H0
  critVal <- qchisq(alpha, df = nparam, lower.tail = FALSE)
  ## return power
  pchisq(critVal, df = nparam, ncp = ncp, lower.tail = FALSE)
}


