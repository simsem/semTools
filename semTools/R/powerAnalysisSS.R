### Alexander M. Schoemann & Terrence D. Jorgensen
### Last updated: 2 June 2022
### Function to apply Satorra & Saris method for chi-squared power analysis


## Steps:
##   1. Specify model (use lavaan syntax based on simulateData)
##   2. get model implied covariance matrix
##   3. Fit model with parameter constrained to 0 (or take a model specification for multiparameter tests?)
##   4. Use chi square from step 3 as non-centrality parameter to get power.
## Alternatively, combine steps 1 and 2 by providing population moments directly


##' Power for model parameters
##'
##' Apply Satorra & Saris (1985) method for chi-squared power analysis.
##'
##' Specify all non-zero parameters in a population model, either by using
##' lavaan syntax (`popModel`) or by submitting a population covariance
##' matrix (`Sigma`) and optional mean vector (`mu`) implied by the
##' population model. Then specify an analysis model that places at least
##' one invalid constraint (note the number in the `nparam` argument).
##'
##' There is also a Shiny app called "power4SEM" that provides a graphical user
##' interface for this functionality (Jak et al., in press).  It can be accessed
##' at <https://sjak.shinyapps.io/power4SEM/>.
##'
##'
##' @importFrom stats qchisq pchisq
##'
##' @param powerModel lavaan [lavaan::model.syntax()] for the model to
##'   be analyzed. This syntax should constrain at least one nonzero parameter
##'   to 0 (or another number).
##' @param n `integer`. Sample size used in power calculation, or a vector
##'   of sample sizes if analyzing a multigroup model. If
##'   `length(n) < length(Sigma)` when `Sigma` is a list, `n` will
##'   be recycled. If `popModel` is used instead of `Sigma`, `n`
##'   must specify a sample size for each group, because that is used to infer
##'   the number of groups.
##' @param nparam `integer`. Number of invalid constraints in `powerModel`.
##' @param popModel lavaan [lavaan::model.syntax()] specifying the
##'   data-generating model. This syntax should specify values for all nonzero
##'   parameters in the model. If `length(n) > 1`, the same population
##'   values will be used for each group, unless different population values are
##'   specified per group, either in the lavaan [lavaan::model.syntax()]
##'   or by utilizing a list of `Sigma` (and optionally `mu`).
##' @param mu `numeric` or `list`. For a single-group model, a vector
##'   of population means. For a multigroup model, a list of vectors (one per
##'   group). If `mu` and `popModel` are missing, mean structure will
##'   be excluded from the analysis.
##' @param Sigma `matrix` or `list`. For a single-group model,
##'   a population covariance matrix. For a multigroup model, a list of matrices
##'   (one per group). If missing, `popModel` will be used to generate a
##'   model-implied Sigma.
##' @param fun character. Name of `lavaan` function used to fit
##'   `powerModel` (i.e., `"cfa"`, `"sem"`, `"growth"`, or
##'   `"lavaan"`).
##' @param alpha Type I error rate used to set a criterion for rejecting H0.
##' @param ... additional arguments to pass to [lavaan::lavaan()].
##'    See also [lavaan::lavOptions()].
##'
##' @author
##' Alexander M. Schoemann (East Carolina University; \email{schoemanna@@ecu.edu})
##'
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##'  Satorra, A., & Saris, W. E. (1985). Power of the likelihood ratio
##'  test in covariance structure analysis. *Psychometrika, 50*(1), 83--90.
##'  \doi{10.1007/BF02294150}
##'
##'  Jak, S., Jorgensen, T. D., Verdam, M. G., Oort, F. J., & Elffers, L.
##'  (2021). Analytical power calculations for structural equation modeling:
##'  A tutorial and Shiny app. *Behavior Research Methods, 53*, 1385--1406.
##'  \doi{10.3758/s13428-020-01479-0}
##'
##' @examples
##' ## Specify population values. Note every parameter has a fixed value.
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
##' Ns <- seq(100, 500, 40)
##' Power <- rep(NA, length(Ns))
##' for(i in 1:length(Ns)) {
##'   Power[i] <- SSpower(powerModel = modelA, popModel = modelP,
##'                       n = Ns[i], nparam = 1, std.lv = TRUE)
##' }
##' plot(x = Ns, y = Power, type = "l", xlab = "Sample Size")
##'
##'
##' ## Optionally specify different values for multiple populations
##'
##' modelP2 <- '
##'   f1 =~ .7*V1 + .7*V2 + .7*V3 + .7*V4
##'   f2 =~ .7*V5 + .7*V6 + .7*V7 + .7*V8
##'   f1 ~~ c(-.3, .3)*f2                  # DIFFERENT ACROSS GROUPS
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
##' modelA2 <- '
##'   f1 =~ V1 + V2 + V3 + V4
##'   f2 =~ V5 + V6 + V7 + V8
##'   f1 ~~ c(psi21, psi21)*f2        # EQUALITY CONSTRAINT ACROSS GROUPS
##' '
##' ## Calculate power
##' SSpower(powerModel = modelA2, popModel = modelP2, n = c(100, 100), nparam = 1,
##'         std.lv = TRUE)
##' ## Get power for a range of sample sizes
##' Ns2 <- cbind(Group1 = seq(10, 100, 10), Group2 = seq(10, 100, 10))
##' Power2 <- apply(Ns2, MARGIN = 1, FUN = function(nn) {
##'   SSpower(powerModel = modelA2, popModel = modelP2, n = nn,
##'           nparam = 1, std.lv = TRUE)
##' })
##' plot(x = rowSums(Ns2), y = Power2, type = "l", xlab = "Total Sample Size",
##'      ylim = 0:1)
##' abline(h = c(.8, .9), lty = c("dotted","dashed"))
##' legend("bottomright", c("80% Power","90% Power"), lty = c("dotted","dashed"))
##'
##' @export
SSpower <- function(powerModel, n, nparam, popModel, mu, Sigma,
                    fun = "sem", alpha = .05, ...) {
  if (missing(Sigma)) {

    ## specify (vector of) sample size(s) for optional multigroup syntax to work
    popMoments <- lavaan::fitted(do.call(fun, list(model = popModel,
                                                   sample.nobs = n),
                                         envir = getNamespace("lavaan")))
    ## without data, can't apply fitted() to multigroup model syntax, so
    ## save the same fitted moments for each group
    if (length(n) > 1L) {
      Sigma <- lapply(popMoments, "[[", i = "cov")
      mu <- if (!is.null(popMoments[[1]]$mean)) {
        lapply(popMoments, "[[", i = "mean")
      } else NULL

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
  fit <- do.call(fun, useArgs, envir = getNamespace("lavaan"))

  ## get NCP from chi square
  ncp <- lavaan::fitmeasures(fit)["chisq"]
  ## critical value under H0
  critVal <- qchisq(alpha, df = nparam, lower.tail = FALSE)
  ## return power
  pchisq(critVal, df = nparam, ncp = ncp, lower.tail = FALSE)
}


