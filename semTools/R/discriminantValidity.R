### Mikko Ronkko
### Last updated: 3 December 2019


##' Calculate discriminant validity statistics
##'
##' Calculate discriminant validity statistics based on a fitted lavaan object
##'
##' Evaluated on the measurement scale level, discriminant validity is commonly
##' evaluated by checking if each pair of latent correlations is sufficiently
##' below one (in absolute value) that the latent variables can be thought of
##' representing two distinct constructs.
##'
##' \code{discriminantValidity} function calculates two sets of statistics that
##' are commonly used in discriminant validity evaluation. The first set are
##' factor correlation estimates and their confidence intervals. The second set
##' is a series of nested model tests, where the baseline model is compared
##' against as set of constrained models that are constructed by constraining
##' each factor correlation to the specified cutoff one at a time.
##'
##' The function assume that the \code{object} is set of confirmatory
##' factor analysis results where the latent variables are scaled by fixing their
##' variances to 1s. If the model is not a CFA model, the function will calculate
##' the statistics for the correlations among exogenous latent variables, but
##' for the \emph{residual} variances with endogenous variables. If the
##' latent variables are scaled in some other way (e.g. fixing the first loadings),
##' the function issues a warning and re-estimates the model by fixing latent
##' variances to 1 (and estimating all loadings) so that factor covariances are
##' already estimated as correlations.
##'
##' The likelihood ratio tests are done by comparing the original baseline model
##' against more constrained alternatives. By default, these alternatives are
##' constructed by fixing each correlation at a time to a cutoff value. The
##' typical purpose of this test is to demonstrate that the estimated factor
##' correlation is well below the cutoff and a significant \eqn{chi^2} statistic
##' thus indicates support for discriminant validity. In some cases, the original
##' correlation estimate may already be greater than the cutoff, making it
##' redundant to fit a "restricted" model. When this happens, the likelihood
##' ratio test will be replaced by comparing the baseline model against itself.
##' For correlations that are estimated to be negative, a negation of the cutoff
##' is used in the constrained model.
##'
##' Another alternative is to do a nested model comparison against a model where
##' two factors are merged as one by setting the \code{merge} argument to
##' \code{TRUE}. In this comparison, the constrained model is constructed by
##' removing one of the correlated factors from the model and assigning its
##' indicators to the factor that remains in the model.
##'
##'
##' @importFrom lavaan lavInspect lavNames parTable
##'
##' @param object The \code{\linkS4class{lavaan}} model object returned by
##'   the \code{\link[lavaan]{cfa}} function.
##' @param cutoff A cutoff to be used in the constrained models in likelihood
##'   ratio tests.
##' @param merge Whether the constrained models should be constructed by merging
##'   two factors as one. Implies \code{cutoff} = 1.
##' @param level The confidence level required.
##'
##' @return A \code{data.frame} of latent variable correlation estimates, their
##' confidence intervals, and a likelihood ratio tests against constrained models.
##' with the following attributes:
##' \describe{
##'  \item{baseline}{The baseline model after possible rescaling.}
##'  \item{constrained}{A \code{list} of the fitted constrained models
##'  used in the likelihood ratio test.}
##' }
##'

## @author ADD MIKKO AFTER ARTICLE IS ACCEPTED, also add to DESCRIPTION:
## person(given = "Mikko", family = "Ronkko", role = "ctb", email = "mikko.ronkko@jyu.fi", comment = c(ORCID = "0000-0001-7988-7609")),

## @references ADD CITATION AFTER ARTICLE HAS A DOI

##'
##' @examples
##'
##' library(lavaan)
##'
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ x7 + x8 + x9 '
##'
##' fit <- cfa(HS.model, data = HolzingerSwineford1939)
##' discriminantValidity(fit)
##' discriminantValidity(fit, merge = TRUE)
##'
##' @export
discriminantValidity <- function(object, cutoff = .9, merge = FALSE, level = .95) {

  free <- lavInspect(object, "free", add.class = FALSE)
  #FIXME: adapt for multiple blocks by looping over groups/levels
  if (lavInspect(object, "ngroups") > 1L | lavInspect(object, "nlevels") > 1L)
    stop("Only implemented for single-group, single-level models so far.")

  # Identify the latent variables that we will use
  lvs <- lavNames(object,"lv")
  if (cutoff <=0 | cutoff >1) stop("The cutoff must be between (0,1]")
  if (merge & ! missing(cutoff) & cutoff != 1)
    message("Merging factors imply constraining factor correlation to 1. ",
            "Cutoff will be ignored.")
  if (length(lvs)==0) stop("The model does not have any exogenous latent variables.")
  if (length(lvs)==1) stop("The model has only one exogenous latent variable. ",
                           "At least two are required for assessing discriminant validity.")
  if (length(lavNames(object, "lv.y")) > 0)
    warning("The model has at least one endogenous latent variable (",
            paste(lavNames(object, "lv.y"), collapse=", "),
            "). The correlations of these variables will be estimated after ",
            "conditioning on their predictors.")

  # Extract the part of psi that contains latent variables
  psi <- free$psi[lvs,lvs]

  # Identify exogenous variances and covariances
  pt <- parTable(object)
  varIndices <- which(pt$lhs == pt$rhs & pt$lhs %in% lvs & pt$op =="~~")
  covIndices <- which(pt$lhs != pt$rhs & pt$lhs %in% lvs & pt$rhs %in% lvs & pt$op =="~~")

  # Check that the diagonal of psi is all zeros
  if (any(diag(psi) != 0)) {

    message("Some of the latent variable variances are estimated instead of ",
            "fixed to 1. The model is re-estimated by scaling the latent ",
            "variables by fixing their variances and freeing all factor loadings.")

    # Identify free exogenous variances
    i <- intersect(varIndices,which(pt$free != 0))

    pt$free[i] <- 0
    pt$ustart[i] <- 1
    pt$user[i] <- 1

    # Free all factor loadings corresponding of lvs where the covariances were just freed
    i <- which(pt$lhs %in% pt$lhs[i] & pt$op =="=~")
    pt$free[i] <- -1
    pt$ustart[i] <- NA

    # Update parameter numbers
    i <- which(pt$free != 0)
    pt$free[i] <- seq_along(i)

    object <- lavaan::update(object, model = pt[,1:12]) # Leave out starting values, estimates and ses from pt

    # Update pt based on the new model
    pt <- parTable(object)
  }

  # At this point we can be sure that all exogenous variances are fixed instead
  # of being estimated. We need to still check that they are fixed to 1s

  est <- lavInspect(object,"est")$psi[lvs,lvs]

  if (any(diag(est) != 1)) {
    message("Some of the latent variable variances are fixed to values other ",
            "than 1. The model is re-estimated by scaling the latent variables",
            " based on the first factor loading.")

    # constrain the exogenous variances to 1
    pt$ustart[varIndices] <- 1
    object <- lavaan::update(object, model = pt[,1:12]) # Leave out starting values, estimates and ses from pt

    # Update pt based on the new estimates
    pt <- parTable(object)
  }

  # At this point we can be sure that all exogenous LVs have their variances
  # fixed to ones and can start constructing the matrix to be returned

  ret <-  lavaan::parameterEstimates(object, ci = TRUE,
                                     level = level)[covIndices,
                                                    c("lhs","op","rhs","est",
                                                      "ci.lower","ci.upper")]
  rownames(ret) <- seq_len(nrow(ret))


  # Add the chi^2 test to all correlation pairs

  constrainedModels <- lapply(covIndices, function(i) {

    thisPt <- pt

    if (merge) {

      lhs <- pt$lhs[i]
      rhs <- pt$rhs[i]

      # Merge the factors by assigning indicator of lhs to rhs
      thisPt$lhs[thisPt$lhs == lhs & thisPt$op == "=~"] <- rhs

      # Then remove all other parameters concering lhs
      thisPt <- thisPt[!(thisPt$lhs == lhs | thisPt$rhs == lhs), ]

      thisPt$id <- seq_len(nrow(thisPt))
    } else {

      # If the correlation is estimated to be greater than the cuttof, constrain it to the estimated alue
      if (abs(pt$est[i]) > cutoff) {
        thisCutoff <- pt$est[i]
      } else {
        thisCutoff <- ifelse(pt$est[i] <0, - cutoff, cutoff)
      }
      thisPt$free[i] <- 0
      thisPt$ustart[i] <- thisCutoff
    }

    # Update parameter numbers
    j <- which(thisPt$free != 0)
    thisPt$free[j] <- seq_along(j)

    lavaan::update(object, model = thisPt[,1:12])
  })

  lrTests <- lapply(constrainedModels, function(constrained) {
    lavaan::lavTestLRT(object,constrained)[2,] # Return the second row of the test
  })

  ret <- cbind(ret,do.call(rbind,lrTests))

  # Store the baseline model
  attr(ret,"baseline") <- object
  attr(ret,"constrained") <- constrainedModels
  ret
}










