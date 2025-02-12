### Leonard Vanbrabant (Roxygen edits by Terrence D. Jorgensen)
### Last updated: 12 February 2025

##' Wrapper for `goric.lavaan()` from the `restriktor` package
##'
##' The `goricaSEM()` function is an interface to [restriktor::goric.lavaan()],
##' allowing users to perform generalized order-restricted information criterion
##' approximation (GORICA) analysis specifically for structural equation
##' models fitted using the \pkg{lavaan} package.
##'
##' @details
##' This function is designed as a wrapper for the [restriktor::goric.lavaan()]
##' function. It calculates GORICA values and weights, which can be used to
##' compare models or hypotheses under inequality constraints.
##'
##' The `hypotheses=` argument allows users to specify constraints in text-based
##' syntax or matrix notation. For text-based syntax, constraints are specified
##' as a string (e.g., `"a1 > a2"`). For matrix notation, a named list with
##' `$constraints`, `$rhs`, and `$neq` elements can be provided.
##'
##' The `comparison=` argument determines whether the specified hypothesis is
##' compared against its `"complement"`, the `"unconstrained"` model, or
##' neither (`"none"`).
##'
##' @param object A [lavaan::lavaan-class] object.
##' @param hypotheses A named `list` of hypotheses to test. See **Details** for
##'   information on how to specify hypotheses.
##' @param comparison A `character` string specifying the type of comparison.
##'   Options are `"unconstrained"`, `"complement"`, or `"none"`.
##'   Default behavior depends on the number of hypotheses.
##' @param type A `character` string indicating the type of analysis, either
##'   `"gorica"` (default) or `"goricac"`.
##' @param standardized `logical` indicating whether standardized estimates are
##'   used in the analysis. Defaults to `FALSE`.
##' @param debug `logical` indicating whether to print debugging information.
##'   Defaults to \code{FALSE}.
##' @param ... Additional arguments passed to [restriktor::goric.lavaan()].
##'
##' @return
##'   A `list` containing the results of the \code{goric.lavaan} function,
##'   including:
##'   \itemize{
##'     \item The log-likelihood.
##'     \item Penalty term.
##'     \item GORIC(A) values and weights.
##'     \item Relative GORIC(A) weights.
##'   }
##'
##' @references
##'   Kuiper, R. M., Hoijtink, H., & Silvapulle, M. J. (2011). An Akaike-type
##'   information criterion for model selection under inequality constraints.
##'   \emph{Biometrika, 98}(2), 495--501. \doi{10.1093/biomet/asr002}
##'
##'   Vanbrabant, L., Van Loey, N., & Kuiper, R. M. (2020). Evaluating a
##'   theory-based hypothesis against its complement using an AIC-type
##'   information criterion with an application to facial burn injury.
##'   \emph{Psychological Methods, 25}(2), 129--142. \doi{10.1037/met0000238}
##'
##' @seealso [restriktor::goric.lavaan()]
##'
##' @author Leonard Vanbrabant and Rebecca Kuiper
##'
##' @examples
##'
##' ## Example: Perform GORICA analysis on a lavaan model
##' library(lavaan)
##' library(restriktor)
##'
##' ## Define the SEM model
##' model <- '
##'   ind60 =~ x1 + x2 + x3
##'   dem60 =~ y1 + a1*y2 + b1*y3 + c1*y4
##'   dem65 =~ y5 + a2*y6 + b2*y7 + c2*y8
##'   dem60 ~ ind60
##'   dem65 ~ ind60 + dem60
##'   y1 ~~ y5
##'   y2 ~~ y4 + y6
##'   y3 ~~ y7
##'   y4 ~~ y8
##'   y6 ~~ y8
##' '
##'
##' ## Fit the model
##' data(PoliticalDemocracy)
##' fit <- sem(model, data = PoliticalDemocracy)
##'
##' ## Define hypotheses
##' myHypothesis <- 'a1 > a2, b1 > b2, c1 > c2'
##'
##' ## Perform GORICA analysis
##' result <- goricaSEM(fit, hypotheses = list(H1 = myHypothesis),
##'                     standardized = FALSE, comparison = "complement",
##'                     type = "gorica")
##'
##' ## Print result
##' print(result)
##'
##' @export
goricaSEM <- function(object, ..., hypotheses = NULL,
                      comparison = NULL,
                      type = "gorica",
                      standardized = FALSE,
                      debug = FALSE) {

  ## Check if the original function is available
  if (!requireNamespace("restriktor", quietly = TRUE)) {
    stop("The 'restriktor' package is required but not installed.")
  }

  ## Call the original function
  restriktor::goric.lavaan(
    object = object,
    hypotheses = hypotheses,
    comparison = comparison,
    type = type,
    standardized = standardized,
    debug = debug,
    ...
  )
}

