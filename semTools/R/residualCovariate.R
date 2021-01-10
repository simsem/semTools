### Sunthud Pornprasertmanit
### Last updated: 10 January 2021


##' Residual-center all target indicators by covariates
##'
##' This function will regress target variables on the covariate and replace the
##' target variables by the residual of the regression analysis. This procedure
##' is useful to control the covariate from the analysis model (Geldhof,
##' Pornprasertmanit, Schoemann, & Little, 2013).
##'
##'
##' @importFrom stats lm
##'
##' @param data The desired data to be transformed.
##' @param targetVar Varible names or the position of indicators that users wish
##' to be residual centered (as dependent variables)
##' @param covVar Covariate names or the position of the covariates using for
##' residual centering (as independent variables) onto target variables
##'
##' @return The data that the target variables replaced by the residuals
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##' @seealso \code{\link{indProd}} For creating the indicator products with no
##' centering, mean centering, double-mean centering, or residual centering.
##'
##' @references Geldhof, G. J., Pornprasertmanit, S., Schoemann, A. M., &
##' Little, T. D. (2013). Orthogonalizing through residual centering:
##' Extended applications and caveats. \emph{Educational and Psychological
##' Measurement, 73}(1), 27--46. \doi{10.1177/0013164412445473}
##'
##' @examples
##'
##' dat <- residualCovariate(attitude, 2:7, 1)
##'
##' @export
residualCovariate <- function(data, targetVar, covVar) {
    x <- as.list(match.call())
    cov <- eval(x$covVar)
    target <- eval(x$targetVar)
    if (all(is.numeric(cov))) cov <- colnames(data)[cov]
    if (all(is.numeric(target))) target <- colnames(data)[target]
    express <- paste("cbind(", paste(target, collapse = ", "), ") ~ ",
                     paste(cov, collapse = " + "), sep = "")
    data[, target] <- lm(express, data = data)$residuals
    return(data)
}


