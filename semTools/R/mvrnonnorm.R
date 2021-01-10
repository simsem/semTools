### Yves Rosseel, Sunthud Pornprasertmanit, & Terrence D. Jorgensen
### Last updated: 10 January 2021


##' Generate Non-normal Data using Vale and Maurelli (1983) method
##'
##' Generate Non-normal Data using Vale and Maurelli (1983) method. The function
##' is designed to be as similar as the popular \code{mvrnorm} function in the
##' \code{MASS} package. The codes are copied from \code{mvrnorm} function in
##' the \code{MASS} package for argument checking and \code{lavaan} package for
##' data generation using Vale and Maurelli (1983) method.
##'
##'
##' @importFrom stats cov2cor
##'
##' @param n Sample size
##' @param mu A mean vector. If elements are named, those will be used as
##' variable names in the returned data matrix.
##' @param Sigma A positive-definite symmetric matrix specifying the covariance
##' matrix of the variables. If rows or columns are named (and \code{mu} is
##' unnamed), those will be used as variable names in the returned data matrix.
##' @param skewness A vector of skewness of the variables
##' @param kurtosis A vector of excessive kurtosis of the variables
##' @param empirical If \code{TRUE}, \code{mu} and \code{Sigma} specify the
##' empirical rather than population mean and covariance matrix
##' @return A data matrix
##' @author The original function is the \code{\link[lavaan]{simulateData}}
##' function written by Yves Rosseel in the \code{lavaan} package. The function
##' is adjusted for a convenient usage by Sunthud Pornprasertmanit
##' (\email{psunthud@@gmail.com}). Terrence D. Jorgensen added the feature to
##' retain variable names from \code{mu} or \code{Sigma}.
##'
##' @references Vale, C. D. & Maurelli, V. A. (1983). Simulating multivariate
##' nonormal distributions. \emph{Psychometrika, 48}(3), 465--471.
##' \doi{10.1007/BF02293687}
##'
##' @examples
##'
##' set.seed(123)
##' mvrnonnorm(20, c(1, 2), matrix(c(10, 2, 2, 5), 2, 2),
##' 	skewness = c(5, 2), kurtosis = c(3, 3))
##' ## again, with variable names specified in mu
##' set.seed(123)
##' mvrnonnorm(20, c(a = 1, b = 2), matrix(c(10, 2, 2, 5), 2, 2),
##' 	skewness = c(5, 2), kurtosis = c(3, 3))
##'
##' @export
mvrnonnorm <- function(n, mu, Sigma, skewness = NULL,
                       kurtosis = NULL, empirical = FALSE) {
  ## number of variables
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) stop("incompatible arguments")
  ## save variable names, if they exist
  varnames <- names(mu)
  if (is.null(varnames)) varnames <- rownames(Sigma)
  if (is.null(varnames)) varnames <- colnames(Sigma)
  ## check for NPD
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -1e-06 * abs(ev[1L])))
    stop("'Sigma' is not positive definite")
  ## simulate X <- NULL
  if (is.null(skewness) && is.null(kurtosis)) {
    X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma, empirical = empirical)
  } else {
    if (empirical) warning(c("The empirical argument does not work when the ",
                             "Vale and Maurelli's method is used."))
    if (is.null(skewness)) skewness <- rep(0, p)
    if (is.null(kurtosis)) kurtosis <- rep(0, p)
    Z <- ValeMaurelli1983copied(n = n, COR = cov2cor(Sigma),
                                skewness = skewness, kurtosis = kurtosis)
    TMP <- scale(Z, center = FALSE, scale = 1/sqrt(diag(Sigma)))[ , , drop = FALSE]
    X <- sweep(TMP, MARGIN = 2, STATS = mu, FUN = "+")
  }
  colnames(X) <- varnames
  X
}


## ----------------
## Hidden Functions
## ----------------

## Copied from lavaan package
##' @importFrom stats nlminb
ValeMaurelli1983copied <- function(n = 100L, COR, skewness, kurtosis,
                                   debug = FALSE) {

    fleishman1978_abcd <- function(skewness, kurtosis) {
        system.function <- function(x, skewness, kurtosis) {
            b.=x[1L]; c.=x[2L]; d.=x[3L]
            eq1 <- b.^2 + 6*b.*d. + 2*c.^2 + 15*d.^2 - 1
            eq2 <- 2*c.*(b.^2 + 24*b.*d. + 105*d.^2 + 2) - skewness
            eq3 <- 24*(b.*d. + c.^2*(1 + b.^2 + 28*b.*d.) +
                       d.^2*(12 + 48*b.*d. + 141*c.^2 + 225*d.^2)) - kurtosis
            eq <- c(eq1,eq2,eq3)
            sum(eq^2) ## SS
        }

        out <- nlminb(start = c(1, 0, 0), objective = system.function,
                      scale = 10, control = list(trace = 0),
                      skewness = skewness, kurtosis = kurtosis)
        if(out$convergence != 0) warning("no convergence")
        b. <- out$par[1L]; c. <- out$par[2L]; d. <- out$par[3L]; a. <- -c.
        c(a.,b.,c.,d.)
    }

    getICOV <- function(b1, c1, d1, b2, c2, d2, R) {
        objectiveFunction <- function(x, b1, c1, d1, b2, c2, d2, R) {
            rho=x[1L]
            eq <- rho*(b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2) +
                  rho^2*(2*c1*c2) + rho^3*(6*d1*d2) - R
            eq^2
        }

        #gradientFunction <- function(x, bcd1, bcd2, R) {
        #
        #}

        out <- nlminb(start=R, objective=objectiveFunction,
                      scale=10, control=list(trace=0),
                      b1=b1, c1=c1, d1=d1, b2=b2, c2=c2, d2=d2, R=R)
        if(out$convergence != 0) warning("no convergence")
        rho <- out$par[1L]
        rho
    }

    # number of variables
    nvar <- ncol(COR)
    # check skewness
    if(is.null(skewness)) {
        SK <- rep(0, nvar)
    } else if(length(skewness) == nvar) {
        SK <- skewness
    } else if(length(skewness) == 1L) {
        SK <- rep(skewness, nvar)
    } else {
        stop("skewness has wrong length")
    }

    if(is.null(kurtosis)) {
        KU <- rep(0, nvar)
    } else if(length(kurtosis) == nvar) {
        KU <- kurtosis
    } else if(length(kurtosis) == 1L) {
        KU <- rep(kurtosis, nvar)
    } else {
        stop("kurtosis has wrong length")
    }

    # create Fleishman table
    FTable <- matrix(0, nvar, 4L)
    for(i in 1:nvar) {
        FTable[i,] <- fleishman1978_abcd(skewness=SK[i], kurtosis=KU[i])
    }

    # compute intermediate correlations between all pairs
    ICOR <- diag(nvar)
    for(j in 1:(nvar-1L)) {
        for(i in (j+1):nvar) {
            if(COR[i,j] == 0) next
            ICOR[i,j] <- ICOR[j,i] <-
                getICOV(FTable[i,2], FTable[i,3], FTable[i,4],
                        FTable[j,2], FTable[j,3], FTable[j,4], R=COR[i,j])
        }
    }

    if(debug) {
         cat("\nOriginal correlations (for Vale-Maurelli):\n")
         print(COR)
         cat("\nIntermediate correlations (for Vale-Maurelli):\n")
         print(ICOR)
         cat("\nEigen values ICOR:\n")
         print( eigen(ICOR)$values )
    }

    # generate Z ## FIXME: replace by rmvnorm once we use that package
    X <- Z <- MASS::mvrnorm(n=n, mu=rep(0,nvar), Sigma=ICOR)

    # transform Z using Fleishman constants
    for(i in 1:nvar) {
        X[,i] <- FTable[i,1L] + FTable[i,2L]*Z[,i] + FTable[i,3L]*Z[,i]^2 +
                 FTable[i,4L]*Z[,i]^3
    }

    X
}

