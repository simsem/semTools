### Edgar Merkle
### Last updated: 10 January 2021
### Kaiser-Dickman (1962) algorithm for generating sample data
### based on the input covmat, which is a covariance matrix.


##' Generate data via the Kaiser-Dickman (1962) algorithm.
##'
##' Given a covariance matrix and sample size, generate raw data that correspond
##' to the covariance matrix.  Data can be generated to match the covariance
##' matrix exactly, or to be a sample from the population covariance matrix.
##'
##' By default, R's `cov()` function divides by `n`-1.  The data
##' generated by this algorithm result in a covariance matrix that matches
##' `covmat`, but you must divide by `n` instead of `n`-1.
##'
##'
##' @importFrom stats cov2cor rnorm
##'
##' @param covmat a symmetric, positive definite covariance matrix
##' @param n the sample size for the data that will be generated
##' @param type type of data generation. `exact` generates data that
##' exactly correspond to `covmat`.  `sample` treats `covmat` as
##' a poulation covariance matrix, generating a sample of size `n`.
##'
##' @return `kd` returns a data matrix of dimension `n` by
##' `nrow(covmat)`.
##'
##' @author Ed Merkle (University of Missouri; \email{merklee@@missouri.edu})
##'
##' @references Kaiser, H. F. and Dickman, K. (1962).  Sample and population
##' score matrices and sample correlation matrices from an arbitrary population
##' correlation matrix.  *Psychometrika, 27*(2), 179--182.
##' \doi{10.1007/BF02289635}
##'
##' @examples
##'
##' #### First Example
##'
##' ## Get data
##' dat <- HolzingerSwineford1939[ , 7:15]
##' hs.n <- nrow(dat)
##'
##' ## Covariance matrix divided by n
##' hscov <- ((hs.n-1)/hs.n) * cov(dat)
##'
##' ## Generate new, raw data corresponding to hscov
##' newdat <- kd(hscov, hs.n)
##'
##' ## Difference between new covariance matrix and hscov is minimal
##' newcov <- (hs.n-1)/hs.n * cov(newdat)
##' summary(as.numeric(hscov - newcov))
##'
##' ## Generate sample data, treating hscov as population matrix
##' newdat2 <- kd(hscov, hs.n, type = "sample")
##'
##' #### Another example
##'
##' ## Define a covariance matrix
##' covmat <- matrix(0, 3, 3)
##' diag(covmat) <- 1.5
##' covmat[2:3,1] <- c(1.3, 1.7)
##' covmat[3,2] <- 2.1
##' covmat <- covmat + t(covmat)
##'
##' ## Generate data of size 300 that have this covariance matrix
##' rawdat <- kd(covmat, 300)
##'
##' ## Covariances are exact if we compute sample covariance matrix by
##' ## dividing by n (vs by n - 1)
##' summary(as.numeric((299/300)*cov(rawdat) - covmat))
##'
##' ## Generate data of size 300 where covmat is the population covariance matrix
##' rawdat2 <- kd(covmat, 300)
##'
##' @export
kd <- function(covmat, n, type=c("exact","sample")) {
  type <- match.arg(type)

  ## Check to ensure that covmat is a valid covariance matrix.
  if (nrow(covmat) != ncol(covmat)) stop("non-square matrix supplied")
  symmetric <- isSymmetric.matrix(covmat)
  if (!symmetric) stop("non-symmetric matrix supplied")
  pd <- all(eigen(covmat, only.values = TRUE)$values > 0)
  if (!pd) stop("covariance matrix is not positive definite")

  p <- nrow(covmat)

  ## Algorithm works on a correlation matrix
  mv.vars <- matrix(0, nrow(covmat), nrow(covmat))
  diag(mv.vars) <- sqrt(diag(covmat))
  cormat <- cov2cor(covmat)

  ## Generate standard normal data and mean center each variable
  Xscore <- matrix(rnorm(p*n), p, n)
  Xsub0 <- t(apply(Xscore, 1, scale, scale = FALSE))

  ## Correlation matrix factored via Cholesky decomposition:
  Fcomp <- t(chol(cormat))

  ## Equation 2 from K&D:
  Zhat <- Fcomp %*% Xscore

  ## Equation 3 from K&D:
  Xsub0.prod <- Xsub0 %*% t(Xsub0)

  ## Get singular value decomp of Xsub0.prod
  Xsub0.svd <- svd(Xsub0.prod)
  M.sqrt <- matrix(0, p, p)
  diag(M.sqrt) <- 1 / sqrt(Xsub0.svd$d)

  ## Equation 5 from K&D:
  Z <- Fcomp %*% M.sqrt %*% t(Xsub0.svd$u) %*% Xsub0
  Z <- Z * sqrt(n)

  dat <- Z
  if (type == "sample") { dat <- Zhat }

  ## Scale data to correspond to covmat
  dat <- t(dat) %*% mv.vars

  ## convert to data.frame, use any existing names from covmat
  dat <- data.frame(dat)
  if(!is.null(colnames(covmat))) names(dat) <- colnames(covmat)

  dat
}


