"kd" <- function(covmat, n, type=c("exact","sample"))
{
  ## Kaiser-Dickman (1962) algorithm for generating sample data 
  ## based on the input covmat, which is a covariance matrix.
  ##
  ## n is desired sample size
  ## type="exact" returns data matrix that yields the exact covmat;
  ## type="sample" returns sample data, treating covmat as population matrix
  ##
  ## Returns the sample data matrix, dat

  ## Code written by Edgar Merkle, University of Missouri

  type <- match.arg(type)
  
  ## Check to ensure that covmat is a valid covariance matrix.
  if(nrow(covmat) != ncol(covmat))
    stop("non-square matrix supplied")
  symmetric <- isSymmetric.matrix(covmat)
  if(!symmetric)
    stop("non-symmetric matrix supplied")
  pd <- all(eigen(covmat, only.values=TRUE)$values > 0)
  if(!pd)
    stop("covariance matrix is not positive definite")
  
  p <- nrow(covmat)
    
  ## Algorithm works on a correlation matrix
  mv.vars <- matrix(0, nrow(covmat), nrow(covmat))
  diag(mv.vars) <- sqrt(diag(covmat))
  cormat <- cov2cor(covmat)

  ## Generate standard normal data and mean center each variable
  Xscore <- matrix(rnorm(p*n), p, n)
  Xsub0 <- t(apply(Xscore, 1, scale, scale=FALSE))

  ## Correlation matrix factored via Cholesky decomposition:
  Fcomp <- t(chol(cormat))

  ## Equation 2 from K&D:
  Zhat <- Fcomp %*% Xscore

  ## Equation 3 from K&D:
  Xsub0.prod <- Xsub0 %*% t(Xsub0)

  ## Get singular value decomp of Xsub0.prod
  Xsub0.svd <- svd(Xsub0.prod)
  M.sqrt <- matrix(0,p,p)
  diag(M.sqrt) <- 1/sqrt(Xsub0.svd$d)

  ## Equation 5 from K&D:
  Z <- Fcomp %*% M.sqrt %*% t(Xsub0.svd$u) %*% Xsub0
  Z <- Z*sqrt(n)

  dat <- Z
  if (type=="sample"){dat <- Zhat}

  ## Scale data to correspond to covmat
  dat <- t(dat) %*% mv.vars

  dat
}
