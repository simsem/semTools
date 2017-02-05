### Ylenio Longo
### Last updated: 20 January 2017

efa.ekc <- function(data = NULL, sample.cov = NULL, sample.nobs = NULL,
                    missing = "default", ordered = NULL, plot = TRUE) {
  ## if data
  if (!is.null(data)) {
    data <- as.data.frame(data)
    R <- lavaan::lavCor(data, missing = missing, ordered = ordered)  #correlations
    j <- dim(data)[2]  #number of variables
    n <- dim(data)[1]  #sample size
  } else {
    ## if covariance matrix
    if (max(diag(sample.cov)) != 1 & min(diag(sample.cov)) != 1) {
      R <- cov2cor(sample.cov)
      j <- dim(R)[2]  #number of variables
      n <- sample.nobs  #sample size
    } else {
      ## if correlation matrix
      R <- sample.cov
      j <- dim(R)[2]  #number of variables
      n <- sample.nobs  #sample size
    }
  }
  g <- j/n  #gamma: var / sample
  l <- (1 + sqrt(g))^2  #1st reference eigenvalue
  e <- eigen(R)$values  #eigenvalues

  v <- cumsum(e)  #Define cumulatively summed eigenvalue vector
  v1 <- v[1:j - 1]  #omit last element
  v2 <- c(0, v1)  #put a zero upfront
  w <- sort(1:j, decreasing = TRUE)  #eigenvalue order vector
  ref <- (((j - v2)/w) * l)  #EKC reference eigenvalues

  # results
  Eigenvalues <- data.frame(Sample = e, Ref = ref)  #sample and reference eigenvalues
  rownames(Eigenvalues) <- 1:j
  class(Eigenvalues) <- c("lavaan.data.frame","data.frame")
  ## add no. factors to extract as attribute, using each criterion
  nfactors_EKC <- which(!(Eigenvalues[, 1] > Eigenvalues[, 2]))[1] - 1 # EKC
  nfactors_KC <- which(!(Eigenvalues[, 1] > 1))[1] - 1  # Kaiser Criterion
  attr(Eigenvalues, "header") <- paste(" Empirical Kaiser Criterion suggests",
                                       nfactors_EKC, "factors.\n",
                                       "Traditional Kaiser Criterion suggests",
                                       nfactors_KC, "factors.")
  attr(Eigenvalues, "nfactors") <- nfactors_EKC
  if (plot) {
    plot(Eigenvalues[, 1], type = "b", pch = 20, cex = 0.9, col = "black",
         main = "Empirical Kaiser Criterion\nScree Plot", ylab = "Eigenvalues",
         ylim = c(min(Eigenvalues), max(ceiling(Eigenvalues))),
         xlab = "Factor Number", xlim = c(1, j))
    lines(Eigenvalues[, 2], lty = "dashed", col = "blue")
    legend("topright", c("  Data", "  Empirical\n  Reference", "  Kaiser Criterion"),
           col = c("black","blue","gray"), bty = "n",
           pch = c(20, NA, NA), lty = c("solid","dashed","solid"), merge = TRUE)
    abline(h = 1, col = "gray")  # Kaiser Criterion
  }
  return(Eigenvalues)
}

