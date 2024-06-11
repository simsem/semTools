### Ylenio Longo
### Last updated: 10 January 2021

##' Empirical Kaiser criterion
##'
##' Identify the number of factors to extract based on the Empirical Kaiser
##' Criterion (EKC). The analysis can be run on a `data.frame` or data
##' `matrix` (`data`), or on a correlation or covariance matrix
##' (`sample.cov`) and the sample size (`sample.nobs`). A
##' `data.frame` is returned with two columns: the eigenvalues from your
##' data or covariance matrix and the reference eigenvalues. The number of
##' factors suggested by the Empirical Kaiser Criterion (i.e. the sample
##' eigenvalues greater than the reference eigenvalues), and the number of
##' factors suggested by the original Kaiser Criterion
##' (i.e. sample eigenvalues > 1) is printed above the output.
##'
##'
##' @importFrom stats cov cov2cor
##'
##' @param data A `data.frame` or data `matrix` containing columns of
##'   variables to be factor-analyzed.
##' @param sample.cov A covariance or correlation matrix can be used, instead of
##'   `data`, to estimate the eigenvalues.
##' @param sample.nobs Number of observations (i.e. sample size) if
##'   `is.null(data)` and `sample.cov=` is used.
##' @param missing If `"listwise"`, incomplete cases are removed listwise from
##'   the `data.frame`. If `"direct"` or `"ml"` or `"fiml"` and the `estimator=`
##'   is maximum likelihood, an EM algorithm is used to estimate an unrestricted
##'   covariance matrix (and mean vector). If `"pairwise"`, pairwise deletion is
##'   used. If `"default"``, the value is set depending on the estimator and the
##'   mimic option (see [lavaan::lavCor()] for details).
##' @param ordered `character` vector. Only used if object is a `data.frame`.
##'   Treat these variables as `ordered=` (ordinal) variables. Importantly, all
##'   other variables will be treated as numeric (unless `is.ordered == TRUE` in
##'   `data`). (see also [lavCor][lavaan::lavCor])
##' @param plot logical. Whether to print a scree plot comparing the sample
##'   eigenvalues with the reference eigenvalues.
##' @return A `data.frame` showing the sample and reference eigenvalues.
##'
##' The number of factors suggested by the Empirical Kaiser Criterion (i.e. the
##' sample eigenvalues greater than the reference eigenvalues) is returned as an
##' attribute (see **Examples**).
##'
##' The number of factors suggested by the original Kaiser Criterion (i.e.
##' sample eigenvalues > 1) is also printed as a header to the `data.frame`
##'
##' @author Ylenio Longo (University of Nottingham;
##' \email{yleniolongo@@gmail.com})
##'
##' Terrence D. Jorgensen (University of Amsterdam;
##' \email{TJorgensen314@@gmail.com})
##'
##' @references Braeken, J., & van Assen, M. A. L. M. (2017). An empirical
##' Kaiser criterion. *Psychological Methods, 22*(3), 450--466.
##' \doi{10.1037/met0000074}
##'
##' @examples
##'
##' ## Simulate data with 3 factors
##' model <- '
##'   f1 =~ .3*x1 + .5*x2 + .4*x3
##'   f2 =~ .3*x4 + .5*x5 + .4*x6
##'   f3 =~ .3*x7 + .5*x8 + .4*x9
##' '
##' dat <- simulateData(model, seed = 123)
##' ## save summary statistics
##' myCovMat <- cov(dat)
##' myCorMat <- cor(dat)
##' N <- nrow(dat)
##'
##' ## Run the EKC function
##' (out <- efa.ekc(dat))
##'
##' ## To extract the recommended number of factors using the EKC:
##' attr(out, "nfactors")
##'
##' ## If you do not have raw data, you can use summary statistics
##' (x1 <- efa.ekc(sample.cov = myCovMat, sample.nobs = N, plot = FALSE))
##' (x2 <- efa.ekc(sample.cov = myCorMat, sample.nobs = N, plot = FALSE))
##'
##' @export
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

