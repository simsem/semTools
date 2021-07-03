### Terrence D. Jorgensen & Andrew R. Johnson
### Last updated: 3 July 2021
### function to derive ordinal-scale moments implied by LRV-scale moments


##' Calculate Population Moments for Ordinal Data Treated as Numeric
##'
##' This function calculates ordinal-scale moments implied by LRV-scale moments
##'
##' Binary and ordinal data are frequently accommodated in SEM by incorporating
##' a threshold model that links each observed categorical response variable to
##' a corresponding latent response variable that is typically assumed to be
##' normally distributed (Kamata & Bauer, 2008; Wirth & Edwards, 2007).
##'
##' @importFrom stats dnorm setNames
##' @importFrom lavaan lavInspect
##' @importFrom pbivnorm pbivnorm
##'
##' @param Sigma Population covariance \code{\link{matrix}}, with variable names
##'   saved in the \code{\link{dimnames}} attribute.
##' @param Mu Optional \code{numeric} vector of population means. If missing,
##'   all means will be set to zero.
##' @param thresholds Either a single \code{numeric} vector of population
##'   thresholds used to discretize each normally distributed variable, or a
##'   named \code{list} of each discretized variable's vector of thresholds.
##'   The discretized variables may be a subset of all variables in \code{Sigma}
##'   if the remaining variables are intended to be observed rather than latent
##'   normally distributed variables.
##' @param cWts Optional (default when missing is to use 0 for the lowest
##'   category, followed by successive integers for each higher category).
##'   Either a single \code{numeric} vector of category weights (if they are
##'   identical across all variables) or a named \code{list} of each
##'   discretized variable's vector of category weights.
##'
##' @return A \code{list} including the LRV-scale population moments (means,
##'   covariance matrix, correlation matrix, and thresholds), the category
##'   weights, a \code{data.frame} of implied univariate moments (means,
##'   \emph{SD}s, skewness, and excess kurtosis (i.e., in excess of 3, which is
##'   the kurtosis of the normal distribution) for discretized data treated as
##'   \code{numeric}, and the implied covariance and correlation matrix of
##'   discretized data treated as \code{numeric}.
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##'   Andrew Johnson (Curtin University; \email{andrew.johnson@@curtin.edu.au})
##'
##' @references
##'
##' Kamata, A., & Bauer, D. J. (2008). A note on the relation between factor
##'   analytic and item response theory models.
##'   \emph{Structural Equation Modeling, 15}(1), 136--153.
##'   \doi{10.1080/10705510701758406}
##'
##' Wirth, R. J., & Edwards, M. C. (2007). Item factor analysis: Current
##'   approaches and future directions. \emph{Psychological Methods, 12}(1),
##'   58--79. \doi{10.1037/1082-989X.12.1.58}
##'
##' @examples
##'
##' ## SCENARIO 1: DIRECTLY SPECIFY POPULATION PARAMETERS
##'
##' ## specify population model in LISREL matrices
##' Nu <- rep(0, 4)
##' Alpha <- c(1, -0.5)
##' Lambda <- matrix(c(1, 1, 0, 0, 0, 0, 1, 1), nrow = 4, ncol = 2,
##'                  dimnames = list(paste0("y", 1:4), paste0("eta", 1:2)))
##' Psi <- diag(c(1, .75))
##' Theta <- diag(4)
##' Beta <- matrix(c(0, .5, 0, 0), nrow = 2, ncol = 2)
##'
##' ## calculate model-implied population means and covariance matrix
##' ## of latent response variables (LRVs)
##' IB <- solve(diag(2) - Beta) # to save time and space
##' Mu_LRV <- Nu + Lambda %*% IB %*% Alpha
##' Sigma_LRV <- Lambda %*% IB %*% Psi %*% t(IB) %*% t(Lambda) + Theta
##'
##' ## Specify (unstandardized) thresholds to discretize normally distributed data
##' ## generated from Mu_LRV and Sigma_LRV, based on marginal probabilities
##' PiList <- list(y1 = c(.25, .5, .25),
##'                y2 = c(.17, .33, .33, .17),
##'                y3 = c(.1, .2, .4, .2, .1),
##'                ## make final variable highly asymmetric
##'                y4 = c(.33, .25, .17, .12, .08, .05))
##' sapply(PiList, sum) # all sum to 100%
##' CumProbs <- sapply(PiList, cumsum)
##' ## unstandardized thresholds
##' TauList <- mapply(qnorm, p = lapply(CumProbs, function(x) x[-length(x)]),
##'                   m = Mu_LRV, sd = sqrt(diag(Sigma_LRV)))
##' for (i in 1:4) names(TauList[[i]]) <- paste0(names(TauList)[i], "|t",
##'                                              1:length(TauList[[i]]))
##'
##' ## assign numeric weights to each category (optional, see default)
##' NumCodes <- list(y1 = c(-0.5, 0, 0.5), y2 = 0:3, y3 = 1:5, y4 = 1:6)
##'
##'
##' ## Calculate Population Moments for Numerically Coded Ordinal Variables
##' lrv2ord(Sigma = Sigma_LRV, Mu = Mu_LRV, thresholds = TauList, cWts = NumCodes)
##'
##'
##' ## SCENARIO 2: USE ESTIMATED PARAMETERS AS POPULATION
##'
##' data(datCat) # already stored as c("ordered","factor")
##' fit <- cfa(' f =~ 1*u1 + 1*u2 + 1*u3 + 1*u4 ', data = datCat)
##' lrv2ord(Sigma = fit, thresholds = fit) # use same fit for both
##' ## or use estimated thresholds with specified parameters, but note that
##' ## lrv2ord() will only extract standardized thresholds
##' dimnames(Sigma_LRV) <- list(paste0("u", 1:4), paste0("u", 1:4))
##' lrv2ord(Sigma = cov2cor(Sigma_LRV), thresholds = fit)
##'
##' @export
lrv2ord <- function(Sigma, Mu, thresholds, cWts) {
  if (inherits(Sigma, "lavaan")) {
    if (lavInspect(Sigma, "ngroups") > 1L || lavInspect(Sigma, "nlevels") > 1L) {
      stop('Sigma= only accepts single-group/level lavaan models')
    }
    fitSigma <- Sigma
    Sigma <- lavInspect(fitSigma, "cov.ov")
  } else stopifnot(is.matrix(Sigma))
  vn <- rownames(Sigma) # variable names
  SDs <- sqrt(diag(Sigma))

  if (missing(Mu)) {
    Mu <- rep(0, nrow(Sigma))
  } else if (inherits(Mu, "lavaan")) {
    if (lavInspect(Mu, "ngroups") > 1L || lavInspect(Mu, "nlevels") > 1L) {
      stop('Mu= only accepts single-group/level lavaan models')
    }
    fitMu <- Mu
    Mu <- lavInspect(fitMu, "mean.ov")
  }
  names(Mu) <- names(SDs) <- vn


  ## If a single vector of thresholds is passed, broadcast to a list
  if (inherits(thresholds, "lavaan")) {
    if (lavInspect(thresholds, "ngroups") > 1L || lavInspect(thresholds, "nlevels") > 1L) {
      stop('thresholds= only accepts single-group/level lavaan models')
    }
    ## check whether diag(Sigma) == 1
    isSTD <- sapply(SDs, function(x) {
      isTRUE(all.equal(x, current = 1, tolerance = .001))
    })
    if (!all(isSTD)) warning('standardized thresholds= extracted from a ',
                             'lavaan object, but Sigma= is not a ',
                             'correlation matrix.')
    fitThr <- thresholds
    thresholds <- lavInspect(fitThr, "th") # STANDARDIZED thresholds

    thresh <- lapply(unique(lavInspect(fitThr, "th.idx")), function(x) {
      thresholds[lavInspect(fitThr, "th.idx") == x]
    })
    names(thresh) <- sapply(thresh, function(x) {
      strsplit(names(x)[1], "|t", fixed = TRUE)[[1]][1]
    })

  } else if (is.atomic(thresholds)) {
    thresh <- sapply(vn, function(x) {thresholds}, simplify = FALSE)
  } else {
    stopifnot(is.list(thresholds)) # must be a list
    stopifnot(length(thresholds) <= nrow(Sigma)) # no more than 1 per variable
    stopifnot(all(names(thresholds) %in% vn)) # names must match
    thresh <- thresholds
  }
  cn <- names(thresh)

  ## If no category weights are passed, default to 0:nCat
  if (missing(cWts)) {
    cWts <- sapply(thresh, function(x) { 0:length(x) }, simplify = FALSE)
  } else if (is.atomic(cWts)) {
    ## If a single vector of category weights is passed, broadcast to a list
    #FIXME: assumes same number of thresholds across variables
    cWts <- sapply(cn, function(x) { cWts }, simplify = FALSE)
  } else {
    stopifnot(is.list(cWts)) # must be a list
    stopifnot(length(cWts) <= nrow(Sigma)) # no more than 1 per variable
    stopifnot(all(names(cWts) %in% vn)) # names must match
    stopifnot(all(cn %in% names(cWts))) # names must match
    cWts <- cWts[cn] # discard any others
  }
  stopifnot(all((sapply(thresh, length) + 1L) == sapply(cWts, length)))

  ## Calculate marginal probabilities implied by thresholds on moments
  get_marg_probs <- function(threshs, m, sd) {
    thr  <- c(-Inf, threshs, Inf)
    sapply(2:length(thr), function(k) {
      pnorm(thr[k], m, sd) - pnorm(thr[k-1], m, sd)
    })
  }
  marginal_probs <- mapply(get_marg_probs, SIMPLIFY = FALSE, threshs = thresh,
                           m = Mu[cn], sd = SDs[cn])

  ## Marginal means
  Mu_ord <- Mu
  Mu_ord[cn] <- mapply(function(p, w) {
    stopifnot(length(p) == length(w))
    sum(p * w)
  }, p = marginal_probs, w = cWts)

  ## marginal variances (fill in covariances below)
  Sigma_ord <- Sigma
  diag(Sigma_ord[cn,cn]) <- mapply(function(p, w, mu) {
    stopifnot(length(p) == length(w))
    sum(p * (w - mu)^2)
  }, p = marginal_probs, w = cWts, mu = Mu_ord[cn])

  ## marginal (standardized) skew
  skew_ord <- setNames(rep(0, nrow(Sigma)), nm = vn)
  skew_ord[cn] <- mapply(function(p, w, mu) {
    stopifnot(length(p) == length(w))
    numerator <- sum(p * (w - mu)^3)
    Moment2 <- sum(p * (w - mu)^2)
    denominator <- sqrt(Moment2)^3
    numerator / denominator
  }, p = marginal_probs, w = cWts, mu = Mu_ord[cn])

  ## marginal (standardized, excess) kurtosis
  kurt_ord <- setNames(rep(0, nrow(Sigma)), nm = vn)
  kurt_ord[cn] <- mapply(function(p, w, mu) {
    stopifnot(length(p) == length(w))
    numerator <- sum(p * (w - mu)^4)
    Moment2 <- sum(p * (w - mu)^2)
    denominator <- sqrt(Moment2)^4
    numerator / denominator
  }, p = marginal_probs, w = cWts, mu = Mu_ord[cn]) - 3 # excess kurtosis

  ## all marginal descriptives
  (margMoments <- data.frame(Mean = Mu_ord, SD = sqrt(diag(Sigma_ord)),
                             Skew = skew_ord, Kurtosis3 = kurt_ord,
                             row.names = vn))
  class(margMoments) <- c("lavaan.data.frame","data.frame") # for printing

  ## save old copies to return with new
  out <- list(Mu_LRV = Mu, Sigma_LRV = Sigma, R_LRV = stats::cov2cor(Sigma),
              Thresholds = thresh, Category_weights = cWts, Uni_ord = margMoments)
  class(out$Mu_LRV)    <- c("lavaan.vector","numeric")
  class(out$Sigma_LRV) <- c("lavaan.matrix.symmetric","matrix")
  class(out$R_LRV)     <- c("lavaan.matrix.symmetric","matrix")
  out$Thresholds       <- lapply(out$Thresholds, "class<-",
                                 c("lavaan.vector","numeric"))
  out$Category_weights <- lapply(out$Category_weights, "class<-",
                                 c("lavaan.vector","numeric"))


  ## function to apply to any pair of indicators (i and j) in Sigma
  getOrdCov <- function(i, j) {
    ## to use apply(), i= can be 2 values indicating the [row, column]
    if (length(i) > 1L) {
      if (!missing(j)) warning("j ignored when i has multiple values")
      if (length(i) > 2L) stop("i had ", length(i), " values. Only the first 2 were used.")
      j <- i[2]
      i <- i[1]
    }
    ## if i/j are numeric, get names
    if (is.numeric(i)) i <- vn[i]
    if (is.numeric(j)) j <- vn[j]
    ## make sure thresholds are standardized
    i.thr <-
    j.thr <-

    ## template for matrices of joint probabilities and cross-products
    JointProbs <- CP <- matrix(0, nrow = length(cWts[[i]]),
                                  ncol = length(cWts[[j]]))
    i.thr <- c(-1e5, (thresh[[i]] - Mu[i]) / SDs[i], 1e5)
    j.thr <- c(-1e5, (thresh[[j]] - Mu[j]) / SDs[j], 1e5)
    tCombos <- cbind(expand.grid(i = i.thr, j = j.thr),
                     expand.grid(cat1 = c(0, seq_along(cWts[[i]])),
                                 cat2 = c(0, seq_along(cWts[[j]]))))
    tCombos$cp <- pbivnorm(x = tCombos$i, y = tCombos$j, rho = out$R_LRV[i,j])
    ## loop over rows & columns
    for (RR in seq_along(cWts[[i]])) for (CC in seq_along(cWts[[j]])) {
      ## calculate joint probabilities
      idx1 <- which(tCombos$cat1 == RR     & tCombos$cat2 == CC    )
      idx2 <- which(tCombos$cat1 == RR - 1 & tCombos$cat2 == CC    )
      idx3 <- which(tCombos$cat1 == RR     & tCombos$cat2 == CC - 1)
      idx4 <- which(tCombos$cat1 == RR - 1 & tCombos$cat2 == CC - 1)
      JointProbs[RR,CC] <- tCombos$cp[idx1] - tCombos$cp[idx2] - tCombos$cp[idx3] + tCombos$cp[idx4]
      ## calculate cross-products
      CP[RR,CC] <- (cWts[[i]][RR] - Mu_ord[i]) * (cWts[[j]][CC] - Mu_ord[j])
    }
    sum(JointProbs * CP) # return covariance
  }

  ## check whether all variables are being discretized
  stayCon <- setdiff(vn, cn)
  if (length(stayCon) == 0) {
    ## all are polychoric
    (ij <- which(lower.tri(Sigma_ord), arr.ind = TRUE))
    Sigma_ord[ij] <- mapply(getOrdCov, i = ij[,1], j = ij[,2])
    Sigma_ord[ ij[,2:1] ] <- Sigma_ord[ij] # copy lower to upper triangle

  } else {
    ## pair by pair, choose polychoric or polyserial
    for (i in vn[-length(vn)]) for (j in vn[(which(vn == i)+1):length(vn)]) {
      if (i %in% stayCon && j %in% stayCon) next
      if (j %in% cn && j %in% cn) {
        ## both discretized, calculate polychoric
        Sigma_ord[i,j] <- Sigma_ord[j,i] <- getOrdCov(i, j)
        next
      }
      ## else, calculate polyserial
      if (i %in% stayCon) {
        CON <- i
        CAT <- j
      } else {
        CAT <- i
        CON <- j
      }
      DENS <- mapply(function(tau, interval, m = 0, sd = 1) {
        dnorm(tau, mean = m, sd = sd) * interval
      }, tau = thresh[[CAT]], interval = diff(cWts[[CAT]]),
         m = Mu[CAT], sd = SDs[CAT])
      ## Note: polyserial correlation divides by sqrt(diag(Sigma_ord)[CAT]),
      ##       but that cancels out when scaling by both SDs to get covariance
      Sigma_ord[CON, CAT] <- Sigma_ord[CAT, CON] <-
        out$R_LRV[CON, CAT] * sum(DENS) * sqrt(diag(out$Sigma_LRV)[CON])

    }
  }

  R_ord <- cov2cor(Sigma_ord)
  class(Sigma_ord) <- c("lavaan.matrix.symmetric","matrix")
  class(R_ord)     <- c("lavaan.matrix.symmetric","matrix")

  c(out, list(Sigma_ord = Sigma_ord, R_ord = R_ord))
}







