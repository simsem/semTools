### Ylenio Longo & Terrence D. Jorgensen
### Last updated: 1 September 2021

##' Assessing Discriminant Validity using Heterotrait--Monotrait Ratio
##'
##' This function assesses discriminant validity through the
##' heterotrait-monotrait ratio (HTMT) of the correlations (Henseler, Ringlet &
##' Sarstedt, 2015). Specifically, it assesses the arithmetic (Henseler et al.,
##' ) or geometric (Roemer et al., 2021) mean correlation
##' among indicators across constructs (i.e. heterotrait--heteromethod
##' correlations) relative to the geometric-mean correlation among indicators
##' within the same construct (i.e. monotrait--heteromethod correlations).
##' The resulting HTMT(2) values are interpreted as estimates of inter-construct
##' correlations. Absolute values of the correlations are recommended to
##' calculate the HTMT matrix, and are required to calculate HTMT2. Correlations
##' are estimated using the \code{\link[lavaan]{lavCor}} function.
##'
##'
##' @importFrom stats cov2cor
##'
##' @param model lavaan \link[lavaan]{model.syntax} of a confirmatory factor
##'   analysis model where at least two factors are required for indicators
##'   measuring the same construct.
##' @param data A \code{data.frame} or data \code{matrix}
##' @param sample.cov A covariance or correlation matrix can be used, instead of
##'   \code{data}, to estimate the HTMT values.
##' @param missing If "listwise", cases with missing values are removed listwise
##'   from the data frame. If "direct" or "ml" or "fiml" and the estimator is
##'   maximum likelihood, an EM algorithm is used to estimate the unrestricted
##'   covariance matrix (and mean vector). If "pairwise", pairwise deletion is
##'   used. If "default", the value is set depending on the estimator and the
##'   mimic option (see details in \link[lavaan]{lavCor}).
##' @param ordered Character vector. Only used if object is a \code{data.frame}.
##'   Treat these variables as ordered (ordinal) variables. Importantly, all
##'   other variables will be treated as numeric (unless \code{is.ordered} in
##'   \code{data}). (see also \link[lavaan]{lavCor})
##' @param absolute \code{logical} indicating whether HTMT values should be
##'   estimated based on absolute correlations (default is \code{TRUE}). This
##'   is recommended for HTMT but required for HTMT2 (so silently ignored).
##' @param htmt2 \code{logical} indicating whether to use the geometric mean
##'   (default, appropriate for congeneric indicators) or arithmetic mean
##'   (which assumes tau-equivalence).
##'
##' @return A matrix showing HTMT(2) values (i.e., discriminant validity)
##'   between each pair of factors.
##'
##' @author
##' Ylenio Longo (University of Nottingham; \email{yleniolongo@@gmail.com})
##'
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Henseler, J., Ringle, C. M., & Sarstedt, M. (2015). A new criterion for
##'   assessing discriminant validity in variance-based structural equation
##'   modeling. \emph{Journal of the Academy of Marketing Science, 43}(1),
##'   115--135. \doi{10.1007/s11747-014-0403-8}
##'
##'   Roemer, E., Schuberth, F., & Henseler, J. (2021). HTMT2-An improved
##'   criterion for assessing discriminant validity in structural equation
##'   modeling. \emph{Industrial Management & Data Systems}.
##'   \doi{10.1108/IMDS-02-2021-0082}
##'
##'   Voorhees, C. M., Brady, M. K., Calantone, R., & Ramirez, E. (2016).
##'   Discriminant validity testing in marketing: An analysis, causes for
##'   concern, and proposed remedies.
##'   \emph{Journal of the Academy of Marketing Science, 44}(1), 119--134.
##'   \doi{10.1007/s11747-015-0455-4}
##'
##' @examples
##'
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ x7 + x8 + x9 '
##'
##' dat <- HolzingerSwineford1939[, paste0("x", 1:9)]
##' htmt(HS.model, dat)
##'
##' ## save covariance matrix
##' HS.cov <- cov(HolzingerSwineford1939[, paste0("x", 1:9)])
##' ## HTMT using arithmetic mean
##' htmt(HS.model, sample.cov = HS.cov, htmt2 = FALSE)
##'
##' @export
htmt <- function(model, data = NULL, sample.cov = NULL, missing = "listwise",
                  ordered = NULL, absolute = TRUE, htmt2 = TRUE) {
  model <- lavaan::lavaanify(model)
  model <- model[model$op %in% "=~", ]
  factors <- unique(model$lhs)
  nf <- length(factors)
  var <- list()
  for (i in 1:nf) {
    var[[i]] <- model$rhs[which(model$lhs %in% factors[i])]
  }
  varnames <- c(unlist(var))
  if(!is.null(data)) { # if data
    if(any(! varnames %in% colnames(data))) {
      absent.vars <- which(! varnames %in% colnames(data))
      stop("Missing observed variables in the dataset: ",
           paste(varnames[absent.vars], collapse = " "))
    }
    data <- data[ , c(varnames)]
    R <- lavaan::lavCor(data, missing = missing, ordered = ordered)
    rownames(R) <- names(data)
    colnames(R) <- names(data)
  } else {
    if (any(! varnames %in% colnames(sample.cov))) {
      absent.vars <- which(! varnames %in% colnames(sample.cov))
      stop("Missing observed variables in the covariance or correlation matrix: ",
           paste(varnames[absent.vars], collapse = " "))
    }
    diagR <- diag(sample.cov)
    if (max(diagR) != 1 & min(diagR) != 1) { #if covariance matrix
      R <- cov2cor(sample.cov[varnames, varnames])
    } else { # if correlation matrix
      R <- sample.cov[varnames, varnames]
    }
  }
  if (absolute || htmt2) {
    R <- abs(R)
  }
  diag(R) <- NA
  m.cor.w <- list()
  for (i in 1:nf) {
    if (htmt2) {
      m.cor.w[[i]] <- exp(mean(log(R[ var[[i]], var[[i]] ]), na.rm = TRUE))
    } else m.cor.w[[i]] <- mean(R[ var[[i]], var[[i]] ], na.rm = TRUE)
  }
  m.cor.w <- as.numeric(m.cor.w)
  comb <- expand.grid(1:nf, 1:nf)
  g <- list()
  for (i in 1:nrow(comb)) {
    g[[i]] <- sqrt(m.cor.w[comb[i, 2]] * m.cor.w[comb[i, 1]])
  }
  g <- as.numeric(g)
  paste(comb[, 2], comb[, 1])
  m.cor.a <- list()
  for (i in 1:nrow(comb)) {
    if (htmt2) {
      m.cor.a[[i]] <- exp(mean(log(R[ var[[comb[i, 2]]],
                                      var[[comb[i, 1]]] ]),
                               na.rm = TRUE))
    } else m.cor.a[[i]] <- mean(R[ var[[ comb[i,2] ]],
                                   var[[ comb[i,1] ]] ], na.rm = TRUE)
  }
  m.cor.a <- as.numeric(m.cor.a)
  outhtmt <- m.cor.a / g
  res <- matrix(outhtmt, nrow = nf, ncol = nf, dimnames = list(factors))
  colnames(res) <- factors
  class(res) <- c("lavaan.matrix.symmetric", "matrix")
  res
}
