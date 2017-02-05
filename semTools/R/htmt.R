### Ylenio Longo
### Last updated:

htmt <- function (model, data = NULL, sample.cov = NULL, missing = "listwise",
                  ordered = NULL, absolute = TRUE) {
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
  if (absolute) {
    R <- abs(R)
  }
  diag(R) <- NA
  m.cor.w <- list()
  for (i in 1:nf) {
    m.cor.w[[i]] <- mean(R[var[[i]], var[[i]]], na.rm = TRUE)
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
    m.cor.a[[i]] <- mean(R[var[[comb[i, 2]]],
                           var[[comb[i, 1]]]], na.rm = TRUE)
  }
  m.cor.a <- as.numeric(m.cor.a)
  outhtmt <- m.cor.a / g
  res <- matrix(outhtmt, nrow = nf, ncol = nf, dimnames = list(factors))
  colnames(res) <- factors
  class(res) <- c("lavaan.matrix.symmetric", "matrix")
  res
}
