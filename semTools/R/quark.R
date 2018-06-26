### Steven R. Chesnut, Danny Squire, Terrence D. Jorgensen
### Last updated: 26 June 2018



#' Quark
#'
#' The \code{quark} function provides researchers with the ability to calculate
#' and include component scores calculated by taking into account the variance
#' in the original dataset and all of the interaction and polynomial effects of
#' the data in the dataset.
#'
#' The \code{quark} function calculates these component scores by first filling
#' in the data via means of multiple imputation methods and then expanding the
#' dataset by aggregating the non-overlapping interaction effects between
#' variables by calculating the mean of the interactions and polynomial
#' effects.  The multiple imputation methods include one of iterative sampling
#' and group mean substitution and multiple imputation using a polytomous
#' regression algorithm (mice). During the expansion process, the dataset is
#' expanded to three times its normal size (in width). The first third of the
#' dataset contains all of the original data post imputation, the second third
#' contains the means of the polynomial effects (squares and cubes), and the
#' final third contains the means of the non-overlapping interaction effects. A
#' full principal componenent analysis is conducted and the individual
#' components are retained. The subsequent \code{\link{combinequark}} function
#' provides researchers the control in determining how many components to
#' extract and retain. The function returns the dataset as submitted (with
#' missing values) and the component scores as requested for a more accurate
#' multiple imputation in subsequent steps.
#'
#' @param data The data frame is a required component for \code{quark}.  In
#' order for \code{quark} to process a data frame, it must not contain any
#' factors or text-based variables.  All variables must be in numeric format.
#' Identifiers and dates can be left in the data; however, they will need to be
#' identified under the \code{id} argument.
#' @param id Identifiers and dates within the dataset will need to be
#' acknowledged as \code{quark} cannot process these.  By acknowledging the
#' identifiers and dates as a vector of column numbers or variable names,
#' \code{quark} will remove them from the data temporarily to complete its main
#' processes.  Among many potential issues of not acknowledging identifiers and
#' dates are issues involved with imputation, product and polynomial effects,
#' and principal component analysis.
#' @param order Order is an optional argument provided by quark that can be
#' used when the imputation procedures in mice fail.  Under some circumstances,
#' mice cannot calculate missing values due to issues with extreme missingness.
#' Should an error present itself stating a failure due to not having any
#' columns selected, set the argument \code{order = 2} in order to reorder the
#' imputation method procedure.  Otherwise, use the default \code{order = 1}.
#' @param silent If \code{FALSE}, the details of the \code{quark} process are
#' printed.
#' @param \dots additional arguments to pass to \code{\link[mice]{mice}}.
#' @return The output value from using the quark function is a list. It will
#' return a list with 7 components.
#'  \item{ID Columns}{Is a vector of the identifier columns entered when
#'   running quark.}
#'  \item{ID Variables}{Is a subset of the dataset that contains the identifiers
#'   as acknowledged when running quark.}
#'  \item{Used Data}{Is a matrix / dataframe of the data provided by user as
#'   the basis for quark to process.}
#'  \item{Imputed Data}{Is a matrix / dataframe of the data after the multiple
#'   method imputation process.}
#'  \item{Big Matrix}{Is the expanded product and polynomial matrix.}
#'  \item{Principal Components}{Is the entire dataframe of principal components
#'   for the dataset.  This dataset will have the same number of rows of the big
#'   matrix, but will have 1 less column (as is the case with principal
#'   component analyses).}
#'  \item{Percent Variance Explained}{Is a vector of the percent variance
#'   explained with each column of principal components.}
#' @author Steven R. Chesnut (University of Southern Mississippi;
#' \email{Steven.Chesnut@@usm.edu})
#'
#' Danny Squire (Texas Tech University)
#'
#' Terrence D. Jorgensen (University of Amsterdam)
#'
#' The PCA code is copied and modified from the \code{FactoMineR} package.
#' @seealso \code{\link{combinequark}}
#' @references Howard, W. J., Rhemtulla, M., & Little, T. D. (2015). Using
#' Principal Components as Auxiliary Variables in Missing Data Estimation.
#' \emph{Multivariate Behavioral Research, 50}(3), 285--299.
#' doi:10.1080/00273171.2014.999267
#' @examples
#'
#' set.seed(123321)
#'
#' dat <- HolzingerSwineford1939[,7:15]
#' misspat <- matrix(runif(nrow(dat) * 9) < 0.3, nrow(dat))
#' dat[misspat] <- NA
#' dat <- cbind(HolzingerSwineford1939[,1:3], dat)
#' \dontrun{
#' quark.list <- quark(data = dat, id = c(1, 2))
#'
#' final.data <- combinequark(quark = quark.list, percent = 80)
#'
#' ## Example to rerun quark after imputation failure:
#' quark.list <- quark(data = dat, id = c(1, 2), order = 2)
#' }
#'
#' @export
quark <- function(data, id, order = 1, silent = FALSE, ...){
  if(!is.data.frame(data) && !is.matrix(data)) {
    stop("Inappropriate data file provided.")
  }
  if(!silent) cat("Data Check Passed.\n")

  if(is.character(id)) id <- match(id, colnames(data))
  for(i in 1:length(id)){
    if(id[i] > ncol(data) || id[i] < 1){
      stop("At least one of the IDs is out of bounds.")
    }
  }
  if(!silent) cat("ID Check Passed.\n")
  if(!(order %in% 1:2)) stop("Currently, the order argument can take either 1 or 2.")

  final.collect <- list()
  final.collect$ID_Columns <- id
  final.collect$ID_Vars <- data[,id]
  final.collect$Used_Data <- data[,-c(id)]
  ##FIXME 26-June-2018: Terrence had to add a logical check for whether mice
  ##                    is installed, otherwise won't pass CRAN checks.
  checkMice <- requireNamespace("mice")
  if (!checkMice) {
    message('The quark function requires the "mice" package to be installed.')
    return(invisible(NULL))
  }
  final.collect$Imputed_Data <- imputequark(data = final.collect$Used_Data,
                                            order = order, silent = silent, ...)
  final.collect$Big_Data_Matrix <- bigquark(data = final.collect$Imputed_Data,
                                            silent = silent)
  cmp <- compquark(data = final.collect$Big_Data_Matrix, silent = silent)
  final.collect$Prin_Components <- cmp[[1]]
  final.collect$Prin_Components_Prcnt <- cmp[[2]]

  return(final.collect)
}



#' Combine the results from the quark function
#'
#' This function builds upon the \code{\link{quark}} function to provide a
#' final dataset comprised of the original dataset provided to
#' \code{\link{quark}} and enough principal components to be able to account
#' for a certain level of variance in the data.
#'
#'
#' @param quark Provide the \code{\link{quark}} object that was returned.  It
#' should be a list of objects.  Make sure to include it in its entirety.
#' @param percent Provide a percentage of variance that you would like to have
#' explained.  That many components (columns) will be extracted and kept with
#' the output dataset.  Enter this variable as a number WITHOUT a percentage
#' sign.
#' @return The output of this function is the original dataset used in quark
#' combined with enough principal component scores to be able to account for
#' the amount of variance that was requested.
#' @author Steven R. Chesnut (University of Southern Mississippi
#' \email{Steven.Chesnut@@usm.edu})
#' @seealso \code{\link{quark}}
#' @examples
#'
#' set.seed(123321)
#' dat <- HolzingerSwineford1939[,7:15]
#' misspat <- matrix(runif(nrow(dat) * 9) < 0.3, nrow(dat))
#' dat[misspat] <- NA
#' dat <- cbind(HolzingerSwineford1939[,1:3], dat)
#'
#' quark.list <- quark(data = dat, id = c(1, 2))
#'
#' final.data <- combinequark(quark = quark.list, percent = 80)
#'
#' @export
combinequark <- function(quark, percent) {
  data <- cbind(quark$ID_Vars, quark$Used_Data)
  pct <- quark$Prin_Components_Prcnt
  comp <- quark$Prin_Components

  for (i in 1:length(pct)) {
    if(pct[i] >= percent) {
      num <- i
      break
    }
  }
  return(cbind(data, comp[,1:num]))
}



## ----------------
## Hidden Functions
## ----------------

imputequark <- function(data, order, silent = FALSE, ...){
  if (order == 1){
    data <- aImp(data = data, silent = silent, ...)
    data <- gImp(data = data, silent = silent)
  } else if(order == 2) {
    data <- gImp(data = data, silent = silent)
    if (length(which(is.na(data > 0)))) {
      data <- aImp(data = data, silent = silent, ...)
    }
  }
  return(data)
}

#' @importFrom stats cor
gImp <- function(data, silent = FALSE) {
  imputed_data <- data
  num_adds <- vector(length = ncol(data)) # number of columns combined into one for averaging.
  data.cor <- cor(data, use = "pairwise", method = "pearson")
  class(data.cor) <- c("lavaan.matrix.symmetric","matrix")
  if (!silent) print(data.cor)
  #populate multiple matrices that can then be utilized to determine if one column should enhance another based upon
  #the correlations they share...
  if (!silent) cat("Imputing Column... \n")

  for (a in 1:ncol(data)) {
    temp_mat <- matrix(ncol = ncol(data), nrow = nrow(data))
    list <- unique(sort(data[,a]))
    if (length(list) > 1 && length(list) <= 10) {
      for (b in 1:nrow(data)) {
        for (c in 1:length(list)) {
          if (data[b, a] == list[c] && !is.na(data[b,a])) {
            temp_mat[b,] <- round(colMeans(subset(data, data[ , a] == list[c]), na.rm = TRUE), digits = 1)
          } else if (is.na(data[b,a])) {
            for (p in 1:ncol(data)) temp_mat[b,p] <- data[b,p]
          }
        }
      }

      # Here I need to determine if the other columns are correlated enough with
      # the reference to ensure accuracy of predictions
      temp_cor <- data.cor[,a]
      # if (countNA(temp_cor)==0) {
      for (i in 1:length(temp_cor)) {
        if (i != a) {
          if (abs(temp_cor[i]) >= .5 && !is.na(temp_cor[i])) { # Using a moderate effect size, column a, will inform other columns.
            for (x in 1:nrow(imputed_data)){
              imputed_data[x,i] <- sum(imputed_data[x,i], temp_mat[x,a], na.rm = TRUE)
            }
            num_adds[i] <- num_adds[i] + 1
          }
        }
      }
      #}
      if (!silent) cat("\t", colnames(data)[a])
    }
  }
  if (!silent) cat("\n")
  imputed_data <- cleanMat(m1 = data, m2 = imputed_data, impact = num_adds)
  imputed_data <- fixData(imputed_data)
  return(imputed_data)
}

cleanMat <- function(m1, m2, impact) {
  #Impact is the number of influences on each column...
  #We need to clean up and then try to determine what final values should be...
  #Go through each of the cells...
  new_mat <- m2
  for (a in 1:ncol(m1)) {
    for (b in 1:nrow(m1)) {
      if (!is.na(m1[b,a])) {
        new_mat[b,a] <- m1[b,a]
      } else if (is.na(m1[b,a])) {
        new_mat[b,a] <- new_mat[b,a] / impact[a]
      }
    }
  }
  return(new_mat)
}

fixData <- function(data) {
  for (a in 1:ncol(data)) {
    for (b in 1:nrow(data)) {
      data[b,a] <- round(data[b,a], digits = 1)
    }
  }

  return(data)
}

aImp <- function(data, silent = FALSE, ...) {
  miceArgs <- list(...)
  miceArgs$data <- data
  miceArgs$maxit <- 1
  miceArgs$m <- 1
  miceArgs$printFlag <- !silent
  requireNamespace("mice")
  if (!("package:mice" %in% search())) attachNamespace("mice")
  if (!silent) cat("Starting Algorithm Imputation...\n")
  impData <- mice::complete(do.call("mice", miceArgs))
  if (!silent) cat("Ending Algorithm Imputation...\n")
  return(impData)
}

bigquark <- function(data, silent = FALSE) {
  if (!silent) cat("Calculating Polynomial Effects.\n")
  poly <- ((data^2)+(data^3))/2
  if (!silent) cat("Creating Matrix for Interaction Effects.\n")
  prod <- matrix(ncol=(ncol(data)-1),nrow=nrow(data))
  if (!silent) cat("Calculating Interaction Effects...0%..")
  for (i in 1:nrow(data)) {
    if (!silent) printpct(percent = i/nrow(data))
    for (j in 1:(ncol(data)-1)) {
      prod[i,j] <- mean(as.numeric(data[i,j])*as.numeric(data[i,(j+1):ncol(data)]))
    }
  }
  cat("\n")
  data <- cbind(data,poly,prod)
  return(data)
}

compquark <- function(data, silent = FALSE) {
  if (!silent) cat("Calculating values for the PCA\n")
  pcam <- pcaquark(data, ncp = ncol(data))
  cmp <- list()
  cmp$pca <- pcam$ind$coord
  cmp$var <- pcam$eig[,3]
  colnames(cmp$pca) <- c(paste0("AuxVar",1:ncol(cmp$pca)))
  return(cmp)
}

printpct <- function(percent) {
  if (round(percent, digits = 10) == 0) cat("0%..")
  if (round(percent, digits = 10) == .10) cat("10%..")
  if (round(percent, digits = 10) == .20) cat("20%..")
  if (round(percent, digits = 10) == .30) cat("30%..")
  if (round(percent, digits = 10) == .40) cat("40%..")
  if (round(percent, digits = 10) == .50) cat("50%..")
  if (round(percent, digits = 10) == .60) cat("60%..")
  if (round(percent, digits = 10) == .70) cat("70%..")
  if (round(percent, digits = 10) == .80) cat("80%..")
  if (round(percent, digits = 10) == .90) cat("90%..")
  if (round(percent, digits = 10) == 1) cat("100%..")
}

# This function is modified from the FactoMinoR package.
pcaquark <- function(X, ncp = 5) {
  moy.p <- function(V, poids) res <- sum(V * poids)/sum(poids)
  ec <- function(V, poids) res <- sqrt(sum(V^2 * poids)/sum(poids))
  X <- as.data.frame(X)
  if (any(is.na(X))) {
      warnings("Missing values are imputed by the mean of the variable: you should use the imputePCA function of the missMDA package")
      X[is.na(X)] <- matrix(apply(X,2,mean,na.rm=TRUE),ncol=ncol(X),nrow=nrow(X),byrow=TRUE)[is.na(X)]
  }
  if (is.null(rownames(X))) rownames(X) <- 1:nrow(X)
  if (is.null(colnames(X))) colnames(X) <- paste("V", 1:ncol(X), sep = "")
  colnames(X)[colnames(X) == ""] <- paste("V", 1:sum(colnames(X)==""),sep="")
  rownames(X)[is.null(rownames(X))] <- paste("row",1:sum(rownames(X)==""),sep="")
  Xtot <- X
  if (any(!sapply(X, is.numeric))) {
      auxi <- NULL
      for (j in 1:ncol(X)) if (!is.numeric(X[, j])) auxi <- c(auxi, colnames(X)[j])
      stop(paste("\nThe following variables are not quantitative: ", auxi))
  }
  ncp <- min(ncp, nrow(X) - 1, ncol(X))
  row.w <- rep(1, nrow(X))
  row.w.init <- row.w
  row.w <- row.w/sum(row.w)
  col.w <- rep(1, ncol(X))
  centre <- apply(X, 2, moy.p, row.w)
  X <- as.matrix(sweep(as.matrix(X), 2, centre, FUN = "-"))
	ecart.type <- apply(X, 2, ec, row.w)
	ecart.type[ecart.type <= 1e-16] <- 1
	X <- sweep(as.matrix(X), 2, ecart.type, FUN = "/")
	dist2.ind <- apply(sweep(X,2,sqrt(col.w),FUN="*")^2,1,sum)
	dist2.var <- apply(sweep(X,1,sqrt(row.w),FUN="*")^2,2,sum)
  tmp <- svd.triplet.quark(X, row.w = row.w, col.w = col.w, ncp = ncp)
  eig <- tmp$vs^2
  vp <- as.data.frame(matrix(NA, length(eig), 3))
  rownames(vp) <- paste("comp", 1:length(eig))
  colnames(vp) <- c("eigenvalue","percentage of variance",
                    "cumulative percentage of variance")
  vp[, "eigenvalue"] <- eig
  vp[, "percentage of variance"] <- (eig/sum(eig)) * 100
  vp[, "cumulative percentage of variance"] <- cumsum(vp[, "percentage of variance"])
  V <- tmp$V
  U <- tmp$U
	eig <- eig[1:ncp]
  coord.ind <- sweep(as.matrix(U), 2, sqrt(eig), FUN = "*")
  coord.var <- sweep(as.matrix(V), 2, sqrt(eig), FUN = "*")
  contrib.var <- sweep(as.matrix(coord.var^2), 2, eig, "/")
  contrib.var <- sweep(as.matrix(contrib.var), 1, col.w, "*")
	dist2 <- dist2.var
  cor.var <- sweep(as.matrix(coord.var), 1, sqrt(dist2), FUN = "/")
  cos2.var <- cor.var^2
  rownames(coord.var) <- rownames(cos2.var) <- rownames(cor.var) <- rownames(contrib.var) <- colnames(X)
  colnames(coord.var) <- colnames(cos2.var) <- colnames(cor.var) <- colnames(contrib.var) <- paste("Dim", c(1:ncol(V)), sep = ".")
  res.var <- list(coord = coord.var[, 1:ncp], cor = cor.var[, 1:ncp],
                  cos2 = cos2.var[, 1:ncp], contrib = contrib.var[, 1:ncp] * 100)
	dist2 <- dist2.ind
  cos2.ind <- sweep(as.matrix(coord.ind^2), 1, dist2, FUN = "/")
  contrib.ind <- sweep(as.matrix(coord.ind^2), 1, row.w/sum(row.w), FUN = "*")
  contrib.ind <- sweep(as.matrix(contrib.ind), 2, eig, FUN = "/")
  rownames(coord.ind) <- rownames(cos2.ind) <- rownames(contrib.ind) <- names(dist2) <- rownames(X)
  colnames(coord.ind) <- colnames(cos2.ind) <- colnames(contrib.ind) <- paste("Dim", c(1:ncol(U)), sep = ".")
  res.ind <- list(coord = coord.ind[, 1:ncp], cos2 = cos2.ind[, 1:ncp],
                  contrib = contrib.ind[, 1:ncp] * 100, dist = sqrt(dist2))
  res <- list(eig = vp, var = res.var, ind = res.ind, svd = tmp)
  class(res) <- c("PCA", "list")
  return(res)
}

# This function is modified from the FactoMinoR package.
svd.triplet.quark <- function (X, row.w = NULL, col.w = NULL, ncp = Inf) {
	tryCatch.W.E <- function(expr) {  ## function proposed by Maechlmr
		W <- NULL
		w.handler <- function(w) { # warning handler
			W <<- w
			invokeRestart("muffleWarning")
		}
		list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
		                                 warning = w.handler), warning = W)
	}
  ncp <- min(ncp,nrow(X)-1,ncol(X))
  row.w <- row.w / sum(row.w)
  X <- sweep(X, 2, sqrt(col.w), FUN = "*")
  X <- sweep(X, 1, sqrt(row.w), FUN = "*")
	if (ncol(X) < nrow(X)) {
		svd.usuelle <- tryCatch.W.E(svd(X, nu = ncp, nv = ncp))$val
		if (names(svd.usuelle)[[1]] == "message") {
		  svd.usuelle <- tryCatch.W.E(svd(t(X), nu = ncp, nv = ncp))$val
		  if (names(svd.usuelle)[[1]] == "d") {
  			aux <- svd.usuelle$u
  			svd.usuelle$u <- svd.usuelle$v
  			svd.usuelle$v <- aux
		  } else {
			  bb <- eigen(t(X) %*% X, symmetric = TRUE)
			  svd.usuelle <- vector(mode = "list", length = 3)
			  svd.usuelle$d[svd.usuelle$d < 0] <- 0
			  svd.usuelle$d <- sqrt(svd.usuelle$d)
			  svd.usuelle$v <- bb$vec[,1:ncp]
			  svd.usuelle$u <- sweep(X %*% svd.usuelle$v, 2, svd.usuelle$d[1:ncp], FUN = "/")
		  }
		}
		U <- svd.usuelle$u
		V <- svd.usuelle$v
		if (ncp > 1) {
		  mult <- sign(apply(V, 2, sum))
		  mult[mult == 0] <- 1
		  U <- sweep(U, 2, mult, FUN = "*")
		  V <- sweep(V, 2, mult, FUN = "*")
		}
		U <- sweep(as.matrix(U), 1, sqrt(row.w), FUN = "/")
		V <- sweep(as.matrix(V), 1, sqrt(col.w), FUN = "/")
	} else {
		svd.usuelle <- tryCatch.W.E(svd(t(X), nu = ncp, nv = ncp))$val
		if (names(svd.usuelle)[[1]] == "message") {
		  svd.usuelle <- tryCatch.W.E(svd(X, nu = ncp, nv = ncp))$val
		  if (names(svd.usuelle)[[1]] == "d") {
  			aux <- svd.usuelle$u
  			svd.usuelle$u <- svd.usuelle$v
  			svd.usuelle$v <- aux
		  } else {
			  bb <- eigen(X%*%t(X),symmetric=TRUE)
			  svd.usuelle <- vector(mode = "list", length = 3)
			  svd.usuelle$d[svd.usuelle$d < 0] <- 0
			  svd.usuelle$d <- sqrt(svd.usuelle$d)
			  svd.usuelle$v <- bb$vec[,1:ncp]
			  svd.usuelle$u <- sweep(t(X) %*% svd.usuelle$v, 2, svd.usuelle$d[1:ncp], FUN = "/")
		  }
		}
		U <- svd.usuelle$v
		V <- svd.usuelle$u
		mult <- sign(apply(V, 2, sum))
		mult[mult == 0] <- 1
		V <- sweep(V, 2, mult, FUN = "*")
		U <- sweep(U, 2, mult, FUN = "*")
		U <- sweep(U, 1, sqrt(row.w), FUN = "/")
		V <- sweep(V, 1, sqrt(col.w), FUN = "/")
	}
  vs <- svd.usuelle$d[1:min(ncol(X), nrow(X) - 1)]
	num <- which(vs[1:ncp] < 1e-15)
  if (length(num)==1) {
	  U[,num] <- U[,num] * vs[num]
    V[,num] <- V[,num] * vs[num]
	}
  if (length(num) > 1) {
	  U[,num] <- sweep(U[,num], 2, vs[num], FUN = "*")
    V[,num] <- sweep(V[,num], 2, vs[num], FUN = "*")
	}
  res <- list(vs = vs, U = U, V = V)
  return(res)
}


