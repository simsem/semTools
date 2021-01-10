### Sunthud Pornprasertmanit and Alexander M. Schoemann
### Last updated: 10 January 2021
### prepare product indicators for 2-way and 3-way interactions in SEM


##' Make products of indicators using no centering, mean centering, double-mean
##' centering, or residual centering
##'
##' The \code{indProd} function will make products of indicators using no
##' centering, mean centering, double-mean centering, or residual centering. The
##' \code{orthogonalize} function is the shortcut of the \code{indProd} function
##' to make the residual-centered indicators products.
##'
##'
##' @aliases indProd orthogonalize
##' @importFrom stats lm
##'
##' @param data The desired data to be transformed.
##' @param var1 Names or indices of the variables loaded on the first factor
##' @param var2 Names or indices of the variables loaded on the second factor
##' @param var3 Names or indices of the variables loaded on the third factor
##' (for three-way interaction)
##' @param match Specify \code{TRUE} to use match-paired approach (Marsh, Wen, &
##' Hau, 2004). If \code{FALSE}, the resulting products are all possible
##' products.
##' @param meanC Specify \code{TRUE} for mean centering the main effect
##' indicator before making the products
##' @param residualC Specify \code{TRUE} for residual centering the products by
##' the main effect indicators (Little, Bovaird, & Widaman, 2006).
##' @param doubleMC Specify \code{TRUE} for centering the resulting products
##' (Lin et. al., 2010)
##' @param namesProd The names of resulting products
##' @return The original data attached with the products.
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com}) Alexander
##' Schoemann (East Carolina University; \email{schoemanna@@ecu.edu})
##' @seealso \itemize{ \item \code{\link{probe2WayMC}} For probing the two-way
##' latent interaction when the results are obtained from mean-centering, or
##' double-mean centering.  \item \code{\link{probe3WayMC}} For probing the
##' three-way latent interaction when the results are obtained from
##' mean-centering, or double-mean centering.  \item \code{\link{probe2WayRC}}
##' For probing the two-way latent interaction when the results are obtained
##' from residual-centering approach.  \item \code{\link{probe3WayRC}} For
##' probing the two-way latent interaction when the results are obtained from
##' residual-centering approach.  \item \code{\link{plotProbe}} Plot the simple
##' intercepts and slopes of the latent interaction.  }
##' @references Marsh, H. W., Wen, Z. & Hau, K. T. (2004). Structural equation
##' models of latent interactions: Evaluation of alternative estimation
##' strategies and indicator construction. \emph{Psychological Methods, 9}(3),
##' 275--300. \doi{10.1037/1082-989X.9.3.275}
##'
##' Lin, G. C., Wen, Z., Marsh, H. W., & Lin, H. S. (2010). Structural equation
##' models of latent interactions: Clarification of orthogonalizing and
##' double-mean-centering strategies. \emph{Structural Equation Modeling, 17}(3),
##' 374--391. \doi{10.1080/10705511.2010.488999}
##'
##' Little, T. D., Bovaird, J. A., & Widaman, K. F. (2006). On the merits of
##' orthogonalizing powered and product terms: Implications for modeling
##' interactions among latent variables. \emph{Structural Equation Modeling,
##' 13}(4), 497--519. \doi{10.1207/s15328007sem1304_1}
##' @examples
##'
##' ## Mean centering / two-way interaction / match-paired
##' dat <- indProd(attitude[ , -1], var1 = 1:3, var2 = 4:6)
##'
##' ## Residual centering / two-way interaction / match-paired
##' dat2 <- indProd(attitude[ , -1], var1 = 1:3, var2 = 4:6, match = FALSE,
##'                 meanC = FALSE, residualC = TRUE, doubleMC = FALSE)
##'
##' ## Double-mean centering / two-way interaction / match-paired
##' dat3 <- indProd(attitude[ , -1], var1 = 1:3, var2 = 4:6, match = FALSE,
##'                 meanC = TRUE, residualC = FALSE, doubleMC = TRUE)
##'
##' ## Mean centering / three-way interaction / match-paired
##' dat4 <- indProd(attitude[ , -1], var1 = 1:2, var2 = 3:4, var3 = 5:6)
##'
##' ## Residual centering / three-way interaction / match-paired
##' dat5 <- orthogonalize(attitude[ , -1], var1 = 1:2, var2 = 3:4, var3 = 5:6,
##'                       match = FALSE)
##'
##' ## Double-mean centering / three-way interaction / match-paired
##' dat6 <- indProd(attitude[ , -1], var1 = 1:2, var2 = 3:4, var3 = 5:6,
##'                 match = FALSE, meanC = TRUE, residualC = TRUE,
##'                 doubleMC = TRUE)
##'
##'
##' ## To add product-indicators to multiple-imputed data sets
##' \dontrun{
##' HSMiss <- HolzingerSwineford1939[ , c(paste0("x", 1:9), "ageyr","agemo")]
##' set.seed(12345)
##' HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
##' age <- HSMiss$ageyr + HSMiss$agemo/12
##' HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
##' library(Amelia)
##' set.seed(12345)
##' HS.amelia <- amelia(HSMiss, m = 3, p2s = FALSE)
##' imps <- HS.amelia$imputations # extract a list of imputations
##' ## apply indProd() to the list of data.frames
##' imps2 <- lapply(imps, indProd,
##'                 var1 = c("x1","x2","x3"), var2 = c("x4","x5","x6"))
##' ## verify:
##' lapply(imps2, head)
##' }
##'
##' @export
indProd <- function(data, var1, var2, var3 = NULL, match = TRUE, meanC = TRUE,
                    residualC = FALSE, doubleMC = TRUE, namesProd = NULL) {
  # Get all variable names
  if (all(is.numeric(var1)))
    var1 <- colnames(data)[var1]
  if (all(is.numeric(var2)))
    var2 <- colnames(data)[var2]
  if (!is.null(var3) && all(is.numeric(var3))) var3 <- colnames(data)[var3]
  dat1 <- data[, var1]
  dat2 <- data[, var2]
  dat3 <- NULL
  if (!is.null(var3)) dat3 <- data[, var3]

  # Mean centering on the original indicators
  if (meanC) {
    dat1 <- scale(dat1, scale = FALSE)
    dat2 <- scale(dat2, scale = FALSE)
    if (!is.null(dat3)) dat3 <- scale(dat3, scale = FALSE)
  }
  if (match) {
    # Check whether the number of variables are equal across variable sets
    if (length(var1) != length(var2))
      stop("If the match-paired approach is used, the number of",
           " variables in all sets must be equal.")
    if (!is.null(var3) && (length(var1) != length(var3)))
      stop("If the match-paired approach is used, the number of",
           " variables in all three sets must be equal.")
    datProd <- NULL
    if (is.null(var3)) {
      # Two-way interaction
      datProd <- dat1 * dat2
      if (residualC) {
        notmissing <- which(!apply(datProd, 1, function(x) any(is.na(x))))
        colnames(datProd) <- paste("interactionProduct", 1:ncol(datProd), sep = "")
        # Write the expression for linear model and residualize the products
        temp <- data.frame(datProd, dat1, dat2)
        express <- paste("cbind(", paste(colnames(datProd), collapse = ", "),
                         ") ~ ", paste(c(colnames(dat1), colnames(dat2)),
                                       collapse = " + "), sep = "")
        datProd[notmissing,] <- lm(express, data = temp)$residuals
      }
    } else {
      # Three-way interaction
      datProd2way <- cbind(dat1 * dat2, dat1 * dat3, dat2 * dat3)
      datProd3way <- dat1 * dat2 * dat3
      if (residualC) {
        notmissing2way <- which(!apply(datProd2way, 1, function(x) any(is.na(x))))
        colnames(datProd2way) <- paste("interaction2Product", 1:ncol(datProd2way), sep = "")

        # Write the expression for linear model and residualize the two-way products
        temp2 <- data.frame(datProd2way, dat1, dat2, dat3)
        express2 <- paste("cbind(", paste(colnames(datProd2way), collapse = ", "),
                          ") ~ ", paste(c(colnames(dat1), colnames(dat2),
                                          colnames(dat3)), collapse = " + "), sep = "")
        datProd2way[notmissing2way,] <- lm(express2, data = temp2)$residuals

        # Making all possible products to residualize the 3-way interaction
        datProd2wayFull <- matrix(0, nrow(data), 1)
        for (i in 1:length(var1)) datProd2wayFull <- data.frame(datProd2wayFull, matrix(rep(dat1[, i], length(var2)), ncol = length(var2)) * dat2)
        for (i in 1:length(var1)) datProd2wayFull <- data.frame(datProd2wayFull, matrix(rep(dat1[, i], length(var3)), ncol = length(var3)) * dat3)
        for (i in 1:length(var2)) datProd2wayFull <- data.frame(datProd2wayFull, matrix(rep(dat2[, i], length(var3)), ncol = length(var3)) * dat3)
        datProd2wayFull <- datProd2wayFull[, -1]
        colnames(datProd2wayFull) <- paste("interaction2Product", 1:ncol(datProd2wayFull), sep = "")

        notmissing3way <- which(!apply(datProd3way, 1, function(x) any(is.na(x))))
        colnames(datProd3way) <- paste("interaction3Product", 1:ncol(datProd3way), sep = "")
        # Write the expression for linear model and residualize the three-way products
        temp3 <- data.frame(datProd3way, dat1, dat2, dat3, datProd2wayFull)
        express3 <- paste("cbind(", paste(colnames(datProd3way), collapse = ", "),
                          ") ~ ", paste(c(colnames(dat1), colnames(dat2), colnames(dat3),
                                          colnames(datProd2wayFull)), collapse = " + "), sep = "")
        datProd3way[notmissing3way,] <- lm(express3, data = temp3)$residuals
      }
      datProd <- cbind(datProd2way, datProd3way)
    }
    ## Mean-centering the final product
    if (doubleMC) datProd <- scale(datProd, scale = FALSE)

    ## Rename the obtained product terms
    if (is.null(namesProd)) {
      if (is.null(var3)) {
        colnames(datProd) <- paste(var1, var2, sep = ".")
      } else {
        colnames(datProd) <- c(paste(var1, var2, sep = "."),
                               paste(var1, var3, sep = "."),
                               paste(var2, var3, sep = "."),
                               paste(var1, var2, var3, sep = "."))
      }
    } else {
      colnames(datProd) <- namesProd
    }
  } else {
    datProd <- NULL
    if (is.null(var3)) {
      # Create all possible combinations of the products of indicators
      datProd <- matrix(0, nrow(data), 1)
      for (i in 1:length(var1)) datProd <- data.frame(datProd, matrix(rep(dat1[, i], length(var2)), ncol = length(var2)) * dat2)
      datProd <- datProd[, -1]
      if (residualC) {
        notmissing <- which(!apply(datProd, 1, function(x) any(is.na(x))))
        colnames(datProd) <- paste("interactionProduct", 1:ncol(datProd), sep = "")
        # Write the expression for linear model and residualize the two-way products
        temp <- data.frame(datProd, dat1, dat2)
        express <- paste("cbind(", paste(colnames(datProd), collapse = ", "),
                         ") ~ ", paste(c(colnames(dat1), colnames(dat2)),
                                       collapse = " + "), sep = "")
        datProd[notmissing,] <- lm(express, data = temp)$residuals
      }
    } else {
      # Create all possible combinations of the products of indicators
      datProd2way <- matrix(0, nrow(data), 1)
      for (i in 1:length(var1)) datProd2way <- data.frame(datProd2way, matrix(rep(dat1[, i], length(var2)), ncol = length(var2)) * dat2)
      for (i in 1:length(var1)) datProd2way <- data.frame(datProd2way, matrix(rep(dat1[, i], length(var3)), ncol = length(var3)) * dat3)
      for (i in 1:length(var2)) datProd2way <- data.frame(datProd2way, matrix(rep(dat2[, i], length(var3)), ncol = length(var3)) * dat3)
      datProd3way <- matrix(0, nrow(data), 1)
      for (i in 1:length(var1)) {
        for(j in 1:length(var2)) {
          datProd3way <- data.frame(datProd3way, matrix(rep(dat1[, i], length(var3)), ncol = length(var3)) * matrix(rep(dat2[, j], length(var3)), ncol = length(var3)) * dat3)
        }
      }
      datProd2way <- datProd2way[, -1]
      datProd3way <- datProd3way[, -1]
      if (residualC) {
        notmissing2way <- which(!apply(datProd2way, 1, function(x) any(is.na(x))))
        colnames(datProd2way) <- paste("interaction2Product", 1:ncol(datProd2way), sep = "")
        # Write the expression for linear model and residualize the two-way products
        temp2 <- data.frame(datProd2way, dat1, dat2, dat3)
        express2 <- paste("cbind(", paste(colnames(datProd2way), collapse = ", "),
                          ") ~ ", paste(c(colnames(dat1), colnames(dat2),
                                          colnames(dat3)), collapse = " + "), sep = "")
        datProd2way[notmissing2way,] <- lm(express2, data = temp2)$residuals
        notmissing3way <- which(!apply(datProd3way, 1, function(x) any(is.na(x))))
        colnames(datProd3way) <- paste("interaction3Product", 1:ncol(datProd3way), sep = "")
        # Write the expression for linear model and residualize the three-way products
        temp3 <- data.frame(datProd3way, dat1, dat2, dat3, datProd2way)
        express3 <- paste("cbind(", paste(colnames(datProd3way), collapse = ", "),
                          ") ~ ", paste(c(colnames(dat1), colnames(dat2),
                                          colnames(dat3), colnames(datProd2way)),
                                        collapse = " + "), sep = "")
        datProd3way[notmissing3way,] <- lm(express3, data = temp3)$residuals
      }
      datProd <- cbind(datProd2way, datProd3way)
    }
    ## Double-mean centering
    if (doubleMC) datProd <- scale(datProd, scale = FALSE)

    ## Name the resulting product terms
    if (is.null(namesProd)) {
      temp <- NULL
      if (is.null(var3)) {
        for (i in 1:length(var1)) temp <- c(temp, paste(var1[i], var2, sep = "."))
      } else {
        for (i in 1:length(var1)) temp <- c(temp, paste(var1[i], var2, sep = "."))
        for (i in 1:length(var1)) temp <- c(temp, paste(var1[i], var3, sep = "."))
        for (i in 1:length(var2)) temp <- c(temp, paste(var2[i], var3, sep = "."))
        for (i in 1:length(var1)) {
          for(j in 1:length(var2)) {
            temp <- c(temp, paste(var1[i], var2[j], var3, sep = "."))
          }
        }
      }
      colnames(datProd) <- temp
    } else {
      colnames(datProd) <- namesProd
    }
  }
  ## Bind the products back to the original data
  data.frame(data, datProd)
}

##' @rdname indProd
##' @export
orthogonalize <- function(data, var1, var2, var3 = NULL,
                          match = TRUE, namesProd = NULL) {
  indProd(data = data, var1 = var1, var2 = var2, var3 = var3,
          match = match, meanC = FALSE, residualC = TRUE, doubleMC = FALSE,
          namesProd = namesProd)
}


