### Sunthud Pornprasertmanit
### Last updated: 27 June 2018
### deprecated because it is obsolete (now available in OpenMx)


#' Find standardized estimates for OpenMx output
#'
#' Find standardized estimates for OpenMx output. This function is applicable
#' for the \code{MxRAMObjective} only.
#'
#'
#' @param object Target OpenMx output using \code{MxRAMObjective}
#' @param free If \code{TRUE}, the function will show only standardized values
#' of free parameters. If \code{FALSE}, the function will show the results for
#' fixed and free parameters.
#'
#' @return A vector of standardized estimates
#'
#' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#'
#' @examples
#'
#' \dontrun{
#' library(OpenMx)
#' data(myFADataRaw)
#' myFADataRaw <- myFADataRaw[,c("x1","x2","x3","x4","x5","x6")]
#' oneFactorModel <- mxModel("Common Factor Model Path Specification",
#' 	type="RAM",
#' 	mxData(
#' 		observed=myFADataRaw,
#' 		type="raw"
#' 	),
#' 	manifestVars=c("x1","x2","x3","x4","x5","x6"),
#' 	latentVars="F1",
#' 	mxPath(from=c("x1","x2","x3","x4","x5","x6"),
#' 		arrows=2,
#' 		free=TRUE,
#' 		values=c(1,1,1,1,1,1),
#' 		labels=c("e1","e2","e3","e4","e5","e6")
#' 	),
#' 	# residual variances
#' 	# -------------------------------------
#' 	mxPath(from="F1",
#' 		arrows=2,
#' 		free=TRUE,
#' 		values=1,
#' 		labels ="varF1"
#' 	),
#' 	# latent variance
#' 	# -------------------------------------
#' 	mxPath(from="F1",
#' 		to=c("x1","x2","x3","x4","x5","x6"),
#' 		arrows=1,
#' 		free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
#' 		values=c(1,1,1,1,1,1),
#' 		labels =c("l1","l2","l3","l4","l5","l6")
#' 	),
#' 	# factor loadings
#' 	# -------------------------------------
#' 	mxPath(from="one",
#' 		to=c("x1","x2","x3","x4","x5","x6","F1"),
#' 		arrows=1,
#' 		free=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE),
#' 		values=c(1,1,1,1,1,1,0),
#' 		labels =c("meanx1","meanx2","meanx3","meanx4","meanx5","meanx6",NA)
#' 	)
#' 	# means
#' 	# -------------------------------------
#' ) # close model
#' # Create an MxModel object
#' # -----------------------------------------------------------------------------
#' oneFactorFit <- mxRun(oneFactorModel)
#' standardizeMx(oneFactorFit)
#'
#' # Compare with lavaan
#' library(lavaan)
#' script <- "f1 =~ x1 + x2 + x3 + x4 + x5 + x6"
#' fit <- cfa(script, data=myFADataRaw, meanstructure=TRUE)
#' standardizedSolution(fit)
#' }
#'
#' @name standardizeMx-deprecated
#' @usage standardizeMx(object, free = TRUE)
#' @seealso \code{\link{semTools-deprecated}}
#' @keywords internal
NULL


#' @rdname semTools-deprecated
#' @section \code{standardizeMx}:
#' The \code{standardizeMx} and \code{fitMeasuresMx} functions will no longer
#' be supported, nor will there be replacement functions.  Their functionality
#' is now available in the \pkg{OpenMx} package, making these functions
#' obsolete. The utility functions \code{nullMx} and \code{saturateMx} will
#' also no longer be supported. These have already been removed from
#' \pkg{semTools}, except that \code{standardizeMx} remains deprecated due to
#' the temporary depndency on it of the \pkg{semPlot} package.  The exception
#' is that \code{\link[OpenMx]{mxStandardizeRAMpaths}} currently only provides
#' standardized estimates of covariance-structure parameters, whereas
#' \code{standardizeMx} also provides standardized means.
#'
#' @export
standardizeMx <- function(object, free = TRUE) {
  # objectOrig <- object
  multigroup <- length(object@submodels) > 0
  if(multigroup) {
    defVars <- lapply(object@submodels, findDefVars)
    defVars <- do.call(c, defVars)
  } else {
    defVars <- findDefVars(object)
  }
  if(length(defVars) > 0) stop("The standardizeMx is not available for the model with definition variable.")
  if(multigroup) {
    object@submodels <- lapply(object@submodels, standardizeMxSingleGroup)
  } else {
    object <- standardizeMxSingleGroup(object)
  }
  vectorizeMx(object, free=free)
}

## ----------------
## Hidden functions
## ----------------

findDefVars <- function(object) {
  ## borrowed from OpenMx::imxIsDefinitionVariable
  imxSeparatorChar <- "."
  imxIsDefinitionVariable <- function (name) {
    if (is.na(name)) {
      return(FALSE)
    }
    components <- unlist(strsplit(name, imxSeparatorChar, fixed = TRUE))
    if (length(components) == 2 && components[[1]] == "data") {
      return(TRUE)
    }
    else if (length(components) > 2 && components[[2]] == "data") {
      return(TRUE)
    }
    else {
      return(FALSE)
    }
  }
  ## end borrowed code
	mat <- lapply(object@matrices, slot, "labels")
	defvars <- sapply(mat, function(x) x[apply(x, c(1,2), imxIsDefinitionVariable)])
	Reduce("c", defvars)
}

vectorizeMx <- function(object, free = TRUE) {
  multigroup <- length(object@submodels) > 0
  if(multigroup) {
    object <- object@submodels
  } else {
    object <- list(object)
  }
  result <- NULL
  for(i in seq_along(object)) {
    name <- ""
    if(multigroup) name <- paste0(object[[i]]@name, ".")
    mat <- object[[i]]@matrices
    for(j in seq_along(mat)) {
      tempname <- paste0(name, mat[[j]]@name)
      lab <- mat[[j]]@labels
      tempfree <- as.vector(mat[[j]]@free)
      madeLab <- paste0(tempname, "[", row(lab), ",", col(lab), "]")
      lab <- as.vector(lab)
      madeLab[!is.na(lab)] <- lab[!is.na(lab)]
      if(!free) tempfree <- rep(TRUE, length(tempfree))
      temp <- mat[[j]]@values[tempfree]
      names(temp) <- madeLab[tempfree]
      result <- c(result, temp)
    }
  }

  result[!duplicated(names(result))]
}

standardizeMxSingleGroup <- function(object) {
  if (!is(object@expectation, "MxExpectationRAM"))
    stop("The standardizeMx function is available for the MxExpectationRAM only.")
  A <- object@matrices$A@values
  I <- diag(nrow(A))
  S <- object@matrices$S@values
  # F <- object@matrices$F@values
  Z <- solve(I - A)
  impliedCov <- Z %*% S %*% t(Z)
  temp <- sqrt(diag(impliedCov))
  if (length(temp) == 1) {
    ImpliedSd <- as.matrix(temp)
  } else {
    ImpliedSd <- diag(temp)
  }
  ImpliedInvSd <- solve(ImpliedSd)
  object@matrices$S@values <- ImpliedInvSd %*% S %*% ImpliedInvSd
  object@matrices$A@values <- ImpliedInvSd %*% A %*% ImpliedSd
  if (!is.null(object@matrices$M)) {
    M <- object@matrices$M@values
    object@matrices$M@values <- M %*% ImpliedInvSd
  }
  object
}


