### Sunthud Pornprasertmanit & Terrence D. Jorgensen
### Last updated: 27 May 2020
### fit and rotate EFA models in lavaan


## -------------------------------------
## Exploratory Factor Analysis in lavaan
## -------------------------------------

##' Analyze Unrotated Exploratory Factor Analysis Model
##'
##' This function will analyze unrotated exploratory factor analysis model. The
##' unrotated solution can be rotated by the \code{\link{orthRotate}} and
##' \code{\link{oblqRotate}} functions.
##'
##' This function will generate a lavaan script for unrotated exploratory factor
##' analysis model such that (1) all factor loadings are estimated, (2) factor
##' variances are fixed to 1, (3) factor covariances are fixed to 0, and (4) the
##' dot products of any pairs of columns in the factor loading matrix are fixed
##' to zero (Johnson & Wichern, 2002). The reason for creating this function
##' to supplement the \code{\link{factanal}} function is that users can enjoy
##' some advanced features from the \code{lavaan} package, such as scaled
##' \eqn{\chi^2}, diagonally weighted least squares estimation for ordinal
##' indicators, or full-information maximum likelihood (FIML) to handle
##' incomplete data.
##'
##' @importFrom lavaan lavInspect parTable
##' @importFrom stats factanal
##'
##' @param data A target \code{data.frame}
##' @param nf The desired number of factors
##' @param varList Target observed variables. If not specified, all variables in
##'   \code{data} will be used (or \code{sample.cov} if \code{is.null(data)};
##'   see \code{\link[lavaan]{cfa}} for argument descriptions).
##' @param start Use starting values in the analysis from the
##'   \code{\link{factanal}} \code{function}. If \code{FALSE}, the starting values
##'   from the \code{lavaan} package will be used. \code{TRUE} is ignored with a
##'   warning if the \code{aux} argument is used.
##' @param aux The list of auxiliary variables. These variables will be included
##'   in the model by the saturated-correlates approach to account for missing
##'   information.
##' @param \dots Other arguments in the \code{\link[lavaan]{cfa}} function in
##'   the \code{lavaan} package, such as \code{ordered}, \code{se},
##'   \code{estimator}, or \code{sample.cov} and \code{sample.nobs}.
##'
##' @return A \code{lavaan} output of unrotated exploratory factor analysis
##' solution.
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @examples
##'
##' unrotated <- efaUnrotate(HolzingerSwineford1939, nf = 3,
##'                          varList = paste0("x", 1:9), estimator = "mlr")
##' summary(unrotated, std = TRUE)
##' inspect(unrotated, "std")
##'
##' dat <- data.frame(HolzingerSwineford1939,
##'                   z = rnorm(nrow(HolzingerSwineford1939), 0, 1))
##' unrotated2 <- efaUnrotate(dat, nf = 2, varList = paste0("x", 1:9), aux = "z")
##'
##' @export
efaUnrotate <- function(data = NULL, nf, varList = NULL,
                        start = TRUE, aux = NULL, ...) {
  efaArgs <- list(...)
  if (is.null(data)) {
    ## check for summary statistics
    sample.cov  <- efaArgs$sample.cov
    sample.mean <- efaArgs$sample.mean
    sample.nobs <- efaArgs$sample.nobs
    sample.th   <- efaArgs$sample.th
    WLS.V       <- efaArgs$WLS.V
    NACOV       <- efaArgs$NACOV

    if (is.null(sample.cov)) stop('User must supply either raw data or ',
                                  'summary statistics to pass to lavaan().')
    if (is.null(varList)) varList <- colnames(sample.cov)

    anyOrdered <- !is.null(sample.th)
    ordNames <- efaArgs$ordered
    if (anyOrdered & (is.logical(ordNames) | is.null(ordNames))) {
      if (is.null(ordNames)) {
        message('Thresholds supplied, but not an ordered= argument. Must ',
                'assume all model variables are ordered.')
      }
      ordNames <- varList
    }

  } else {
    sample.cov  <- NULL
    sample.mean <- NULL
    sample.nobs <- NULL
    sample.th   <- NULL
    WLS.V       <- NULL
    NACOV       <- NULL

    efaArgs$data <- data
    if (is.null(varList)) varList <- colnames(data)

    anyOrdered <- checkOrdered(data, varList, ...)
    ordNames <- checkOrdered(data, varList, ..., return.names = TRUE)

    if (!is.null(efaArgs$group)) stop("Multi-group EFA is not currently supported.")
  }

  if (!is.null(aux)) {
    if (anyOrdered) {
      stop("The analysis model or the analysis data have ordered categorical",
           " variables. The auxiliary variable feature is not available for",
           " the models for categorical variables with the weighted least",
           " square approach.")
    }
    efaArgs$fixed.x <- FALSE
    efaArgs$missing <- "fiml"
    efaArgs$aux <- aux
    lavaancfa <- function(...) { cfa.auxiliary(...)}
  } else lavaancfa <- function(...) { lavaan::cfa(...)}
  nvar <- length(varList)
  facnames <- paste0("factor", 1:nf)
  loading <- outer(1:nvar, 1:nf, function(x, y) paste0("load", x, "_", y))
  syntax <- ""
  for (i in 1:nf) {
    variablesyntax <- paste(paste0(loading[,i], "*", varList), collapse = " + ")
    factorsyntax <- paste0(facnames[i], " =~ NA*", varList[1], " + ", variablesyntax, "\n")
    syntax <- paste(syntax, factorsyntax)
  }
  syntax <- paste(syntax, paste(paste0(facnames, " ~~ 1*", facnames),
                                collapse = "\n"), "\n")

  dataSupportsMeans <- length(setdiff(varList, ordNames)) && !(is.null(data) && is.null(sample.mean))
  meanstructure <- efaArgs$meanstructure
  if (is.null(meanstructure)) meanstructure <- anyOrdered #FIXME: wise default for EFA?
  stopifnot(is.logical(meanstructure))
  if (meanstructure && dataSupportsMeans) {
    syntax <- paste(syntax, paste(paste0(setdiff(varList, ordNames),
                                         " ~ 1"), collapse = "\n"), "\n")
  }

  if (nf > 1) {
    covsyntax <- outer(facnames, facnames,
                       function(x, y) paste0(x, " ~~ 0*", y, "\n"))[lower.tri(diag(nf), diag = FALSE)]
    syntax <- paste(syntax, paste(covsyntax, collapse = " "))
    for (i in 2:nf) {
      for (j in 1:(i - 1)) {
        loadconstraint <- paste(paste0(loading[,i], "*", loading[,j]), collapse=" + ")
        syntax <- paste(syntax, paste0("0 == ", loadconstraint), "\n")
      }
    }
  }
  if (start) {
    if (is.null(aux)) {
      List <- c(list(model = syntax, data = data), list(...))
      List$do.fit <- FALSE
      outtemp <- do.call(lavaancfa, List)
      covtemp <- lavInspect(outtemp, "sampstat")$cov
      partemp <- parTable(outtemp)
      err <- try(startload <- factanal(factors = nf, covmat = covtemp)$loadings[],
                 silent = TRUE)
      if (is(err, "try-error")) stop("The starting values from the factanal",
                                     " function cannot be calculated. Please",
                                     " use start = FALSE instead.")
      startval <- sqrt(diag(diag(covtemp))) %*% startload
      partemp$ustart[match(as.vector(loading), partemp$label)] <- as.vector(startval)
      partemp$est <- partemp$se <- NULL
      syntax <- partemp
    } else warning("The 'start' argument was ignored because factanal() does",
                   " not support auxiliary variables.  When using auxiliary",
                   " variables, set 'start = FALSE' ")
  } else {
    ## FIXME: optimizer can't escape loadings == 0 without starting values from
    ##        factanal() or by using theta parameterization
    ##        https://groups.google.com/d/msg/lavaan/ujkHmCVirEY/-LGut4ewAwAJ
    parameterization <- efaArgs$parameterization
    if (is.null(parameterization)) parameterization <- lavaan::lavOptions("parameterization")
    if (anyOrdered && parameterization != "theta")
      warning('If the default parameterization = "delta" returns results with ',
              'all factor loadings equal to zero, try either setting start ',
              '= TRUE or setting parameterization = "theta" instead.')
  }
  efaArgs$model <- syntax
  do.call(lavaancfa, efaArgs)
}



## -----------------
## Class and Methods
## -----------------

##' Class For Rotated Results from EFA
##'
##' This class contains the results of rotated exploratory factor analysis
##'
##'
##' @name EFA-class
##' @aliases EFA-class show,EFA-method summary,EFA-method
##' @docType class
##' @slot loading Rotated standardized factor loading matrix
##' @slot rotate Rotation matrix
##' @slot gradRotate gradient of the objective function at the rotated loadings
##' @slot convergence Convergence status
##' @slot phi: Factor correlation matrix. Will be an identity matrix if
##'  orthogonal rotation is used.
##' @slot se Standard errors of the rotated standardized factor loading matrix
##' @slot method Method of rotation
##' @slot call The command used to generate this object
##' @section Objects from the Class: Objects can be created via the
##' \code{\link{orthRotate}} or \code{\link{oblqRotate}} function.
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##' @seealso \code{\link{efaUnrotate}}; \code{\link{orthRotate}};
##' \code{\link{oblqRotate}}
##' @examples
##'
##' unrotated <- efaUnrotate(HolzingerSwineford1939, nf = 3,
##'                          varList = paste0("x", 1:9), estimator = "mlr")
##' summary(unrotated, std = TRUE)
##' lavInspect(unrotated, "std")
##'
##' # Rotated by Quartimin
##' rotated <- oblqRotate(unrotated, method = "quartimin")
##' summary(rotated)
##'
setClass("EFA", representation(loading = "matrix",
                               rotate = "matrix",
                               gradRotate = "matrix",
                               convergence = "logical",
                               phi = "matrix",
                               se = "matrix",
                               method = "character",
                               call = "call"))

##' @rdname EFA-class
##' @aliases show,EFA-method
##' @export
setMethod("show", signature(object = "EFA"), function(object) {
	cat("Standardized Rotated Factor Loadings\n")
	print(printLoadings(object))
	cat("\nFactor Correlation\n")
	print(object@phi)
	cat("\nMethod of rotation:\t")
	cat(object@method, "\n")
	message("The standard errors are close but do not match with other packages.",
	        " Be mindful when using it.")
})

##' @rdname EFA-class
##' @aliases summary,EFA-method
##' @param object object of class \code{EFA}
##' @param suppress any standardized loadings less than the specified value
##'  will not be printed to the screen
##' @param sort \code{logical}. If \code{TRUE} (default), factor loadings will
##'  be sorted by their size in the console output
##' @export
setMethod("summary", signature(object = "EFA"),
          function(object, suppress = 0.1, sort = TRUE) {
	cat("Standardized Rotated Factor Loadings\n")
	print(printLoadings(object, suppress = suppress, sort = sort))
	cat("\nFactor Correlation\n")
	print(object@phi)
	cat("\nMethod of rotation:\t")
	cat(object@method, "\n")
	cat("\nTest Statistics for Standardized Rotated Factor Loadings\n")
	print(testLoadings(object))
})



## ------------------------------
## Rotation Constructor Functions
## ------------------------------

##' Implement orthogonal or oblique rotation
##'
##' These functions will implement orthogonal or oblique rotation on
##' standardized factor loadings from a lavaan output.
##'
##' These functions will rotate the unrotated standardized factor loadings by
##' orthogonal rotation using the \code{\link[GPArotation]{GPForth}} function or
##' oblique rotation using the \code{\link[GPArotation]{GPFoblq}} function the
##' \code{GPArotation} package. The resulting rotation matrix will be used to
##' calculate standard errors of the rotated standardized factor loading by
##' delta method by numerically computing the Jacobian matrix by the
##' \code{\link[lavaan]{lav_func_jacobian_simple}} function.
##'
##' @aliases orthRotate oblqRotate funRotate
##' @rdname rotate
##' @param object A lavaan output
##' @param method The method of rotations, such as \code{"varimax"},
##' \code{"quartimax"}, \code{"geomin"}, \code{"oblimin"}, or any gradient
##' projection algorithms listed in the \code{\link[GPArotation]{GPA}} function
##' in the \code{GPArotation} package.
##' @param fun The name of the function that users wish to rotate the
##' standardized solution. The functions must take the first argument as the
##' standardized loading matrix and return the \code{GPArotation} object. Check
##' this page for available functions: \code{\link[GPArotation]{rotations}}.
##' @param \dots Additional arguments for the \code{\link[GPArotation]{GPForth}}
##' function (for \code{orthRotate}), the \code{\link[GPArotation]{GPFoblq}}
##' function (for \code{oblqRotate}), or the function that users provide in the
##' \code{fun} argument.
##' @return An \code{linkS4class{EFA}} object that saves the rotated EFA solution
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##' @examples
##'
##' \dontrun{
##'
##' unrotated <- efaUnrotate(HolzingerSwineford1939, nf = 3,
##'                          varList = paste0("x", 1:9), estimator = "mlr")
##'
##' # Orthogonal varimax
##' out.varimax <- orthRotate(unrotated, method = "varimax")
##' summary(out.varimax, sort = FALSE, suppress = 0.3)
##'
##' # Orthogonal Quartimin
##' orthRotate(unrotated, method = "quartimin")
##'
##' # Oblique Quartimin
##' oblqRotate(unrotated, method = "quartimin")
##'
##' # Geomin
##' oblqRotate(unrotated, method = "geomin")
##'
##' # Target rotation
##' library(GPArotation)
##' target <- matrix(0, 9, 3)
##' target[1:3, 1] <- NA
##' target[4:6, 2] <- NA
##' target[7:9, 3] <- NA
##' colnames(target) <- c("factor1", "factor2", "factor3")
##' ## This function works with GPArotation version 2012.3-1
##' funRotate(unrotated, fun = "targetQ", Target = target)
##' }
##'
##' @export
orthRotate <- function(object, method = "varimax", ...) {
	requireNamespace("GPArotation")
	if (!("package:GPArotation" %in% search())) attachNamespace("GPArotation")
	mc <- match.call()
	initL <- getLoad(object)
	rotated <- GPArotation::GPForth(initL, method=method, ...)
	# rotateMat <- t(solve(rotated$Th)) # defined but never used
	LIST <- seStdLoadings(rotated, object, fun = GPArotation::GPForth,
	                      MoreArgs = c(method = method, list(...)))
	# orthogonal <- rotated$orthogonal # define but never used
	loading <- rotated$loadings
	rotate <- rotated$Th
	gradRotate <- rotated$Gq
	convergence <- rotated$convergence
	method <- rotated$method
	phi <- diag(ncol(loading))
	lv.names <- colnames(loading)
	dimnames(phi) <- list(lv.names, lv.names)
	new("EFA", loading = loading, rotate = rotate, gradRotate = gradRotate,
	    convergence = convergence, phi = phi, se = LIST, method = method, call = mc)
}



##' @rdname rotate
##' @export
oblqRotate <- function(object, method = "quartimin", ...) {
	requireNamespace("GPArotation")
	if (!("package:GPArotation" %in% search())) attachNamespace("GPArotation")
	mc <- match.call()
	initL <- getLoad(object)
	rotated <- GPArotation::GPFoblq(initL, method = method, ...)
	# rotateMat <- t(solve(rotated$Th)) # defined but never used
	LIST <- seStdLoadings(rotated, object, fun = GPArotation::GPFoblq,
	                      MoreArgs = c(method = method, list(...)))
	# orthogonal <- rotated$orthogonal # defined but never used
	loading <- rotated$loadings
	rotate <- rotated$Th
	gradRotate <- rotated$Gq
	convergence <- rotated$convergence
	method <- rotated$method
	phi <- rotated$Phi
	lv.names <- colnames(loading)
	dimnames(phi) <- list(lv.names, lv.names)
	new("EFA", loading=loading, rotate = rotate, gradRotate = gradRotate,
	    convergence = convergence, phi = phi, se = LIST, method = method, call = mc)
}



##' @rdname rotate
##' @export
funRotate <- function(object, fun, ...) {
	stopifnot(is.character(fun))
	requireNamespace("GPArotation")
	if (!("package:GPArotation" %in% search())) attachNamespace("GPArotation")
	mc <- match.call()
	initL <- getLoad(object)
	rotated <- do.call(fun, c(list(L = initL), list(...)))
	# rotateMat <- t(solve(rotated$Th)) # defined but never used
	gradRotate <- rotated$Gq
	LIST <- seStdLoadings(rotated, object, fun = fun, MoreArgs = list(...))
	# orthogonal <- rotated$orthogonal # defined but never used
	loading <- rotated$loadings
	rotate <- rotated$Th
	convergence <- rotated$convergence
	method <- rotated$method
	phi <- rotated$Phi
	if (is.null(phi)) phi <- diag(ncol(loading))
	lv.names <- colnames(loading)
	dimnames(phi) <- list(lv.names, lv.names)
	new("EFA", loading = loading, rotate = rotate, gradRotate = gradRotate,
	    convergence = convergence, phi = phi, se = LIST, method = method, call = mc)
}



## ----------------
## Hidden Functions
## ----------------

printLoadings <- function(object, suppress = 0.1, sort = TRUE) {
  loading <- object@loading
  nf <- ncol(loading)
  loadingText <- sprintf("%.3f", object@loading)
  sig <- ifelse(testLoadings(object)$p < 0.05, "*", " ")
  loadingText <- paste0(loadingText, sig)
  loadingText[abs(loading) < suppress] <- ""
  loadingText <- matrix(loadingText, ncol = nf, dimnames = dimnames(loading))
  lead <- apply(abs(loading), 1, which.max)
  ord <- NULL
  if (sort) {
    for (i in 1:nf) {
      ord <- c(ord, intersect(order(abs(loading[,i]), decreasing = TRUE), which(lead == i)))
    }
    loadingText <- loadingText[ord,]
  }
  as.data.frame(loadingText)
}

##' @importFrom stats pnorm qnorm
testLoadings <- function(object, level = 0.95) {
  se <- object@se
  loading <- object@loading
  lv.names <- colnames(loading)
  z <- loading / se
  p <- 2 * (1 - pnorm( abs(z) ))
  crit <- qnorm(1 - (1 - level)/2)
  est <- as.vector(loading)
  se <- as.vector(se)
  warning("The standard error is currently invalid because it does not account",
          " for the variance of the rotation function. It is simply based on",
          " the delta method.")
  out <- data.frame(lhs=lv.names[col(loading)], op = "=~",
                    rhs = rownames(loading)[row(loading)], std.loading = est,
                    se = se, z = as.vector(z), p = as.vector(p),
                    ci.lower = (est - crit*se), ci.upper = (est + crit*se))
  class(out) <- c("lavaan.data.frame", "data.frame")
  out
}

##' @importFrom lavaan lavInspect
getLoad <- function(object, std = TRUE) {
	out <- lavInspect(object, "est")$lambda #FIXME: check for multiple groups
	if (std) {
		impcov <- lavaan::fitted.values(object)$cov
		impsd <- sqrt(diag(diag(impcov)))
		out <- solve(impsd) %*% out
	}
	rownames(out) <- lavaan::lavNames(object@ParTable, "ov", group = 1)
	if (!is.null(object@external$aux)) {
		out <- out[!(rownames(out) %in% object@external$aux),]
	}
	class(out) <- c("loadings", out)
	out
}

fillMult <- function(X, Y, fillrowx = 0, fillrowy = 0, fillcolx = 0, fillcoly = 0) {
	tempX <- matrix(0, nrow = nrow(X) + fillrowx, ncol = ncol(X) + fillcolx)
	tempY <- matrix(0, nrow = nrow(Y) + fillrowy, ncol = ncol(Y) + fillcoly)
	tempX[1:nrow(X), 1:ncol(X)] <- X
	tempY[1:nrow(Y), 1:ncol(Y)] <- Y
	result <- tempX %*% tempY
	result[1:nrow(X), 1:ncol(Y)]
}

##' @importFrom lavaan lavInspect
stdRotatedLoadings <- function(est, object, fun, aux = NULL, rotate = NULL, MoreArgs = NULL) {
	ov.names <- lavaan::lavNames(object@ParTable, "ov", group = 1)
  lv.names <- lavaan::lavNames(object@ParTable, "lv", group = 1)
	ind.names <- setdiff(ov.names, aux)

	# Compute model-implied covariance matrix
	partable <- object@ParTable
	# LY
	load.idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names))
	loading <- matrix(est[load.idx], ncol = length(lv.names))
	loading <- rbind(loading, matrix(0, length(aux), ncol(loading)))
	# Nu
	if (lavInspect(object, "options")$meanstructure) {
	  int.idx <- which(partable$op == "~1" & (partable$rhs == "") & (partable$lhs %in% ov.names))
	  intcept <- matrix(est[int.idx], ncol = 1)
	}
	# Theta
	th.idx <- which(partable$op == "~~" & (partable$rhs %in% ov.names) & (partable$lhs %in% ov.names))
	theta <- matrix(0, nrow = length(ov.names), ncol = length(ov.names),
	                dimnames = list(ov.names, ov.names))
	for (i in th.idx) {
		theta[partable$lhs[i], partable$rhs[i]] <- theta[partable$rhs[i], partable$lhs[i]] <- est[i]
	}
	OV <- loading %*% t(loading) + theta
	invsd <- solve(sqrt(diag(diag(OV))))

	requireNamespace("GPArotation")
	if (!("package:GPArotation" %in% search())) attachNamespace("GPArotation")
	# Compute standardized results
	loading <- invsd %*% loading

	obj <- do.call(fun, c(list(loading), MoreArgs))

	# GPArotation::GPFoblq(loading, method = "geomin")
	loading <- obj$loadings
	rotMat <- t(solve(obj$Th))

	# %*% rotate
	est[load.idx] <- as.vector(loading[seq_along(ind.names),])
	if (lavInspect(object, "options")$meanstructure) {
	  est[int.idx] <- as.vector(invsd %*% intcept)
	}
	theta <- invsd %*% theta %*% invsd
	rownames(theta) <- colnames(theta) <- ov.names
	for(i in th.idx)  est[i] <- theta[partable$lhs[i], partable$rhs[i]]

	# Put phi
	rv.idx <- which(partable$op == "~~" & partable$rhs %in% lv.names)
	templhs <- match(partable$lhs[rv.idx], lv.names)
	temprhs <- match(partable$rhs[rv.idx], lv.names)
	# rotate2 <- t(solve(rotate))
	# phi <- t(rotate2) %*% rotate2
	phi <- obj$Phi
	if (!is.null(phi)) {
		for (i in seq_along(templhs)) {
			est[rv.idx[i]] <- phi[templhs[i], temprhs[i]]
		}
	}
	est
}

##' @importFrom lavaan lavInspect parTable
seStdLoadings <- function(rotate, object, fun, MoreArgs) {
	# object <- efaUnrotate(HolzingerSwineford1939, nf=3, varList=paste0("x", 1:9), estimator="mlr")
	# initL <- getLoad(object)
	# rotate <- GPArotation::GPFoblq(initL, method="oblimin")

	rotMat <- t(solve(rotate$Th))
	gradient <- rotate$Gq
	loading <- rotate$loadings
	phi <- rotate$Phi
	if (is.null(phi)) phi <- diag(ncol(loading))

	est <- lavaan::parameterEstimates(object)$est
	aux <- object@external$aux

	# Standardized results
	JAC1 <- lavaan::lav_func_jacobian_simple(func = stdRotatedLoadings,
	                                         x = object@Fit@est, object = object,
	                                         aux = aux, rotate = rotMat,
	                                         fun = fun, MoreArgs = MoreArgs)

	LIST <- lavInspect(object, "list")
	free.idx <- which(LIST$free > 0L)
	m <- ncol(phi)
	phi.idx <- which(LIST$op == "~~" & LIST$lhs != LIST$rhs & (LIST$lhs %in% paste0("factor", 1:m)))
	JAC1 <- JAC1[c(free.idx, phi.idx), free.idx]
	VCOV <- as.matrix(lavaan::vcov(object, labels = FALSE))
	if(object@Model@eq.constraints) {
		JAC1 <- JAC1 %*% object@Model@eq.constraints.K
	}
	COV1 <- JAC1 %*% VCOV %*% t(JAC1)
	# I1 <- MASS::ginv(COV1)
	# I1p <- matrix(0, nrow(I1) + length(phi.idx), ncol(I1) + length(phi.idx))
	# I1p[1:nrow(I1), 1:ncol(I1)] <- I1
	# phi.idx2 <- nrow(I1) + 1:length(phi.idx)


	# p <- nrow(loading)
	# dconlambda <- matrix(0, m^2 - m, p*m)
	# gradphi <- gradient %*% solve(phi)
	# lambgradphi <- t(loading) %*% gradphi
	# lambphi <- loading %*% solve(phi)
	# lamblamb <- t(loading) %*% loading
	# runrow <- 1
	# descript <- NULL
	# for(u in 1:m) {
		# for(v in setdiff(1:m, u)) {
			# runcol <- 1
			# for(r in 1:m) {
				# for(i in 1:p) {
					# mir <- (1 - 1/p) * sum(loading[i,]^2) + sum(loading[,r]^2)/p - loading[i,r]^2
					# dur <- 0
					# if(u == r) dur <- 1
					# dconlambda[runrow, runcol] <- dur * gradphi[i, v] + 4 * mir * loading[i, u] * phi[r, v] + (8 - 8/p)*loading[i,r]*loading[i,u]*lambphi[i,v] + 8*loading[i,r]*lamblamb[u,r]*phi[r,v]/p - 8*loading[i,r]^2*loading[i,u]*phi[r,v]
					# descript <- rbind(descript, c(runrow, runcol, u, v, i, r))
					# runcol <- runcol + 1
				# }
			# }
			# runrow <- runrow + 1
		# }
	# }

	# dconphi <- matrix(0, m^2 - m, m*(m-1)/2)
	# runrow <- 1
	# descript2 <- NULL
	# for(u in 1:m) {
		# for(v in setdiff(1:m, u)) {
			# runcol <- 1
			# for(x in 2:m) {
				# for(y in 1:(x - 1)) {
					# dux <- 0
					# if(u == x) dux <- 1
					# duy <- 0
					# if(u == y) duy <- 1
					# dconphi[runrow, runcol] <- -(dux * phi[y, v] + duy * phi[x, v]) * lambgradphi[u, u]
					# descript2 <- rbind(descript2, c(runrow, runcol, u, v, x, y))
					# runcol <- runcol + 1
				# }
			# }
			# runrow <- runrow + 1
		# }
	# }
	# I2 <- matrix(0, nrow(I1p) + m^2 - m, ncol(I1p) + m^2 - m)
	# I2[1:nrow(I1p), 1:ncol(I1p)] <- I1p
	# I2[lamb.idx, 1:(m^2 - m) + nrow(I1p)] <- t(dconlambda)
	# I2[1:(m^2 - m) + nrow(I1p), lamb.idx] <- dconlambda
	# I2[phi.idx2, 1:(m^2 - m) + nrow(I1p)] <- t(dconphi)
	# I2[1:(m^2 - m) + nrow(I1p), phi.idx2] <- dconphi
	# COV2 <- MASS::ginv(I2)[1:nrow(I1p), 1:ncol(I1p)]

	COV2 <- COV1
	LIST <- LIST[,c("lhs", "op", "rhs", "group")]
	LIST$se <- rep(NA, length(LIST$lhs))
	LIST$se[c(free.idx, phi.idx)] <- sqrt(diag(COV2))
	tmp.se <- ifelse( LIST$se == 0.0, NA, LIST$se)
  lv.names <- lavaan::lavNames(object@ParTable, "lv", group = 1)
	partable <- parTable(object)
	idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names))
	matrix(LIST$se[idx], ncol = length(lv.names))
}

checkOrdered <- function(dat, varnames, ..., return.names = FALSE) {
  ord <- list(...)$ordered
  if (is.null(ord)) {
    ord <- FALSE
  } else {
    ord <- TRUE
  }
  if (is.null(dat)) {
    orderedVar <- FALSE
  } else {
    orderedVar <- sapply(dat[,varnames], function(x) "ordered" %in% is(x))
  }

  if (return.names) {
    ## added by TDJ 4-Feb-2020
    return(unique(c(list(...)$ordered, names(orderedVar[orderedVar]))))
  }
  any(c(ord, orderedVar))
}


