### Sunthud Pornprasertmanit
### Last updated: 14 October 2016
### run EFA model in lavaan

setClass("EFA", representation(loading = "matrix", rotate="matrix", gradRotate="matrix", convergence="logical", phi="matrix", se = "matrix", method = "character", call="call"))

printLoadings <- function(object, suppress = 0.1, sort=TRUE) {
	loading <- object@loading
	nf <- ncol(loading)
	loadingText <- sprintf("%.3f", object@loading)
	sig <- ifelse(testLoadings(object)$p < 0.05, "*", " ")
	loadingText <- paste0(loadingText, sig)
	loadingText[abs(loading) < suppress] <- ""
	loadingText <- matrix(loadingText, ncol=nf, dimnames=dimnames(loading))
	lead <- apply(abs(loading), 1, which.max)
	ord <- NULL
	if(sort) {
		for(i in 1:nf) {
			ord <- c(ord, intersect(order(abs(loading[,i]), decreasing=TRUE), which(lead==i)))
		}
		loadingText <- loadingText[ord,]
	}
	as.data.frame(loadingText)
}

testLoadings <- function(object, level=0.95) {
	se <- object@se
	loading <- object@loading
	lv.names <- colnames(loading)
	z <- loading/se
	p <- 2 * (1 - pnorm( abs(z) ))
	crit <- qnorm(1 - (1 - level)/2)
	est <- as.vector(loading)
	se <- as.vector(se)
	warnings("The standard error is currently invalid because it does not account for the variance of rotation function. It is simply based on delta method.")
	out <- data.frame(lhs=lv.names[col(loading)], op="=~", rhs=rownames(loading)[row(loading)], std.loading=est, se=se, z=as.vector(z), p=as.vector(p), ci.lower=(est - crit*se), ci.upper=(est + crit*se))
	class(out) <- c("lavaan.data.frame", "data.frame")
	out
}

setMethod("show", signature(object = "EFA"), function(object) {	
	cat("Standardized Rotated Factor Loadings\n")
	print(printLoadings(object))
	cat("\nFactor Correlation\n")
	print(object@phi)
	cat("\nMethod of rotation:\t")
	cat(object@method, "\n")
	print("The standard errors are close but do not match with other packages. Be mindful when using it.")
}) 

setMethod("summary", signature(object = "EFA"), function(object, suppress = 0.1, sort = TRUE) {
	cat("Standardized Rotated Factor Loadings\n")
	print(printLoadings(object, suppress = suppress, sort = sort))
	cat("\nFactor Correlation\n")
	print(object@phi)	
	cat("\nMethod of rotation:\t")
	cat(object@method, "\n")
	cat("\nTest Statistics for Standardized Rotated Factor Loadings\n")
	print(testLoadings(object))
}) 

efaUnrotate <- function(data, nf, varList=NULL, start=TRUE, aux=NULL, ...) {
	if(is.null(varList)) varList <- colnames(data)
	lavaancfa <- function(...) { lavaan::cfa(...)}
	nvar <- length(varList)
	facnames <- paste0("factor", 1:nf)
	loading <- outer(1:nvar, 1:nf, function(x, y) paste0("load", x, "_", y))
	syntax <- ""
	for(i in 1:nf) {
		variablesyntax <- paste(paste0(loading[,i], "*", varList), collapse=" + ")
		factorsyntax <- paste0(facnames[i], " =~ NA*", varList[1], " + ", variablesyntax, "\n")
		syntax <- paste(syntax, factorsyntax)
	}
	syntax <- paste(syntax, paste(paste0(facnames, " ~~ 1*", facnames), collapse="\n"), "\n")

	isOrdered <- checkOrdered(data, varList, ...)
	if(!isOrdered) {
		syntax <- paste(syntax, paste(paste0(varList, " ~ 1"), collapse="\n"), "\n")
	}
	
	if(nf > 1) {
		covsyntax <- outer(facnames, facnames, function(x, y) paste0(x, " ~~ 0*", y, "\n"))[lower.tri(diag(nf), diag=FALSE)]
		syntax <- paste(syntax, paste(covsyntax, collapse = " "))	
		for(i in 2:nf) {
			for(j in 1:(i - 1)) {
				loadconstraint <- paste(paste0(loading[,i], "*", loading[,j]), collapse=" + ")
				syntax <- paste(syntax, paste0("0 == ", loadconstraint), "\n")
			}
		}
	}
	if(start) {
		List <- c(list(model=syntax, data=data), list(...))
		List$do.fit <- FALSE
		outtemp <- do.call(lavaancfa, List)
		covtemp <- lavaan::lavInspect(outtemp, "sampstat")$cov
		partemp <- lavaan::parTable(outtemp)
		err <- try(startload <- factanal(factors=nf, covmat=covtemp)$loadings[], silent = TRUE)
		if(is(err, "try-error")) stop("The starting values from the factanal function cannot be calculated. Please use start=FALSE instead.")
		startval <- sqrt(diag(diag(covtemp))) %*% startload
		partemp$ustart[match(as.vector(loading), partemp$label)] <- as.vector(startval)
		partemp$est <- partemp$se <- partemp$start <- NULL
		syntax <- partemp
	} 
	args <- list(...)
	args$model <- syntax
	args$data <- data
	if(!is.null(aux)) {
		if(isOrdered) {
			stop("The analysis model or the analysis data have ordered categorical variables. The auxiliary variable feature is not available for the models for categorical variables with the weighted least square approach.")
		}
		auxResult <- craftAuxParTable(syntax, aux = aux)
		args$model <- auxResult$model
		args$fixed.x <- FALSE
		args$missing <- "fiml"		
		result <- do.call(lavaancfa, args)
		codeNull <- nullAuxiliary(aux, auxResult$indName, NULL, any(syntax$op == "~1"), max(syntax$group))
		resultNull <- lavaan::lavaan(codeNull, data=data, ...)
		result <- as(result, "lavaanStar")
		fit <- lavaan::fitMeasures(resultNull)
		name <- names(fit)
		fit <- as.vector(fit)
		names(fit) <- name
		result@nullfit <- fit
		result@auxNames <- aux
		return(result)
	} else {
		return(do.call(lavaancfa, args))
	}
}

getLoad <- function(object, std = TRUE) {
	out <- lavaan::inspect(object, "coef")$lambda
	if(std) {
		impcov <- lavaan::fitted.values(object)$cov
		impsd <- sqrt(diag(diag(impcov)))
		out <- solve(impsd) %*% out
	}
	rownames(out) <- lavaan::lavNames(object@ParTable, "ov", group = 1)
	if(is(object, "lavaanStar")) {
		out <- out[!(rownames(out) %in% object@auxNames),]
	}
	class(out) <- c("loadings", out)
	out
}

orthRotate <- function(object, method="varimax", ...) {
	requireNamespace("GPArotation")
	if(!("package:GPArotation" %in% search())) attachNamespace("GPArotation")
	mc <- match.call()
	initL <- getLoad(object)
	rotated <- GPArotation::GPForth(initL, method=method, ...)
	rotateMat <- t(solve(rotated$Th))
	LIST <- seStdLoadings(rotated, object, fun = GPArotation::GPForth, MoreArgs = c(method = method, list(...)))
	orthogonal <- rotated$orthogonal
	loading <- rotated$loadings
	rotate <- rotated$Th
	gradRotate <- rotated$Gq
	convergence <- rotated$convergence
	method <- rotated$method
	phi <- diag(ncol(loading))
	lv.names <- colnames(loading)
	dimnames(phi) <- list(lv.names, lv.names)
	new("EFA", loading=loading, rotate=rotate, gradRotate=gradRotate, convergence=convergence, phi=phi, se=LIST, method=method, call=mc) 
}

oblqRotate <- function(object, method="quartimin", ...) {
	requireNamespace("GPArotation")
	if(!("package:GPArotation" %in% search())) attachNamespace("GPArotation")
	mc <- match.call()
	initL <- getLoad(object)
	rotated <- GPArotation::GPFoblq(initL, method=method, ...)
	rotateMat <- t(solve(rotated$Th))
	LIST <- seStdLoadings(rotated, object, fun = GPArotation::GPFoblq, MoreArgs = c(method = method, list(...)))
	orthogonal <- rotated$orthogonal
	loading <- rotated$loadings
	rotate <- rotated$Th
	gradRotate <- rotated$Gq
	convergence <- rotated$convergence
	method <- rotated$method
	phi <- rotated$Phi
	lv.names <- colnames(loading)
	dimnames(phi) <- list(lv.names, lv.names)
	new("EFA", loading=loading, rotate=rotate, gradRotate=gradRotate, convergence=convergence, phi=phi, se=LIST, method=method, call=mc) 
}

funRotate <- function(object, fun, ...) {
	stopifnot(is.character(fun))
	requireNamespace("GPArotation")
	if(!("package:GPArotation" %in% search())) attachNamespace("GPArotation")
	mc <- match.call()
	initL <- getLoad(object)
	rotated <- do.call(fun, c(list(L = initL), list(...)))
	rotateMat <- t(solve(rotated$Th))
	gradRotate <- rotated$Gq
	LIST <- seStdLoadings(rotated, object, fun = fun, MoreArgs = list(...))
	orthogonal <- rotated$orthogonal
	loading <- rotated$loadings
	rotate <- rotated$Th
	convergence <- rotated$convergence
	method <- rotated$method
	phi <- rotated$Phi
	if(is.null(phi)) phi <- diag(ncol(loading))
	lv.names <- colnames(loading)
	dimnames(phi) <- list(lv.names, lv.names)
	new("EFA", loading=loading, rotate=rotate, gradRotate=gradRotate, convergence=convergence, phi=phi, se=LIST, method=method, call=mc) 
}

fillMult <- function(X, Y, fillrowx = 0, fillrowy = 0, fillcolx = 0, fillcoly = 0) {
	tempX <- matrix(0, nrow = nrow(X) + fillrowx, ncol = ncol(X) + fillcolx)
	tempY <- matrix(0, nrow = nrow(Y) + fillrowy, ncol = ncol(Y) + fillcoly)
	tempX[1:nrow(X), 1:ncol(X)] <- X
	tempY[1:nrow(Y), 1:ncol(Y)] <- Y
	result <- tempX %*% tempY
	result[1:nrow(X), 1:ncol(Y)]
}

stdRotatedLoadings <- function(est, object, fun, aux=NULL, rotate=NULL, MoreArgs = NULL) {
	ov.names <- lavaan::lavNames(object@ParTable, "ov", group = 1)
    lv.names <- lavaan::lavNames(object@ParTable, "lv", group = 1)
	ind.names <- setdiff(ov.names, aux)

	# Compute model-implied covariance matrix
	partable <- object@ParTable
	# LY
	load.idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names))
	loading <- matrix(est[load.idx], ncol=length(lv.names))
	loading <- rbind(loading, matrix(0, length(aux), ncol(loading)))
	# Nu
	int.idx <- which(partable$op == "~1" & (partable$rhs == "") & (partable$lhs %in% ov.names))
	intcept <- matrix(est[int.idx], ncol = 1)
	
	# Theta
	th.idx <- which(partable$op == "~~" & (partable$rhs %in% ov.names) & (partable$lhs %in% ov.names))
	theta <- matrix(0, length(ov.names), length(ov.names), dimnames = list(ov.names, ov.names))
	for(i in th.idx) {
		theta[partable$lhs[i], partable$rhs[i]] <- theta[partable$rhs[i], partable$lhs[i]] <- est[i]
	}
	OV <- loading %*% t(loading) + theta
	invsd <- solve(sqrt(diag(diag(OV))))

	requireNamespace("GPArotation")
	if(!("package:GPArotation" %in% search())) attachNamespace("GPArotation")
	# Compute standardized results
	loading <- invsd %*% loading 
	
	obj <- do.call(fun, c(list(loading), MoreArgs))
	
	# GPArotation::GPFoblq(loading, method="geomin")
	loading <- obj$loadings
	rotMat <- t(solve(obj$Th))

	# %*% rotate
	est[load.idx] <- as.vector(loading[seq_along(ind.names),])
	intcept <- invsd %*% intcept
	est[int.idx] <- as.vector(intcept)
	theta <- invsd %*% theta %*% invsd
	rownames(theta) <- colnames(theta) <- ov.names
	for(i in th.idx) {
		est[i] <- theta[partable$lhs[i], partable$rhs[i]]
	}
	
	# Put phi
	rv.idx <- which(partable$op == "~~" & partable$rhs %in% lv.names)
	templhs <- match(partable$lhs[rv.idx], lv.names)
	temprhs <- match(partable$rhs[rv.idx], lv.names)
	# rotate2 <- t(solve(rotate))
	# phi <- t(rotate2) %*% rotate2
	phi <- obj$Phi
	if(!is.null(phi)) {
		for(i in seq_along(templhs)) {
			est[rv.idx[i]] <- phi[templhs[i], temprhs[i]]
		}
	}
	est
}

seStdLoadings <- function(rotate, object, fun, MoreArgs) {
	# object <- efaUnrotate(HolzingerSwineford1939, nf=3, varList=paste0("x", 1:9), estimator="mlr")
	# initL <- getLoad(object)
	# rotate <- GPArotation::GPFoblq(initL, method="oblimin")
	
	rotMat <- t(solve(rotate$Th))
	gradient <- rotate$Gq
	loading <- rotate$loadings
	phi <- rotate$Phi
	if(is.null(phi)) phi <- diag(ncol(loading))

	est <- lavaan::parameterEstimates(object)$est
	aux <- NULL
	if(is(object, "lavaanStar")) {
		aux <- object@auxNames
	}
	# Standardized results
	JAC1 <- lavaan::lav_func_jacobian_simple(func=stdRotatedLoadings, x=object@Fit@est, object=object, aux=aux, rotate = rotMat, fun = fun, MoreArgs = MoreArgs)
		
	LIST <- lavaan::lavInspect(object, "list")
	free.idx <- which(LIST$free > 0L)
	m <- ncol(phi)
	phi.idx <- which(LIST$op == "~~" & LIST$lhs != LIST$rhs & (LIST$lhs %in% paste0("factor", 1:m)))
	JAC1 <- JAC1[c(free.idx, phi.idx), free.idx]
	VCOV <- as.matrix(lavaan::vcov(object, labels=FALSE))
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
	partable <- lavaan::parTable(object)
	idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names))
	matrix(LIST$se[idx], ncol=length(lv.names))
}
