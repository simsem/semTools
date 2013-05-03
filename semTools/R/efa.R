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

efaUnrotate <- function(data, nf, varList=NULL, start=TRUE, ...) {
	if(is.null(varList)) varList <- colnames(data)
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
		outtemp <- do.call("cfa", List)
		covtemp <- outtemp@SampleStats@cov[[1]]
		partemp <- parTable(outtemp)
		startval <- sqrt(diag(diag(covtemp))) %*% factanal(factors=nf, covmat=covtemp)$loadings[] 
		partemp$ustart[match(as.vector(loading), partemp$label)] <- as.vector(startval)
		return(cfa(partemp, data, ...))
	} else {
		return(cfa(model=syntax, data=data, ...))
	}
}

stdLoad <- function(object) {
	out <- solve(sqrt(diag(diag(fitted.values(object)$cov)))) %*% inspect(object, "coef")$lambda
	rownames(out) <- lavaan:::vnames(object@ParTable, "ov", group = 1)
	class(out) <- c("loadings", out)
	out
}

orthRotate <- function(object, method="varimax", ...) {
	mc <- match.call()
	library(GPArotation)
	initL <- stdLoad(object)
	rotated <- GPForth(initL, method=method, ...)
	rotateMat <- t(solve(rotated$Th))
	LIST <- seStdLoadings(rotateMat, object)
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
	mc <- match.call()
	library(GPArotation)
	initL <- stdLoad(object)
	rotated <- GPFoblq(initL, method=method, ...)
	rotateMat <- t(solve(rotated$Th))
	LIST <- seStdLoadings(rotateMat, object)
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
	mc <- match.call()
	library(GPArotation)
	initL <- stdLoad(object)
	rotated <- do.call(fun, c(list(L = initL), list(...)))
	rotateMat <- t(solve(rotated$Th))
	LIST <- seStdLoadings(rotateMat, object)
	orthogonal <- rotated$orthogonal
	loading <- rotated$loadings
	rotate <- rotated$Th
	gradRotate <- rotated$Gq
	convergence <- rotated$convergence
	method <- rotated$method
	phi <- rotated$Phi
	if(!is.null(phi)) phi <- diag(ncol(loading))
	lv.names <- colnames(loading)
	dimnames(phi) <- list(lv.names, lv.names)
	new("EFA", loading=loading, rotate=rotate, gradRotate=gradRotate, convergence=convergence, phi=phi, se=LIST, method=method, call=mc) 
}

rotateStdLoadings <- function(est, object, rotate=NULL) {
	ov.names <- lavaan:::vnames(object@ParTable, "ov", group = 1)
    lv.names <- lavaan:::vnames(object@ParTable, "lv", group = 1)
	OV <- sqrt(lavaan:::computeVY(object@Model, samplestats = object@SampleStats)[[1]])
	partable <- object@ParTable
	idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names))
	loading <- (solve(diag(OV)) %*% matrix(est[idx], ncol=length(lv.names))) %*% rotate
	est[idx] <- as.vector(loading)
	# Cannot find standard error for Phi because the values are fixed	
	rv.idx <- which(partable$op == "~~" & partable$rhs %in% 
		lv.names)
	# phi <- diag(length(lv.names))
	templhs <- match(partable$lhs[rv.idx], lv.names)
	temprhs <- match(partable$rhs[rv.idx], lv.names)
	# tempval <- est[rv.idx]
	# for(i in seq_along(tempval)) {
		# phi[templhs[i], temprhs[i]] <- tempval[i]
	# }
	rotate2 <- t(solve(rotate))
	phi <- t(rotate2) %*% rotate2
	for(i in seq_along(templhs)) {
		est[rv.idx[i]] <- phi[templhs[i], temprhs[i]]
	}
	est
}

seStdLoadings <- function(rotate, object) {
	est <- object@Fit@est
	JAC <- lavaan:::lavJacobianD(func=rotateStdLoadings, x=object@Fit@est, object=object, rotate=rotate)
	LIST <- inspect(object, "list")
	unco.idx <- which(LIST$unco > 0L)
	LIST <- LIST[,c("lhs", "op", "rhs", "group")]
	JAC <- JAC[unco.idx,unco.idx]
	VCOV <- as.matrix(vcov(object, labels=FALSE))
	if(object@Model@eq.constraints) {
		JAC <- JAC %*% object@Model@eq.constraints.K
	}
	COV <- JAC %*% VCOV %*% t(JAC)
	LIST$se <- rep(NA, length(LIST$lhs))
	LIST$se[unco.idx] <- sqrt(diag(COV))
	tmp.se <- ifelse( LIST$se == 0.0, NA, LIST$se)
    lv.names <- lavaan:::vnames(object@ParTable, "lv", group = 1)
	partable <- object@ParTable
	idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names))
	matrix(LIST$se[idx], ncol=length(lv.names))
}
