## Title: Data Diagnosis
## Author: Sunthud Pornprasertmanit
# Description: Diagnose data for its distribution
# Remark: initial version from the simsem package

# centralMoment
# Calculate central moments of a variable
# Argument:
#	x: vector of a variable
# 	ord: order of the moment
# 	weight: weight variable

centralMoment <- function(x, ord) {
	if(ord < 2) stop("Central moment can be calculated for order 2 or more in an integer.")
	wm <- mean(x)
	result <- sum((x - wm)^(ord))/length(x)
	return(result)
}
# Example
# centralMoment(1:5, 2)

# kStat
# Calculate the k-statistic (i.e., unbiased estimator of a cumulant) of a variable
# Argument:
#	x: vector of a variable
# 	ord: order of the k-statistics

kStat <- function(x, ord) {
	# Formula from mathworld wolfram
	n <- length(x)
	if(ord == 1) {
		return(mean(x))
	} else if (ord == 2) {
		return(centralMoment(x, 2) * n / (n - 1))
	} else if (ord == 3) {
		return(centralMoment(x, 3) * n^2 / ((n - 1) * (n - 2)))
	} else if (ord == 4) {
		num1 <- n^2
		num2 <- (n + 1) * centralMoment(x, 4)
		num3 <- 3 * (n - 1) * centralMoment(x, 2)^2
		denom <- (n - 1) * (n - 2) * (n - 3)
		return((num1 * (num2 - num3))/denom)
	} else {
		stop("Order can be 1, 2, 3, or 4 only.")
	}
}
# Example
# kStat(1:5, 4)

# skew
# Calculate the skewness of a vector
# Argument:
#	object:	The target vector
#	population: The vector represents population values or sample values
skew <- function(object, population=FALSE) {
	if(any(is.na(object))) {
		object <- object[!is.na(object)]
		warning("Missing observations are removed from a vector.")
	}
	if(population) {
		return(centralMoment(object, 3)/(centralMoment(object, 2)^(3/2)))
	} else {
		est <- kStat(object, 3)/(kStat(object, 2)^(3/2))
		se <- sqrt(6/length(object))
		z <- est/se
		p <- (1 - pnorm(abs(z)))*2
		return(c("skew (g1)"=est, se=se, z=z, p=p))
	}
}

# kurtosis
# Calculate the (excessive) kurtosis of a vector
# Argument:
#	object:	The target vector
#	population: The vector represents population values or sample values
kurtosis <- function(object, population=FALSE) {
	if(any(is.na(object))) {
		object <- object[!is.na(object)]
		warning("Missing observations are removed from a vector.")
	}	
	if(population) {
		return((centralMoment(object, 4)/(centralMoment(object, 2)^2)) - 3)
	} else {
		est <- kStat(object, 4)/(kStat(object, 2)^(2))
		se <- sqrt(24/length(object))
		z <- est/se
		p <- (1 - pnorm(abs(z)))*2
		return(c("Excess Kur (g2)"=est, se=se, z=z, p=p))
	}
}

# mardiaSkew
# Calculate the Mardia's skewness 
# Argument:
#	dat: Datasets with multiple variables
mardiaSkew <- function(dat, use = "everything") {
	centeredDat <- scale(dat, center=TRUE, scale=FALSE)
	invS <- solve(cov(dat, use = use))
	FUN <- function(vec1, vec2, invS) {
		as.numeric(t(as.matrix(vec1)) %*% invS %*% as.matrix(vec2))
	}
	FUN2 <- function(vec1, listVec2, invS) {
		sapply(listVec2, FUN, vec1=vec1, invS=invS)
	}
	indivTerm <- sapply(as.list(data.frame(t(centeredDat))), FUN2, listVec2=as.list(data.frame(t(centeredDat))), invS=invS)
	b1d <- sum(indivTerm^3)/(nrow(dat)^2)
	d <- ncol(dat)
	chi <- nrow(dat) * b1d / 6
	df <- d * (d + 1) * (d + 2) / 6
	p <- pchisq(chi, df = df, lower.tail = FALSE)
	return(c(b1d = b1d, chi = chi, df=df, p=p))
}

# mardiaKurtosis
# Calculate the Mardia's Kurtosis 
# Argument:
#	dat: Datasets with multiple variables
mardiaKurtosis <- function(dat, use = "everything") {
	centeredDat <- scale(dat, center=TRUE, scale=FALSE)
	invS <- solve(cov(dat, use = use))
	FUN <- function(vec, invS) {
		as.numeric(t(as.matrix(vec)) %*% invS %*% as.matrix(vec))
	}
	indivTerm <- sapply(as.list(data.frame(t(centeredDat))), FUN, invS=invS)
	b2d <- sum(indivTerm^2)/nrow(dat)
	d <- ncol(dat)
	m <- d * (d + 2)
	v <- 8 * d * (d + 2) / nrow(dat)
	z <- (b2d - m)/sqrt(v)
	p <- pnorm(-abs(z)) * 2
	return(c(b2d = b2d, z = z, p=p))
}
