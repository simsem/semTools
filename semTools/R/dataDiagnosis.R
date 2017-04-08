### Sunthud Pornprasertmanit
### Last updated: 2 April 2017
### Higher-order moments. Initial version from the simsem package.



#' Finding skewness
#'
#' Finding skewness (\eqn{g_{1}}) of an object
#'
#' The skewness computed is \eqn{g_{1}}. The parameter skewness \eqn{\gamma_{2}}
#' formula is
#'
#' \deqn{\gamma_{2} = \frac{\mu_{3}}{\mu^{3/2}_{2}},}
#'
#' where \eqn{\mu_{i}} denotes the \eqn{i} order central moment.
#'
#' The excessive kurtosis formula for sample statistic \eqn{g_{2}} is
#'
#' \deqn{g_{2} = \frac{k_{3}}{k^{2}_{2}},}
#'
#' where \eqn{k_{i}} are the \eqn{i} order \emph{k}-statistic.
#'
#' The standard error of the skewness is
#'
#' \deqn{Var(\hat{g}_2) = \frac{6}{N}}
#'
#' where \eqn{N} is the sample size.
#'
#'
#' @importFrom stats pnorm
#' 
#' @param object A vector used to find a skewness
#' @param population \code{TRUE} to compute the parameter formula. \code{FALSE}
#' to compute the sample statistic formula.
#' @return A value of a skewness with a test statistic if the population is
#' specified as \code{FALSE}
#' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#' @seealso \itemize{
#'   \item \code{\link{kurtosis}} Find the univariate excessive kurtosis
#'    of a variable
#'   \item \code{\link{mardiaSkew}} Find Mardia's multivariate skewness
#'    of a set of variables
#'   \item \code{\link{mardiaKurtosis}} Find the Mardia's multivariate
#'    kurtosis of a set of variables
#'  }
#' @references Weisstein, Eric W. (n.d.). \emph{Skewness}. Retrived from
#'  \emph{MathWorld}--A Wolfram Web Resource:
#'  \url{http://mathworld.wolfram.com/Skewness.html}
#' @examples
#'
#' skew(1:5)
#'
#' @export
skew <- function(object, population = FALSE) {
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



#' Finding excessive kurtosis
#'
#' Finding excessive kurtosis (\eqn{g_{2}}) of an object
#'
#' The excessive kurtosis computed is \eqn{g_{2}}. The parameter excessive
#' kurtosis \eqn{\gamma_{2}} formula is
#'
#' \deqn{\gamma_{2} = \frac{\mu_{4}}{\mu^{2}_{2}} - 3,}
#'
#' where \eqn{\mu_{i}} denotes the \eqn{i} order central moment.
#'
#' The excessive kurtosis formula for sample statistic \eqn{g_{2}} is
#'
#' \deqn{g_{2} = \frac{k_{4}}{k^{2}_{2}},}
#'
#' where \eqn{k_{i}} are the \eqn{i} order \emph{k}-statistic.
#'
#' The standard error of the excessive kurtosis is
#'
#' \deqn{Var(\hat{g}_{2}) = \frac{24}{N}}
#'
#' where \eqn{N} is the sample size.
#'
#'
#' @importFrom stats pnorm
#' 
#' @param object A vector used to find a excessive kurtosis
#' @param population \code{TRUE} to compute the parameter formula. \code{FALSE}
#' to compute the sample statistic formula.
#' @return A value of an excessive kurtosis with a test statistic if the
#' population is specified as \code{FALSE}
#' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#' @seealso \itemize{
#'   \item \code{\link{skew}} Find the univariate skewness of a variable
#'   \item \code{\link{mardiaSkew}} Find the Mardia's multivariate
#'    skewness of a set of variables
#'   \item \code{\link{mardiaKurtosis}} Find the Mardia's multivariate kurtosis
#'    of a set of variables
#' }
#' @references Weisstein, Eric W. (n.d.). \emph{Kurtosis.} Retrived from
#' \emph{MathWorld}--A Wolfram Web Resource:
#' \url{http://mathworld.wolfram.com/Kurtosis.html}
#' @examples
#'
#' kurtosis(1:5)
#'
#' @export
kurtosis <- function(object, population = FALSE) {
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



#' Finding Mardia's multivariate skewness
#'
#' Finding Mardia's multivariate skewness of multiple variables
#'
#' The Mardia's multivariate skewness formula (Mardia, 1970) is
#'  \deqn{ b_{1, d} = \frac{1}{n^2}\sum^n_{i=1}\sum^n_{j=1}\left[
#'  \left(\bold{X}_i - \bold{\bar{X}} \right)^{'} \bold{S}^{-1}
#'  \left(\bold{X}_j - \bold{\bar{X}} \right) \right]^3, }
#' where \eqn{d} is the number of variables, \eqn{X} is the target dataset
#' with multiple variables, \eqn{n} is the sample size, \eqn{\bold{S}} is
#' the sample covariance matrix of the target dataset, and \eqn{\bold{\bar{X}}}
#' is the mean vectors of the target dataset binded in \eqn{n} rows.
#' When the population multivariate skewness is normal, the
#' \eqn{\frac{n}{6}b_{1,d}} is asymptotically distributed as \eqn{\chi^2}
#' distribution with \eqn{d(d + 1)(d + 2)/6} degrees of freedom.
#'
#'
#' @importFrom stats cov pchisq
#' 
#' @param dat The target matrix or data frame with multiple variables
#' @param use Missing data handling method from the \code{\link[stats]{cov}}
#' function.
#' @return A value of a Mardia's multivariate skewness with a test statistic
#' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#' @seealso \itemize{
#'   \item \code{\link{skew}} Find the univariate skewness of a variable
#'   \item \code{\link{kurtosis}} Find the univariate excessive
#'     kurtosis of a variable
#'   \item \code{\link{mardiaKurtosis}} Find the Mardia's multivariate
#'     kurtosis of a set of variables
#' }
#' @references Mardia, K. V. (1970). Measures of multivariate skewness and
#'   kurtosis with applications. \emph{Biometrika, 57}(3), 519-530.
#'   doi:10.2307/2334770
#' @examples
#'
#' library(lavaan)
#' mardiaSkew(HolzingerSwineford1939[ , paste0("x", 1:9)])
#'
#' @export
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



#' Finding Mardia's multivariate kurtosis
#'
#' Finding Mardia's multivariate kurtosis of multiple variables
#'
#' The Mardia's multivariate kurtosis formula (Mardia, 1970) is
#'  \deqn{ b_{2, d} = \frac{1}{n}\sum^n_{i=1}\left[ \left(\bold{X}_i -
#'  \bold{\bar{X}} \right)^{'} \bold{S}^{-1} \left(\bold{X}_i -
#'  \bold{\bar{X}} \right) \right]^2, }
#' where \eqn{d} is the number of variables, \eqn{X} is the target
#' dataset with multiple variables, \eqn{n} is the sample size, \eqn{\bold{S}}
#' is the sample covariance matrix of the target dataset, and
#' \eqn{\bold{\bar{X}}} is the mean vectors of the target dataset binded in
#' \eqn{n} rows. When the population multivariate kurtosis is normal, the
#' \eqn{b_{2,d}} is asymptotically distributed as normal distribution with the
#' mean of \eqn{d(d + 2)} and variance of \eqn{8d(d + 2)/n}.
#'
#'
#' @importFrom stats cov pnorm
#' 
#' @param dat The target matrix or data frame with multiple variables
#' @param use Missing data handling method from the \code{\link[stats]{cov}}
#' function.
#' @return A value of a Mardia's multivariate kurtosis with a test statistic
#' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#' @seealso \itemize{
#'  \item \code{\link{skew}} Find the univariate skewness of a variable
#'  \item \code{\link{kurtosis}} Find the univariate excessive kurtosis
#'   of a variable
#'  \item \code{\link{mardiaSkew}} Find the Mardia's multivariate skewness
#'   of a set of variables
#' }
#' @references Mardia, K. V. (1970). Measures of multivariate skewness and
#'  kurtosis with applications. \emph{Biometrika, 57}(3), 519-530.
#'  doi:10.2307/2334770
#' @examples
#'
#' library(lavaan)
#' mardiaKurtosis(HolzingerSwineford1939[ , paste0("x", 1:9)])
#'
#' @export
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



## ----------------
## Hidden Functions
## ----------------

## centralMoment
## Calculate central moments of a variable
## Arguments:
##	 x: vector of a variable
## 	 ord: order of the moment
## 	 weight: weight variable
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
# Arguments:
#	  x: vector of a variable
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

