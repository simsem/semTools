### Sunthud Pornprasertmanit, Alexander M. Schoemann, Kristopher J. Preacher, Donna Coffman
### Last updated: 22 April 2022


##' Plot power curves for RMSEA
##'
##' Plots power of RMSEA over a range of sample sizes
##'
##' This function creates plot of power for RMSEA against a range of sample
##' sizes. The plot places sample size on the horizontal axis and power on the
##' vertical axis. The user should indicate the lower and upper values for
##' sample size and the sample size between each estimate ("step size") We
##' strongly urge the user to read the sources below (see References) before
##' proceeding.  A web version of this function is available at:
##' \url{http://quantpsy.org/rmsea/rmseaplot.htm}. This function is also
##' implemented in the web application "power4SEM":
##' \url{https://sjak.shinyapps.io/power4SEM/}
##'
##'
##' @importFrom stats qchisq pchisq
##'
##' @param rmsea0 Null RMSEA
##' @param rmseaA Alternative RMSEA
##' @param df Model degrees of freedom
##' @param nlow Lower sample size
##' @param nhigh Upper sample size
##' @param steps Increase in sample size for each iteration. Smaller values of
##' steps will lead to more precise plots. However, smaller step sizes means a
##' longer run time.
##' @param alpha Alpha level used in power calculations
##' @param group The number of group that is used to calculate RMSEA.
##' @param \dots The additional arguments for the plot function.
##'
##' @return Plot of power for RMSEA against a range of sample sizes
##'
##' @author
##' Alexander M. Schoemann (East Carolina University; \email{schoemanna@@ecu.edu})
##'
##' Kristopher J. Preacher (Vanderbilt University; \email{kris.preacher@@vanderbilt.edu})
##'
##' Donna L. Coffman (Pennsylvania State University; \email{dlc30@@psu.edu})
##'
##' @seealso \itemize{
##' \item \code{\link{plotRMSEAdist}} to visualize the RMSEA distributions
##' \item \code{\link{findRMSEApower}} to find the statistical power based on
##'   population RMSEA given a sample size
##' \item \code{\link{findRMSEAsamplesize}} to find the minium sample size for
##'   a given statistical power based on population RMSEA
##' }
##'
##' @references
##' MacCallum, R. C., Browne, M. W., & Cai, L. (2006). Testing
##' differences between nested covariance structure models: Power analysis and
##' null hypotheses. \emph{Psychological Methods, 11}(1), 19--35.
##' \doi{10.1037/1082-989X.11.1.19}
##'
##' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis
##' and determination of sample size for covariance structure modeling.
##' \emph{Psychological Methods, 1}(2), 130--149. \doi{10.1037/1082-989X.1.2.130}
##'
##' MacCallum, R. C., Lee, T., & Browne, M. W. (2010). The issue of isopower in
##' power analysis for tests of structural equation models. \emph{Structural
##' Equation Modeling, 17}(1), 23--41. \doi{10.1080/10705510903438906}
##'
##' Preacher, K. J., Cai, L., & MacCallum, R. C. (2007). Alternatives to
##' traditional model comparison strategies for covariance structure models. In
##' T. D. Little, J. A. Bovaird, & N. A. Card (Eds.), \emph{Modeling contextual
##' effects in longitudinal studies} (pp. 33--62). Mahwah, NJ: Lawrence Erlbaum
##' Associates.
##'
##' Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA fit
##' index. \emph{Structural Equation Modeling, 5}(4), 411--419.
##' \doi{10.1080/10705519809540115}
##'
##' Steiger, J. H., & Lind, J. C. (1980, June). \emph{Statistically based tests
##' for the number of factors.} Paper presented at the annual meeting of the
##' Psychometric Society, Iowa City, IA.
##'
##' Jak, S., Jorgensen, T. D., Verdam, M. G., Oort, F. J., & Elffers, L.
##' (2021). Analytical power calculations for structural equation modeling:
##' A tutorial and Shiny app. \emph{Behavior Research Methods, 53}, 1385--1406.
##' \doi{10.3758/s13428-020-01479-0}
##'
##' @examples
##'
##' plotRMSEApower(rmsea0 = .025, rmseaA = .075, df = 23,
##'                nlow = 100, nhigh = 500, steps = 10)
##'
##' @export
plotRMSEApower <- function(rmsea0, rmseaA, df, nlow, nhigh, steps = 1,
                           alpha = .05, group = 1, ...) {
	pow1 <- 0
	nseq <- seq(nlow,nhigh, by=steps)
	for(i in nseq){
		ncp0 <- ((i-1)*df*rmsea0^2)/group
		ncpa <- ((i-1)*df*rmseaA^2)/group
		#Compute power
		if(rmsea0 < rmseaA) {
		cval <- qchisq(alpha,df,ncp=ncp0,lower.tail=FALSE)
		pow <- pchisq(cval,df,ncp=ncpa,lower.tail=FALSE)
		}
		if(rmsea0 > rmseaA) {
		cval <- qchisq(1-alpha, df, ncp=ncp0, lower.tail=FALSE)
		pow <- 1-pchisq(cval,df,ncp=ncpa,lower.tail=FALSE)
		}
		pow1<-c(pow1, pow)
	}
	pow1 <- pow1[-1]

	plot(nseq,pow1,xlab="Sample Size",ylab="Power",main="Compute Power for RMSEA",type="l", ...)
}



##' Plot the sampling distributions of RMSEA
##'
##' Plots the sampling distributions of RMSEA based on the noncentral chi-square
##' distributions
##'
##' This function creates overlappling plots of the sampling distribution of
##' RMSEA based on noncentral \eqn{\chi^2} distribution (MacCallum, Browne, &
##' Suguwara, 1996). First, the noncentrality parameter (\eqn{\lambda}) is
##' calculated from RMSEA (Steiger, 1998; Dudgeon, 2004) by \deqn{\lambda = (N -
##' 1)d\varepsilon^2 / K,} where \eqn{N} is sample size, \eqn{d} is the model
##' degree of freedom, \eqn{K} is the number of group, and \eqn{\varepsilon} is
##' the population RMSEA. Next, the noncentral \eqn{\chi^2} distribution with a
##' specified \emph{df} and noncentrality parameter is plotted. Thus,
##' the x-axis represents the sample \eqn{\chi^2} value. The sample \eqn{\chi^2}
##' value can be transformed to the sample RMSEA scale (\eqn{\hat{\varepsilon}})
##' by \deqn{\hat{\varepsilon} = \sqrt{K}\sqrt{\frac{\chi^2 - d}{(N - 1)d}},}
##' where \eqn{\chi^2} is the \eqn{\chi^2} value obtained from the noncentral
##' \eqn{\chi^2} distribution.
##'
##'
##' @importFrom stats qchisq
##'
##' @param rmsea The vector of RMSEA values to be plotted
##' @param n Sample size of a dataset
##' @param df Model degrees of freedom
##' @param ptile The percentile rank of the distribution of the first RMSEA that
##' users wish to plot a vertical line in the resulting graph
##' @param caption The name vector of each element of \code{rmsea}
##' @param rmseaScale If \code{TRUE}, the RMSEA scale is used in the x-axis. If
##' \code{FALSE}, the chi-square scale is used in the x-axis.
##' @param group The number of group that is used to calculate RMSEA.
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' @seealso \itemize{
##'  \item \code{\link{plotRMSEApower}} to plot the statistical power
##'    based on population RMSEA given the sample size
##'  \item \code{\link{findRMSEApower}} to find the statistical power based on
##'    population RMSEA given a sample size
##'  \item \code{\link{findRMSEAsamplesize}} to find the minium sample size for
##'   a given statistical power based on population RMSEA
##' }
##'
##' @references
##' Dudgeon, P. (2004). A note on extending Steiger's (1998)
##' multiple sample RMSEA adjustment to other noncentrality parameter-based
##' statistic. \emph{Structural Equation Modeling, 11}(3), 305--319.
##' \doi{10.1207/s15328007sem1103_1}
##'
##' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis
##' and determination of sample size for covariance structure modeling.
##' \emph{Psychological Methods, 1}(2), 130--149. \doi{10.1037/1082-989X.1.2.130}
##'
##' Steiger, J. H. (1998). A note on multiple sample extensions of the RMSEA fit
##' index. \emph{Structural Equation Modeling, 5}(4), 411--419.
##' \doi{10.1080/10705519809540115}
##'
##' @examples
##'
##' plotRMSEAdist(c(.05, .08), n = 200, df = 20, ptile = .95, rmseaScale = TRUE)
##' plotRMSEAdist(c(.05, .01), n = 200, df = 20, ptile = .05, rmseaScale = FALSE)
##'
##' @export
plotRMSEAdist <- function(rmsea, n, df, ptile = NULL, caption = NULL,
                          rmseaScale = TRUE, group = 1) {
	graph <- cbind(rmsea, df)
	ncp <- apply(graph, MARGIN = 1,
	             FUN = function(x, n, group) ((n - 1) * x[2] * (x[1]^2))/group,
	             n = n, group = group)
	graph <- cbind(graph, ncp)
	dens <- lapply(as.list(data.frame(t(graph))), function(x) findDensity("chisq", df = x[2], ncp=x[3]))
	if (rmseaScale) dens <- lapply(dens, function(x, df, n, group) { x[,1] <- (x[,1] - df)/(n-1); x[(x[,1] < 0),1] <- 0; x[,1] <- sqrt(group) * sqrt(x[,1]/df); return(x) }, df=df, n=n, group=group)
	cutoff <- NULL
	if (!is.null(ptile)) {
		cutoff <- qchisq(ptile, df = graph[1, 2], ncp = graph[1, 3])
		if (rmseaScale) cutoff <- sqrt(group) * sqrt((cutoff - df)/(df * (n - 1)))
	}
	if (is.null(caption)) caption <- sapply(graph[,1],
	                                        function(x) paste("Population RMSEA = ",
	                                                          format(x, digits = 3),
	                                                          sep = ""))
	plotOverlapDensity(dens, cutoff, caption, ylab = "Density",
	                   xlab = ifelse(rmseaScale, "RMSEA", "Chi-Square"))
	equal0 <- sapply(dens, function(x) x[,1] == 0)
	if (any(equal0)) warning("The density at RMSEA = 0 cannot be trusted",
	                         " because the plots are truncated.")
}



##' Find the statistical power based on population RMSEA
##'
##' Find the proportion of the samples from the sampling distribution of RMSEA
##' in the alternative hypothesis rejected by the cutoff dervied from the
##' sampling distribution of RMSEA in the null hypothesis. This function can be
##' applied for both test of close fit and test of not-close fit (MacCallum,
##' Browne, & Suguwara, 1996)
##'
##' This function find the proportion of sampling distribution derived from the
##' alternative RMSEA that is in the critical region derived from the sampling
##' distribution of the null RMSEA. If \code{rmseaA} is greater than
##' \code{rmsea0}, the test of close fit is used and the critical region is in
##' the right hand side of the null sampling distribution. On the other hand, if
##' \code{rmseaA} is less than \code{rmsea0}, the test of not-close fit is used
##' and the critical region is in the left hand side of the null sampling
##' distribution (MacCallum, Browne, & Suguwara, 1996).
##'
##' There is also a Shiny app called "power4SEM" that provides a graphical user
##' interface for this functionality (Jak et al., in press).  It can be accessed
##' at \url{https://sjak.shinyapps.io/power4SEM/}.
##'
##'
##' @importFrom stats qchisq pchisq
##'
##' @param rmsea0 Null RMSEA
##' @param rmseaA Alternative RMSEA
##' @param df Model degrees of freedom
##' @param n Sample size of a dataset
##' @param alpha Alpha level used in power calculations
##' @param group The number of group that is used to calculate RMSEA.
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' @seealso \itemize{
##'  \item \code{\link{plotRMSEApower}} to plot the statistical power based on
##'   population RMSEA given the sample size
##'  \item \code{\link{plotRMSEAdist}} to visualize the RMSEA distributions
##'  \item \code{\link{findRMSEAsamplesize}} to find the minium sample size for
##'   a given statistical power based on population RMSEA
##' }
##'
##' @references
##' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis
##' and determination of sample size for covariance structure modeling.
##' \emph{Psychological Methods, 1}(2), 130--149. \doi{10.1037/1082-989X.1.2.130}
##'
##'  Jak, S., Jorgensen, T. D., Verdam, M. G., Oort, F. J., & Elffers, L.
##'  (2021). Analytical power calculations for structural equation modeling:
##'  A tutorial and Shiny app. \emph{Behavior Research Methods, 53}, 1385--1406.
##'  \doi{10.3758/s13428-020-01479-0}
##'
##' @examples
##'
##' findRMSEApower(rmsea0 = .05, rmseaA = .08, df = 20, n = 200)
##'
##' @export
findRMSEApower <- function(rmsea0, rmseaA, df, n, alpha = .05, group = 1) {
	ncp0 <- ((n-1)*df*rmsea0^2)/group
	ncpa <- ((n-1)*df*rmseaA^2)/group
	if (rmsea0<rmseaA) {
		cval <- qchisq(alpha,df,ncp=ncp0,lower.tail=FALSE)
		pow <- pchisq(cval,df,ncp=ncpa,lower.tail=FALSE)
	} else {
		cval <- qchisq(1-alpha,df,ncp=ncp0,lower.tail=FALSE)
		pow <- 1-pchisq(cval,df,ncp=ncpa,lower.tail=FALSE)
	}
	return(pow)
}



##' Find the minimum sample size for a given statistical power based on
##' population RMSEA
##'
##' Find the minimum sample size for a specified statistical power based on
##' population RMSEA. This function can be applied for both test of close fit
##' and test of not-close fit (MacCallum, Browne, & Suguwara, 1996)
##'
##' This function find the minimum sample size for a specified power based on an
##' iterative routine. The sample size keep increasing until the calculated
##' power from \code{\link{findRMSEApower}} function is just over the specified
##' power. If \code{group} is greater than 1, the resulting sample size is the
##' sample size per group.
##'
##' @param rmsea0 Null RMSEA
##' @param rmseaA Alternative RMSEA
##' @param df Model degrees of freedom
##' @param power Desired statistical power to reject misspecified model (test of
##' close fit) or retain good model (test of not-close fit)
##' @param alpha Alpha level used in power calculations
##' @param group The number of group that is used to calculate RMSEA.
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' @seealso \itemize{
##'  \item \code{\link{plotRMSEApower}} to plot the statistical power based on
##'   population RMSEA given the sample size
##'  \item \code{\link{plotRMSEAdist}} to visualize the RMSEA distributions
##'  \item \code{\link{findRMSEApower}} to find the statistical power based on
##'   population RMSEA given a sample size
##' }
##'
##' @references
##' MacCallum, R. C., Browne, M. W., & Sugawara, H. M. (1996). Power analysis
##' and determination of sample size for covariance structure modeling.
##' \emph{Psychological Methods, 1}(2), 130--149. \doi{10.1037/1082-989X.1.2.130}
##'
##'  Jak, S., Jorgensen, T. D., Verdam, M. G., Oort, F. J., & Elffers, L.
##'  (2021). Analytical power calculations for structural equation modeling:
##'  A tutorial and Shiny app. \emph{Behavior Research Methods, 53}, 1385--1406.
##'  \doi{10.3758/s13428-020-01479-0}
##'
##' @examples
##'
##' findRMSEAsamplesize(rmsea0 = .05, rmseaA = .08, df = 20, power = 0.80)
##'
##' @export
findRMSEAsamplesize <- function(rmsea0, rmseaA, df, power = 0.80,
                                alpha = .05, group = 1) {
	n <- 5:500
	pow <- findRMSEApower(rmsea0, rmseaA, df, n, alpha, group=group)
	if(all(pow > power)) {
		return("Sample Size <= 5")
	} else if (all(power > pow)) {
		repeat {
			n <- n + 500
			pow <- findRMSEApower(rmsea0, rmseaA, df, n, alpha, group=group)
			if(any(pow > power)) {
				index <- which(pow > power)[1]
				return(n[index]/group)
			}
		}
	} else {
		index <- which(pow > power)[1]
		return(n[index]/group)
	}
}



## ----------------
## Hidden Functions
## ----------------

## findDensity
## Find the x and y coordinate of a distribution in order to plot a density of a distribution
## dist: target distribution in text, such as "chisq"
## ...: Additional argument of the distribution
## Return the data frame with x and y coordinates for plotting density
findDensity <- function(dist, ...) {
  FUN <- list()
  FUN[[1]] <- get(paste("q", dist, sep=""))
  FUN[[2]] <- c(0.001, 0.999)
  FUN <- c(FUN, ...)
  bound <- eval(as.call(FUN))
  val <- seq(bound[1], bound[2], length.out=1000)
  FUN[[1]] <- get(paste("d", dist, sep=""))
  FUN[[2]] <- val
  height <- eval(as.call(FUN))
  return(cbind(val, height))
}
##Example Code
##findDensity("chisq", df=10)

## plotOverlapDensity
## Plot the overlapping distributions using density
## dat: A list of data frame where each data frame has the x coordinate as the variable 1 and y coordinate as the variable 2
## vline: Vertical line in the graph
## caption: The name of each density line
## ...: Additional argument of the plot function
plotOverlapDensity <- function(dat, vline = NULL, caption = NULL, ...) {
  if (!is.list(dat)) {
    temp <- list()
    temp[[1]] <- dat
    dat <- temp
  }
  stack <- do.call(rbind, dat)
  lim <- apply(stack, 2, function(x) c(min(x), max(x)))
  plot(stack, xlim = lim[,1], ylim = lim[,2], type = "n", ...)
  for (i in 1:length(dat)) lines(dat[[i]], col = i, lwd = 1.5)
  for (i in 1:length(vline)) abline(v = vline[i], lwd = 1.5)
  if (!is.null(caption))
    legend(0.50 * (lim[2,1] - lim[1,1]) + lim[1,1], 0.99 * (lim[2,2] - lim[1,2]) + lim[1,2], caption, col=1:length(dat), lty=1)
}



