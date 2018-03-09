### Sunthud Pornprasertmanit, Bell Clinton, Pavel Panko
### Last updated: 9 March 2018


#' Find power given a sample size in nested model comparison
#'
#' Find the sample size that the power in rejection the samples from the
#' alternative pair of RMSEA is just over the specified power.
#'
#'
#' @importFrom stats qchisq pchisq
#'
#' @param rmsea0A The \eqn{H_0} baseline RMSEA
#' @param rmsea0B The \eqn{H_0} alternative RMSEA (trivial misfit)
#' @param rmsea1A The \eqn{H_1} baseline RMSEA
#' @param rmsea1B The \eqn{H_1} alternative RMSEA (target misfit to be rejected)
#' @param dfA degree of freedom of the more-restricted model
#' @param dfB degree of freedom of the less-restricted model
#' @param n Sample size
#' @param alpha The alpha level
#' @param group The number of group in calculating RMSEA
#' @author Bell Clinton
#'
#' Pavel Panko (Texas Tech University; \email{pavel.panko@@ttu.edu})
#'
#' Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#' @seealso \itemize{
#'  \item \code{\link{plotRMSEApowernested}} to plot the statistical power for
#'   nested model comparison based on population RMSEA given the sample size
#'  \item \code{\link{findRMSEAsamplesizenested}} to find the minium sample
#'   size for a given statistical power in nested model comparison based on
#'   population RMSEA
#' }
#' @references
#' MacCallum, R. C., Browne, M. W., & Cai, L. (2006). Testing
#' differences between nested covariance structure models: Power analysis and
#' null hypotheses. \emph{Psychological Methods, 11}(1), 19--35.
#' doi:10.1037/1082-989X.11.1.19
#' @examples
#'
#' findRMSEApowernested(rmsea0A = 0.06, rmsea0B = 0.05, rmsea1A = 0.08,
#'                      rmsea1B = 0.05, dfA = 22, dfB = 20, n = 200,
#'                      alpha = 0.05, group = 1)
#'
#' @export
findRMSEApowernested <- function(rmsea0A = NULL, rmsea0B = NULL, rmsea1A, rmsea1B = NULL, dfA, dfB, n, alpha = 0.05, group = 1) {
  if(is.null(rmsea0A)) rmsea0A <- 0
  if(is.null(rmsea0B)) rmsea0B <- 0
  if(is.null(rmsea1B)) rmsea1B <- rmsea0B
  if(dfA <= dfB) stop("The degree of freedom of the more-restricted model (dfA) should be greater than the degree of freedom of the less-restricted model (dfB)")
  if(rmsea0A < rmsea0B) stop("In the null-hypothesis models, the RMSEA of the more-restricted model (rmsea0A) should have a higher than or equal to the RMSEA of the less-restricted model (rmsea0B).")
  if(rmsea1A < rmsea1B) stop("In the alternative-hypothesis models, the RMSEA of the more-restricted model (rmsea1A) should have a higher than or equal to the RMSEA of the less-restricted model (rmsea1B).")
  ddiff <- dfA-dfB
  f0a <- (dfA*rmsea0A^2)/group
  f0b <- (dfB*rmsea0B^2)/group
  f1a <- (dfA*rmsea1A^2)/group
  f1b <- (dfB*rmsea1B^2)/group
  ncp0 <- (n-1)*(f0a-f0b)
  ncp1 <- (n-1)*(f1a-f1b)
  cval <- qchisq(1-alpha,ddiff,ncp0)
  Power <- 1-pchisq(cval,ddiff,ncp1)
  Power
}



#' Find sample size given a power in nested model comparison
#'
#' Find the sample size that the power in rejection the samples from the
#' alternative pair of RMSEA is just over the specified power.
#'
#'
#' @param rmsea0A The \eqn{H_0} baseline RMSEA
#' @param rmsea0B The \eqn{H_0} alternative RMSEA (trivial misfit)
#' @param rmsea1A The \eqn{H_1} baseline RMSEA
#' @param rmsea1B The \eqn{H_1} alternative RMSEA (target misfit to be rejected)
#' @param dfA degree of freedom of the more-restricted model.
#' @param dfB degree of freedom of the less-restricted model.
#' @param power The desired statistical power.
#' @param alpha The alpha level.
#' @param group The number of group in calculating RMSEA.
#' @author Bell Clinton
#'
#' Pavel Panko (Texas Tech University; \email{pavel.panko@@ttu.edu})
#'
#' Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#'
#' @seealso \itemize{
#'  \item \code{\link{plotRMSEApowernested}} to plot the statistical power for
#'   nested model comparison based on population RMSEA given the sample size
#'  \item \code{\link{findRMSEApowernested}} to find the power for a given
#'   sample size in nested model comparison based on population RMSEA
#' }
#' @references
#' MacCallum, R. C., Browne, M. W., & Cai, L. (2006). Testing
#' differences between nested covariance structure models: Power analysis and
#' null hypotheses. \emph{Psychological Methods, 11}(1), 19-35.
#' doi:10.1037/1082-989X.11.1.19
#' @examples
#'
#' findRMSEAsamplesizenested(rmsea0A = 0, rmsea0B = 0, rmsea1A = 0.06,
#'                           rmsea1B = 0.05, dfA = 22, dfB = 20, power = 0.80,
#'                           alpha = .05, group = 1)
#'
#' @export
findRMSEAsamplesizenested <- function(rmsea0A = NULL, rmsea0B = NULL, rmsea1A,
                                      rmsea1B = NULL, dfA, dfB, power = 0.80,
                                      alpha = .05, group = 1) {
	if(is.null(rmsea0A)) rmsea0A <- 0
	if(is.null(rmsea0B)) rmsea0B <- 0
	if(is.null(rmsea1B)) rmsea1B <- rmsea0B
	if(dfA <= dfB) stop("The degree of freedom of the more-restricted model (dfA) should be greater than the degree of freedom of the less-restricted model (dfB)")
	if(rmsea0A < rmsea0B) stop("In the null-hypothesis models, the RMSEA of the more-restricted model (rmsea0A) should have a higher than or equal to the RMSEA of the less-restricted model (rmsea0B).")
	if(rmsea1A < rmsea1B) stop("In the alternative-hypothesis models, the RMSEA of the more-restricted model (rmsea1A) should have a higher than or equal to the RMSEA of the less-restricted model (rmsea1B).")
	n <- 5:500
	pow <- findRMSEApowernested(rmsea0A, rmsea0B, rmsea1A, rmsea1B, dfA, dfB, n, alpha, group = group)
	if(all(pow > power)) {
		return("Sample Size <= 5")
	} else if (all(power > pow)) {
		repeat {
			n <- n + 500
			pow <- findRMSEApowernested(rmsea0A, rmsea0B, rmsea1A, rmsea1B, dfA, dfB, n, alpha, group = group)
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



#' Plot power of nested model RMSEA
#'
#' Plot power of nested model RMSEA over a range of possible sample sizes.
#'
#'
#' @param rmsea0A The \eqn{H_0} baseline RMSEA
#' @param rmsea0B The \eqn{H_0} alternative RMSEA (trivial misfit)
#' @param rmsea1A The \eqn{H_1} baseline RMSEA
#' @param rmsea1B The \eqn{H_1} alternative RMSEA (target misfit to be rejected)
#' @param dfA degree of freedom of the more-restricted model
#' @param dfB degree of freedom of the less-restricted model
#' @param nlow Lower bound of sample size
#' @param nhigh Upper bound of sample size
#' @param steps Step size
#' @param alpha The alpha level
#' @param group The number of group in calculating RMSEA
#' @param \dots The additional arguments for the plot function.
#' @author Bell Clinton
#'
#' Pavel Panko (Texas Tech University; \email{pavel.panko@@ttu.edu})
#'
#' Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#'
#' @seealso \itemize{
#'  \item \code{\link{findRMSEApowernested}} to find the power for a given
#'   sample size in nested model comparison based on population RMSEA
#'  \item \code{\link{findRMSEAsamplesizenested}} to find the minium sample
#'   size for a given statistical power in nested model comparison based on
#'   population RMSEA
#' }
#' @references
#' MacCallum, R. C., Browne, M. W., & Cai, L. (2006). Testing
#' differences between nested covariance structure models: Power analysis and
#' null hypotheses. \emph{Psychological Methods, 11}(1), 19-35.
#' doi:10.1037/1082-989X.11.1.19
#' @examples
#'
#' plotRMSEApowernested(rmsea0A = 0, rmsea0B = 0, rmsea1A = 0.06,
#'                      rmsea1B = 0.05, dfA = 22, dfB = 20, nlow = 50,
#'                      nhigh = 500, steps = 1, alpha = .05, group = 1)
#'
#' @export
plotRMSEApowernested <- function(rmsea0A = NULL, rmsea0B = NULL, rmsea1A,
                                 rmsea1B = NULL, dfA, dfB, nlow, nhigh,
                                 steps = 1, alpha = .05, group = 1, ...){
	nseq <- seq(nlow,nhigh, by=steps)
	pow1 <- findRMSEApowernested(rmsea0A = rmsea0A, rmsea0B = rmsea0B, rmsea1A = rmsea1A, rmsea1B = rmsea1B, dfA = dfA, dfB = dfB, n = nseq, alpha = alpha, group = group)
	plot(nseq, pow1, xlab="Sample Size", ylab="Power", main="Compute Power for Nested RMSEA", type="l", ...)
}


