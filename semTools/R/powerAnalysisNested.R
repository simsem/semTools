# Power analysis for nested model comparison

# Note: Model0 = Null hypothesis
#		Model1 = Alternative hypothesis
#		ModelA = More-restricted models (higher df; higher RMSEA)
#		ModelB = Less-restricted models (lower df; lower RMSEA)

# findRMSEApowernested
# Find the proportion of the samples from the alternative pair of RMSEAs in nested model comparison rejected by the cutoff dervied from the null pair of RMSEAs in nested model comparison
# rmsea0A: The H0 baseline RMSEA
# rmsea0B: The H0 alternative RMSEA (trivial misfit)
# rmsea1A: The H1 baseline RMSEA
# rmsea1B: The H1 alternative RMSEA (target misfit to be rejected)
# n: sample size
# dfA: degree of freedom of the more-restricted model
# dfB: degree of freedom of the less-restricted model
# alpha: The alpha level
# group: The number of group in calculating RMSEA
# Return power

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

test.findRMSEApowernested <- function() {
	alpha <- 0.05
	rmsea0A <- 0.06
	rmsea0B <- 0.05
	rmsea1A <- 0.08
	rmsea1B <- 0.05
	dfA <- 22
	dfB <- 20
	n <- 200
	group <- 1
	findRMSEApowernested(rmsea0A, rmsea0B, rmsea1A, rmsea1B, dfA, dfB, n, alpha, group)
}

# findRMSEAsamplesizenested
# Find the sample size that the power in rejection the samples from the alternative pair of RMSEA is just over the specified power
# rmsea0A: The H0 baseline RMSEA
# rmsea0B: The H0 alternative RMSEA (trivial misfit)
# rmsea1A: The H1 baseline RMSEA
# rmsea1B: The H1 alternative RMSEA (target misfit to be rejected)
# dfA: degree of freedom of the more-restricted model
# dfB: degree of freedom of the less-restricted model
# power: The desired statistical power
# alpha: The alpha level
# group: The number of group in calculating RMSEA
# Return The estimated sample size

findRMSEAsamplesizenested <- function(rmsea0A = NULL, rmsea0B = NULL, rmsea1A, rmsea1B = NULL, dfA, dfB, power=0.80, alpha=.05, group=1) {
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

test.findRMSEAsamplesizenested <- function() {
	alpha <- 0.05
	rmseaA <- 0.06
	rmseaB <- 0.05
	da <- 22
	db <- 20
	powd <- 0.8
	G <- 1
	findRMSEAsamplesizenested(rmsea0A = 0, rmsea0B = 0, rmsea1A = rmseaA, rmsea1B = rmseaB, da, db, power=0.80, alpha=.05, group=1) 
}

# plotRMSEApowernested
#Plot power of nested model RMSEA over a range of possible sample sizes
# rmsea0A: The H0 baseline RMSEA
# rmsea0B: The H0 alternative RMSEA (trivial misfit)
# rmsea1A: The H1 baseline RMSEA
# rmsea1B: The H1 alternative RMSEA (target misfit to be rejected)
# dfA: degree of freedom of the more-restricted model
# dfB: degree of freedom of the less-restricted model
# nlow: Lower bound of sample size
# nhigh: Upper bound of sample size
# steps: Step size
# alpha: The alpha level
# group: The number of group in calculating RMSEA
# ...: Additional parameters for graphs
# Return plot of power

plotRMSEApowernested <- function(rmsea0A = NULL, rmsea0B = NULL, rmsea1A, rmsea1B = NULL, dfA, dfB, nlow, nhigh, steps=1, alpha=.05, group=1, ...){ 
	nseq <- seq(nlow,nhigh, by=steps)
	pow1 <- findRMSEApowernested(rmsea0A = rmsea0A, rmsea0B = rmsea0B, rmsea1A = rmsea1A, rmsea1B = rmsea1B, dfA = dfA, dfB = dfB, n = nseq, alpha = alpha, group = group)
	plot(nseq, pow1, xlab="Sample Size", ylab="Power", main="Compute Power for Nested RMSEA", type="l", ...)
}

test.plotRMSEApowernested <- function() {
	alpha <- 0.05
	rmseaA <- 0.06
	rmseaB <- 0.05
	da <- 22
	db <- 20
	plotRMSEApowernested(rmsea0A = 0, rmsea0B = 0, rmsea1A = rmseaA, rmsea1B = rmseaB, da, db, nlow=50, nhigh=500, steps=1, alpha=.05, group=1)  
}
