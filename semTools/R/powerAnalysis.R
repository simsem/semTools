# plotRMSEApower
#Plot power of RMSEA over a range of possible sample sizes
#input: rmsea of null and alternative model, degress of freedom, lower sampel size, upper sample sample, sample size steps, alpha, the number of group in calculating RMSEA
#Output: plot of power
#Alexander M. Schoemann, Kristopher J. Preacher, Donna Coffman
#5/30/2012

plotRMSEApower <- function(rmsea0, rmseaA, df, nlow, nhigh, steps=1, alpha=.05, group=1, ...) { 
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

#Example Code
#plotRMSEApower(.025, .075, 23, 100, 500, 10)

# findDensity
# Find the x and y coordinate of a distribution in order to plot a density of a distribution
# dist: target distribution in text, such as "chisq"
# ...: Additional argument of the distribution
# Return the data frame with x and y coordinates for plotting density
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
#Example Code
#findDensity("chisq", df=10)

# plotOverlapDensity
# Plot the overlapping distributions using density
# dat: A list of data frame where each data frame has the x coordinate as the variable 1 and y coordinate as the variable 2
# vline: Vertical line in the graph
# caption: The name of each density line 
# ...: Additional argument of the plot function
plotOverlapDensity <- function(dat, vline = NULL, caption=NULL, ...) {
	if(!is.list(dat)) {
		temp <- list()
		temp[[1]] <- dat
		dat <- temp
	}
	stack <- do.call(rbind, dat)
	lim <- apply(stack, 2, function(x) c(min(x), max(x)))
	plot(stack, xlim=lim[,1], ylim=lim[,2], type="n", ...)
	for(i in 1:length(dat)) lines(dat[[i]], col = i, lwd=1.5)
	for(i in 1:length(vline)) abline(v = vline[i], lwd=1.5)
	if(!is.null(caption))
	legend(0.50 * (lim[2,1] - lim[1,1]) + lim[1,1], 0.99 * (lim[2,2] - lim[1,2]) + lim[1,2], caption, col=1:length(dat), lty=1)
}

# plotRMSEAdist
# Plot the overlapping distributions of RMSEA based on noncentral chi-square distribution
# rmsea: A vector of RMSEA
# n: sample size
# df: degree of freedom of the chi-square distribution
# ptile: The percentile rank of the first specified rmsea to put the vertical line
# caption: The description of each rmsea
# rmseaScale: If TRUE, use RMSEA as the scale in x-axis. If FALSE, use chi-square as the scale in x-axis.
# group: The number of group in calculating RMSEA
plotRMSEAdist <- function(rmsea, n, df, ptile=NULL, caption=NULL, rmseaScale = TRUE, group=1) {
	graph <- cbind(rmsea, df)
	ncp <- apply(graph, 1, function(x, n, group) ((n - 1) * x[2] * (x[1]^2))/group, n=n, group=group)
	graph <- cbind(graph, ncp)
	dens <- lapply(as.list(data.frame(t(graph))), function(x) findDensity("chisq", df = x[2], ncp=x[3])) 
	if(rmseaScale) dens <- lapply(dens, function(x, df, n, group) { x[,1] <- (x[,1] - df)/(n-1); x[(x[,1] < 0),1] <- 0; x[,1] <- sqrt(group) * sqrt(x[,1]/df); return(x) }, df=df, n=n, group=group)
	cutoff <- NULL
	if(!is.null(ptile)) { 
		cutoff <- qchisq(ptile, df=graph[1, 2], ncp=graph[1, 3])
		if(rmseaScale) cutoff <- sqrt(group) * sqrt((cutoff - df)/(df * (n - 1)))
	}
	if(is.null(caption)) caption <- sapply(graph[,1], function(x) paste("Population RMSEA = ", format(x, digits=3), sep="")) 
	plotOverlapDensity(dens, cutoff, caption, xlab=ifelse(rmseaScale, "RMSEA", "Chi-Square"), ylab="Density")
	equal0 <- sapply(dens, function(x) x[,1] == 0) 
	if(any(equal0)) warning("The density at RMSEA = 0 cannot be trusted because the plots are truncated.")
}

# findRMSEApower
# Find the proportion of the samples from the alternative RMSEA rejected by the cutoff dervied from the null RMSEA
# rmsea0: The null RMSEA
# rmseaA: The alternative RMSEA
# n: sample size
# df: degree of freedom of the chi-square distribution
# alpha: The alpha level
# group: The number of group in calculating RMSEA
# Return power
findRMSEApower <- function(rmsea0, rmseaA, df, n, alpha=.05, group=1) {
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

# findRMSEAsamplesize
# Find the sample size that the power in rejection the samples from the alternative RMSEA is just over the specified power 
# rmsea0: The null RMSEA
# rmseaA: The alternative RMSEA
# power: The desired statistical power
# df: degree of freedom of the chi-square distribution
# alpha: The alpha level
# group: The number of group in calculating RMSEA
# Return The estimated sample size
findRMSEAsamplesize <- function(rmsea0, rmseaA, df, power=0.80, alpha=.05, group=1) {
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
