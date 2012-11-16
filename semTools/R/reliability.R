## Title: Reliability of factors
## Author: Sunthud Pornprasertmanit <psunthud@ku.edu>; Yves Rosseel <Yves.Rosseel@UGent.be>
## Description: Find the relability values of each factor
##----------------------------------------------------------------------------##

reliability <- function(object) {
	param <- inspect(object, "coef")
	ngroup <- object@Data@ngroups
	name <- names(param)
	ly <- param[name == "lambda"]
	ps <- impliedFactorCov(object)
	if(ngroup == 1) ps <- list(ps)
	te <- param[name == "theta"]
	S <- object@SampleStats@cov
	SigmaHat <- object@Fit@Sigma.hat
	threshold <- object@SampleStats@th
	flag <- FALSE
	result <- list()
	for(i in 1:ngroup) {
		if(!is.null(threshold[[i]])) flag <- TRUE
		common <- (apply(ly[[i]], 2, sum)^2) * diag(ps[[i]])
		error <- rep(NA, length(common))
		alpha <- rep(NA, length(common))
		total <- rep(NA, length(common))
		impliedTotal <- rep(NA, length(common))
		for(j in 1:length(error)) {
			index <- which(ly[[i]][,j] != 0)
			error[j] <- sum(te[[i]][index, index])
			sigma <- S[[i]][index, index]
			alpha[j] <- length(index)/(length(index) - 1) * (1.0 - sum(diag(sigma))/sum(sigma))
			total[j] <- sum(sigma)
			impliedTotal[j] <- sum(SigmaHat[[i]][index, index])
		}
		omega <- common/(common + error)
		singleIndicator <- apply(ly[[i]], 2, function(x) sum(x != 0)) %in% 0:1
		result[[i]] <- rbind(alpha=alpha, omega=omega, omega2=1.0 - error/impliedTotal,omega3=common/total)[,!singleIndicator]
	}
	if(flag) warning("The alpha is calculated from polychoric (polyserial) correlation not from Pearson's correlation.\n")
	if(ngroup == 1) {
		result <- result[[1]]
	} else {
		names(result) <- object@Data@group.label
	}
	result
}

reliabilityL2 <- function(object, secondFactor) {
	param <- inspect(object, "coef")
	ngroup <- object@Data@ngroups
	name <- names(param)
	ly <- param[name == "lambda"]
	ve <- impliedFactorCov(object)
	if(ngroup == 1) ve <- list(ve)
	ps <- param[name == "psi"]
	te <- param[name == "theta"]
	be <- param[name == "beta"]
	SigmaHat <- object@Fit@Sigma.hat
	S <- object@SampleStats@cov
	threshold <- object@SampleStats@th
	result <- list()
	for(i in 1:ngroup) {
		
		# Prepare for higher-order reliability
		l2var <- ve[[i]][secondFactor, secondFactor]
		l2load <- be[[1]][,secondFactor]
		indexl2 <- which(l2load != 0)
		commonl2 <- (sum(l2load)^2) * l2var
		errorl2 <- sum(ps[[i]][indexl2, indexl2])

		# Prepare for lower-order reliability
		indexl1 <- which(apply(ly[[i]][,indexl2], 1, function(x) sum(x != 0)) > 0)
		l1load <- ly[[i]][,indexl2] %*% as.matrix(be[[1]][indexl2,secondFactor])
		commonl1 <- (sum(l1load)^2) * l2var
		errorl1 <- sum(te[[i]][indexl1, indexl1])
		uniquel1 <- 0
		for (j in seq_along(indexl2)) {
			uniquel1 <- uniquel1 + (sum(ly[[i]][,indexl2[j]])^2) * ps[[i]][indexl2[j], indexl2[j]]
		}
		
		# Adjustment for direct loading from L2 to observed variables
		if(any(ly[[i]][,secondFactor] != 0)) {
			indexind <- which(ly[[i]][,secondFactor] != 0)
			if(length(intersect(indexind, indexl1)) > 0) stop("Direct and indirect loadings of higher-order factor to observed variables are specified at the same time.")
			commonl2 <- sum(c(ly[[i]][,secondFactor], l2load))^2 * l2var
			errorl2 <- errorl2 + sum(te[[i]][indexind, indexind])
			commonl1 <- sum(c(ly[[i]][,secondFactor], l1load))^2 * l2var
			errorl1 <- errorl1 + sum(te[[i]][indexind, indexind])
		}
		
		# Calculate Reliability
		omegaL1 <- commonl1 / (commonl1 + uniquel1 + errorl1)
		omegaL2 <- commonl2 / (commonl2 + errorl2)
		partialOmegaL1 <- commonl1 / (commonl1 + errorl1)
		result[[i]] <- c(omegaL1=omegaL1, omegaL2=omegaL2, partialOmegaL1=partialOmegaL1)
	}
	if(ngroup == 1) {
		result <- result[[1]]
	} else {
		names(result) <- object@Data@group.label
	}
	result
}
