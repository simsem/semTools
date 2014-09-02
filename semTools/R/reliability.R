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
	SigmaHat <- object@Fit@Sigma.hat
	tau <- param[name = "tau"]
	implied <- fitted.values(object)[name = "cov"]
	categorical <- !is.null(tau[[1]])
	threshold <- NULL
	S <- object@SampleStats@cov
	if(categorical) {
		polycor <- polycorLavaan(object) # Check the order of variable names!!!!!!!!!!!!!!!!!!!!!!!!!!
		# Check when the ordered argument is not specified!!!!!!!!!!!!!!!!!
		if(ngroup == 1) polycor <- list(polycor)
		S <- polycor
		threshold <- getThreshold(object)
	}
	flag <- FALSE
	result <- list()
	for(i in 1:ngroup) {
		common <- (apply(ly[[i]], 2, sum)^2) * diag(ps[[i]])
		truevar <- ly[[i]]%*%ps[[i]]%*%t(ly[[i]])
		error <- rep(NA, length(common))
		alpha <- rep(NA, length(common))
		total <- rep(NA, length(common))
		omega1 <- omega2 <- omega3 <- rep(NA, length(common))
		impliedTotal <- rep(NA, length(common))
		avevar <- rep(NA, length(common))
		for(j in 1:length(error)) {
			index <- which(ly[[i]][,j] != 0)
			error[j] <- sum(te[[i]][index, index])
			sigma <- S[[i]][index, index]
			alpha[j] <- computeAlpha(sigma, length(index))
			total[j] <- sum(sigma)
			impliedTotal[j] <- sum(SigmaHat[[i]][index, index])
			trueitem <- diag(truevar[index, index])
			erritem <- diag(te[[i]][index, index])
			avevar[j] <- mean(trueitem / (trueitem + erritem))
			if(categorical) {
				omega1[j] <- omegaCat(truevar[index, index], SigmaHat[[i]][index, index], threshold[[i]][index], truevar[index, index] + te[[i]][index, index])
				omega2[j] <- omegaCat(truevar[index, index], SigmaHat[[i]][index, index], threshold[[i]][index], SigmaHat[[i]][index, index])
				omega3[j] <- omegaCat(truevar[index, index], SigmaHat[[i]][index, index], threshold[[i]][index], sigma)
			} else {
				omega1[j] <- common[j] / (common[j] + error[j])
				omega2[j] <- common[j] / impliedTotal[j]
				omega3[j] <- common[j] / total[j]
			}
		}
		alpha <- c(alpha, total = computeAlpha(S[[i]], nrow(S[[i]])))
		names(alpha) <- c(names(common), "total")
		if(categorical) {
			omega1 <- c(omega1, total = omegaCat(truevar, SigmaHat[[i]], threshold[[i]], truevar + te[[i]]))
			omega2 <- c(omega2, total = omegaCat(truevar, SigmaHat[[i]], threshold[[i]], SigmaHat[[i]]))
			omega3 <- c(omega3, total = omegaCat(truevar, SigmaHat[[i]], threshold[[i]], S[[i]]))
		} else {
			omega1 <- c(omega1, total = sum(truevar) / (sum(truevar) + sum(te[[i]])))
			omega2 <- c(omega2, total = sum(truevar) / (sum(SigmaHat[[i]])))
			omega3 <- c(omega3, total = sum(truevar) / (sum(S[[i]])))
		}
		avevar <- c(avevar, total = mean(diag(truevar) / (diag(truevar) + diag(te[[i]]))))
		singleIndicator <- apply(ly[[i]], 2, function(x) sum(x != 0)) %in% 0:1
		result[[i]] <- rbind(alpha=alpha, omega=omega1, omega2=omega2,omega3=omega3, avevar = avevar)[,!singleIndicator]
	}
	if(flag) warning("The alpha and the average variance extracted are calculated from polychoric (polyserial) correlation not from Pearson's correlation.\n")
	if(ngroup == 1) {
		result <- result[[1]]
	} else {
		names(result) <- object@Data@group.label
	}
	result
}

computeAlpha <- function(S, k) k/(k - 1) * (1.0 - sum(diag(S))/sum(S))

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

omegaCat <- function(truevar, implied, threshold, denom) {
	# denom could be polychoric correlation, model-implied correlation, or model-implied without error correlation
	library(mnormt)
	polyc <- truevar
	invstdvar <- 1 / sqrt(diag(implied))
	polyr <- diag(invstdvar) %*% polyc %*% diag(invstdvar)
	nitem <- ncol(implied)
	sumnum <- 0
	addden <- 0
	for(j in 1:nitem) {
	for(jp in 1:nitem) {
		sumprobn2 <- 0
		addprobn2 <- 0
		t1 <- threshold[[j]]
		t2 <- threshold[[jp]]
		for(c in 1:length(t1)) {
		for(cp in 1:length(t2)) {
			sumprobn2 <- sumprobn2 + p2(t1[c], t2[cp], polyr[j, jp])
			addprobn2 <- addprobn2 + p2(t1[c], t2[cp], denom[j, jp])
		}
		}
		sumprobn1 <- sum(pnorm(t1))
		sumprobn1p <- sum(pnorm(t2))
		sumnum <- sumnum + (sumprobn2 - sumprobn1 * sumprobn1p)
		addden <- addden + (addprobn2 - sumprobn1 * sumprobn1p)
	}
	}
	reliab <- sumnum / addden
	reliab
}


p2 <- function(t1, t2, r) {
	pmnorm(c(t1, t2), c(0,0), matrix(c(1, r, r, 1), 2, 2))
}


polycorLavaan <- function(object) {
	ngroups <- object@Data@ngroups
	coef <- inspect(object, "coef")
	targettaunames <- NULL
	if(ngroups == 1) {
		targettaunames <- rownames(coef$tau)
	} else {
		targettaunames <- rownames(coef[[1]]$tau)
	}
	barpos <- sapply(strsplit(targettaunames, ""), function(x) which(x == "|"))
	varnames <- unique(apply(data.frame(targettaunames, barpos - 1), 1, function(x) substr(x[1], 1, x[2])))
	script <- ""
	for(i in 2:length(varnames)) {
		temp <- paste0(varnames[1:(i - 1)], collapse = " + ")
		temp <- paste0(varnames[i], "~~", temp, "\n")
		script <- paste(script, temp)
	}
	newobject <- refit(script, object)
	if(ngroups == 1) {
		return(inspect(newobject, "coef")$theta)
	} else {
		return(lapply(inspect(newobject, "coef"), "[[", "theta"))
	}
}

getThreshold <- function(object) {
	ngroups <- object@Data@ngroups
	coef <- inspect(object, "coef")
	result <- NULL
	if(ngroups == 1) {
		targettaunames <- rownames(coef$tau)
		barpos <- sapply(strsplit(targettaunames, ""), function(x) which(x == "|"))
		varthres <- apply(data.frame(targettaunames, barpos - 1), 1, function(x) substr(x[1], 1, x[2]))
		result <- list(split(coef$tau, varthres))
	} else {
		result <- list()
		for(g in 1:ngroups) {
			targettaunames <- rownames(coef[[g]]$tau)
			barpos <- sapply(strsplit(targettaunames, ""), function(x) which(x == "|"))
			varthres <- apply(data.frame(targettaunames, barpos - 1), 1, function(x) substr(x[1], 1, x[2]))
			result[[g]] <- split(coef[[g]]$tau, varthres)
		}
	}
	return(result)
}

invGeneralRelia <- function(w, truevar, totalvar) {
	1-(t(w) %*% truevar %*% w) / (t(w) %*% totalvar %*% w) 
}

invGeneralReliaCat <- function(w, polyr, threshold, denom, nitem) {
	# denom could be polychoric correlation, model-implied correlation, or model-implied without error correlation
	upper <- matrix(NA, nitem, nitem)
	lower <- matrix(NA, nitem, nitem)
	for(j in 1:nitem) {
	for(jp in 1:nitem) {
		sumprobn2 <- 0
		addprobn2 <- 0
		t1 <- threshold[[j]]
		t2 <- threshold[[jp]]
		for(c in 1:length(t1)) {
		for(cp in 1:length(t2)) {
			sumprobn2 <- sumprobn2 + p2(t1[c], t2[cp], polyr[j, jp])
			addprobn2 <- addprobn2 + p2(t1[c], t2[cp], denom[j, jp])
		}
		}
		sumprobn1 <- sum(pnorm(t1))
		sumprobn1p <- sum(pnorm(t2))
		upper[j, jp] <- (sumprobn2 - sumprobn1 * sumprobn1p)
		lower[j, jp] <- (addprobn2 - sumprobn1 * sumprobn1p)
	}
	}
	1 - (t(w) %*% upper %*% w) / (t(w) %*% lower %*% w) 
}


calcMaximalRelia <- function(truevar, totalvar, varnames) {
	start <- rep(1, nrow(truevar))
	out <- nlminb(start, invGeneralRelia, truevar = truevar, totalvar = totalvar)
	if(out$convergence != 0) stop("The numerical method for finding the maximal reliability was not converged.")
	result <- 1 - out$objective
	weight <- out$par
	weight <- weight/mean(weight)
	names(weight) <- varnames
	attr(result, "weight") <- weight
	result
}

calcMaximalReliaCat <- function(polyr, threshold, denom, nitem, varnames) {
	start <- rep(1, nrow(polyr))
	out <- nlminb(start, invGeneralReliaCat, polyr = polyr, threshold = threshold, denom = denom, nitem = nitem)
	if(out$convergence != 0) stop("The numerical method for finding the maximal reliability was not converged.")
	result <- 1 - out$objective
	weight <- out$par
	weight <- weight/mean(weight)
	names(weight) <- varnames
	attr(result, "weight") <- weight
	result
}

maximalRelia <- function(object) {
	param <- inspect(object, "coef")
	ngroup <- object@Data@ngroups
	name <- names(param)
	ly <- param[name == "lambda"]
	ps <- impliedFactorCov(object)
	if(ngroup == 1) ps <- list(ps)
	te <- param[name == "theta"]
	SigmaHat <- object@Fit@Sigma.hat
	tau <- param[name = "tau"]
	implied <- fitted.values(object)[name = "cov"]
	categorical <- !is.null(tau[[1]])
	threshold <- NULL
	S <- object@SampleStats@cov
	result <- list()
	if(categorical) {
		polycor <- polycorLavaan(object)
		if(ngroup == 1) polycor <- list(polycor)
		S <- polycor
		threshold <- getThreshold(object)
	}
	for(i in 1:ngroup) {
		truevar <- ly[[i]]%*%ps[[i]]%*%t(ly[[i]])
		varnames <- colnames(truevar)
		if(categorical) {
			invstdvar <- 1 / sqrt(diag(SigmaHat[[i]]))
			polyr <- diag(invstdvar) %*% truevar %*% diag(invstdvar)
			nitem <- ncol(SigmaHat[[i]])
			result[[i]] <- calcMaximalReliaCat(polyr, threshold[[i]], S[[i]], nitem, varnames)
		} else {
			result[[i]] <- calcMaximalRelia(truevar, S[[i]], varnames)
		}
	}
	if(ngroup == 1) {
		result <- result[[1]]
	} else {
		names(result) <- object@Data@group.label
	}
	result
}
