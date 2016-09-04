# miPowerFit: Evaluate model fit by Satorra, Saris, & van der Weld (2009) method

miPowerFit <- function(lavaanObj, stdLoad=0.4, cor=0.1, stdBeta=0.1, intcept=0.2, stdDelta=NULL, delta=NULL, cilevel=0.90) {
	mi <- lavaan::lavInspect(lavaanObj, "mi")
	mi <- mi[mi$op != "==",]
	sigma <- mi[,"epc"] / sqrt(mi[,"mi"])
	if(is.null(delta)) {
		if(is.null(stdDelta)) stdDelta <- getTrivialEpc(mi, stdLoad=stdLoad, cor=cor, stdBeta=stdBeta, intcept=intcept)
		if(length(stdDelta) == 1) stdDelta <- rep(stdDelta, nrow(mi)) 
		delta <- unstandardizeEpc(mi, stdDelta, findTotalVar(lavaanObj))
	}
	if(length(delta) == 1) delta <- rep(delta, nrow(mi))
	ncp <- (delta / sigma)^2
	alpha <- 0.05
	desiredPow <- 0.80
	cutoff <- qchisq(1 - alpha, df = 1)
	pow <- 1 - pchisq(cutoff, df = 1, ncp=ncp)
	sigMI <- mi[,"mi"] > cutoff
	highPow <- pow > desiredPow
	group <- rep(1, nrow(mi))
	if("group" %in% colnames(mi)) group <- mi[,"group"]
	decision <- mapply(decisionMIPow, sigMI=sigMI, highPow=highPow, epc=mi[,"epc"], trivialEpc=delta)
	if(is.null(stdDelta)) stdDelta <- standardizeEpc(mi, findTotalVar(lavaanObj), delta=delta)
	result <- cbind(mi[,1:3], group, as.numeric(mi[,"mi"]), mi[,"epc"], delta, standardizeEpc(mi, findTotalVar(lavaanObj)), stdDelta, sigMI, highPow, decision)
	# New method
	crit <- abs(qnorm((1 - cilevel)/2))
	seepc <- abs(result[,6]) / sqrt(abs(result[,5]))
	lowerepc <- result[,6] - crit * seepc
	upperepc <- result[,6] + crit * seepc
	stdlowerepc <- standardizeEpc(mi, findTotalVar(lavaanObj), delta = lowerepc)
	stdupperepc <- standardizeEpc(mi, findTotalVar(lavaanObj), delta = upperepc)
	isVar <- mi[,"op"] == "~~" & mi[,"lhs"] == mi[,"rhs"]
	decisionci <- mapply(decisionCIEpc, targetval=as.numeric(stdDelta), lower=stdlowerepc, upper=stdupperepc, positiveonly=isVar)
	
	result <- cbind(result, seepc, lowerepc, upperepc, stdlowerepc, stdupperepc, decisionci)
	result <- result[!is.na(decision),]
	colnames(result) <- c("lhs", "op", "rhs", "group", "mi", "epc", "target.epc", "std.epc", "std.target.epc", "significant.mi", "high.power", "decision.pow", "se.epc", "lower.epc", "upper.epc", "lower.std.epc", "upper.std.epc", "decision.ci")
	result <- format(result, scientific=FALSE, digits=4)
	return(result)
}

# totalFacVar: Find total factor variances when regression coeffient matrix and factor residual covariance matrix are specified

totalFacVar <- function(beta, psi) {
	ID <- diag(nrow(psi))
    total <- solve(ID - beta) %*% psi %*% t(solve(ID - beta))
    return(diag(total))
}

# findTotalVar: find the total indicator and factor variances

findTotalVar <- function(lavaanObj) {
	result <- list()
	nGroups <- lavaan::lavInspect(lavaanObj, "ngroups")
	cov.all <- lavaan::lavInspect(lavaanObj, "cov.all")
	if(nGroups == 1) cov.all <- list(cov.all)
	for(i in 1:nGroups) {
		temp <- diag(cov.all[[i]])
		names(temp) <- rownames(cov.all[[i]])
		result[[i]] <- temp
	}
	return(result)
}

# getTrivialEpc: find the trivial misspecified expected parameter changes given the type of parameters in each row of modification indices

getTrivialEpc <- function(mi, stdLoad=0.4, cor=0.1, stdBeta=0.1, intcept=0.2) {
	op <- mi[,"op"]
	result <- gsub("=~", stdLoad, op)
	result <- gsub("~~", cor, result)
	result <- gsub("~1", intcept, result)
	result <- gsub("~", stdBeta, result)
	return(result)
}

# unstandardizeEpc: Transform from standardized EPC to unstandardized EPC

unstandardizeEpc <- function(mi, delta, totalVar) {
	name <- names(totalVar[[1]])
	lhsPos <- match(mi[,"lhs"], name)
	rhsPos <- match(mi[,"rhs"], name)
	group <- rep(1, nrow(mi))
	if("group" %in% colnames(mi)) group <- mi[,"group"]
	getVar <- function(pos, group) totalVar[[group]][pos]
	lhsVar <- mapply(getVar, pos=lhsPos, group=group)
	rhsVar <- mapply(getVar, pos=rhsPos, group=group)
	FUN <- function(op, lhsVar, rhsVar, delta) {
		if(op == "|") return(NA)
		lhsSD <- sqrt(lhsVar)
		rhsSD <- sqrt(rhsVar)
		if(!is.numeric(delta)) delta <- as.numeric(delta)
		if(op == "=~") {
			return((rhsSD * delta) / lhsSD)
		} else if (op == "~~") {
			return(lhsSD * delta * rhsSD)
		} else if (op == "~1") {
			return(lhsSD * delta)
		} else if (op == "~") {
			return((lhsSD * delta) / rhsSD)
		} else {
			return(NA)
		}
	}
	unstdDelta <- mapply(FUN, op=mi[,"op"], lhsVar=lhsVar, rhsVar=rhsVar, delta=delta)
	return(unstdDelta)
}

# unstandardizeEpc: Transform from unstandardized EPC to standardized EPC. If delta is null, the unstandardized epc from the modification indices data.frame are used

standardizeEpc <- function(mi, totalVar, delta=NULL) {
	if(is.null(delta)) delta <- mi[,"epc"]
	name <- names(totalVar[[1]])
	lhsPos <- match(mi[,"lhs"], name)
	rhsPos <- match(mi[,"rhs"], name)
	group <- rep(1, nrow(mi))
	if("group" %in% colnames(mi)) group <- mi[,"group"]
	getVar <- function(pos, group) totalVar[[group]][pos]
	lhsVar <- mapply(getVar, pos=lhsPos, group=group)
	rhsVar <- mapply(getVar, pos=rhsPos, group=group)
	FUN <- function(op, lhsVar, rhsVar, delta) {
		lhsSD <- sqrt(lhsVar)
		rhsSD <- sqrt(rhsVar)
		if(!is.numeric(delta)) delta <- as.numeric(delta)
		if(op == "=~") {
			#stdload = beta * sdlatent / sdindicator = beta * lhs / rhs
			return((delta / rhsSD) * lhsSD)
		} else if (op == "~~") {
			#r = cov / (sd1 * sd2)
			return(delta / (lhsSD * rhsSD))
		} else if (op == "~1") {
			#d = meanDiff/sd
			return(delta / lhsSD)
		} else if (op == "~") {
			#beta = b * sdX / sdY = b * rhs / lhs
			return((delta / lhsSD) * rhsSD)
		} else {
			return(NA)
		}
	}
	stdDelta <- mapply(FUN, op=mi[,"op"], lhsVar=lhsVar, rhsVar=rhsVar, delta=delta)
	return(stdDelta)
}

# decisionMIPow: provide the decision given the significance of modification indices and power to detect trivial misspecification

decisionMIPow <- function(sigMI, highPow, epc, trivialEpc) {
	if(is.na(sigMI) | is.na(highPow)) return(NA)
	if(sigMI & highPow) {
		if(abs(epc) > abs(trivialEpc)) {
			return("EPC:M")
		} else {
			return("EPC:NM")
		}
	} else if (sigMI & !highPow) {
		return("M")
	} else if (!sigMI & highPow) {
		return("NM")
	} else if (!sigMI & !highPow) {
		return("I")
	} else {
		return(NA)
	}
}

decisionCIEpc <- function(targetval, lower, upper, positiveonly = FALSE) {
	if(is.na(lower) | is.na(upper)) return(NA)
	if(positiveonly) {
		if(lower > targetval) {
			return("M")
		} else if (upper < targetval) {
			return("NM")
		} else {
			return("I")
		}
	} else {
		negtargetval <- -targetval
		if(lower > targetval | upper < negtargetval) {
			return("M")
		} else if (upper < targetval & negtargetval < lower) {
			return("NM")
		} else {
			return("I")
		}	
	}
}
