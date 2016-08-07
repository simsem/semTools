## Title: Probing Interaction
## Author: Sunthud Pornprasertmanit <psunthud@ku.edu>
## Description: Probing Interaction with Residual Centering
##----------------------------------------------------------------------------##

probe2WayMC <- function(fit, nameX, nameY, modVar, valProbe) {
	# Check whether modVar is correct
	if(is.character(modVar)) {
		modVar <- match(modVar, nameX)
	} 
	if(is.na(modVar) || !(modVar %in% 1:2)) stop("The moderator name is not in the name of independent factors or not 1 or 2.")

	# Check whether the fit object does not use mlm, mlr, or mlf estimator (because the variance-covariance matrix of parameter estimates cannot be computed
	estSpec <- lavaan::lavInspect(fit, "call")$estimator
	if(!is.null(estSpec) && (estSpec %in% c("mlr", "mlm", "mlf"))) stop("This function does not work when 'mlr', 'mlm', or 'mlf' is used as the estimator because the covariance matrix of the parameter estimates cannot be computed.")
	
	# Get the parameter estimate values from the lavaan object
	est <- lavaan::lavInspect(fit, "coef")

	# Compute the intercept of no-centering
	betaNC <- as.matrix(est$beta[nameY, nameX]); colnames(betaNC) <- nameY

	# Extract all varEst
	varEst <- lavaan::vcov(fit) 
	
	# Check whether intercept are estimated
	targetcol <- paste(nameY, "~", 1, sep="") 
	estimateIntcept <- targetcol %in% rownames(varEst) 
	
	pvalue <- function(x) (1 - pnorm(abs(x))) * 2 
	
	resultIntcept <- NULL 
	resultSlope <- NULL 
	if(estimateIntcept) {
		# Extract SE from residual centering
		targetcol <- c(targetcol, paste(nameY, "~", nameX, sep="")) 

		# Transform it to non-centering SE
		usedVar <-  varEst[targetcol, targetcol] 
		usedBeta <- rbind(est$alpha[nameY,], betaNC) 
		
		# Change the order of usedVar and usedBeta if the moderator variable is listed first 
		if(modVar == 1) {
			usedVar <- usedVar[c(1, 3, 2, 4), c(1, 3, 2, 4)]
			usedBeta <- usedBeta[c(1, 3, 2, 4)]
		}
		
		# Find simple intercept  
		simpleIntcept <- usedBeta[1] + usedBeta[3] * valProbe 
		varIntcept <- usedVar[1, 1] + 2 * valProbe * usedVar[1, 3] + (valProbe^2) * usedVar[3, 3]
		zIntcept <- simpleIntcept/sqrt(varIntcept)
		pIntcept <- round(pvalue(zIntcept),6) #JG: rounded values to make them more readable
		resultIntcept <- cbind(valProbe, simpleIntcept, sqrt(varIntcept), zIntcept, pIntcept)
		colnames(resultIntcept) <- c(nameX[modVar], "Intcept", "SE", "Wald", "p")
		
		# Find simple slope 
		simpleSlope <- usedBeta[2] + usedBeta[4] * valProbe
		varSlope <- usedVar[2, 2] + 2 * valProbe * usedVar[2, 4] + (valProbe^2) * usedVar[4, 4]
		zSlope <- simpleSlope/sqrt(varSlope)
		pSlope <- pvalue(zSlope)
		resultSlope <- cbind(valProbe, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "Slope", "SE", "Wald", "p")
	} else {
		targetcol <- paste(nameY, "~", nameX, sep="") 

		# Transform it to non-centering SE 
		usedVar <-  varEst[targetcol, targetcol]
		usedBeta <- betaNC
		
		# Change the order of usedVar and usedBeta if the moderator variable is listed first 
		if(modVar == 2) {
			usedVar <- usedVar[c(2, 1, 3), c(2, 1, 3)]
			# usedBeta <- usedBeta[c(2, 1, 3)]
		}
		
		# Find simple slope 
		simpleSlope <- usedBeta[1] + usedBeta[3] * valProbe
		varSlope <- usedVar[1, 1] + 2 * valProbe * usedVar[1, 3] + (valProbe^2) * usedVar[3, 3]
		zSlope <- simpleSlope/sqrt(varSlope)
		pSlope <- round(pvalue(zSlope),6) #JG: rounded values to make them more readable
		resultSlope <- cbind(valProbe, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "Slope", "SE", "Wald", "p")
	}
	
	return(list(SimpleIntcept = resultIntcept, SimpleSlope = resultSlope))
}

probe2WayRC <- function(fit, nameX, nameY, modVar, valProbe) {
	# Check whether modVar is correct
	if(is.character(modVar)) {
		modVar <- match(modVar, nameX)
	} 
	if(is.na(modVar) || !(modVar %in% 1:2)) stop("The moderator name is not in the name of independent factors or not 1 or 2.")

	# Check whether the fit object does not use mlm, mlr, or mlf estimator (because the variance-covariance matrix of parameter estimates cannot be computed
	estSpec <- lavaan::lavInspect(fit, "call")$estimator
	if(!is.null(estSpec) && (estSpec %in% c("mlr", "mlm", "mlf"))) stop("This function does not work when 'mlr', 'mlm', or 'mlf' is used as the estimator because the covariance matrix of the parameter estimates cannot be computed.")
	
	# Get the parameter estimate values from the lavaan object
	est <- lavaan::lavInspect(fit, "coef")
	
	# Find the mean and covariance matrix of independent factors
	varX <- est$psi[nameX, nameX]
	meanX <- as.matrix(est$alpha[nameX,]); colnames(meanX) <- "intcept"

	# Find the intercept, regression coefficients, and residual variance of residual-centered regression
	intceptRC <- est$alpha[nameY,]
	resVarRC <- est$psi[nameY, nameY]
	betaRC <- as.matrix(est$beta[nameY, nameX]); colnames(betaRC) <- nameY

	# Find the number of observations
	numobs <- lavaan::lavInspect(fit, "nobs")
	
	# Compute SSRC 
	meanXwith1 <- rbind(1, meanX) 
	varXwith0 <- cbind(0, rbind(0, varX)) 
	SSRC <- numobs * (varXwith0 + (meanXwith1 %*% t(meanXwith1))) 

	# Compute Mean(Y) and Var(Y)
	betaRCWithIntcept <- rbind(intceptRC, betaRC)
	meanY <- t(meanXwith1) %*% betaRCWithIntcept 
	varY <- (t(betaRCWithIntcept) %*% SSRC %*% betaRCWithIntcept)/numobs - meanY^2 + resVarRC 

	# Compute Cov(Y, X)
	covY <- as.matrix((varX %*% betaRC)[1:2,]) 

	# Compute E(XZ)
	meanX[3] <- meanX[1] * meanX[2] + varX[1, 2] 
	
	# Compute Var(XZ)
	varX[3, 3] <- meanX[1]^2 * varX[2, 2] + meanX[2]^2 * varX[1, 1] + 2 * meanX[1] * meanX[2] * varX[1, 2] + varX[1, 1] * varX[2, 2] + varX[1, 2]^2  

	# Compute Cov(X, XZ), Cov(Z, XZ)
	varX[1, 3] <- varX[3, 1] <- meanX[1] * varX[1, 2] + meanX[2] * varX[1, 1] 
	varX[2, 3] <- varX[3, 2] <- meanX[1] * varX[2, 2] + meanX[2] * varX[1, 2] 

	# Compute Cov(Y, XZ) and regression coefficients of no-centering
	betaNC <- solve(varX[1:2,1:2], covY - rbind(varX[1,3] * betaRC[3,1], varX[2, 3] * betaRC[3,1])) 
	betaNC <- rbind(betaNC, betaRC[3, 1]) 
	covY <- rbind(covY, (varX %*% betaNC)[3, 1]) 

	# Aggregate the non-centering sufficient statistics (Just show how to do but not necessary)
	fullCov <- rbind(cbind(varX, covY), c(covY, varY)) 
	fullMean <- rbind(meanX, meanY) 

	# Compute the intercept of no-centering
	intceptNC <- meanY - t(betaNC) %*% meanX 

	# Compute SSNC
	betaNCWithIntcept <- rbind(intceptNC, betaNC)
	meanXwith1 <- rbind(1, meanX) 
	varXwith0 <- rbind(0, cbind(0, varX)) 
	SSNC <- numobs * (varXwith0 + (meanXwith1 %*% t(meanXwith1))) 

	# Compute residual variance on non-centering
	resVarNC <- varY - (t(betaNCWithIntcept) %*% SSNC %*% betaNCWithIntcept)/numobs + meanY^2 

	# Extract all varEst
	varEst <- lavaan::vcov(fit) 
	
	# Check whether intercept are estimated
	targetcol <- paste(nameY, "~", 1, sep="") 
	estimateIntcept <- targetcol %in% rownames(varEst) 
	
	pvalue <- function(x) (1 - pnorm(abs(x))) * 2 
	
	resultIntcept <- NULL 
	resultSlope <- NULL 
	if(estimateIntcept) {
		# Extract SE from residual centering
		targetcol <- c(targetcol, paste(nameY, "~", nameX, sep="")) 
		varEstSlopeRC <- varEst[targetcol, targetcol] 

		# Transform it to non-centering SE
		usedVar <-  as.numeric(resVarNC/resVarRC) * (varEstSlopeRC %*% SSRC %*% solve(SSNC)) 
		usedBeta <- betaNCWithIntcept 
		
		# Change the order of usedVar and usedBeta if the moderator variable is listed first 
		if(modVar == 1) {
			usedVar <- usedVar[c(1, 3, 2, 4), c(1, 3, 2, 4)]
			usedBeta <- usedBeta[c(1, 3, 2, 4)]
		}
		
		# Find simple intercept  
		simpleIntcept <- usedBeta[1] + usedBeta[3] * valProbe 
		varIntcept <- usedVar[1, 1] + 2 * valProbe * usedVar[1, 3] + (valProbe^2) * usedVar[3, 3]
		zIntcept <- simpleIntcept/sqrt(varIntcept)
		pIntcept <- round(pvalue(zIntcept),6) #JG: rounded values to make them more readable
		resultIntcept <- cbind(valProbe, simpleIntcept, sqrt(varIntcept), zIntcept, pIntcept)
		colnames(resultIntcept) <- c(nameX[modVar], "Intcept", "SE", "Wald", "p")
		
		# Find simple slope 
		simpleSlope <- usedBeta[2] + usedBeta[4] * valProbe
		varSlope <- usedVar[2, 2] + 2 * valProbe * usedVar[2, 4] + (valProbe^2) * usedVar[4, 4]
		zSlope <- simpleSlope/sqrt(varSlope)
		pSlope <- pvalue(zSlope)
		resultSlope <- cbind(valProbe, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "Slope", "SE", "Wald", "p")
	} else {
		targetcol <- paste(nameY, "~", nameX, sep="") 
		varEstSlopeRC <- varEst[targetcol, targetcol] 

		# Transform it to non-centering SE 
		usedVar <-  as.numeric(resVarNC/resVarRC) * (varEstSlopeRC %*% SSRC[2:4, 2:4] %*% solve(SSNC[2:4, 2:4]))
		usedBeta <- betaNC
		
		# Change the order of usedVar and usedBeta if the moderator variable is listed first 
		if(modVar == 2) {
			usedVar <- usedVar[c(2, 1, 3), c(2, 1, 3)]
			# usedBeta <- usedBeta[c(2, 1, 3)]
		}
		
		# Find simple slope 
		simpleSlope <- usedBeta[1] + usedBeta[3] * valProbe
		varSlope <- usedVar[1, 1] + 2 * valProbe * usedVar[1, 3] + (valProbe^2) * usedVar[3, 3]
		zSlope <- simpleSlope/sqrt(varSlope)
		pSlope <- round(pvalue(zSlope),6) #JG: rounded values to make them more readable
		resultSlope <- cbind(valProbe, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "Slope", "SE", "Wald", "p")
	}
	
	return(list(SimpleIntcept = resultIntcept, SimpleSlope = resultSlope))
}


probe3WayMC <- function(fit, nameX, nameY, modVar, valProbe1, valProbe2) {
	# Check whether modVar is correct
	if(is.character(modVar)) {
		modVar <- match(modVar, nameX)
	} 
	if((NA %in% modVar) || !(do.call("&", as.list(modVar %in% 1:3)))) stop("The moderator name is not in the list of independent factors and is not 1, 2 or 3.") 

	# Check whether the fit object does not use mlm, mlr, or mlf estimator (because the variance-covariance matrix of parameter estimates cannot be computed
	estSpec <- lavaan::lavInspect(fit, "call")$estimator
	if(!is.null(estSpec) && (estSpec %in% c("mlr", "mlm", "mlf"))) stop("This function does not work when 'mlr', 'mlm', or 'mlf' is used as the estimator because the covariance matrix of the parameter estimates cannot be computed.")

	# Get the parameter estimate values from the lavaan object
	est <- lavaan::lavInspect(fit, "coef")

	# Compute the intercept of no-centering
	betaNC <- as.matrix(est$beta[nameY, nameX]); colnames(betaNC) <- nameY

	# Extract all varEst
	varEst <- lavaan::vcov(fit)
	
	# Check whether intercept are estimated
	targetcol <- paste(nameY, "~", 1, sep="")
	estimateIntcept <- targetcol %in% rownames(varEst)
	
	pvalue <- function(x) (1 - pnorm(abs(x))) * 2
	
	# Find the order to rearrange
	ord <- c(setdiff(1:3, modVar), modVar)
	ord <- c(ord, 7 - rev(ord))
	
	resultIntcept <- NULL
	resultSlope <- NULL
	if(estimateIntcept) {
		# Extract SE from residual centering
		targetcol <- c(targetcol, paste(nameY, "~", nameX, sep=""))

		# Transform it to non-centering SE
		usedVar <-  varEst[targetcol, targetcol] 
		usedBeta <- rbind(est$alpha[nameY,], betaNC) 
		if(sum(diag(usedVar) < 0) > 0) stop("This method does not work. The resulting calculation provided negative standard errors.") # JG: edited this error
		
		# Change the order of usedVar and usedBeta if the moderator variable is listed first
		usedVar <- usedVar[c(1, ord+1, 8), c(1, ord+1, 8)]
		usedBeta <- usedBeta[c(1, ord+1, 8)]
				
		# Find probe value
		val <- expand.grid(valProbe1, valProbe2)
		
		# Find simple intercept
		simpleIntcept <- usedBeta[1] + usedBeta[3] * val[,1] + usedBeta[4] * val[,2] + usedBeta[7] * val[,1] * val[,2]
		varIntcept <- usedVar[1, 1] + val[,1]^2 * usedVar[3, 3] + val[,2]^2 * usedVar[4, 4] + val[,1]^2 * val[,2]^2 * usedVar[7, 7] + 2 * val[,1] * usedVar[1, 3] + 2 * val[,2] * usedVar[1, 4] + 2 * val[,1] * val[,2] * usedVar[1, 7] + 2 * val[,1] * val[,2] * usedVar[3, 4] + 2 * val[,1]^2 * val[,2] * usedVar[3, 7] + 2* val[,1] * val[,2]^2 * usedVar[4, 7]
		zIntcept <- simpleIntcept/sqrt(varIntcept)
		pIntcept <- pvalue(zIntcept)
		resultIntcept <- cbind(val, simpleIntcept, sqrt(varIntcept), zIntcept, pIntcept)
		colnames(resultIntcept) <- c(nameX[modVar], "Intcept", "SE", "Wald", "p")
		
		# Find simple slope
		simpleSlope <- usedBeta[2] + usedBeta[5] * val[,1] + usedBeta[6] * val[,2] + usedBeta[8] * val[,1] * val[,2]
		varSlope <- usedVar[2, 2] + val[,1]^2 * usedVar[5, 5] + val[,2]^2 * usedVar[6, 6] + val[,1]^2 * val[,2]^2 * usedVar[8, 8] + 2 * val[,1] * usedVar[2, 5] + 2 * val[,2] * usedVar[2, 6] + 2 * val[,1] * val[,2] * usedVar[2, 8] + 2 * val[,1] * val[,2] * usedVar[5, 6] + 2 * val[,1]^2 * val[,2] * usedVar[5, 8] + 2 * val[,1] * val[,2]^2 * usedVar[6, 8]
		zSlope <- simpleSlope/sqrt(varSlope)
		pSlope <- round(pvalue(zSlope),6) # JG: rounded values
		resultSlope <- cbind(val, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "Slope", "SE", "Wald", "p")
	} else {
		targetcol <- paste(nameY, "~", nameX, sep="")

		# Transform it to non-centering SE
		usedVar <-  varEst[targetcol, targetcol]
		usedBeta <- betaNC
		if(sum(diag(usedVar) < 0) > 0) stop("This method does not work. The resulting calculation provided negative standard errors.") # JG: edited this error
		
		# Change the order of usedVar and usedBeta if the moderator variable is listed first
		usedVar <- usedVar[c(ord, 7), c(ord, 7)]
		usedBeta <- usedBeta[c(ord, 7)]

		# Find probe value
		val <- expand.grid(valProbe1, valProbe2)
		
		# Find simple slope
		simpleSlope <- usedBeta[1] + usedBeta[4] * val[,1] + usedBeta[5] * val[,2] + usedBeta[7] * val[,1] * val[,2]
		varSlope <- usedVar[1, 1] + val[,1]^2 * usedVar[4, 4] + val[,2]^2 * usedVar[5, 5] + val[,1]^2 * val[,2]^2 * usedVar[7, 7] + 2 * val[,1] * usedVar[1, 4] + 2 * val[,2] * usedVar[1, 5] + 2 * val[,1] * val[,2] * usedVar[1, 7] + 2 * val[,1] * val[,2] * usedVar[4, 5] + 2 * val[,1]^2 * val[,2] * usedVar[4, 7] + 2 * val[,1] * val[,2]^2 * usedVar[5, 7]
		zSlope <- simpleSlope/sqrt(varSlope)
		pSlope <- round(pvalue(zSlope),6) # JG: rounded values
		resultSlope <- cbind(val, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "Slope", "SE", "Wald", "p")
	}
	
	return(list(SimpleIntcept = resultIntcept, SimpleSlope = resultSlope))
}

probe3WayRC <- function(fit, nameX, nameY, modVar, valProbe1, valProbe2) {
	# Check whether modVar is correct
	if(is.character(modVar)) {
		modVar <- match(modVar, nameX)
	} 
	if((NA %in% modVar) || !(do.call("&", as.list(modVar %in% 1:3)))) stop("The moderator name is not in the list of independent factors and is not 1, 2 or 3.") # JG: Changed error

	# Check whether the fit object does not use mlm, mlr, or mlf estimator (because the variance-covariance matrix of parameter estimates cannot be computed
	estSpec <- lavaan::lavInspect(fit, "call")$estimator
	if(!is.null(estSpec) && (estSpec %in% c("mlr", "mlm", "mlf"))) stop("This function does not work when 'mlr', 'mlm', or 'mlf' is used as the estimator because the covariance matrix of the parameter estimates cannot be computed.")

	# Get the parameter estimate values from the lavaan object
	est <- lavaan::lavInspect(fit, "coef")
	
	# Find the mean and covariance matrix of independent factors
	varX <- est$psi[nameX, nameX]
	meanX <- as.matrix(est$alpha[nameX,]); colnames(meanX) <- "intcept"

	# Find the intercept, regression coefficients, and residual variance of residual-centered regression
	intceptRC <- est$alpha[nameY,]
	resVarRC <- est$psi[nameY, nameY]
	if(resVarRC < 0) stop("The residual variance is negative. The model did not converge!") # JG: Changed error
	betaRC <- as.matrix(est$beta[nameY, nameX]); colnames(betaRC) <- nameY

	# Find the number of observations
	numobs <- lavaan::lavInspect(fit, "nobs")

	# Compute SSRC
	meanXwith1 <- rbind(1, meanX)
	varXwith0 <- cbind(0, rbind(0, varX))
	SSRC <- numobs * (varXwith0 + (meanXwith1 %*% t(meanXwith1)))

	# Compute Mean(Y) and Var(Y)
	betaRCWithIntcept <- rbind(intceptRC, betaRC)
	meanY <- t(meanXwith1) %*% betaRCWithIntcept
	varY <- (t(betaRCWithIntcept) %*% SSRC %*% betaRCWithIntcept)/numobs - meanY^2 + resVarRC

	# Compute Cov(Y, X)
	covY <- as.matrix((varX %*% betaRC)[1:3,])

	# Compute E(XZ), E(XW), E(ZW), E(XZW)
	meanX[4] <- expect2NormProd(meanX[c(1,2)], varX[c(1,2), c(1,2)])
	meanX[5] <- expect2NormProd(meanX[c(1,3)], varX[c(1,3), c(1,3)])
	meanX[6] <- expect2NormProd(meanX[c(2,3)], varX[c(2,3), c(2,3)])
	meanX[7] <- expect3NormProd(meanX[1:3], varX[1:3, 1:3])
	
	# Compute Var(XZ), Var(XW), Var(ZW), Var(XZW)
	varX[4, 4] <- var2NormProd(meanX[c(1,2)], varX[c(1,2), c(1,2)])
	varX[5, 5] <- var2NormProd(meanX[c(1,3)], varX[c(1,3), c(1,3)])
	varX[6, 6] <- var2NormProd(meanX[c(2,3)], varX[c(2,3), c(2,3)])
	varX[7, 7] <- var3NormProd(meanX[1:3], varX[1:3, 1:3])
	
	# Compute All covariances
	varX[4, 1] <- varX[1, 4] <- expect3NormProd(meanX[c(1, 2, 1)], varX[c(1, 2, 1),c(1, 2, 1)]) - expect2NormProd(meanX[c(1,2)], varX[c(1,2),c(1,2)]) * meanX[1]
	varX[5, 1] <- varX[1, 5] <- expect3NormProd(meanX[c(1, 3, 1)], varX[c(1, 3, 1),c(1, 3, 1)]) - expect2NormProd(meanX[c(1,3)], varX[c(1,3),c(1,3)]) * meanX[1]
	varX[6, 1] <- varX[1, 6] <- expect3NormProd(meanX[c(2, 3, 1)], varX[c(2, 3, 1),c(2, 3, 1)]) - expect2NormProd(meanX[c(2,3)], varX[c(2,3),c(2,3)]) * meanX[1]
	varX[7, 1] <- varX[1, 7] <- expect4NormProd(meanX[c(1,2,3,1)], varX[c(1,2,3,1),c(1,2,3,1)]) - expect3NormProd(meanX[c(1,2,3)], varX[c(1,2,3),c(1,2,3)]) * meanX[1]
	
	varX[4, 2] <- varX[2, 4] <- expect3NormProd(meanX[c(1, 2, 2)], varX[c(1, 2, 2),c(1, 2, 2)]) - expect2NormProd(meanX[c(1,2)], varX[c(1,2),c(1,2)]) * meanX[2]
	varX[5, 2] <- varX[2, 5] <- expect3NormProd(meanX[c(1, 3, 2)], varX[c(1, 3, 2),c(1, 3, 2)]) - expect2NormProd(meanX[c(1,3)], varX[c(1,3),c(1,3)]) * meanX[2]
	varX[6, 2] <- varX[2, 6] <- expect3NormProd(meanX[c(2, 3, 2)], varX[c(2, 3, 2),c(2, 3, 2)]) - expect2NormProd(meanX[c(2,3)], varX[c(2,3),c(2,3)]) * meanX[2]
	varX[7, 2] <- varX[2, 7] <- expect4NormProd(meanX[c(1,2,3,2)], varX[c(1,2,3,2),c(1,2,3,2)]) - expect3NormProd(meanX[c(1,2,3)], varX[c(1,2,3),c(1,2,3)]) * meanX[2]
	
	varX[4, 3] <- varX[3, 4] <- expect3NormProd(meanX[c(1, 2, 3)], varX[c(1, 2, 3),c(1, 2, 3)]) - expect2NormProd(meanX[c(1,2)], varX[c(1,2),c(1,2)]) * meanX[3]
	varX[5, 3] <- varX[3, 5] <- expect3NormProd(meanX[c(1, 3, 3)], varX[c(1, 3, 3),c(1, 3, 3)]) - expect2NormProd(meanX[c(1,3)], varX[c(1,3),c(1,3)]) * meanX[3]
	varX[6, 3] <- varX[3, 6] <- expect3NormProd(meanX[c(2, 3, 3)], varX[c(2, 3, 3),c(2, 3, 3)]) - expect2NormProd(meanX[c(2,3)], varX[c(2,3),c(2,3)]) * meanX[3]
	varX[7, 3] <- varX[3, 7] <- expect4NormProd(meanX[c(1,2,3,3)], varX[c(1,2,3,3),c(1,2,3,3)]) - expect3NormProd(meanX[c(1,2,3)], varX[c(1,2,3),c(1,2,3)]) * meanX[3]
	
	varX[5, 4] <- varX[4, 5] <- expect4NormProd(meanX[c(1,3,1,2)], varX[c(1,3,1,2),c(1,3,1,2)]) - expect2NormProd(meanX[c(1,3)], varX[c(1,3),c(1,3)]) * expect2NormProd(meanX[c(1,2)], varX[c(1,2),c(1,2)])
	varX[6, 4] <- varX[4, 6] <- expect4NormProd(meanX[c(2,3,1,2)], varX[c(2,3,1,2),c(2,3,1,2)]) - expect2NormProd(meanX[c(2,3)], varX[c(2,3),c(2,3)]) * expect2NormProd(meanX[c(1,2)], varX[c(1,2),c(1,2)])
	varX[7, 4] <- varX[4, 7] <- expect5NormProd(meanX[c(1,2,3,1,2)], varX[c(1,2,3,1,2),c(1,2,3,1,2)]) - expect3NormProd(meanX[c(1, 2, 3)], varX[c(1, 2, 3),c(1, 2, 3)]) * expect2NormProd(meanX[c(1,2)], varX[c(1,2),c(1,2)])

	varX[6, 5] <- varX[5, 6] <- expect4NormProd(meanX[c(2,3,1,3)], varX[c(2,3,1,3),c(2,3,1,3)]) - expect2NormProd(meanX[c(2,3)], varX[c(2,3),c(2,3)]) * expect2NormProd(meanX[c(1,3)], varX[c(1,3),c(1,3)])
	varX[7, 5] <- varX[5, 7] <- expect5NormProd(meanX[c(1,2,3,1,3)], varX[c(1,2,3,1,3),c(1,2,3,1,3)]) - expect3NormProd(meanX[c(1, 2, 3)], varX[c(1, 2, 3),c(1, 2, 3)]) * expect2NormProd(meanX[c(1,3)], varX[c(1,3),c(1,3)])
	varX[7, 6] <- varX[6, 7] <- expect5NormProd(meanX[c(1,2,3,2,3)], varX[c(1,2,3,2,3),c(1,2,3,2,3)]) - expect3NormProd(meanX[c(1, 2, 3)], varX[c(1, 2, 3),c(1, 2, 3)]) * expect2NormProd(meanX[c(2,3)], varX[c(2,3),c(2,3)])
	
	# Find the meanX and varX without XZW
	meanXReducedWith1 <- rbind(1, as.matrix(meanX[1:6]))
	varXReducedWith0 <- cbind(0, rbind(0, varX[1:6, 1:6]))
	SSMCReduced <- numobs * (varXReducedWith0 + (meanXReducedWith1 %*% t(meanXReducedWith1)))
	
	# Find product of main and two-way onto three-way
	covXZWwith0 <- rbind(0, as.matrix(varX[7, 1:6]))
	meanXZWwith1 <- meanX[7] * meanXReducedWith1
	SSXZW <- numobs * (covXZWwith0 + meanXZWwith1) # should the mean vector be squared (postmultiplied by its transpose)?
	
	# Compute a vector and b4, b5, b6
	a <- solve(SSMCReduced) %*% as.matrix(SSXZW)
	betaTemp <- betaRC[4:6] - (as.numeric(betaRC[7]) * a[5:7])
	betaTemp <- c(betaTemp, betaRC[7])
		
	# Compute Cov(Y, XZ) and regression coefficients of no-centering
	betaNC <- solve(varX[1:3,1:3], as.matrix(covY) - (t(varX[4:7, 1:3]) %*% as.matrix(betaTemp)))
	betaNC <- rbind(as.matrix(betaNC), as.matrix(betaTemp))
	covY <- rbind(covY, as.matrix((varX %*% betaNC)[4:7, 1]))

	# Aggregate the non-centering sufficient statistics (Just show how to do but not necessary)
	fullCov <- rbind(cbind(varX, covY), c(covY, varY))
	fullMean <- rbind(meanX, meanY)

	# Compute the intercept of no-centering
	intceptNC <- meanY - t(betaNC) %*% meanX

	# Compute SSNC
	betaNCWithIntcept <- rbind(intceptNC, betaNC)
	meanXwith1 <- rbind(1, meanX) #JG: redundant
	varXwith0 <- rbind(0, cbind(0, varX)) #JG: redundant
	SSNC <- numobs * (varXwith0 + (meanXwith1 %*% t(meanXwith1)))

	# Compute residual variance on non-centering
	resVarNC <- varY - (t(betaNCWithIntcept) %*% SSNC %*% betaNCWithIntcept)/numobs + meanY^2

	# Extract all varEst
	varEst <- lavaan::vcov(fit)
	
	# Check whether intercept are estimated
	targetcol <- paste(nameY, "~", 1, sep="")
	estimateIntcept <- targetcol %in% rownames(varEst)
	
	pvalue <- function(x) (1 - pnorm(abs(x))) * 2
	
	# Find the order to rearrange
	ord <- c(setdiff(1:3, modVar), modVar)
	ord <- c(ord, 7 - rev(ord))
	
	resultIntcept <- NULL
	resultSlope <- NULL
	if(estimateIntcept) {
		# Extract SE from residual centering
		targetcol <- c(targetcol, paste(nameY, "~", nameX, sep=""))
		varEstSlopeRC <- varEst[targetcol, targetcol]

		# Transform it to non-centering SE
		usedVar <-  as.numeric(resVarNC/resVarRC) * (varEstSlopeRC %*% SSRC %*% solve(SSNC))
		usedBeta <- betaNCWithIntcept
		if(sum(diag(usedVar) < 0) > 0) stop("This method does not work. The resulting calculation provided negative standard errors.") # JG: edited this error
		
		# Change the order of usedVar and usedBeta if the moderator variable is listed first
		usedVar <- usedVar[c(1, ord+1, 8), c(1, ord+1, 8)]
		usedBeta <- usedBeta[c(1, ord+1, 8)]
				
		# Find probe value
		val <- expand.grid(valProbe1, valProbe2)
		
		# Find simple intercept
		simpleIntcept <- usedBeta[1] + usedBeta[3] * val[,1] + usedBeta[4] * val[,2] + usedBeta[7] * val[,1] * val[,2]
		varIntcept <- usedVar[1, 1] + val[,1]^2 * usedVar[3, 3] + val[,2]^2 * usedVar[4, 4] + val[,1]^2 * val[,2]^2 * usedVar[7, 7] + 2 * val[,1] * usedVar[1, 3] + 2 * val[,2] * usedVar[1, 4] + 2 * val[,1] * val[,2] * usedVar[1, 7] + 2 * val[,1] * val[,2] * usedVar[3, 4] + 2 * val[,1]^2 * val[,2] * usedVar[3, 7] + 2* val[,1] * val[,2]^2 * usedVar[4, 7]
		zIntcept <- simpleIntcept/sqrt(varIntcept)
		pIntcept <- pvalue(zIntcept)
		resultIntcept <- cbind(val, simpleIntcept, sqrt(varIntcept), zIntcept, pIntcept)
		colnames(resultIntcept) <- c(nameX[modVar], "Intcept", "SE", "Wald", "p")
		
		# Find simple slope
		simpleSlope <- usedBeta[2] + usedBeta[5] * val[,1] + usedBeta[6] * val[,2] + usedBeta[8] * val[,1] * val[,2]
		varSlope <- usedVar[2, 2] + val[,1]^2 * usedVar[5, 5] + val[,2]^2 * usedVar[6, 6] + val[,1]^2 * val[,2]^2 * usedVar[8, 8] + 2 * val[,1] * usedVar[2, 5] + 2 * val[,2] * usedVar[2, 6] + 2 * val[,1] * val[,2] * usedVar[2, 8] + 2 * val[,1] * val[,2] * usedVar[5, 6] + 2 * val[,1]^2 * val[,2] * usedVar[5, 8] + 2 * val[,1] * val[,2]^2 * usedVar[6, 8]
		zSlope <- simpleSlope/sqrt(varSlope)
		pSlope <- round(pvalue(zSlope),6) # JG: rounded values
		resultSlope <- cbind(val, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "Slope", "SE", "Wald", "p")
	} else {
		targetcol <- paste(nameY, "~", nameX, sep="")
		varEstSlopeRC <- varEst[targetcol, targetcol]

		# Transform it to non-centering SE
		usedVar <-  as.numeric(resVarNC/resVarRC) * (varEstSlopeRC %*% SSRC[2:8, 2:8] %*% solve(SSNC[2:8, 2:8]))
		usedBeta <- betaNC
		if(sum(diag(usedVar) < 0) > 0) stop("This method does not work. The resulting calculation provided negative standard errors.") # JG: edited this error
		
		# Change the order of usedVar and usedBeta if the moderator variable is listed first
		usedVar <- usedVar[c(ord, 7), c(ord, 7)]
		usedBeta <- usedBeta[c(ord, 7)]

		# Find probe value
		val <- expand.grid(valProbe1, valProbe2)
		
		# Find simple slope
		simpleSlope <- usedBeta[1] + usedBeta[4] * val[,1] + usedBeta[5] * val[,2] + usedBeta[7] * val[,1] * val[,2]
		varSlope <- usedVar[1, 1] + val[,1]^2 * usedVar[4, 4] + val[,2]^2 * usedVar[5, 5] + val[,1]^2 * val[,2]^2 * usedVar[7, 7] + 2 * val[,1] * usedVar[1, 4] + 2 * val[,2] * usedVar[1, 5] + 2 * val[,1] * val[,2] * usedVar[1, 7] + 2 * val[,1] * val[,2] * usedVar[4, 5] + 2 * val[,1]^2 * val[,2] * usedVar[4, 7] + 2 * val[,1] * val[,2]^2 * usedVar[5, 7]
		zSlope <- simpleSlope/sqrt(varSlope)
		pSlope <- round(pvalue(zSlope),6) # JG: rounded values
		resultSlope <- cbind(val, simpleSlope, sqrt(varSlope), zSlope, pSlope)
		colnames(resultSlope) <- c(nameX[modVar], "Slope", "SE", "Wald", "p")
	}
	
	return(list(SimpleIntcept = resultIntcept, SimpleSlope = resultSlope))
}

# Find the expected value of the product of two normal variates
# m = the mean of each normal variate
# s = the covariance matrix of all variates
expect2NormProd <- function(m, s) {
	return(prod(m) + s[1, 2])
}

# Find the expected value of the product of three normal variates
# m = the mean of each normal variate
# s = the covariance matrix of all variates
expect3NormProd <- function(m, s) {
	return(prod(m) + m[3] * s[1, 2] + m[2] * s[1, 3] + m[1] * s[2, 3])
}

# Find the expected value of the product of four normal variates
# m = the mean of each normal variate
# s = the covariance matrix of all variates
expect4NormProd <- function(m, s) {
	first <- prod(m)
	com <- combn(1:4, 2)
	forSecond <- function(draw, meanval, covval, index) {
		draw2 <- setdiff(index, draw)
		prod(meanval[draw2]) * covval[draw[1], draw[2]]
	}
	second <- sum(apply(com, 2, forSecond, meanval=m, covval=s, index=1:4))

	com2 <- com[,1:3] #select only first three terms containing the first element only
	forThird <- function(draw, covval, index) {
		draw2 <- setdiff(index, draw)
		covval[draw[1], draw[2]] * covval[draw2[1], draw2[2]]
	}
	third <- sum(apply(com2, 2, forThird, covval=s, index=1:4))
	return(first + second + third)
}

# Find the expected value of the product of five normal variates
# m = the mean of each normal variate
# s = the covariance matrix of all variates
expect5NormProd <- function(m, s) {  
	first <- prod(m)
	com <- combn(1:5, 2) 
	forSecond <- function(draw, meanval, covval, index) {
		draw2 <- setdiff(index, draw)
		prod(meanval[draw2]) * covval[draw[1], draw[2]]
	}
	second <- sum(apply(com, 2, forSecond, meanval=m, covval=s, index=1:5))

	com2 <- combn(1:5, 4)
	forThirdOuter <- function(index, m, s, indexall) {
		targetMean <- m[setdiff(indexall, index)]
		cominner <- combn(index, 2)[,1:3] #select only first three terms containing the first element only
		forThirdInner <- function(draw, covval, index) {
			draw2 <- setdiff(index, draw)
			covval[draw[1], draw[2]] * covval[draw2[1], draw2[2]]
		}
		thirdInner <- targetMean * sum(apply(cominner, 2, forThirdInner, covval=s, index=index))
		return(thirdInner)
	}
	third <- sum(apply(com2, 2, forThirdOuter, m=m, s=s, indexall=1:5))
	return(first + second + third)
}

# Find the variance of the product of two normal variates
# m = the mean of each normal variate
# s = the covariance matrix of all variates
var2NormProd <- function(m, s) {
	first <- m[2]^2 * s[1, 1] + m[1]^2 * s[2, 2]
	second <- 2 * m[1] * m[2] * s[1, 2]
	third <- s[1, 1] * s[2, 2]
	fourth <- s[1, 2]^2
	return(first + second + third + fourth)
}

# Find the variance of the product of three normal variates
# m = the mean of each normal variate
# s = the covariance matrix of all variates
var3NormProd <- function(m, s) {
	com <- combn(1:3, 2)
	forFirst <- function(draw, meanval, covval, index) {
		# draw = 2, 3; draw2 = 1
		draw2 <- setdiff(index, draw)
		term1 <- meanval[draw[1]]^2 * meanval[draw[2]]^2 * covval[draw2, draw2]
		term2 <- 2 * meanval[draw2]^2 * meanval[draw[1]] * meanval[draw[2]] * covval[draw[1], draw[2]]
		term3 <- (meanval[draw2]^2 * covval[draw[1], draw[1]] * covval[draw[2], draw[2]]) + (meanval[draw2]^2 * covval[draw[1], draw[2]]^2)
		term4 <- 4 * meanval[draw[1]] * meanval[draw[2]] * covval[draw2, draw2] * covval[draw[1], draw[2]]
		term5 <- 6 * meanval[draw[1]] * meanval[draw[2]] * covval[draw2, draw[1]] * covval[draw2, draw[2]]
		term1 + term2 + term3 + term4 + term5
	}
	first <- sum(apply(com, 2, forFirst, meanval=m, covval=s, index=1:3))
	second <- prod(diag(s))
	third <- 2 * s[3, 3] * s[1, 2]^2 + 2 * s[2, 2] * s[1, 3]^2 + 2 * s[1, 1] * s[2, 3]^2
	fourth <- 8 * s[1, 2] * s[1, 3] * s[2, 3]
	return(first + second + third + fourth)
}

# plotProbe: plot the probing interaction result 
plotProbe <- function(object, xlim, xlab="Indepedent Variable", ylab="Dependent Variable", ...) {
	if(length(xlim) != 2) stop("The x-limit should be specified as a numeric vector with the length of 2.")
	
	# Extract simple slope
	slope <- object$SimpleSlope
	
	# Check whether the object is the two-way or three-way interaction result
	numInt <- 2
	if(ncol(slope) == 6) numInt <- 3
	estSlope <- slope[,ncol(slope) - 3]
	
	# Get whether the simple slope is significant. If so, the resulting lines will be shown as red. If not, the line will be black.
	estSlopeSig <- (slope[,ncol(slope)] < 0.05) + 1
	
	# Extract simple intercept. If the simple intercept is not provided, the intercept will be fixed as 0.
	estIntercept <- NULL
	if(!is.null(object$SimpleIntcept)) estIntercept <- object$SimpleIntcept[,ncol(slope) - 3]
	if(numInt == 2) {
		plotSingleProbe(estSlope, estIntercept, xlim=xlim, xlab=xlab, ylab=ylab, colLine=estSlopeSig, legendMain=colnames(slope)[1], legendVal=slope[,1], ...)
	} else if (numInt == 3) {
		# Three-way interaction; separate lines for the first moderator, separate graphs for the second moderator
		mod2 <- unique(slope[,2])
		mod1 <- unique(slope[,1])
		
		# Use multiple graphs in a figure
		if (length(mod2) == 2) {
			obj <- par(mfrow = c(1, 2))
		} else if (length(mod2) == 3) {
			obj <- par(mfrow = c(1, 3))
		} else if (length(mod2) > 3) {
			obj <- par(mfrow = c(2, ceiling(length(mod2)/2)))
		} else if (length(mod2) == 1) {
			# Intentionally leaving as blank
		} else {
			stop("Some errors occur")
		}
		for(i in 1:length(mod2)) {
			select <- slope[,2] == mod2[i]
			plotSingleProbe(estSlope[select], estIntercept[select], xlim=xlim, xlab=xlab, ylab=ylab, colLine=estSlopeSig[select], main=paste(colnames(slope)[2], "=", mod2[i]), legendMain=colnames(slope)[1], legendVal=mod1, ...)
		}
		if (length(mod2) > 1) 
			par(obj)
	} else {
		stop("Please make sure that the object argument is obtained from 'probe2wayMC', 'probe2wayRC', 'probe3wayMC', or 'probe3wayRC'.")
	}
}

# plotSingleProbe : plot the probing interaction result specific for only one moderator
# estSlope = slope of each line
# estIntercept = intercept of each line
# xlim = the minimum and maximum values of the independent variable (x-axis)
# xlab = the label for the independent variable
# ylab = the lable for the dependent variable
# main = the title of the graph
# colLine = the color of each line
# legendMain = the title of the legend
# legendVal = the description of each line representing in the plot
plotSingleProbe <- function(estSlope, estIntercept=NULL, xlim, xlab="Indepedent Variable", ylab="Dependent Variable", main=NULL, colLine="black", legendMain=NULL, legendVal=NULL, ...) {
	if(is.null(estIntercept)) estIntercept <- rep(0, length(estSlope))
	if(length(colLine) == 1) colLine <- rep(colLine, length(estSlope))
	lower <- estIntercept + (xlim[1] * estSlope)
	upper <- estIntercept + (xlim[2] * estSlope)
	ylim <- c(min(c(lower, upper)), max(c(lower, upper)))
	plot(cbind(xlim, ylim), xlim=xlim, ylim=ylim, type="n", xlab=xlab, ylab=ylab, main=main, ...)
	for(i in 1:length(estSlope)) {
		lines(cbind(xlim, c(lower[i], upper[i])), col = colLine[i], lwd=1.5, lty=i)
	}
	if(!is.null(legendVal)) {
		positionX <- 0.25
		if(all(estSlope > 0)) positionX <- 0.01
		if(all(estSlope < 0)) positionX <- 0.50
		legend(positionX * (xlim[2] - xlim[1]) + xlim[1], 0.99 * (ylim[2] - ylim[1]) + ylim[1], legendVal, col=colLine, lty=1:length(estSlope), title=legendMain)
	}
}
