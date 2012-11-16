## Title: Model-implied statistics
## Author: Sunthud Pornprasertmanit <psunthud@ku.edu>
## Description: Find the model-implied statistics
##----------------------------------------------------------------------------##

impliedFactorCov <- function(object) {
	param <- inspect(object, "coef")
	ngroup <- object@Data@ngroups
	name <- names(param)
	be <- param[name == "beta"]
	ps <- param[name == "psi"]
	result <- list()
	if (length(be) == ngroup) {
		for(i in 1:ngroup) {
			nf <- nrow(be[[i]])
			result[[i]] <- solve(diag(nf) - be[[i]]) %*% ps[[i]] %*% t(solve(diag(nf) - be[[i]]))
		}
	} else if (length(be) == 0) {
		result <- ps
	} else {
		stop("The number of beta matrix and the number of groups are not equal.")
	}
	if(ngroup == 1) {
		result <- result[[1]]
	} else {
		names(result) <- object@Data@group.label
	}
	result
}

impliedFactorMean <- function(object) {
	if(!object@Options$meanstructure) stop("This model does not estimate the mean structure.")
	param <- inspect(object, "coef")
	ngroup <- object@Data@ngroups
	name <- names(param)
	be <- param[name == "beta"]
	al <- param[name == "alpha"]
	result <- list()
	if (length(be) == ngroup) {
		for(i in 1:ngroup) {
			nf <- nrow(be[[i]])
			result[[i]] <- solve(diag(nf) - be[[i]]) %*% al[[i]]
		}
	} else if (length(be) == 0) {
		result <- al
	} else {
		stop("The number of beta matrix and the number of groups are not equal.")
	}
	if(ngroup == 1) {
		result <- result[[1]]
	} else {
		names(result) <- object@Data@group.label
	}
	result
}

impliedFactorStat <- function(object) {
	result <- list()
	if(object@Options$meanstructure) {
		result$mean <- impliedFactorMean(object)
	}
	result$cov <- impliedFactorCov(object)
	result
}

