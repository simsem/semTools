## Title: Coefficient omega calculation
## Author: Sunthud Pornprasertmanit <psunthud@ku.edu>
## Description: Find the coefficient omega of each factor
##----------------------------------------------------------------------------##

reliability <- function(object) {
	param <- inspect(object, "coef")
	ngroup <- object@Data@ngroups
	name <- names(param)
	if(!is.null(param$beta)) warning("Some factors are endogenous. The reliability values of endogenous variables are not trustworthy because the residual latent varaince (instead of total variance) is used.")
	ly <- param[name == "lambda"]
	ps <- param[name == "psi"]
	te <- param[name == "theta"]
	result <- list()
	for(i in 1:ngroup) {
		common <- (apply(ly[[i]], 2, sum)^2) * diag(ps[[i]])
		error <- rep(NA, length(common))
		for(j in 1:length(error)) {
			index <- which(ly[[i]][,j] != 0)
			error[j] <- sum(te[[i]][index, index])
		}
		result[[i]] <- common/(common + error)
	}
	if(ngroup == 1) {
		result <- result[[1]]
	} else {
		names(result) <- object@Data@group.label
	}
	result
}
