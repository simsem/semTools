imposeStart <- function(out, expr, silent = TRUE) {
	if(!is(out, "lavaan")) stop("The first argument of the function must be a lavaan output.")
	template2 <- template <- substitute(expr)
	template2$do.fit <- FALSE
	model <- eval(expr = template2, enclos = parent.frame())
	ptmodel <- lavaan::parTable(model)
	coefmodel <- lavaan::coef(model)
	coefout <- lavaan::coef(out)
	start <- coefout[match(names(coefmodel), names(coefout))]
	ptmodel$start[ptmodel$free != 0] <- start[ptmodel$free[ptmodel$free != 0]]
	ptmodel$est <- NULL
	ptmodel$se <- NULL
	if(!silent) {
		cat("########## Model with imposed starting values #########\n")
		print(ptmodel)
	}
	if("model" %in% names(template)) {
		template$model <- ptmodel
	} else {
		template[[2]] <- ptmodel
	}
	eval(expr = template, enclos = parent.frame())
}
