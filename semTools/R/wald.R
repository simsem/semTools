
wald <- function(object, syntax) {
	model <- unlist( strsplit(syntax, "\n") )
	
    # remove comments starting with '#' or '!'
    model <- gsub("#.*","", model); model <- gsub("!.*","", model)

    # replace semicolons by newlines and split in lines again
    model <- gsub(";","\n", model); model <- unlist( strsplit(model, "\n") )

    # strip all white space
    model <- gsub("[[:space:]]+", "", model)

    # keep non-empty lines only
    idx <- which(nzchar(model))
    model <- model[idx]

	beta <- lavaan::coef(object)
	contrast <- matrix(0, length(model), length(beta))
	
    for(i in 1:length(model)) {
        rhs <- model[i]
        out <- NULL
		sign <- NULL
		if(substr(rhs, 1, 1) == "-") {
			sign <- "-"
			rhs <- substr(rhs, 2, nchar(rhs))
		} else {
			sign <- "+"
		}
		cont <- TRUE
		while(cont) {
			pos <- regexpr("[+-]", rhs)
			if(pos == -1) {
				out <- c(out, rhs)
				cont <- FALSE
			} else {
				out <- c(out, substr(rhs, 1, pos - 1))
				sign <- c(sign, substr(rhs, pos, pos))
				rhs <- substr(rhs, pos + 1, nchar(rhs))
			}
		}
		
		num <- rep(NA, length(out))
		vname <- rep(NA, length(out))
		for(j in seq_along(out)) {
			pos <- regexpr("[*]", out[j])
			tmp <- 1
			if(pos == -1) {
				vname[j] <- out[j]
			} else {
				tmp <- substr(out[j], 1, pos-1)
				vname[j] <- substr(out[j], pos + 1, nchar(out[j]))
			}
			if(is.character(tmp) && regexpr("[/^]", tmp) != -1) tmp <- eval(parse(text = tmp))
			if(is.character(tmp)) tmp <- as.numeric(tmp)
			num[j] <- tmp
			if(sign[j] == "-") num[j] <- -num[j]
		}
		posmatch <- match(vname, names(beta))
		if(any(is.na(posmatch))) {
			stop(paste("Unknown parameters:", paste(vname[is.na(posmatch)], collapse = ", ")))
		}
		contrast[i,posmatch] <- num
	}
	result <- waldContrast(object, contrast)
	print(round(result, 6))
	invisible(result)
}

waldContrast <- function(object, contrast) {
	beta <- lavaan::coef(object)
	acov <- lavaan::vcov(object)
	chisq <- t(contrast %*% beta) %*%  solve(contrast %*% as.matrix(acov) %*% t(contrast)) %*% (contrast %*% beta)
	df <- nrow(contrast)
	p <- pchisq(chisq, df, lower.tail=FALSE)
	c(chisq = chisq, df = df, p = p)
}