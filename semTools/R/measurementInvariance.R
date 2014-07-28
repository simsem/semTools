measurementInvariance <- measurementinvariance <- function(...,
    strict=FALSE, quiet=FALSE) {

    # check for a group.equal argument in ...
    dotdotdot <- list(...)
    if(!is.null(dotdotdot$group.equal))
        stop("lavaan ERROR: group.equal argument should not be used")

    res <- list()
    # base-line model: configural invariance
	
	configural <- dotdotdot
	configural$group.equal <- ""
    res$fit.configural <- do.call("cfa", configural)
	
    # fix loadings across groups
	loadings <- dotdotdot
	loadings$group.equal <- c("loadings")
    res$fit.loadings <- do.call("cfa", loadings)

    # fix loadings + intercepts across groups
	intercepts <- dotdotdot
	intercepts$group.equal <- c("loadings", "intercepts")
    res$fit.intercepts <- do.call("cfa", intercepts)

    if(strict) {
        # fix loadings + intercepts + residuals
		residuals <- dotdotdot
		residuals$group.equal <- c("loadings", "intercepts", "residuals")
		res$fit.residuals <- do.call("cfa", residuals)

        # fix loadings + residuals + intercepts + means
		means <- dotdotdot
		means$group.equal <- c("loadings", "intercepts", "residuals", "means")
		res$fit.means <- do.call("cfa", means)
	
    } else {
        # fix loadings + intercepts + means
		means <- dotdotdot
		means$group.equal <- c("loadings", "intercepts", "means")
		res$fit.means <- do.call("cfa", means)
    }

    if(!quiet) {
        cat("\nMeasurement invariance tests:\n")
        cat("\nModel 1: configural invariance:\n")
        printFitLine(res$fit.configural)

        cat("\nModel 2: weak invariance (equal loadings):\n")
        printFitLine(res$fit.loadings)

        cat("\n[Model 1 versus model 2]\n")
        difftest(res$fit.configural, res$fit.loadings)

        cat("\nModel 3: strong invariance (equal loadings + intercepts):\n")
        printFitLine(res$fit.intercepts)
        cat("\n[Model 1 versus model 3]\n")
        difftest(res$fit.configural, res$fit.intercepts)
        cat("\n[Model 2 versus model 3]\n")
        difftest(res$fit.loadings, res$fit.intercepts)


        if(strict) {
            cat("\nModel 4: strict invariance (equal loadings + intercepts + residuals):\n")
            printFitLine(res$fit.residuals)
            cat("\n[Model 1 versus model 4]\n")
            difftest(res$fit.configural, res$fit.residuals)
            cat("\n[Model 3 versus model 4]\n")
            difftest(res$fit.intercepts, res$fit.residuals)
  
            cat("\nModel 5: equal loadings + intercepts + residuals + means:\n")
            printFitLine(res$fit.means,horizontal=TRUE)
            cat("\n[Model 1 versus model 5]\n")
            difftest(res$fit.configural, res$fit.means)
            cat("\n[Model 4 versus model 5]\n")
            difftest(res$fit.residuals, res$fit.means)
        } else {
            cat("\nModel 4: equal loadings + intercepts + means:\n")
            printFitLine(res$fit.means)
            cat("\n[Model 1 versus model 4]\n")
            difftest(res$fit.configural, res$fit.means)
            cat("\n[Model 3 versus model 4]\n")
            difftest(res$fit.intercepts, res$fit.means)
        }
    }
    invisible(res)
}


printFitLine <- function(object, horizontal=TRUE) {

    # which `tests' do we have?
    scaled <- FALSE
    TESTS <- unlist(lapply(object@Fit@test, "[", "test"))
    if(any(c("satorra.bentler", "yuan.bentler", "scaled.shifted") %in% TESTS)) {
        scaled <- TRUE
    }

    if(!scaled) {
        out <- fitMeasures(object, c("chisq", "df", "pvalue",
                                      "cfi", "rmsea", "bic"))
    } else {
        out <- fitMeasures(object, c("chisq.scaled", "df.scaled",
                                      "pvalue.scaled",
                                      "cfi.scaled", "rmsea.scaled", "bic"))
        names(out) <- c("chisq.scaled", "df", "pvalue",
                        "cfi.scaled", "rmsea.scaled", "bic")
    }
    
    print(out)
}

difftest <- function(model1, model2) {
	if(model1@Fit@test[[1]]$df > model2@Fit@test[[1]]$df) {
        fit0 <- model1
        fit1 <- model2
    } else {
        fit0 <- model2
        fit1 <- model1
    }
	d <- unlist(lavTestLRT(fit1, fit0)[2,5:7])
	i0 <- fitMeasures(fit0)
	i1 <- fitMeasures(fit1)
	names(d) <- c("delta.chisq", "delta.df", "delta.p.value")
	if("cfi" %in% names(i0) & "cfi" %in% names(i1)) {
		temp <- i1["cfi"] - i0["cfi"]
		names(temp) <- NULL
		d <- c(d, "delta.cfi" = temp)
	}
	if("cfi.scaled" %in% names(i0) & "cfi.scaled" %in% names(i1)) {
		temp <- i1["cfi.scaled"] - i0["cfi.scaled"]
		names(temp) <- NULL
		d <- c(d, "delta.cfi.scaled" = temp)
	}
	print.lavaan.vector(d)
    invisible(d)
}

print.lavaan.vector <- function(x, ..., nd=3) {
    y <- unclass(x)
    #if(!is.null(names(x))) {
    #    names(y) <- abbreviate(names(x), minlength = nd + 3)
    #}
    print( round(y, nd), ... )
    invisible(x)
}
