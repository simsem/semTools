measurementInvariance <- measurementinvariance <- function(..., std.lv = FALSE,
    strict=FALSE, quiet=FALSE) {

    # check for a group.equal argument in ...
    dotdotdot <- list(...)
    if(!is.null(dotdotdot$group.equal))
        stop("lavaan ERROR: group.equal argument should not be used")

    res <- list()
    # base-line model: configural invariance
	
	configural <- dotdotdot
	configural$group.equal <- ""
	template <- do.call("cfa", configural)
	pttemplate <- partable(template)
	varnames <- unique(pttemplate$rhs[pttemplate$op == "=~"])
	facnames <- unique(pttemplate$lhs[(pttemplate$op == "=~") & (pttemplate$rhs %in% varnames)])
	ngroups <- max(pttemplate$group)
	if(ngroups <= 1) stop("Well, the number of groups is 1. Measurement invariance across 'groups' cannot be done.")

	if(std.lv) {
		for(i in facnames) {
			pttemplate <- fixParTable(pttemplate, i, "~~", i, 1:ngroups, 1)
		}
		fixloadings <- which(pttemplate$op == "=~" & pttemplate$free == 0)
		for(i in fixloadings) {
			pttemplate <- freeParTable(pttemplate, pttemplate$lhs[i], "=~", pttemplate$rhs[i], pttemplate$group[i])
		}
		res$configural <- refit(pttemplate, template)
	} else {
		res$configural <- template
	}
	
    # fix loadings across groups
	if(std.lv) {
		findloadings <- which(pttemplate$op == "=~" & pttemplate$free != 0 & pttemplate$group == 1)
		for(i in findloadings) {
			pttemplate <- constrainParTable(pttemplate, pttemplate$lhs[i], "=~", pttemplate$rhs[i], 1:ngroups)
		}
		for(i in facnames) {
			pttemplate <- freeParTable(pttemplate, i, "~~", i, 2:ngroups)
		}		
		res$metric <- refit(pttemplate, template)
	} else {
		loadings <- dotdotdot
		loadings$group.equal <- c("loadings")
		res$metric <- do.call("cfa", loadings)
	}
	
    # fix loadings + intercepts across groups
	if(std.lv) {
		findintcepts <- which(pttemplate$op == "~1" & pttemplate$lhs %in% varnames & pttemplate$free != 0 & pttemplate$group == 1)
		for(i in findintcepts) {
			pttemplate <- constrainParTable(pttemplate, pttemplate$lhs[i], "~1", "", 1:ngroups)
		}
		for(i in facnames) {
			pttemplate <- freeParTable(pttemplate, i, "~1", "", 2:ngroups)
		}	
		res$scalar <- refit(pttemplate, template)
	} else {
		intercepts <- dotdotdot
		intercepts$group.equal <- c("loadings", "intercepts")
		res$scalar <- do.call("cfa", intercepts)
	}
	
    if(strict) {
		if(std.lv) {
			findresiduals <- which(pttemplate$op == "~~" & pttemplate$lhs %in% varnames & pttemplate$rhs == pttemplate$lhs & pttemplate$free != 0 & pttemplate$group == 1)
			for(i in findresiduals) {
				pttemplate <- constrainParTable(pttemplate, pttemplate$lhs[i], "~~", pttemplate$rhs[i], 1:ngroups)
			}
			res$strict <- refit(pttemplate, template)
			for(i in facnames) {
				pttemplate <- fixParTable(pttemplate, i, "~1", "", 1:ngroups, 0)
			}
			res$means <- refit(pttemplate, template)
		} else {
			# fix loadings + intercepts + residuals
			residuals <- dotdotdot
			residuals$group.equal <- c("loadings", "intercepts", "residuals")
			res$strict <- do.call("cfa", residuals)

			# fix loadings + residuals + intercepts + means
			means <- dotdotdot
			means$group.equal <- c("loadings", "intercepts", "residuals", "means")
			res$means <- do.call("cfa", means)
		}
    } else {
		if(std.lv) {
			for(i in facnames) {
				pttemplate <- fixParTable(pttemplate, i, "~1", "", 1:ngroups, 0)
			}
			res$means <- refit(pttemplate, template)
		} else {
			# fix loadings + intercepts + means
			means <- dotdotdot
			means$group.equal <- c("loadings", "intercepts", "means")
			res$means <- do.call("cfa", means)
		}
    }

    if(!quiet) {
        cat("\nMeasurement invariance tests:\n")
        cat("\nModel 1: configural invariance:\n")
        printFitLine(res$configural)

        cat("\nModel 2: weak invariance (equal loadings):\n")
        printFitLine(res$metric)

        cat("\n[Model 1 versus model 2]\n")
        difftest(res$configural, res$metric)

        cat("\nModel 3: strong invariance (equal loadings + intercepts):\n")
        printFitLine(res$scalar)
        cat("\n[Model 1 versus model 3]\n")
        difftest(res$configural, res$scalar)
        cat("\n[Model 2 versus model 3]\n")
        difftest(res$metric, res$scalar)


        if(strict) {
            cat("\nModel 4: strict invariance (equal loadings + intercepts + residuals):\n")
            printFitLine(res$strict)
            cat("\n[Model 1 versus model 4]\n")
            difftest(res$configural, res$strict)
            cat("\n[Model 3 versus model 4]\n")
            difftest(res$scalar, res$strict)
  
            cat("\nModel 5: equal loadings + intercepts + residuals + means:\n")
            printFitLine(res$means,horizontal=TRUE)
            cat("\n[Model 1 versus model 5]\n")
            difftest(res$configural, res$means)
            cat("\n[Model 4 versus model 5]\n")
            difftest(res$strict, res$means)
        } else {
            cat("\nModel 4: equal loadings + intercepts + means:\n")
            printFitLine(res$means)
            cat("\n[Model 1 versus model 4]\n")
            difftest(res$configural, res$means)
            cat("\n[Model 3 versus model 4]\n")
            difftest(res$scalar, res$means)
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
