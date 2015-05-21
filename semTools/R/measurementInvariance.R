measurementInvariance <- measurementinvariance <- function(..., std.lv = FALSE,
    strict=FALSE, quiet=FALSE, fit.measures = "default", method = "satorra.bentler.2001") {

    # check for a group.equal argument in ...
    dotdotdot <- list(...)
    if(!is.null(dotdotdot$group.equal))
        stop("lavaan ERROR: group.equal argument should not be used")

    res <- list()
    # base-line model: configural invariance
	
	configural <- dotdotdot
	configural$group.equal <- ""
	template <- do.call("cfa", configural)
	pttemplate <- lavaan::partable(template)
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
		res$fit.configural <- refit(pttemplate, template)
	} else {
		res$fit.configural <- template
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
		res$fit.loadings <- refit(pttemplate, template)
	} else {
		loadings <- dotdotdot
		loadings$group.equal <- c("loadings")
		res$fit.loadings <- do.call("cfa", loadings)
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
		res$fit.intercepts <- refit(pttemplate, template)
	} else {
		intercepts <- dotdotdot
		intercepts$group.equal <- c("loadings", "intercepts")
		res$fit.intercepts <- do.call("cfa", intercepts)
	}
	
    if(strict) {
		if(std.lv) {
			findresiduals <- which(pttemplate$op == "~~" & pttemplate$lhs %in% varnames & pttemplate$rhs == pttemplate$lhs & pttemplate$free != 0 & pttemplate$group == 1)
			for(i in findresiduals) {
				pttemplate <- constrainParTable(pttemplate, pttemplate$lhs[i], "~~", pttemplate$rhs[i], 1:ngroups)
			}
			res$fit.residuals <- refit(pttemplate, template)
			for(i in facnames) {
				pttemplate <- fixParTable(pttemplate, i, "~1", "", 1:ngroups, 0)
			}
			res$fit.means <- refit(pttemplate, template)
		} else {
			# fix loadings + intercepts + residuals
			residuals <- dotdotdot
			residuals$group.equal <- c("loadings", "intercepts", "residuals")
			res$fit.residuals <- do.call("cfa", residuals)

			# fix loadings + residuals + intercepts + means
			means <- dotdotdot
			means$group.equal <- c("loadings", "intercepts", "residuals", "means")
			res$fit.means <- do.call("cfa", means)
		}
    } else {
		if(std.lv) {
			for(i in facnames) {
				pttemplate <- fixParTable(pttemplate, i, "~1", "", 1:ngroups, 0)
			}
			res$fit.means <- refit(pttemplate, template)
		} else {
			# fix loadings + intercepts + means
			means <- dotdotdot
			means$group.equal <- c("loadings", "intercepts", "means")
			res$fit.means <- do.call("cfa", means)
		}
    }

	if(!quiet) {
        printInvarianceResult(res, fit.measures, method)
    }
	
    invisible(res)
}

printInvarianceResult <- function(FIT, fit.measures, method) {
	# compare models
	NAMES <- names(FIT); names(FIT) <- NULL
	TABLE <- do.call("lavTestLRT", c(FIT, list(model.names = NAMES,
									 method = method)))

	if(length(fit.measures) == 1L && fit.measures == "default") {
		# scaled test statistic?
		if(length(FIT[[1]]@Fit@test) > 1L) {
			fit.measures <- c("cfi.scaled", "rmsea.scaled")
		} else {
			fit.measures <- c("cfi", "rmsea")
		}
	}

	# add some fit measures
	if(length(fit.measures)) {

		FM <- lapply(FIT, lavaan::fitMeasures, fit.measures)
		FM.table1 <- sapply(fit.measures, function(x) sapply(FM, "[[", x))
		if(length(FM) == 1L) {
			FM.table1 <- rbind( rep(as.numeric(NA), length(fit.measures)),
								FM.table1 )
		}
		if(length(FM) > 1L) {
			FM.table2 <- rbind(as.numeric(NA),
							   abs(apply(FM.table1, 2, diff)))
			colnames(FM.table2) <- paste(colnames(FM.table2), ".delta", sep="")
			FM.TABLE <- as.data.frame(cbind(FM.table1, FM.table2))
		} else {
			FM.TABLE <- as.data.frame(FM.table1)
		} 
		rownames(FM.TABLE) <- rownames(TABLE)
		class(FM.TABLE) <- c("lavaan.data.frame", "data.frame")
	}
	cat("\n")
	cat("Measurement invariance models:\n\n")
	cat(paste(paste("Model", seq_along(FIT), ":", NAMES), collapse = "\n"))
	cat("\n\n")

	print(TABLE)
	if(length(fit.measures)) {
		cat("\n\n")
		cat("Fit measures:\n\n")
		print(FM.TABLE)
		cat("\n")
	}
}
