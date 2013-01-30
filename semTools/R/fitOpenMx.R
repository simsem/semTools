
saturateMx <- function(data, groupLab = NULL) {
	library(OpenMx)

	multipleGroup <- FALSE
	if(is.data.frame(data) && !is.null(groupLab) && groupLab %in% colnames(data)) multipleGroup <- TRUE
	if(is.list(data) && !is.data.frame(data)) multipleGroup <- TRUE
	if(multipleGroup) {
		if(is.data.frame(data)) {
			data.l <- split(data, data[,groupLab])
			data.l <- lapply(data.l, function(x) x[-ncol(x)])
			ngroups <- length(data.l)
		} else if(is.list(data)) {
			data.l <- data
			ngroups <- length(data.l)
		} else {
			stop("The data argument must be a data frame or a list of MxData objects")
		}
		temp <- mapply(saturateMxSingleGroup, data = data.l, title = paste0("group", 1:ngroups), groupnum = 1:ngroups, SIMPLIFY=FALSE)
		title <- "Multiple group Saturate Model"
		algebra <- mxAlgebra("", name="allobjective")
		groupnames <- paste0("group", 1:ngroups)
		groupnames <- paste0(groupnames, ".objective")
		groupnames <- lapply(groupnames, as.name)
		algebra@formula <- as.call(c(list(as.name("sum")), groupnames))
		objective <- mxAlgebraObjective("allobjective")
		Saturate <- mxModel(title, unlist(temp), algebra, objective)
	} else {
		Saturate <- saturateMxSingleGroup(data, title = "Saturate Model")
	}
	capture.output(fit <- mxRun(Saturate, suppressWarnings = FALSE))
	fit
}

saturateMxSingleGroup <- function(data, title = "Saturate Model", groupnum = NULL) {
	if(!is(data, "MxData")) {
		data <- mxData(
			observed=data,
			type="raw")
	}
	p <- ncol(data@observed)
	if(data@type == "raw") {
		categorical <- rep(FALSE, p)
		for(i in seq_len(p)) {
			categorical[i] <- "ordered" %in% class(data@observed[,i])
		}
		startMeans <- apply(data@observed, 2, function(x) mean(as.numeric(x), na.rm=TRUE))
		startVar <- apply(data@observed, 2, var, na.rm=TRUE)
	} else {
		categorical <- rep(FALSE, p)
		if(!all(is.na(data@means))) {
			startMeans <- data@means
		} else {
			startMeans <- rep(0, p)
		}
		startVar <- diag(data@observed)
	}
	startCor <- diag(p)

	# library(polycor)
	# for(i in seq_len(p)) {
		# for(j in seq_len(i - 1)) {
			# if(categorical[i] & categorical[j]) {
				# startCor[i, j] <- polychor(data@observed[,i], data@observed[,j])
			# } else if (categorical[i] & !categorical[j]) {
				# startCor[i, j] <- polyserial(data@observed[,j], data@observed[,i])
			# } else if (!categorical[i] & categorical[j]) {
				# startCor[i, j] <- polyserial(data@observed[,i], data@observed[,j])
			# } else {
				# startCor[i, j] <- cor(sapply(data@observed[,c(i,j)], as.numeric), use="pairwise.complete.obs")[2,1]
			# }
			# startCor[j, i] <- startCor[i, j]
		# }
	# }
	
	startVar[categorical] <- 1
	startMeans[categorical] <- 0
	startCov <- lavaan::cor2cov(startCor, sqrt(startVar))
	lab <- outer(1:p, 1:p, function(x, y) paste0("cov", x, y, "_", groupnum))
	lab2 <- outer(1:p, 1:p, function(x, y) paste0("cov", y, x, "_", groupnum))
	lab[upper.tri(lab)] <- lab2[upper.tri(lab2)]
	freeMean <- !categorical
	freeCov <- matrix(TRUE, p, p)
	diag(freeCov) <- !categorical
	if(any(categorical)) {
		labCategorical <- colnames(data@observed)[categorical]
		datCategorical <- data@observed[,categorical, drop=FALSE]
		numCat <- apply(datCategorical, 2, function(x) length(unique(x)))
		maxCat <- max(numCat)
		FUN <- function(x, tot) c(rep(TRUE, x), rep(FALSE, tot-x))
		freeThreshold <- sapply(numCat - 1, FUN, maxCat - 1)
		FUN2 <- function(x, tot) {
			x <- x[!is.na(x)]
			f <- table(x)/length(x)
			f <- cumsum(f)[-length(f)]
			f <- qnorm(f)
			c(f, rep(NA, tot - length(f)))
		}
		valueThreshold <- sapply(datCategorical, FUN2, maxCat - 1)
		T <- mxMatrix(
				type="Full",
				nrow=maxCat - 1,
				ncol=length(labCategorical),
				free=freeThreshold,
				values=valueThreshold,
				dimnames=list(c(), labCategorical),
				byrow=TRUE,
				name="thresh"
		)
		Saturate <- mxModel(title,
			data,
			# means
			mxMatrix(
				type="Full",
				nrow=1,
				ncol=p,
				values=startMeans,
				free=freeMean,
				labels=paste0("mean", 1:p, "_", groupnum),
				name="M"
			),
			# symmetric paths
			mxMatrix(
				type="Symm",
				nrow=p,
				ncol=p,
				values=startCov,
				free=freeCov,
				labels=lab,
				byrow=TRUE,
				name="S"
			),
			T,
			mxFIMLObjective(
				covariance="S",
				means="M",
				dimnames=colnames(data@observed),
				thresholds = "thresh"				
			)
		)
	} else {
		if(data@type == "raw") {
			obj <- mxFIMLObjective(
				covariance="S",
				means="M",
				dimnames=colnames(data@observed)
			)
			modelMean <- mxMatrix(
				type="Full",
				nrow=1,
				ncol=p,
				values=startMeans,
				free=freeMean,
				labels=paste0("mean", 1:p, "_", groupnum),
				name="M"
			)
		} else {
			if(!all(is.na(data@means))) {
				modelMean <- mxMatrix(
					type="Full",
					nrow=1,
					ncol=p,
					values=startMeans,
					free=freeMean,
					labels=paste0("mean", 1:p, "_", groupnum),
					name="M"
				)
				obj <- mxMLObjective(
					covariance="S",
					means="M",
					dimnames=colnames(data@observed)
				)
			} else {
				modelMean <- NULL
				obj <- mxMLObjective(
					covariance="S",
					dimnames=colnames(data@observed)
				)
			}
			
		}
		Saturate <- mxModel(title,
			data,
			# means
			modelMean,
			# symmetric paths
			mxMatrix(
				type="Symm",
				nrow=p,
				ncol=p,
				values=startCov,
				free=freeCov,
				labels=lab,
				byrow=TRUE,
				name="S"
			),
			obj
		)
	}
	Saturate
}

nullMx <- function(data, groupLab = NULL) {
	library(OpenMx)
	multipleGroup <- FALSE
	if(is.data.frame(data) && !is.null(groupLab) && groupLab %in% colnames(data)) multipleGroup <- TRUE
	if(is.list(data) && !is.data.frame(data)) multipleGroup <- TRUE
	if(multipleGroup) {
		if(is.data.frame(data)) {
			data.l <- split(data, data[,groupLab])
			data.l <- lapply(data.l, function(x) x[-ncol(x)])
			ngroups <- length(data.l)
		} else if(is.list(data)) {
			data.l <- data
			ngroups <- length(data.l)
		} else {
			stop("The data argument must be a data frame or a list of MxData objects")
		}
		temp <- mapply(nullMxSingleGroup, data = data.l, title = paste0("group", 1:ngroups), groupnum = 1:ngroups, SIMPLIFY=FALSE)
		title <- "Multiple group Null Model"
		algebra <- mxAlgebra("", name="allobjective")
		groupnames <- paste0("group", 1:ngroups)
		groupnames <- paste0(groupnames, ".objective")
		groupnames <- lapply(groupnames, as.name)
		algebra@formula <- as.call(c(list(as.name("sum")), groupnames))
		objective <- mxAlgebraObjective("allobjective")
		Null <- mxModel(title, unlist(temp), algebra, objective)
	} else {
		Null <- nullMxSingleGroup(data, title = "Null Model")
	}
	capture.output(fit <- mxRun(Null, suppressWarnings = FALSE))
	fit
}

nullMxSingleGroup <- function(data, title = "Null Model", groupnum = NULL) {
	if(!is(data, "MxData")) {
		data <- mxData(
			observed=data,
			type="raw")
	}
	p <- ncol(data@observed)
	if(data@type == "raw") {
		categorical <- rep(FALSE, p)
		for(i in seq_len(p)) {
			categorical[i] <- "ordered" %in% class(data@observed[,i])
		}
		startMeans <- apply(data@observed, 2, function(x) mean(as.numeric(x), na.rm=TRUE))
		startVar <- apply(data@observed, 2, var, na.rm=TRUE)
	} else {
		categorical <- rep(FALSE, p)
		if(!all(is.na(data@means))) {
			startMeans <- data@means
		} else {
			startMeans <- rep(0, p)
		}
		startVar <- diag(data@observed)
	}
	startVar[categorical] <- 1
	startMeans[categorical] <- 0
	lab <- paste0("var", 1:p, "_", groupnum)
	freeMean <- !categorical
	
	if(any(categorical)) {
		labCategorical <- colnames(data@observed)[categorical]
		datCategorical <- data@observed[,categorical, drop=FALSE]
		numCat <- apply(datCategorical, 2, function(x) length(unique(x)))
		maxCat <- max(numCat)
		FUN <- function(x, tot) c(rep(TRUE, x), rep(FALSE, tot-x))
		freeThreshold <- sapply(numCat - 1, FUN, maxCat - 1)
		FUN2 <- function(x, tot) {
			f <- table(x)/length(x)
			f <- cumsum(f)[-length(f)]
			f <- qnorm(f)
			c(f, rep(NA, tot - length(f)))
		}
		valueThreshold <- sapply(datCategorical, FUN2, maxCat - 1)
		T <- mxMatrix(
				type="Full",
				nrow=maxCat - 1,
				ncol=length(labCategorical),
				free=freeThreshold,
				values=valueThreshold,
				dimnames=list(c(), labCategorical),
				byrow=TRUE,
				name="thresh"
		)
		NullModel <- mxModel(title,
			data,
			# means
			mxMatrix(
				type="Full",
				nrow=1,
				ncol=p,
				values=startMeans,
				free=freeMean,
				labels=paste0("mean", 1:p, "_", groupnum),
				name="M"
			),
			# symmetric paths
			mxMatrix(
				type="Diag",
				nrow=p,
				ncol=p,
				values=startVar,
				free=freeMean,
				labels=lab,
				byrow=TRUE,
				name="S"
			),
			T,
			mxFIMLObjective(
				covariance="S",
				means="M",
				dimnames=colnames(data@observed),
				thresholds = "thresh"				
			)
		)
	} else {
		if(data@type == "raw") {
			obj <- mxFIMLObjective(
				covariance="S",
				means="M",
				dimnames=colnames(data@observed)
			)
			modelMean <- mxMatrix(
				type="Full",
				nrow=1,
				ncol=p,
				values=startMeans,
				free=freeMean,
				labels=paste0("mean", 1:p, "_", groupnum),
				name="M"
			)
		} else {
			if(!all(is.na(data@means))) {
				modelMean <- mxMatrix(
					type="Full",
					nrow=1,
					ncol=p,
					values=startMeans,
					free=freeMean,
					labels=paste0("mean", 1:p, "_", groupnum),
					name="M"
				)
				obj <- mxMLObjective(
					covariance="S",
					means="M",
					dimnames=colnames(data@observed)
				)
			} else {
				modelMean <- NULL
				obj <- mxMLObjective(
					covariance="S",
					dimnames=colnames(data@observed)
				)
			}
			
		}
		NullModel <- mxModel(title,
			data,
			# means
			modelMean,
			# symmetric paths
			mxMatrix(
				type="Diag",
				nrow=p,
				ncol=p,
				values=startVar,
				free=freeMean,
				labels=lab,
				byrow=TRUE,
				name="S"
			),
			obj
		)
	}
	NullModel
}

checkConvergence <- function(object) {
	(object@output$status[[1]] %in% c(0,1)) & (object@output$status[[2]] == 0)
}

fitMeasuresMx <- function(object, fit.measures="all") {
	library(OpenMx)
	
	if(length(object@submodels) > 1) {
		varnames <- lapply(object@submodels, function(x) {
			out <- x@objective@dims
			if(any(is.na(out))) out <- x@manifestVars
			out
		})
		dat <- lapply(object@submodels, slot, "data")
		FUN <- function(x, var) {
			if(x@type == "raw") {
				x@observed <- x@observed[,intersect(var, colnames(x@observed))]	
			} 
			x
		}
		dat <- mapply(FUN, x=dat, var=varnames)
	} else {
		varnames <- object@objective@dims
		if(any(is.na(varnames))) varnames <- object@manifestVars
		dat <- object@data
		dat@observed <- dat@observed[,intersect(varnames, colnames(dat@observed))]
	}

	if(length(object@output) == 0) {
		stop("The target model has not been estimated yet.")
	}

    if(!checkConvergence(object)) {
        warning("The target model may be not convergent.")
    }
	
    if("all" %in% fit.measures) {
       class.flag <- TRUE
    } else {
       class.flag <- FALSE
    }

    # reference: Muthen technical Appendix 5

    # collect info from the lavaan slots
	if(length(object@submodels) > 1) {
		N <- sum(sapply(dat, slot, "numObs"))
	} else {
		N <- dat@numObs
	}
    #q <- length(vnames(object@ParTable, "ov.x"))
    #p <- nvar - q
    npar <- length(object@output$estimate)
    #fx <- object@Fit@fx
    #fx.group <- object@Fit@fx.group
	
    multigroup    <- length(object@submodels) > 1
    G <- length(object@submodels) # number of groups
	if(G == 0) G <- 1 # Correct when there is no submodel

	if(multigroup) {
		if(is(object@submodels[[1]]@objective, "MxRAMObjective")) {
			meanstructure <- !all(is.na(object@submodels[[1]]@objective@M)) # Only ML objective
			categorical   <- !all(is.na(object@submodels[[1]]@objective@thresholds)) # Only ML objective	
		} else {
			meanstructure <- !all(is.na(object@submodels[[1]]@objective@means)) # Only ML objective
			categorical   <- !all(is.na(object@submodels[[1]]@objective@thresholds)) # Only ML objective	
		}
	} else {
		if(is(object@objective, "MxRAMObjective")) {
			meanstructure <- !all(is.na(object@objective@M)) # Only ML objective
			categorical   <- !all(is.na(object@objective@thresholds)) # Only ML objective		
		} else {
			meanstructure <- !all(is.na(object@objective@means)) # Only ML objective
			categorical   <- !all(is.na(object@objective@thresholds)) # Only ML objective
		}
	}
	# define 'sets' of fit measures:

    # basic chi-square test
    fit.chisq <- c("chisq", "df", "pvalue")

    # basline model
    fit.baseline <- c("baseline.chisq", "baseline.df", "baseline.pvalue")

    # cfi/tli
    fit.cfi.tli <- c("cfi", "tli")

    # more incremental fit indices
    fit.incremental <- c("cfi", "tli", "nnfi", "rfi", "nfi", "pnfi",
                         "ifi", "rni")
    
    # likelihood based measures
    fit.logl <- c("logl", "unrestricted.logl", "npar", "aic", "bic",
                  "ntotal", "bic2")

    # rmsea
    fit.rmsea <- c("rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue")

    # srmr
    if(categorical) {
        fit.srmr <- character(0L)
        fit.srmr2 <- character(0L)
    } else {
        fit.srmr <- c("srmr")
        fit.srmr2 <- c("rmr", "rmr_nomean", "srmr", "srmr_nomean")
    }

    # various
    fit.other <- c("cn_05","cn_01","mfi")
    if(!categorical && G == 1) {
        fit.other <- c(fit.other, "ecvi")
    }

    # lower case
    fit.measures <- tolower(fit.measures)

    # select 'default' fit measures
	
	# Check ML Categorical in OpenMx
    if(length(fit.measures) == 1L) {
        if(fit.measures == "default") {
			fit.measures <- c(fit.chisq, fit.baseline, 
							  fit.cfi.tli, fit.logl, 
							  fit.rmsea, fit.srmr, "saturate.status", "null.status")
        } else if(fit.measures == "all") {
			fit.measures <- c(fit.chisq, fit.baseline, fit.incremental, 
							  fit.logl, fit.rmsea, fit.srmr2, fit.other, "saturate.status", "null.status")
        }
    }
    
    # main container
    indices <- list()

	if("logl" %in% fit.measures ||
	   "npar" %in% fit.measures ||
	   "aic" %in% fit.measures ||
	   "bic" %in% fit.measures) {
		
		logl.H0 <- object@output$Minus2LogLikelihood * (-1/2)

		if("logl" %in% fit.measures) {
			indices["logl"] <- logl.H0
		}

		# Number of free parameters
		if("npar" %in% fit.measures) {
			indices["npar"] <- npar
		}

		# AIC
		AIC <-  -2*logl.H0 + 2*npar
		if("aic" %in% fit.measures) {
			indices["aic"] <- AIC
		}

		# BIC
		if("bic" %in% fit.measures) {
			BIC <- -2*logl.H0 + npar*log(N)
			indices["bic"] <- BIC

			# add sample-size adjusted bic
			N.star <- (N + 2) / 24
			BIC2 <- -2*logl.H0 + npar*log(N.star)
			indices["bic2"] <- BIC2
		}
	}
	if(multigroup) {
		defVars <- lapply(object@submodels, findDefVars)
		defVars <- do.call(c, defVars)
	} else {
		defVars <- findDefVars(object)	
	}
	if(length(defVars) > 0) {
		out <- unlist(indices[intersect(fit.measures, names(indices))])
        return(out)
	}
	
	objectSat <- saturateMx(dat)
	objectNull <- nullMx(dat)
	
	if("unrestricted.logl" %in% fit.measures) {
		logl.H1 <- objectSat@output$Minus2LogLikelihood * (-1/2)
		indices["unrestricted.logl"] <- logl.H1
    }


	if(length(objectSat@output) == 0) {
		stop("The saturated model has not been estimated yet.")
	}

	if(length(objectNull@output) == 0) {
		stop("The null model has not been estimated yet.")
	}
	
	if(length(objectNull@output) == 0) {
		stop("The null model has not been estimated yet.")
	}

	if(length(object@output) == 0) {
		stop("The model has not been estimated yet.")
	}
	
    # has the model converged?

	if(!checkConvergence(objectSat)) {
        warning("The saturated model may be not convergent.")
    }

	if(!checkConvergence(objectNull)) {
        warning("The null model may be not convergent.")
    }
	X2 <- object@output$Minus2LogLikelihood - objectSat@output$Minus2LogLikelihood
	df <- length(objectSat@output$estimate) - length(object@output$estimate)

	indices["saturate.status"] <- objectSat@output$status[[1]]
	indices["null.status"] <- objectNull@output$status[[1]]
	
	if(objectSat@output$status[[2]] != 0) indices["saturate.status"] <- -1
	if(objectNull@output$status[[2]] != 0) indices["null.status"] <- -1
	
    # Chi-square value estimated model (H0)
    if(any("chisq" %in% fit.measures)) {
	indices["chisq"] <- X2
    }
    if(any("df" %in% fit.measures)) {
        indices["df"] <- df
    }
	
    if(any(c("pvalue") %in% fit.measures)) {
        indices["pvalue"] <- pchisq(X2, df, lower.tail = FALSE)
    }


    if(any(c("cfi", "tli", 
             "nnfi", "pnfi", 
             "rfi", "nfi", 
             "ifi", "rni", 
             "baseline.chisq", 
             "baseline.pvalue") %in% fit.measures)) {
        
        # call explicitly independence model
        # this is not strictly needed for ML, but it is for
        # GLS and WLS
        # and MLM and MLR to get the scaling factor(s)!
        #if(estimator == "ML") {
        #    if(object@SampleStats@missing.flag) {
        #        do.fit <- TRUE
        #    } else {
        #        do.fit <- FALSE
        #    }
        #} else {
        #    do.fit <- TRUE
        #}

        #OV.X <- character(0L)
        #if(object@Options$mimic == "Mplus") 
        #    OV.X <- vnames(object@ParTable, type="ov.x")

        #indep.syntax <- 
        #    syntax.independence.model(ov.names   = object@Data@ov.names,
        #                              ov.names.x = OV.X,
        #                              sample.cov = object@SampleStats@cov)
        #fit.indep <- update(object, model = indep.syntax, 
        #                    se = "none", do.fit=TRUE, 
        #                    constraints = "",
        #                    verbose = FALSE, warn = FALSE)
        #OCALL <- as.list(object@call); OCALL$env <- NULL; OCALL[[1]] <- NULL
        #NCALL <- list(model = indep.syntax, se = "none", do.fit = TRUE, 
        #              constraints = "", verbose = FALSE, warn = FALSE)
        #CALL  <- modifyList(OCALL, NCALL)
        #cat("DEBUG!\n"); print(as.list(object@call$env)); cat("*******\n")
        #fit.indep <- do.call("lavaan", args=CALL, envir=object@call$env)
        #fit.indep <- do.call("lavaan", args=CALL)
        
		X2.null <- objectNull@output$Minus2LogLikelihood - objectSat@output$Minus2LogLikelihood
		df.null <- length(objectSat@output$estimate) - length(objectNull@output$estimate)
	
        # check for NAs
        if(is.na(X2) || is.na(df) || is.na(X2.null) || is.na(df.null)) {
            indices[fit.incremental] <- as.numeric(NA)
        } else {
            if("baseline.chisq" %in% fit.measures) {
                indices["baseline.chisq"] <- X2.null
            }
            if("baseline.df" %in% fit.measures) {
                indices["baseline.df"] <- df.null
            }
            if("baseline.pvalue" %in% fit.measures) {
                indices["baseline.pvalue"] <- pchisq(X2.null, df.null, lower.tail = FALSE)
            }

            # CFI - comparative fit index (Bentler, 1990) 
            if("cfi" %in% fit.measures) {
                t1 <- max( c(X2 - df, 0) )
                t2 <- max( c(X2 - df, X2.null - df.null, 0) )
                if(t1 == 0 && t2 == 0) {
                    indices["cfi"] <- 1
                } else {
                    indices["cfi"] <- 1 - t1/t2
                }
            }

            # TLI - Tucker-Lewis index (Tucker & Lewis, 1973)
            # same as
            # NNFI - nonnormed fit index (NNFI, Bentler & Bonett, 1980)
            if("tli" %in% fit.measures || "nnfi" %in% fit.measures) {
                if(df > 0) {
                    t1 <- X2.null/df.null - X2/df
                    t2 <- X2.null/df.null - 1 
                    # note: TLI original formula was in terms of fx/df, not X2/df
                    # then, t1 <- fx_0/df.null - fx/df
                    #       t2 <- fx_0/df.null - 1/N (or N-1 for wishart)
                    if(t1 < 0 && t2 < 0) {
                        TLI <- 1
                    } else {
                        TLI <- t1/t2
                    }
                } else {
                   TLI <- 1
                }
                indices["tli"] <- indices["nnfi"] <- TLI
            }
    
            # RFI - relative fit index (Bollen, 1986; Joreskog & Sorbom 1993)
            if("rfi" %in% fit.measures) {
                if(df > 0) {
                    t1 <- X2.null/df.null - X2/df
                    t2 <- X2.null/df.null
                    if(t1 < 0 || t2 < 0) {
                        RLI <- 1
                    } else {
                        RLI <- t1/t2
                    }
                } else {
                   RLI <- 1
                }
                indices["rfi"] <- RLI
            }

            # NFI - normed fit index (Bentler & Bonett, 1980)
            if("nfi" %in% fit.measures) {
                t1 <- X2.null - X2
                t2 <- X2.null
                NFI <- t1/t2
                indices["nfi"] <- NFI
            }

            # PNFI - Parsimony normed fit index (James, Mulaik & Brett, 1982)
            if("pnfi" %in% fit.measures) {
                t1 <- X2.null - X2
                t2 <- X2.null
                PNFI <- (df/df.null) * t1/t2
                indices["pnfi"] <- PNFI
            }

            # IFI - incremental fit index (Bollen, 1989; Joreskog & Sorbom, 1993)
            if("ifi" %in% fit.measures) {
                t1 <- X2.null - X2
                t2 <- X2.null - df
                if(t2 < 0) {
                    IFI <- 1
                } else {
                    IFI <- t1/t2
                }
                indices["ifi"] <- IFI
            }

            # RNI - relative noncentrality index (McDonald & Marsh, 1990)
            if("rni" %in% fit.measures) {
                t1 <- X2 - df
                t2 <- X2.null - df.null
                if(t1 < 0 || t2 < 0) {
                    RNI <- 1
                } else {
                    RNI <- 1 - t1/t2
                }
                indices["rni"] <- RNI
            }
        }
    }

    N.RMSEA <- max(N, X2*4) # FIXME: good strategy??
    if(any("rmsea" %in% fit.measures)) {
        # RMSEA
        if(is.na(X2) || is.na(df)) {
            RMSEA <- as.numeric(NA)
        } else if(df > 0) {
			GG <- 0
			RMSEA <- sqrt( max( c((X2/N)/df - 1/(N-GG), 0) ) ) * sqrt(G)
        } else {
            RMSEA <- 0
        }
        indices["rmsea"] <- RMSEA
    }

    if("rmsea.ci.lower" %in% fit.measures) {
        lower.lambda <- function(lambda) {
            (pchisq(X2, df=df, ncp=lambda) - 0.95)
        }
        if(is.na(X2) || is.na(df)) {
            indices["rmsea.ci.lower"] <- NA
        } else if(df < 1 || lower.lambda(0) < 0.0) {
            indices["rmsea.ci.lower"] <- 0
        } else {
            lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=X2)$root)
            if(inherits(lambda.l, "try-error")) { lambda.l <- NA }
			GG <- 0
			indices["rmsea.ci.lower"] <- 
				sqrt( lambda.l/((N-GG)*df) ) * sqrt(G)
        }
    }

    if("rmsea.ci.upper" %in% fit.measures) {
        upper.lambda <- function(lambda) {
            (pchisq(X2, df=df, ncp=lambda) - 0.05)
        }
        if(is.na(X2) || is.na(df)) {
            indices["rmsea.ci.upper"] <- NA
        } else if(df < 1 || upper.lambda(N.RMSEA) > 0 || upper.lambda(0) < 0) {
            indices["rmsea.ci.upper"] <- 0
        } else {
            lambda.u <- try(uniroot(f=upper.lambda, lower=0, upper=N.RMSEA)$root)
            if(inherits(lambda.u, "try-error")) { lambda.u <- NA }
			GG <- 0
			indices["rmsea.ci.upper"] <- 
				sqrt( lambda.u/((N-GG)*df) ) * sqrt(G)
        }
    }

    if("rmsea.pvalue" %in% fit.measures) {
        if(is.na(X2) || is.na(df)) {
            indices["rmsea.pvalue"] <- as.numeric(NA)
        } else if(df > 0) {
			GG <- 0
			ncp <- (N-GG)*df*0.05^2/G
			indices["rmsea.pvalue"] <- 
				1 - pchisq(X2, df=df, ncp=ncp)
        } else {
            indices["rmsea.pvalue"] <- 1
        }
    }

    if(any(c("rmr","srmr") %in% fit.measures)) {
        # RMR and SRMR
        rmr.group <- numeric(G)
        rmr_nomean.group <- numeric(G)
        srmr.group <- numeric(G)
        srmr_nomean.group <- numeric(G)
        for(g in 1:G) {
			if(G > 1) {
				if(is(objectSat@submodels[[g]]@objective, "MxRAMObjective")) {
					impliedSat <- getImpliedStatRAM(objectSat@submodels[[g]])
				} else {
					impliedSat <- getImpliedStatML(objectSat@submodels[[g]])
				}
			} else {
				if(is(objectSat@objective, "MxRAMObjective")) {
					impliedSat <- getImpliedStatRAM(objectSat)
				} else {
					impliedSat <- getImpliedStatML(objectSat)
				}
			}
			S <- impliedSat[[2]]
			M <- matrix(impliedSat[[1]], ncol=1)

            nvar <- ncol(S)

            # estimated
			if(G > 1) {
				if(is(object@submodels[[g]]@objective, "MxRAMObjective")) {
					implied <- getImpliedStatRAM(object@submodels[[g]])
				} else {
					implied <- getImpliedStatML(object@submodels[[g]])
				}
			} else {
				if(is(object@objective, "MxRAMObjective")) {
					implied <- getImpliedStatRAM(object)
				} else {
					implied <- getImpliedStatML(object)
				}
			}
			Sigma.hat <- implied[[2]]
			Mu.hat <- matrix(implied[[1]], ncol=1)

            # standardized residual covariance matrix
            # this is the Hu and Bentler definition, not the Bollen one!
            sqrt.d <- 1/sqrt(diag(S))
            D <- diag(sqrt.d, ncol=length(sqrt.d))
            R <- D %*% (S - Sigma.hat) %*% D
            RR <- (S - Sigma.hat) # not standardized

            # this is what the Mplus documentation suggest, 
            # but is not what is used!
            #sqrt.d2 <- 1/sqrt(diag(Sigma.hat))
            #D2 <- diag(sqrt.d2, ncol=length(sqrt.d2))
            #R <- D %*% S %*% D   - D2 %*% Sigma.hat %*% D2

            if(meanstructure) {
                # standardized residual mean vector
                R.mean <- D %*% (M - Mu.hat)
                RR.mean <- (M - Mu.hat) # not standardized
                #R.mean <-  D %*% M - D2 %*% Mu.hat
                e <- nvar*(nvar+1)/2 + nvar
                srmr.group[g] <- sqrt( (sum(R[lower.tri(R, diag=TRUE)]^2) +
                                        sum(R.mean^2))/ e )
                rmr.group[g] <- sqrt( (sum(RR[lower.tri(RR, diag=TRUE)]^2) +
                                       sum(RR.mean^2))/ e )
                e <- nvar*(nvar+1)/2
                srmr_nomean.group[g] <- sqrt( sum(R[lower.tri(R, diag=TRUE)]^2) / e )
                rmr_nomean.group[g] <- sqrt( sum(RR[lower.tri(RR, diag=TRUE)]^2) / e )
            } else {
                e <- nvar*(nvar+1)/2
                srmr_nomean.group[g] <- srmr.group[g] <- sqrt( sum(R[lower.tri(R, diag=TRUE)]^2) / e )
                rmr_nomean.group[g] <- rmr.group[g] <- sqrt( sum(RR[lower.tri(RR, diag=TRUE)]^2) / e )
                
            }
        }
        
        if(G > 1) {
            ## FIXME: get the scaling right
			nEachGroup <- as.matrix(t(sapply(dat, slot, "numObs")))
            SRMR <- as.numeric( (nEachGroup %*% srmr.group) / N )
            SRMR_NOMEAN <- as.numeric( (nEachGroup %*% srmr_nomean.group) / N )
            RMR <- as.numeric( (nEachGroup %*% rmr.group) / N )
            RMR_NOMEAN <- as.numeric( (nEachGroup %*% rmr_nomean.group) / N )
        } else {
            SRMR <- srmr.group[1]
            SRMR_NOMEAN <- srmr_nomean.group[1]
            RMR <- rmr.group[1]
            RMR_NOMEAN <- rmr_nomean.group[1]
        }

        indices["srmr"] <- SRMR
        indices["srmr_nomean"] <- SRMR_NOMEAN
        indices["rmr"] <- RMR
        indices["rmr_nomean"] <- RMR_NOMEAN
    }

    if(any(c("cn_05", "cn_01") %in% fit.measures)) {
        CN_05 <- qchisq(p=0.95, df=df)/(X2/N) + 1
        CN_01 <- qchisq(p=0.99, df=df)/(X2/N) + 1
        indices["cn_05"] <- CN_05
        indices["cn_01"] <- CN_01
    }

    # MFI - McDonald Fit Index (McDonald, 1989)
    if("mfi" %in% fit.measures) { 
        #MFI <- exp(-0.5 * (X2 - df)/(N-1)) # Hu & Bentler 1998 Table 1
        MFI <- exp(-0.5 * (X2 - df)/N)
        indices["mfi"] <- MFI
    }

    # ECVI - cross-validation index (Brown & Cudeck, 1989)
    # not defined for multiple groups and/or models with meanstructures
    if("ecvi" %in% fit.measures) {
        if(G > 1 || meanstructure) {
            ECVI <- as.numeric(NA)
        } else {
            ECVI <- X2/N + (2*npar)/N
        }
        indices["ecvi"] <- ECVI
    }

    if("ntotal" %in% fit.measures) {
        indices["ntotal"] <- N
    }

    # do we have everything that we requested?
    idx.missing <- which(is.na(match(fit.measures, names(indices))))
    if(length(idx.missing) > 0L) {
        cat("lavaan WARNING: some requested fit measure(s) are not available for this model:\n")
        print( fit.measures[ idx.missing ] )
        cat("\n")
    }
    
    out <- unlist(indices[fit.measures])
	
    if(length(out) > 0L) {
        return(out)
    } else {
        return( invisible(numeric(0)) )
    }
}

findDefVars <- function(object) {
	mat <- lapply(object@matrices, slot, "labels")
	defvars <- sapply(mat, function(x) x[apply(x, c(1,2), OpenMx::imxIsDefinitionVariable)])
	Reduce("c", defvars)
}

getImpliedStatML <- function(xxxobjectxxx, xxxcovdatatxxx = NULL) {
	xxxmatxxx <- xxxobjectxxx@matrices
	xxxmatnamexxx <- names(xxxmatxxx)
	xxxmatvalxxx <- lapply(xxxmatxxx, slot, "values")
	for(i in seq_along(xxxmatnamexxx)) {
		assign(xxxmatnamexxx[i], xxxmatvalxxx[[i]])
	}
	if(!is.null(xxxcovdatatxxx)) {
		xxxmatlabxxx <- lapply(xxxmatxxx, slot, "labels")
		xxxdefvarsxxx <- lapply(xxxmatlabxxx, function(x) apply(x, c(1,2), OpenMx::imxIsDefinitionVariable))
		for(i in seq_along(xxxmatnamexxx)) {
			if(any(xxxdefvarsxxx[[i]])) {
				xxxtempxxx <- get(xxxmatnamexxx[i])
				for(j in seq_len(length(xxxdefvarsxxx[[i]]))) {
					if(xxxdefvarsxxx[[i]][j]) {
						xxxtempnamexxx <- gsub("data.", "", xxxmatlabxxx[[i]][j])
						xxxtempxxx[j] <- xxxcovdatatxxx[xxxtempnamexxx]
					}
				}
				assign(xxxmatnamexxx[i], xxxtempxxx)
			}
		}
	}
	xxxalgebraxxx <- xxxobjectxxx@algebras
	xxxalgebranamexxx <- names(xxxalgebraxxx)
	xxxalgebraformulaxxx <- lapply(xxxalgebraxxx, slot, "formula")
	for(i in seq_along(xxxalgebranamexxx)) {
		assign(xxxalgebranamexxx[i], eval(xxxalgebraformulaxxx[[i]]))
	}
	
	xxximpliedCovxxx <- get(xxxobjectxxx@objective@covariance)
	
	if(is.na(xxxobjectxxx@objective@means)) {
		xxximpliedMeanxxx <- rep(0, nrow(xxximpliedCovxxx))
	} else {
		xxximpliedMeanxxx <- get(xxxobjectxxx@objective@means)
	}
	
	if(is.na(xxxobjectxxx@objective@thresholds)) {
		xxximpliedThresholdxxx <- NA
	} else {
		xxximpliedThresholdxxx <- get(xxxobjectxxx@objective@thresholds)
	}
	list(xxximpliedMeanxxx, xxximpliedCovxxx, xxximpliedThresholdxxx)
}

getImpliedStatRAM <- function(object) {
	A <- object@matrices$A@values
	I <- diag(nrow(A))
	S <- object@matrices$S@values
	F <- object@matrices$F@values
	Z <- solve(I - A)
	impliedCov <- F %*% Z %*% S %*% t(Z) %*% t(F)
	if(!is.null(object@matrices$M)) {
		M <- object@matrices$M@values
		impliedMean <- t(F %*% Z %*% t(M))
	} else {
		impliedMean <- rep(0, nrow(impliedCov))
	}
	list(impliedMean, impliedCov)
}

standardizeMx <- function(object, free = TRUE) {
	objectOrig <- object
	multigroup <- length(object@submodels) > 0
	if(multigroup) {
		object@submodels <- lapply(object@submodels, standardizeMxSingleGroup)
	} else {
		object <- standardizeMxSingleGroup(object)	
	}
	vectorizeMx(object, free=free)
}

standardizeMxSingleGroup <- function(object) {
	if(!is(object@objective, "MxRAMObjective")) stop("The standardizeMx function is available for the MxRAMObjective only.")
	A <- object@matrices$A@values
	I <- diag(nrow(A))
	S <- object@matrices$S@values
	F <- object@matrices$F@values
	Z <- solve(I - A)
	impliedCov <- Z %*% S %*% t(Z)
	ImpliedSd <- diag(sqrt(diag(impliedCov)))
	ImpliedInvSd <- solve(ImpliedSd)
	object@matrices$S@values <- ImpliedInvSd %*% S %*% ImpliedInvSd
	object@matrices$A@values <- ImpliedInvSd %*% A %*% ImpliedSd
	if(!is.null(object@matrices$M)) {
		M <- object@matrices$M@values
		object@matrices$M@values <- M %*% ImpliedInvSd
	} 
	return(object)
}


vectorizeMx <- function(object, free = TRUE) {
	multigroup <- length(object@submodels) > 0
	if(multigroup) {
		object <- object@submodels
	} else {
		object <- list(object)	
	}
	result <- NULL
	for(i in seq_along(object)) {
		name <- ""
		if(multigroup) name <- paste0(object[[i]]@name, ".")
		mat <- object[[i]]@matrices
		for(j in seq_along(mat)) {
			tempname <- paste0(name, mat[[j]]@name)
			lab <- mat[[j]]@labels
			tempfree <- as.vector(mat[[j]]@free)
			madeLab <- paste0(tempname, "[", row(lab), ",", col(lab), "]")
			lab <- as.vector(lab)
			madeLab[!is.na(lab)] <- lab[!is.na(lab)]
			if(!free) tempfree <- rep(TRUE, length(tempfree))
			temp <- mat[[j]]@values[tempfree]
			names(temp) <- madeLab[tempfree]
			result <- c(result, temp)
		}
	}
	
	result[!duplicated(names(result))]
}
