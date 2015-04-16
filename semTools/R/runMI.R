##  Functon to impute missing data, run Lavaan on each one 
##  input: data frames of raw data with missing data, model specification (lavaan script), number of imputations wanted.)
##  Output: lavaanStar object which filled with the appropriate information 
##  Alexander Schoemann, Patrick Miller, Mijke Rhemtulla, Sunthud Pornprasertmanit, Alexander Robitzsch, Mauricio Garnier Villarreal
##  Last modified 5/25/2012

##Currently outputs a list of parameter estimates, standard errors, fit indices and fraction missing information

cfa.mi <- function(model, data, m, miArgs=list(), miPackage="Amelia", chi="all", seed=12345, ...) {
	runMI(model=model, data=data, m=m, miArgs=miArgs, chi=chi, miPackage=miPackage, seed=seed, fun="cfa", ...)
}

sem.mi <- function(model, data, m, miArgs=list(), miPackage="Amelia", chi="all", seed=12345, ...) {
	runMI(model=model, data=data, m=m, miArgs=miArgs, chi=chi, miPackage=miPackage, seed=seed, fun="sem", ...)
}

growth.mi <- function(model, data, m, miArgs=list(), miPackage="Amelia", chi="all", seed=12345, ...) {
	runMI(model=model, data=data, m=m, miArgs=miArgs, chi=chi, miPackage=miPackage, seed=seed, fun="growth", ...)
}

lavaan.mi <- function(model, data, m, miArgs=list(), miPackage="Amelia", chi="all", seed=12345, ...) {
	runMI(model=model, data=data, m=m, miArgs=miArgs, chi=chi, miPackage=miPackage, seed=seed, fun="lavaan", ...)
}


runMI <- function(model, data, m, miArgs=list(), chi="all", miPackage="Amelia", seed=12345, fun, ...) 
{
	set.seed(seed)
	chi <- tolower(chi)
	if(!(chi %in% c("none", "mplus", "mr", "lmrr", "all"))) {
		stop("The chi argument should be one of the followings only: 'none, 'mr', 'lmrr', 'mplus', or 'all'.")
	}
	
	imputed.data <- is.list(data) & (!is.data.frame(data))
	imputed.l <- NULL
	if (!imputed.data){		
		if( ( miPackage!="Amelia" )  &  ( miPackage !="mice")  ) { 
			stop("Currently runMI only supports imputation by Amelia or mice") 
		}
		if(miPackage=="Amelia"){
			imputed.l<-imputeMissingAmelia(data,m, miArgs)
		} else if(miPackage=="mice"){
			imputed.l<-imputeMissingMice(data,m, miArgs)
		}
	} else { 
		imputed.l <- data 
		m <- length( data )
		data <- data[[1]]
	}
	out <- list(model=model, data=imputed.l[[1]], se="none", do.fit=FALSE)
	out <- c(out, list(...))
	template <- do.call(fun, out)

    imputed.results.l <- suppressWarnings(lapply(imputed.l, runlavaanMI, syntax=model, fun=fun, ...))
	converged.l <- sapply(imputed.results.l, function(x) x@Fit@converged)
	
	coefAll <- sapply(imputed.results.l, function(x) x@Fit@est)
	seAll <- sapply(imputed.results.l, function(x) x@Fit@se)
	partableImp <- lavaan::partable(imputed.results.l[[1]])
	posVar <- (partableImp$op == "~~") & (partableImp$lhs == partableImp$rhs)
	convergedtemp <- converged.l
	properSE <- apply(seAll, 2, function(x) all(x >= 0) & !(all(x == 0)))
	properVariance <- apply(coefAll[posVar, ,drop=FALSE], 2, function(x) all(x >= 0))
	converged.l <- converged.l & properSE & properVariance
	if(sum(converged.l) < 2) {
		tab <- cbind(convergedtemp, properSE, properVariance, converged.l)
		colnames(tab) <- c("1. Convergence", "2. Proper SE", "3. Proper Variance Estimate", "Used for pooling")
		print(tab)
		stop("Please increase the number of imputations. The number of convergent replications is less than or equal to 1. See the details above.")
	}
	
	mOriginal <- m
	m <- sum(converged.l)
	convergenceRate <- m/mOriginal
	imputed.results.l <- imputed.results.l[converged.l]
	
	coefs <- sapply(imputed.results.l, function(x) x@Fit@est)
	se <- sapply(imputed.results.l, function(x) x@Fit@se)
	Sigma.hat <- lapply(imputed.results.l, function(object) object@Fit@Sigma.hat)
    Mu.hat <- lapply(imputed.results.l, function(object) object@Fit@Mu.hat)
	meanSigmaHat <- list()
	meanMuHat <- list()
	for(g in seq_len(template@SampleStats@ngroups)) {
		tempSigma <- lapply(Sigma.hat, "[[", g)
		meanSigmaHat[[g]] <- Reduce("+", tempSigma)/m
		tempMu <- lapply(Mu.hat, "[[", g)
		meanMuHat[[g]] <- Reduce("+", tempMu)/m		
	}
	template@Fit@Sigma.hat <- meanSigmaHat
	template@Fit@Mu.hat <- meanMuHat
	comb.results <- miPoolVector(t(coefs),t(se), m)
	template@Fit@est <- comb.results$coef
	template@Fit@se <- comb.results$se
	template@Fit@x <- comb.results$coef[comb.results$se != 0]
	template@Model <- imposeGLIST(template@Model, comb.results$coef, template@ParTable)
	
	fmi.results <- cbind(lavaan::parameterEstimates(template)[,1:3], group=template@ParTable$group, fmi1 = comb.results[[3]], fmi2 = comb.results[[4]])

	fit <- imputed.results.l[[1]]@Fit@test
	df <- fit[[1]]$df
  if (df == 0) chi <- "none" # for saturated models, no model fit available
	chi1 <- sapply(imputed.results.l, function(x) x@Fit@test[[1]]$stat)

	if(any(template@Data@ov$type == "ordered") | (length(fit) > 1)) {
		if(chi=="all") chi <- "lmrr"
		if(chi %in% c("mplus", "mr")) {
			stop("The 'mplus' or 'mr' method for pooling chi-square values is not available with categorical variables.")
		}
	}
	
	chiScaled1 <- NULL
	dfScaled <- NULL
	if(length(fit) > 1) {
		chiScaled1 <- sapply(imputed.results.l, function(x) x@Fit@test[[2]]$stat)
		dfScaled <- fit[[2]]$df
	}
	
	if(template@SampleStats@ngroups == 1) {
		fit[[1]]$stat.group <- mean(sapply(imputed.results.l, function(x) x@Fit@test[[1]]$stat.group))
	} else {
		fit[[1]]$stat.group <- rowMeans(sapply(imputed.results.l, function(x) x@Fit@test[[1]]$stat.group))
	}
	nullModel <- lavaan::lav_partable_independence(template)
    null.results <- suppressWarnings(lapply(imputed.l, runlavaanMI, syntax=nullModel, fun=fun, ...))

	convergedNull.l <- sapply(null.results, function(x) x@Fit@converged)
	seNullAll <- sapply(null.results, function(x) x@Fit@se)
	convergedNull.l <- convergedNull.l & apply(seNullAll, 2, function(x) all(x >= 0))
	
	mNull <- sum(convergedNull.l)
	convergenceNullRate <- mNull/mOriginal
	null.results <- null.results[convergedNull.l]
	
	chiNull <- sapply(null.results, function(x) x@Fit@test[[1]]$stat)
	dfNull <- null.results[[1]]@Fit@test[[1]]$df
	
	chiNullScaled1 <- NULL
	dfNullScaled <- NULL
	if(length(fit) > 1) {
		chiNullScaled1 <- sapply(null.results, function(x) x@Fit@test[[2]]$stat)
		dfNullScaled <- null.results[[1]]@Fit@test[[2]]$df
	}
	outNull <- list(model=nullModel, data=imputed.l[[1]], se="none", do.fit=FALSE)
	outNull <- c(outNull, list(...))
	templateNull <- suppressWarnings(do.call(fun, outNull))
	
	coefsNull <- sapply(null.results, function(x) x@Fit@est)
	seNull <- sapply(null.results, function(x) x@Fit@se)
	
	comb.results.null <- miPoolVector(t(coefsNull),t(seNull), mNull)
	fitNull <- null.results[[1]]@Fit@test

	
	lmrr <- NULL
	lmrrNull <- NULL
	mr <- NULL
	mrNull <- NULL
	mplus <- NULL
	mplusNull <- NULL
	lmrrScaled <- NULL
	lmrrScaledNull <- NULL
	
	if(chi %in% c("lmrr", "all")){
		lmrr <- lmrrPooledChi(chi1, df)
		lmrrNull <- lmrrPooledChi(chiNull, dfNull)
		fit[[1]]$stat <- as.numeric(lmrr[1] * lmrr[2])
		fit[[1]]$pvalue <- as.numeric(lmrr[4])
		fitNull[[1]]$stat <- as.numeric(lmrrNull[1] * lmrrNull[2])
		fitNull[[1]]$pvalue <- as.numeric(lmrrNull[4])
		
		if(!is.null(chiScaled1)) {
			lmrrScaled <- lmrrPooledChi(chiScaled1, dfScaled)
			lmrrScaledNull <- lmrrPooledChi(chiNullScaled1, dfNullScaled)
			fit[[2]] <- imputed.results.l[[1]]@Fit@test[[2]]
			fit[[2]]$stat <- as.numeric(lmrrScaled[1] * lmrrScaled[2])
			fit[[2]]$pvalue <- as.numeric(lmrrScaled[4])
			fitNull[[2]] <- null.results[[1]]@Fit@test[[2]]
			fitNull[[2]]$stat <- as.numeric(lmrrScaledNull[1] * lmrrScaledNull[2])
			fitNull[[2]]$pvalue <- as.numeric(lmrrScaledNull[4])
			template@Options$estimator <- imputed.results.l[[1]]@Options$estimator 
			template@Options$test <- imputed.results.l[[1]]@Options$test 
			templateNull@Options$estimator <- null.results[[1]]@Options$estimator 
			templateNull@Options$test <- null.results[[1]]@Options$test 
		}
	}
	
	if(chi %in% c("mplus", "mr", "all")){
		mrplus <- mrplusPooledChi(template, imputed.l[converged.l], chi1, df, coef=comb.results$coef, m=m, fun=fun, ...)
		mrplusNull <- mrplusPooledChi(templateNull, imputed.l[convergedNull.l], chiNull, dfNull, coef=comb.results.null$coef, m=mNull, fun=fun, par.sat=satPartable(template), ...)
		
		if(chi %in% c("mr", "all")){
			mr <- mrPooledChi(mrplus[1], mrplus[2], mrplus[3], mrplus[4])
			mrNull <- mrPooledChi(mrplusNull[1], mrplusNull[2], mrplusNull[3], mrplusNull[4])
			
			fit[[1]]$stat <- as.numeric(mr[1] * mr[2])
			fit[[1]]$pvalue <- as.numeric(mr[4])
			fitNull[[1]]$stat <- as.numeric(mrNull[1] * mrNull[2])
			fitNull[[1]]$pvalue <- as.numeric(mrNull[4])
		}
		if(chi %in% c("mplus", "all")){
			mplus <- mplusPooledChi(mrplus[1], mrplus[3], mrplus[4])
			mplusNull <- mplusPooledChi(mrplusNull[1], mrplusNull[3], mrplusNull[4])
			fit[[1]]$stat <- as.numeric(mplus[1])
			fit[[1]]$pvalue <- as.numeric(mplus[3])
			fitNull[[1]]$stat <- as.numeric(mplusNull[1])
			fitNull[[1]]$pvalue <- as.numeric(mplusNull[3])
		}
	}
	template@Fit@test <- fit
	template@Fit@npar <- imputed.results.l[[1]]@Fit@npar
	template@Options <- imputed.results.l[[1]]@Options
	templateNull@Fit@test <- fitNull
	result <- as(template, "lavaanStar")
    ## HACK! YR
    templateNull@Fit@converged <- TRUE ### ! to trick fitMeasures
    ## 
	notused <- capture.output(fitVec <- suppressWarnings(lavaan::fitMeasures(templateNull)))
	name <- names(fitVec)
	fitVec <- as.vector(fitVec)
	names(fitVec) <- name
	result@nullfit <- fitVec
	
	result@Fit@iterations <- as.integer(m)
	result@Fit@converged <- TRUE
	
	summaryImputed <- list()
	summaryImputed[[1]] <- c("target model" = convergenceRate, "null model" = convergenceNullRate)
	summaryImputed[[2]] <- fmi.results
	summaryImputed[[3]] <- list(lmrr = lmrr, mr = mr, mplus = mplus)
	summaryImputed[[4]] <- list(lmrr = lmrrNull, mr = mrNull, mplus = mplusNull)
	if(!is.null(lmrrScaled)) {
		summaryImputed[[5]] <- list(lmrr = lmrrScaled)
		summaryImputed[[6]] <- list(lmrr = lmrrScaledNull)
		names(summaryImputed) <- c("convergenceRate", "fractionMissing", "targetFit", "nullFit", "targetFit.scaled", "nullFit.scaled")
	} else {
		names(summaryImputed) <- c("convergenceRate", "fractionMissing", "targetFit", "nullFit")
	}
	result@imputed <- summaryImputed

	return(result)
}
  
#Convenient function to run lavaan models and get results out. For easy use with lapply
runlavaanMI <- function(MIdata, syntax, fun, ...) {
	out <- list(model=syntax, data=MIdata)
	out <- c(out, list(...))
	fit <- NULL
	try(fit <- do.call(fun, out), silent=TRUE)
    return(fit)
}
	
#Conveniance function to run impuations on data and only return list of data
imputeMissingAmelia <- function(data,m, miArgs){
  # pull out only the imputations
  library(Amelia)
  out <- c(list(Amelia::amelia, x = data, m = m, p2s=0), miArgs)
  temp.am <- eval(as.call(out))
  return(temp.am$imputations)

} # end imputeMissingAmelia

imputeMissingMice <- function(data,m, miArgs){
  # pull out only the imputations
  library(mice)
  out <- c(list(mice::mice, data=data, m = m, diagnostics=FALSE, printFlag=FALSE), miArgs)
  temp.mice <- eval(as.call(out))
  temp.mice.imp <- NULL
  for(i in 1:m) {
	temp.mice.imp[[i]] <- mice::complete(x=temp.mice, action=i, include=FALSE) 
  }
  return(temp.mice.imp)
} # end imputeMissingAmelia



 
# miPoolVector
# Function -- simsem package
# Pool MI results that providing in matrix or vector formats
# Argument:
#	MI.param: 	Coefficients matrix (row = imputation, col = parameters)
#	MI.se:		Standard errors matrix (row = imputation, col = parameters)
#	imps:		Number of imputations
# Return:
#	coef:	Parameter estimates
#	se:		Standard error combining the between and within variances
#	FMI.1:	Fraction missing?
#	FMI.2:	Fraction missing?
# Author: 	Mijke Rhumtella
#			Alex Schoemann
#			Sunthud Pornprasertmanit (University of Kansas; psunthud@ku.edu)
# Date Modified: February 8, 2012

miPoolVector <- function(MI.param, MI.se, imps) {
   #compute parameter estimates
  Estimates <- colMeans(MI.param)

#compute between-imputation variance: variance of parameter estimates
  Bm <- apply(MI.param,2,var)

#compute within-imputation variance: average of squared estimated SEs 
#Um <- colSums(MI.se^2/m)
  Um <- apply(MI.se^2,2,mean)

#Total variance
#Tm <- Um + (Bm)*((1+m)/m+1)
#compute total variance: sum of between- and within- variance with correction
  TV <- Um + ((imps+1)/imps)*Bm

#compute correction factor for fraction of missing info
  nu <- (imps-1)*((((1+1/imps)*Bm)/TV)^-2)

#compute 2 estimates of fraction of missing information
  FMI.1 <- 1-(Um/TV)
  FMI.2 <- 1- ((nu+1)*Um)/((nu+3)*TV)
  FMI.2[is.nan(FMI.2)] <- 0
  FMI<-rbind(FMI.1,FMI.2)

#Get rid of estimates from fixed variables
#fixedParam <- Bm==0

#Estimates <- subset(Estimates, !fixedParam)
#TV <- subset(TV, !fixedParam)
#FMI.1 <- subset(FMI.1, !fixedParam)
#FMI.2 <- subset(FMI.2, !fixedParam)
SE <- sqrt(TV)
MI.res<-list(Estimates,SE,FMI.1,FMI.2)
names(MI.res)<-c('coef','se','FMI.1','FMI.2')
#compute chi-square proportion (is this useful?)
#(MI.fit.mat$chisq.p is a placeholder for however we'll index the p-value of chi square)
#chisq <- sum(MI.fit.mat$chisq.pval<.05)/m
  return(MI.res)
}
#Examples:
#param <- matrix(c(0.7, 0.1, 0.5,
#					0.75, 0.12, 0.54,
#					0.66, 0.11, 0.56,
#					0.74, 0.09, 0.55), nrow=4, byrow=T)
#SE <- matrix(c(0.1, 0.01, 0.05,
#				0.11, 0.023, 0.055,
#				0.10, 0.005, 0.04,
#				0.14, 0.012, 0.039), nrow=4, byrow=T)
#nimps <- 4
#miPoolVector(param, SE, nimps)

# lmrrPooledChi
# Function -- simsem package
# Pool Chi-square statistic based on Li, Meng, Raghunathan, & Rubin (1991) adapted from http://psychology.clas.asu.edu/files/CombiningLikelihoodRatioChi-SquareStatisticsFromaMIAnalysis.sas
# Argument:
#	chis: 	vector of chi-square values
#	df:		degree of freedom
# Author: 	Craig Enders
#			Sunthud Pornprasertmanit (University of Kansas; psunthud@ku.edu)
# Date Modified: March 31, 2012

lmrrPooledChi <- function(chis, df) {
	# From Li, Meng, Raghunathan, & Rubin (1991)
	if(is.matrix(chis)) {
		ifelse(ncol(chis) == 1 | nrow(chis) == 1, chis <- as.vector(chis), stop("Please put a vector of chi-square values"))
	}
	m <- length(chis)
	dbar <- mean(chis)
	sqrtd <- sqrt(chis)
	xbarsqrtd <- mean(sqrtd)
	# Equation 2.2
	r <- (1 + 1/m) * (sum((sqrtd - xbarsqrtd)^2)/(m - 1))
	# Equation 2.1
	D <- (dbar/df - ((m + 1) * r /(m - 1)))/(1 + r)
	if(D < 0) D <- 0
	# Equation 2.16 and 2.17
	aw <- df^(-(3/m)) * (m - 1) * (1 + (1/r))^2
	p <- 1 - pf(D, df, aw)
	result <- c(D, df, aw, p)
	names(result) <- c("F", "df1", "df2", "p.F")
	return(result)
}
#Examples:
#lmrrPooledChi(c(89.864, 81.116,71.500,49.022,61.986,64.422,55.256,57.890,79.416,63.944), 2)


## function that builds a lavaan parameter table of the saturate model
## using the information in a lavaan object from the cfa function
satPartable <- function(fit.alt){
  
	par.alt <- lavaan::partable(fit.alt) #get the parameter table form the original model
	ngroups <- fit.alt@Data@ngroups # get the number of groups 
	# gets the parameter table from the null model
	
	par.null <- lavaan::lav_partable_independence(fit.alt)
	lhs.diag <- par.null$lhs
	op.diag <- par.null$op
	rhs.diag <- par.null$rhs
	gnull <- par.null$group
	#combine the variable names to set al the covariances
	pairs <- t(combn(lavaan::lavaanNames(fit.alt, type="ov"), 2))
	lhs.up <- rep(pairs[, 1],times=ngroups)
	op.up <- rep("~~", length(lhs.up))
	rhs.up <- rep(pairs[, 2],times=ngroups)
	galt <- sort(rep(1:ngroups,times=length(lhs.up)/ngroups))
	#put together the null table and the covariances
	lhs.all <- c(lhs.up, lhs.diag)
	id <- seq(1:length(lhs.all))
	op.all <- c(op.up, op.diag)
	rhs.all <- c(rhs.up, rhs.diag)
	user <- rep(1,length(lhs.all))
	group <- as.integer(c(galt,gnull))
	free <- as.integer(id)
	ustart <- rep(NA, length(lhs.all))
	exo <- rep(0, length(lhs.all))
	label <- rep("", length(lhs.all))
	plabel <- paste0(".p", id, ".")
	par.sat <- list(id, lhs.all, op.all, rhs.all, user, group,
				  free, ustart, exo, label, plabel)
	names(par.sat) <- c("id", "lhs", "op", "rhs", "user", "group", "free", "ustart", "exo", "label", "plabel")
	return(par.sat)
}

##### function that does the part of the MR and Mplus combination methods are equal 
mrplusPooledChi <- function(template, imputed.l, chi1, df, coef, m, fun, par.sat=NULL, ...) {
	
	if(is.null(par.sat)) par.sat <- satPartable(template)
	comb.sat <- suppressWarnings(lapply(imputed.l, runlavaanMI, syntax=par.sat, fun=fun, ...))
	converged.sat1 <- sapply(comb.sat, function(x) x@Fit@converged)
	
	coefs.sat1 <- sapply(comb.sat, function(x) x@Fit@est)
	est.sat1 <- rowMeans(coefs.sat1[,converged.sat1])
	par.sat2 <- par.sat
	par.sat2$free <- as.integer(rep(0, length(par.sat2$free)))
	par.sat2$ustart <- est.sat1
	comb.sat2 <- suppressWarnings(lapply(imputed.l, runlavaanMI, syntax=par.sat2, fun=fun, ...))
	
    comb.sat2 <- lapply(comb.sat2, forceTest)
	fit.sat2 <- sapply(comb.sat2, function(x) inspect(x, "fit")["logl"])
	
	par.alt2 <- lavaan::partable(template)
	par.alt2$free <- as.integer(rep(0, length(par.alt2$free)))
	par.alt2$ustart <- coef
	comb.alt2 <- suppressWarnings(lapply(imputed.l, runlavaanMI, syntax=par.alt2, fun=fun, ...))
    comb.alt2 <- lapply(comb.alt2, forceTest)
	fit.alt2 <- sapply(comb.alt2, function(x) inspect(x, "fit")["logl"])
  	
	chinew <- cbind(fit.sat2, fit.alt2, (fit.sat2-fit.alt2)*2)
	chimean <- mean(chinew[,3])
	
	fit.altcc <- mean(chi1)
	ariv <- ((m+1)/((m-1)*df))*(fit.altcc-chimean)
	resmrCHI <- c(chimean, m, df, ariv)
	return(resmrCHI)
}

##### function that does the calculations for the Mplus chi combination
mplusPooledChi <- function(chimean, k, ariv){
  comb.chi.mplus <- matrix(NA, nrow=1, ncol=3)
  comb.chi.mplus[1] <- chimean/(1+ariv)
  comb.chi.mplus[2] <- k
  comb.chi.mplus[3] <- 1 - pchisq(comb.chi.mplus[1], comb.chi.mplus[2])
  colnames(comb.chi.mplus) <- c("chisq", "df", "pvalue")
  comb.chi.mplus <- as.data.frame(comb.chi.mplus)
  rownames(comb.chi.mplus) <- ""
  return(comb.chi.mplus)
}

##### function that does the calculations for the MR chi combination
mrPooledChi <-function(chimean, m, k, ariv){
  km <- m*k
  kmtest <- km-k
  
  if(kmtest<=4){
    v4 <- 4+(km-k-4)*(1+(1-(2/kmtest))*(1/ariv))^2
  }
  else{
    v4 <- (kmtest*(1+k^-1)*(1+(1/ariv))^2)/2          
  }       
  comb.chi.mr <- matrix(NA, nrow=1, ncol=4)
  comb.chi.mr[1] <- chimean/((1+ariv)*k)
  comb.chi.mr[2] <- k
  comb.chi.mr[3] <- v4
  comb.chi.mr[4] <- 1 - pf(comb.chi.mr[1], comb.chi.mr[2], comb.chi.mr[3])
  colnames(comb.chi.mr) <- c("F", "df1", "df2", "pvalue")
  comb.chi.mr <- as.data.frame(comb.chi.mr)
  rownames(comb.chi.mr) <- ""
  return(comb.chi.mr)
}

forceTest <- function(object) {
	previousCall <- object@call
	args <- as.list(previousCall[-1])
	args$model <- lavaan::partable(object)
	args$control <- list(optim.method="none", optim.force.converged=TRUE)
	funcall <- as.character(previousCall[[1]])
	lav <- do.call(funcall[length(funcall)], args)
	lav
}

imposeGLIST <- function(object, coef, partable) {
	GLIST <- object@GLIST
	for(i in 1:length(GLIST)) {
		if(!is.matrix(GLIST[[i]])) GLIST[[i]] <- as.matrix(GLIST[[i]])
		dimnames(GLIST[[i]]) <- object@dimNames[[i]]
	}
	for(i in 1:length(coef)) {
		group <- partable$group[i]
		lhs <- partable$lhs[i]
		rhs <- partable$rhs[i]
		if(partable$op[i] == "=~") {
			targetName <- "lambda"
			if(!(rhs %in% rownames(GLIST[names(GLIST) == "lambda"][[group]]))) targetName <- "beta"
			GLIST[names(GLIST) == targetName][[group]][rhs, lhs] <- coef[i]
		} else if (partable$op[i] == "~~") {
			if(lhs %in% rownames(GLIST[names(GLIST) == "psi"][[group]])) {
				GLIST[names(GLIST) == "psi"][[group]][rhs, lhs] <- coef[i]
				GLIST[names(GLIST) == "psi"][[group]][lhs, rhs] <- coef[i]
			} else {
				GLIST[names(GLIST) == "theta"][[group]][rhs, lhs] <- coef[i]
				GLIST[names(GLIST) == "theta"][[group]][lhs, rhs] <- coef[i]
			}
		} else if (partable$op[i] == "~") {
			targetName <- "beta"
			if(!(rhs %in% colnames(GLIST[names(GLIST) == "beta"][[group]]))) targetName <- "gamma"
			GLIST[names(GLIST) == targetName][[group]][lhs, rhs] <- coef[i]
		} else if (partable$op[i] == "~1") {
			if(lhs %in% rownames(GLIST[names(GLIST) == "alpha"][[group]])) {
				GLIST[names(GLIST) == "alpha"][[group]][lhs, 1] <- coef[i]
			} else {
				GLIST[names(GLIST) == "nu"][[group]][lhs, 1] <- coef[i]
			}		
		} else if (partable$op[i] == "|") {
			GLIST[names(GLIST) == "tau"][[group]][paste0(lhs, "|", rhs), 1] <- coef[i]
		}
	}
	object@GLIST <- GLIST
	object
}
