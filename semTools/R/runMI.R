##  Functon to impute missing data, run Lavaan on each one 
##  input: data frames of raw data with missing data, model specification (lavaan script), number of imputations wanted, number of digits to print in output.)
##  Output: list of results with:  fit parameter estimates, standard errors fit indices, and two types of fraction of missing information
##  Alexander Schoemann, Patrick Miller, Mijke Rhemtulla, Sunthud Pornprasertmanit, Alexander Robitzsch, Mauricio Garnier Villarreal
##  Last modified 5/25/2012

##Currently outputs a list of parameter estimates, standard errors, fit indices and fraction missing information

cfa.mi <- function(model, data, m, miArgs=list(), miPackage="Amelia", chi="all", digits=3, seed=12345, std.lv = FALSE, estimator = "ML", group = NULL, group.equal = "", ...) {
	runMI(data.mat=data, data.model=model, m=m, miArgs=miArgs, chi=chi, miPackage=miPackage, digits=digits, seed=seed, std.lv = std.lv, estimator = estimator, group = group, group.equal = group.equal, fun="cfa", ...)
}

sem.mi <- function(model, data, m, miArgs=list(), miPackage="Amelia", chi="all", digits=3, seed=12345, std.lv = FALSE, estimator = "ML", group = NULL, group.equal = "", ...) {
	runMI(data.mat=data, data.model=model, m=m, miArgs=miArgs, chi=chi, miPackage=miPackage, digits=digits, seed=seed, std.lv = std.lv, estimator = estimator, group = group, group.equal = group.equal, fun="sem", ...)
}

growth.mi <- function(model, data, m, miArgs=list(), miPackage="Amelia", chi="all", digits=3, seed=12345, std.lv = FALSE, estimator = "ML", group = NULL, group.equal = "", ...) {
	runMI(data.mat=data, data.model=model, m=m, miArgs=miArgs, chi=chi, miPackage=miPackage, digits=digits, seed=seed, std.lv = std.lv, estimator = estimator, group = group, group.equal = group.equal, fun="growth", ...)
}

lavaan.mi <- function(model, data, m, miArgs=list(), miPackage="Amelia", chi="all", digits=3, seed=12345, std.lv = FALSE, estimator = "ML", group = NULL, group.equal = "", ...) {
	runMI(data.mat=data, data.model=model, m=m, miArgs=miArgs, chi=chi, miPackage=miPackage, digits=digits, seed=seed, std.lv = std.lv, estimator = estimator, group = group, group.equal = group.equal, fun="lavaan", ...)
}


runMI <- function(data.mat,data.model, m, miArgs=list(), chi="all", miPackage="Amelia", digits=3, seed=12345, std.lv = FALSE, estimator = "ML", group = NULL, group.equal = "", fun, ...) 
{
	test <- HolzingerSwineford1939[-301,-5]
	data.model <- ' visual  =~ NA*x1 + x2 + x3
               textual =~ NA*x4 + x5 + x6
               speed   =~ NA*x7 + x8 + x9 
			   visual ~~1*visual
			   textual ~~ 1*textual
			   speed ~~ 1*speed
			   visual ~ speed'

	##Impose missing data to test
	log.mat1 <- matrix(FALSE, nrow=dim(test)[1], ncol=dim(test)[2])
	log.mat1[,9] <- test$x6>3
	test[log.mat1] <- NA
	data.mat <- test
	m <- 20
	miArgs=list(); chi="all"; miPackage="Amelia"; digits=3; seed=12345; std.lv = FALSE; estimator = "ML"; group = "grade"; group.equal = ""
	fun <- "sem"
	
	set.seed(seed)
	imputed.data <- is.list(data.mat) & (!is.data.frame(data.mat))
	imputed.l <- NULL
	if (!imputed.data){		
		if( ( miPackage!="Amelia" )  &  ( miPackage !="mice")  ) { 
			stop("Currently runMI only supports imputation by Amelia or mice") 
		}
		if(miPackage=="Amelia"){
			imputed.l<-imputeMissingAmelia(data.mat,m, miArgs)
		} else if(miPackage=="mice"){
			imputed.l<-imputeMissingMice(data.mat,m, miArgs)
		}
	} else { 
		imputed.l <- data.mat 
		m <- length( data.mat )
	}
	out <- list(model=data.model, data=data.mat, std.lv = std.lv, estimator = estimator, group = group, group.equal = group.equal, meanstructure = TRUE, se="none", do.fit=FALSE)
	out <- c(out, list(...))
	template <- do.call(fun, out)

    imputed.results.l <- lapply(imputed.l, runlavaanMI, syntax=data.model, std.lv = std.lv, 
	estimator = estimator, group = group, group.equal = group.equal, fun=fun)
    
	coefs <- sapply(imputed.results.l, function(x) x@Fit@est)
	se <- sapply(imputed.results.l, function(x) x@Fit@se)
	
	comb.results <- miPoolVector(t(coefs),t(se), m)
	template@Fit@est <- comb.results$coef
	template@Fit@se <- comb.results$se
	template@Fit@x <- comb.results$coef[comb.results$se != 0]
	# Do not need to change imputed.results.l@Model@GLIST because the methods from the lavaan object does not use that information
	
	fmi.results <- cbind(parameterEstimates(template)[,1:3], group=length(template@ParTable$group), fmi1 = comb.results[[3]], fmi2 = comb.results[[4]])

	fit <- imputed.results.l[[1]]@Fit@test
	df <- fit[[1]]$df
	chi1 <- sapply(imputed.results.l, function(x) x@Fit@test[[1]]$stat)
	fit[[1]]$stat.group <- rowMeans(sapply(imputed.results.l, function(x) x@Fit@test[[1]]$stat.group))

	nullModel <- partable(lavaan:::independence.model.fit(template))
    null.results <- lapply(imputed.l, runlavaanMI, syntax=nullModel, std.lv = std.lv, 
	estimator = estimator, group = group, group.equal = group.equal, fun=fun)
	chiNull <- sapply(null.results, function(x) x@Fit@test[[1]]$stat)
	dfNull <- null.results[[1]]@Fit@test[[1]]$df
	
	outNull <- list(model=nullModel, data=data.mat, std.lv = std.lv, estimator = estimator, group = group, group.equal = group.equal, meanstructure = TRUE, se="none", do.fit=FALSE)
	outNull <- c(outNull, list(...))
	templateNull <- do.call(fun, outNull)
	
	coefsNull <- sapply(null.results, function(x) x@Fit@est)
	seNull <- sapply(null.results, function(x) x@Fit@se)
	
	comb.results.null <- miPoolVector(t(coefsNull),t(seNull), m)
	fitNull <- null.results[[1]]@Fit@test
	
	if(chi %in% c("LMRR", "all")){
		lmrr <- lmrrPooledChi(chi1, df)
		lmrrNull <- lmrrPooledChi(chiNull, df)
		fit[[1]]$stat <- lmrr[1] * lmrr[2]
		fit[[1]]$pvalue <- lmrr[4]
		fitNull[[1]]$stat <- lmrrNull[1] * lmrrNull[2]
		fitNull[[1]]$pvalue <- lmrrNull[4]
	}
	
	if(chi %in% c("Mplus", "MR", "all")){
		mrplus <- mrplusPooledChi(template, chi1, df, coef=comb.results$coef, estimator = estimator, std.lv = std.lv, 
		group = group, group.equal = group.equal, fun)
		mrplusNull <- mrplusPooledChi(templateNull, chiNull, dfNull, coef=comb.results.null$coef, estimator = estimator, std.lv = std.lv, group = group, group.equal = group.equal, fun, par.sat=satPartable(template))
		
		if(chi %in% c("MR", "all")){
			mr <- mrPooledChi(mrplus[1], mrplus[2], mrplus[3], mrplus[4], digits = digits)
			mrNull <- mrPooledChi(mrplusNull[1], mrplusNull[2], mrplusNull[3], mrplusNull[4], digits = digits)
			fit[[1]]$stat <- mr[1] * mr[2]
			fit[[1]]$pvalue <- mr[4]
			fitNull[[1]]$stat <- mrNull[1] * mrNull[2]
			fitNull[[1]]$pvalue <- mrNull[4]
		}
		if(chi %in% c("Mplus", "all")){
			mplus <- mplusPooledChi(mrplus[1], mrplus[3], mrplus[4], digits = digits)
			mplusNull <- mplusPooledChi(mrplusNull[1], mrplusNull[3], mrplusNull[4], digits = digits)
			fit[[1]]$stat <- as.numeric(mplus[1])
			fit[[1]]$pvalue <- as.numeric(mplus[3])
			fitNull[[1]]$stat <- as.numeric(mplusNull[1])
			fitNull[[1]]$pvalue <- as.numeric(mplusNull[3])
		}
	}
	template@Fit@test <- fit
	templateNull@Fit@test <- fitNull
	result <- as(template, "lavaanStar")
	fitVec <- fitMeasures(templateNull)
	name <- names(fitVec)
	fitVec <- as.vector(fitVec)
	names(fitVec) <- name
	result@nullfit <- fitVec
	
	# Fix the SRMR by average the model implied stats
	# Fix the shown warning sign
	
	return(result)
}
  
#Convenient function to run lavaan models and get results out. For easy use with lapply
runlavaanMI <- function(MIdata, syntax, std.lv = FALSE, estimator = "ML", group = NULL, group.equal = "", fun, ...) {
	out <- list(model=syntax, data=MIdata, std.lv = std.lv, estimator = estimator, group = group, group.equal = group.equal, meanstructure = TRUE)
	out <- c(out, list(...))
	fit <- NULL
	try(fit <- do.call(fun, out), silent=TRUE)
    # FitIndices <- inspect(fit, 'fit')
	# Converged = TRUE
	# if(sum(unlist(lapply(inspect(fit, "se"), sum))) == 0) Converged = FALSE
	# params <- parameterEstimates(fit,standardized=T)
    return(fit)
}
	
testMI <- function() {
##Shamelessly using the example in lavaan

test<-HolzingerSwineford1939[,-5]
HS.model <- ' visual  =~ NA*x1 + x2 + x3
               textual =~ NA*x4 + x5 + x6
               speed   =~ NA*x7 + x8 + x9 
			   visual ~~1*visual
			   textual ~~ 1*textual
			   speed ~~ 1*speed
			   visual ~ speed'
summary(cfa(HS.model,data=test), fit.measures=TRUE)

test[(test$x6>3 && test$x6<4)]
##Impose missing data to test
log.mat1 <- matrix(FALSE, nrow=dim(test)[1], ncol=dim(test)[2])
log.mat1[,9] <- test$x6>3
test[log.mat1] <- NA

runMI(test,HS.model,3, idvars='id')



# fit@Options (May need to be changed some options)
# fit@Model@GLIST
# fit@Fit@est; fit@Fit@sd; fit@Fit@test 

}


#Conveniance function to run impuations on data and only return list of data
imputeMissingAmelia <- function(data.mat,m, miArgs){
  # pull out only the imputations
  require(Amelia)
  out <- c(list(amelia, x = data.mat, m = m, p2s=0), miArgs)
  temp.am <- eval(as.call(out))
  return(temp.am$imputations)

} # end imputeMissingAmelia

imputeMissingMice <- function(data.mat,m, miArgs){
  # pull out only the imputations
  require(mice)
  out <- c(list(mice, data=data.mat, m = m, diagnostics=FALSE, printFlag=FALSE), miArgs)
  temp.mice <- eval(as.call(out))
  temp.mice.imp <- NULL
  for(i in 1:m) {
	temp.mice.imp[[i]] <- complete(x=temp.mice, action=i, include=FALSE) 
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
  
	par.alt<-partable(fit.alt) #get the parameter table form the original model
	ngroups <- fit.alt@Data@ngroups # get the number of groups 
	# gets the parameter table from the null model
	par.null <- partable(lavaan:::independence.model.fit(fit.alt))
	lhs.diag <- par.null$lhs
	op.diag <- par.null$op
	rhs.diag <- par.null$rhs
	gnull <- par.null$group
	#combine the variable names to set al the covariances
	pairs <- t(combn(lavaanNames(fit.alt, type="ov"), 2))
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
	eq.id <- exo
	unco <- id
	par.sat <- list(id, lhs.all, op.all, rhs.all, user, group,
				  free, ustart, exo, label, eq.id, unco)
	names(par.sat)<-colnames(par.alt)
	return(par.sat)
}

##### function that does the part of the MR and Mplus combination methods are equal 
mrplusPooledChi <- function(template, chi1, df, coef, estimator, std.lv, group, group.equal, fun, par.sat=NULL) {
	if(is.null(par.sat)) par.sat <- satPartable(template)
	par.alt<-partable(template)
	comb.sat <- lapply(imputed.l, runlavaanMI, syntax=par.sat, std.lv = std.lv, 
	estimator = estimator, group = group, group.equal = group.equal, fun=fun)

	coefs.sat1 <- sapply(comb.sat, function(x) x@Fit@est)
	
	fit.altcc <- mean(chi1)
	est.sat1 <- rowMeans(coefs.sat1)
	est.alt1 <- coef

	par.sat2 <- par.sat
	par.sat2$free <- as.integer(rep(0, length(par.sat2$free)))
	par.sat2$ustart <- est.sat1

	par.alt2 <- par.alt
	par.alt2$free <- as.integer(rep(0, length(par.alt2$free)))
	par.alt2$ustart <- est.alt1

	comb.sat2 <- lapply(imputed.l, runlavaanMI, syntax=par.sat2, std.lv = std.lv, 
	estimator = estimator, group = group, group.equal = group.equal, fun=fun)

	comb.alt2 <- lapply(imputed.l, runlavaanMI, syntax=par.alt2, std.lv = std.lv, 
	estimator = estimator, group = group, group.equal = group.equal, fun=fun)
					  
	fit.sat2 <- sapply(comb.sat2, function(x) inspect(x, "fit")["logl"])
	fit.alt2 <- sapply(comb.alt2, function(x) inspect(x, "fit")["logl"])
  
	chinew <- cbind(fit.sat2, fit.alt2, (fit.sat2-fit.alt2)*2)
	chimean <- mean(chinew[,3])
	
	ariv <- ((m+1)/((m-1)*df))*(fit.altcc-chimean)
	resmrCHI <- c(chimean, m, df, ariv)
	return(resmrCHI)
}

##### function that does the calculations for the Mplus chi combination
mplusPooledChi <- function(chimean, k, ariv, digits = digits){
  comb.chi.mplus <- matrix(NA, nrow=1, ncol=3)
  comb.chi.mplus[1] <- chimean/(1+ariv)
  comb.chi.mplus[2] <- k
  comb.chi.mplus[3] <- 1-pchisq(comb.chi.mplus[1], comb.chi.mplus[2])
  colnames(comb.chi.mplus) <- c("chisq", "df", "pvalue")
  comb.chi.mplus <- as.data.frame(round(comb.chi.mplus, digits=digits))
  rownames(comb.chi.mplus) <- ""
  return(comb.chi.mplus)
}

##### function that does the calculations for the MR chi combination
mrPooledChi <-function(chimean, m, k, ariv, digits = digits){
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
  comb.chi.mr[4] <- 1-pf(comb.chi.mr[1], comb.chi.mr[2], comb.chi.mr[3])
  colnames(comb.chi.mr) <- c("F", "df1", "df2", "pvalue")
  comb.chi.mr <- as.data.frame(round(comb.chi.mr, digits=digits))
  rownames(comb.chi.mr) <- ""
  return(comb.chi.mr)
}