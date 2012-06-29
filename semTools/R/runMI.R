##  Functon to impute missing data, run Lavaan on each one 
##  input: data frames of raw data with missing data, model specification (lavaan script), number of imputations wanted, number of digits to print in output.)
##  Output: list of results with:  fit parameter estimates, standard errors fit indices, and two types of fraction of missing information
##  Alexander Schoemann, Patrick Miller, Mijke Rhemtulla, Sunthud Pornprasertmanit, Alexander Robitzsch, Mauricio Garnier Villarreal
##  Last modified 5/25/2012

##Currently outputs a list of parameter estimates, standard errors, fit indices and fraction missing information

runMI <- function(data.mat,data.model, m, miPackage="Amelia", digits=3, seed=12345, std.lv = FALSE, estimator = "ML", group = NULL, group.equal = "", ...) 
{

set.seed(seed)
imputed.data <- is.list(data.mat) & (!is.data.frame(data.mat))
imputed.l <- NULL
if (!imputed.data){		
  if( ( miPackage!="Amelia" )  &  ( miPackage !="mice")  )
		{ stop("Currently runMI only supports imputation by Amelia or mice") }
		
  args <- list(...)

  if(miPackage=="Amelia"){
  imputed.l<-imputeMissingAmelia(data.mat,m, ...)
  }
  
  if(miPackage=="mice"){
  imputed.l<-imputeMissingMice(data.mat,m, ...)
  }

		} else { 
				imputed.l <- data.mat 
				m <- length( data.mat )
					}
    imputed.results.l <- lapply(imputed.l, runlavaanMI, syntax=data.model, std.lv = std.lv, 
	estimator = estimator, group = group, group.equal = group.equal)
    
	coefs <- matrix(NA, nrow = m, ncol = length(imputed.results.l[[1]][[1]]$est))
	stdlv <- matrix(NA, nrow = m, ncol = length(imputed.results.l[[1]][[1]]$std.lv))
	stdall <- matrix(NA, nrow = m, ncol = length(imputed.results.l[[1]][[1]]$std.all))
	stdnox <- matrix(NA, nrow = m, ncol = length(imputed.results.l[[1]][[1]]$std.nox))
	se <- coefs
	fit <- matrix(NA, nrow = m, ncol = length(imputed.results.l[[1]][[2]]))
	
	for(i in 1:length(imputed.results.l)){
		coefs[i,] <- imputed.results.l[[i]][[1]]$est
		stdlv[i,] <- imputed.results.l[[i]][[1]]$std.lv
		stdall[i,] <- imputed.results.l[[i]][[1]]$std.all
		stdnox[i,] <- imputed.results.l[[i]][[1]]$std.nox
		se[i,] <- imputed.results.l[[i]][[1]]$se
		fit[i,] <- imputed.results.l[[i]][[2]]
		}

  comb.results <- miPoolVector(coefs,se, m)
  comb.stdlv <- miPoolVector(stdlv,se, m)
  comb.stdall <- miPoolVector(stdall,se, m)
  comb.stdnox <- miPoolVector(stdnox,se, m)
  Wald <- comb.results[[1]]/comb.results[[2]]
  p <- 2*pnorm(-abs(Wald))

  comb.results <- cbind(comb.results[[1]],comb.results[[2]], Wald, p, comb.stdlv[[1]], comb.stdall[[1]], comb.stdnox[[1]], comb.results[[3]], comb.results[[4]])
  comb.results <- as.data.frame(comb.results)
  comb.results <- round(comb.results, digits=digits)
  pval <- comb.results[,"p"] == 0
  pval[is.na(pval)] <- FALSE
  comb.results[pval,"p"] <- paste("<.", paste(rep(0, (digits-1)),collapse=""), 1, sep="")
  colnames(comb.results) <- c('coef', 'se', 'Wald', 'p', 'std.lv', 'std.all', 'std.nox',
  'FMI.1', 'FMI.2')
  
    fixedParam <- is.nan(comb.results[,8])
  comb.results[,2][fixedParam] <- ""
  comb.results[,3][fixedParam] <- ""
  comb.results[,4][fixedParam] <- ""
  comb.results[,5][fixedParam] <- ""
  comb.results[,6][fixedParam] <- ""
  comb.results[,7][fixedParam] <- ""
  comb.results[,8][fixedParam] <- ""
  comb.results[,9][fixedParam] <- ""
  
  headings <- cbind(lhs = imputed.results.l[[1]][[1]]$lhs, op = imputed.results.l[[1]][[1]]$op, rhs = imputed.results.l[[1]][[1]]$rhs, group = imputed.results.l[[1]][[1]]$group)
  comb.results <- data.frame(headings, comb.results)
  

  
  comb.fit <- colMeans(fit)
  names(comb.fit) <- names(imputed.results.l[[1]][[2]])
  comb.fit <- data.frame(comb.fit)
  comb.fit <- round(comb.fit, digits=digits)
  colnames(comb.fit) <- ""
  comb.chi <- miPoolChi(fit[,1], fit[1,2])
  comb.chi <- data.frame(comb.chi)
  comb.chi <- round(comb.chi, digits=digits)
  colnames(comb.chi) <- "" 
  
  fit.results <- list(comb.fit, comb.chi)
  names(fit.results) <- c('Average fit statistics. These may not be trustworthy!', 'Pooled Chi-square statistic')
  
  results <- list(fit.results, comb.results)
  names(results) <- c('fit', 'parameters')
  
  return(results)
}
  
#Conveniance function to run lavaan models and get results out. For easy use with lapply
runlavaanMI <- function(MIdata,syntax, std.lv = FALSE, estimator = "ML", group = NULL, group.equal = "") {
     fit <- cfa(syntax, data=MIdata, std.lv = std.lv, estimator = estimator, 
	 group = group, group.equal = group.equal, meanstructure = TRUE)
     FitIndices <- inspect(fit, 'fit')
	Converged = TRUE
	if(sum(unlist(lapply(inspect(fit, "se"), sum))) == 0) Converged = FALSE
	params <- parameterEstimates(fit,standardized=T)
    return(list(params, FitIndices, Converged))
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

}


#Conveniance function to run impuations on data and only return list of data
imputeMissingAmelia <- function(data.mat,m, ...){
  # pull out only the imputations
  require(Amelia)
  temp.am <- amelia(data.mat,m, p2s=0, ...)
  return(temp.am$imputations)

} # end imputeMissingAmelia

imputeMissingMice <- function(data.mat,m, ...){
  # pull out only the imputations
  require(mice)
  temp.mice <- mice(data.mat,m, diagnostics=FALSE, printFlag=FALSE, ...)
  temp.mice.imp <- NULL
  for (i in 1:m) {
  temp.mice.imp[[i]] <- complete(temp.mice, action=i)
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

# miPoolChi
# Function -- simsem package
# Pool Chi-square statistic based on Li, Meng, Raghunathan, & Rubin (1991) adapted from http://psychology.clas.asu.edu/files/CombiningLikelihoodRatioChi-SquareStatisticsFromaMIAnalysis.sas
# Argument:
#	chis: 	vector of chi-square values
#	df:		degree of freedom
# Author: 	Craig Enders
#			Sunthud Pornprasertmanit (University of Kansas; psunthud@ku.edu)
# Date Modified: March 31, 2012

miPoolChi <- function(chis, df) {
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
#miPoolChi(c(89.864, 81.116,71.500,49.022,61.986,64.422,55.256,57.890,79.416,63.944), 2)


