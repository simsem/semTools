##  Functon to impute missing data, run Lavaan on each one 
##  input: data frames of raw data with missing data, model specification (lavaan script), number of imputations wanted, number of digits to print in output.)
##  Output: list of results with:  fit parameter estimates, standard errors fit indices, and two types of fraction of missing information
##  Alexander Schoemann, Patrick Miller, Mijke Rhemtulla, Sunthud Pornprasertmanit, Alexander Robitzsch, Mauricio Garnier Villarreal
##  Last modified 5/25/2012

##Currently outputs a list of parameter estimates, standard errors, fit indices and fraction missing information

runMI <- function(data.mat,data.model, m, chi="all", miPackage="Amelia", digits=3, seed=12345, std.lv = FALSE, estimator = "ML", group = NULL, group.equal = "", ...) 
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

############# here we start the 3 chi square combination methods

if(chi=="LMRR" | chi=="all"){ ## here is the code that get the LMRR chi square combination
  comb.chi.lmrr <- data.frame(t(miPoolChi(fit[,1], fit[1,2])))
  comb.chi.lmrr <- round(comb.chi.lmrr, digits=digits)
  colnames(comb.chi.lmrr) <- c("F","df1","df2","pvalue")
  rownames(comb.chi.lmrr)<-""
}
#### here we do the mplus and mr methods
if(chi == "Mplus" | chi == "MR" | chi == "all"){
  fit.alt <- cfa(data.model, data=imputed.l[[1]], std.lv = std.lv, meanstructure=T,
                 estimator = estimator, group = group, group.equal = group.equal)
    
  mrplus <- mrplusCHI(fit.alt, imputed.l, imputed.results.l, comb.results, estimator = estimator, std.lv = std.lv, group = group, group.equal = group.equal)

  if(chi == "Mplus" | chi == "all"){
    comb.chi.mplus <- mplusCHI(mrplus[1], mrplus[3], mrplus[4], digits = digits)
    }
  
if(chi == "MR" | chi == "all"){  
  comb.chi.mr <- mrCHI(mrplus[1], mrplus[2], mrplus[3], mrplus[4], digits = digits)
  }
}
### here we put all the results together depending on the method to combine the chi square
if(chi == "LMRR") {
  results <- list(comb.fit, comb.chi.lmrr, comb.results)
  names(results) <- c('Average fit statistics. These may not be trustworthy!', 
                      'Pooled Chi-square statistic LMRR', 
                      'parameters')
}
if(chi == "Mplus") {
  results <- list(comb.fit, comb.chi.mplus, comb.results)
  names(results) <- c('Average fit statistics. These may not be trustworthy!', 
                      'Pooled Chi-square statistic Mplus',
                      'parameters')
}
if(chi == "MR") {
  results <- list(comb.fit, comb.chi.mr, comb.results)
  names(results) <- c('Average fit statistics. These may not be trustworthy!', 
                      'Pooled Chi-square statistic MR', 
                      'parameters')
}
if(chi == "all"){
  results <- list(comb.fit, comb.chi.mr, comb.chi.mplus,
                  comb.chi.lmrr, comb.results)
  names(results) <- c('Average fit statistics. These may not be trustworthy!', 
                      'Pooled Chi-square statistic MR',   
                      'Pooled Chi-square statistic Mplus',
                      'Pooled Chi-square statistic LMRR', 
                      'parameters')
}

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
mrplusCHI<-function(fit.alt, imputed.l, imputed.results.l, comb.results, estimator = estimator, std.lv = std.lv, group = group, group.equal = group.equal){
  
  par.sat <- satPartable(fit.alt)
  
  par.alt<-partable(fit.alt)
  
  comb.sat <- lapply(imputed.l, runlavaanMI, par.sat, std.lv = std.lv, estimator = estimator, group = group, group.equal = group.equal)
  
  coefs.sat1 <- matrix(NA, nrow = length(imputed.l), ncol = length(comb.sat[[1]][[1]]$est))
  fit.alt1 <- matrix(NA, nrow = length(imputed.l), ncol = length(imputed.results.l[[1]][[2]]["chisq"]))
  
  for(i in 1:length(comb.sat)){
    coefs.sat1[i,] <- comb.sat[[i]][[1]]$est
    fit.alt1[i,] <- imputed.results.l[[i]][[2]]["chisq"]
  }
  
  fit.altcc <- colMeans(fit.alt1)
  est.sat1 <- colMeans(coefs.sat1)
  est.alt1 <- comb.results$coef
  
  par.sat2 <- par.sat
  par.sat2$free <- as.integer(rep(0, length(par.sat2$free)))
  par.sat2$ustart <- est.sat1
  
  par.alt2 <- par.alt
  par.alt2$free <- as.integer(rep(0, length(par.alt2$free)))
  par.alt2$ustart <- est.alt1
  
  comb.sat2 <- lapply(imputed.l, runlavaanMI, par.sat2, std.lv = std.lv, 
                      estimator = estimator, group = group, group.equal = group.equal)
  
  comb.alt2 <- lapply(imputed.l, runlavaanMI, par.alt2, std.lv = std.lv, 
                      estimator = estimator, group = group, group.equal = group.equal)
  
  
  fit.sat2 <- matrix(NA, nrow = length(imputed.l), ncol = length(comb.sat2[[1]][[2]]["logl"]))
  fit.alt2 <- matrix(NA, nrow = length(imputed.l), ncol = length(comb.alt2[[1]][[2]]["logl"]))
  
  for(i in 1:length(comb.alt2)){
    fit.alt2[i,] <- comb.alt2[[i]][[2]]["logl"]
    fit.sat2[i,] <- comb.sat2[[i]][[2]]["logl"]
  }
  chinew <- matrix(NA,ncol=3,nrow=length(imputed.l))
  chinew[,1] <- fit.sat2
  chinew[,2] <- fit.alt2
  chinew[,3] <- (chinew[,1]-chinew[,2])*2
  chimean <- mean(chinew[,3])
  
  m <- length(imputed.l)
  k <- fitMeasures(fit.alt)["df"]
  ariv <- ((m+1)/((m-1)*k))*(fit.altcc-chimean)
  resmrCHI <- c(chimean, m, k, ariv)
  return(resmrCHI)
}

##### function that does the calculations for the Mplus chi combination
mplusCHI <- function(chimean, k, ariv, digits = digits){
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
mrCHI <-function(chimean, m, k, ariv, digits = digits){
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