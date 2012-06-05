## Monte Carlo test of mediation for complex mediation cases
## Corbin Quick, Alex Schoemann, Kris Preacher, James Selig
## Function that takes an expression for an indirect effect, related parameter estimates and SEs and outputs a Monte Carlo SE
##Output: matrix of LL and UL, optional plot of indirect effect, or values of indirect effect.

monteCarloMed<-function(expression, ..., rep=20000, CI=95, plot=FALSE, outputValues=FALSE, ACOV=NULL){
  input<- c(...)
 # input<-c(1,.5,2,.7)
 
 #Get names and the number of unique variables in the expression
  uniquepar<-function(var){
    var<-gsub(" ","",var) 
    var<-strsplit(var,'+',fixed=TRUE)
    var<-strsplit(var[[1]],'*',fixed=TRUE)
    varb<-var[[1]]
    if(length(var)>1){
      for(i in 2:length(var)){
        varb<-c(varb,var[[i]])}
        var<-unique(varb)}
    if(is.list(var)){var<-var[[1]]}
    return(var)}
	
  #Get coefficents and SE from user input
  if(is.null(ACOV)){
  means<-as.list(input[seq(1,length(input),2)])
  sds<-as.list(input[seq(2,length(input),2)])
   names(means)<-uniquepar(expression)
  names(sds)<-paste(names(means),"sd",sep="")
  }
  if(is.matrix(ACOV)){
  means<-as.list(input)
  names(means)<-uniquepar(expression)
  }
  
  
  #names(means)<-uniquepar(expression)
  #names(sds)<-paste(names(means),"sd",sep="")
  vecs<-list()
    #Use ACOV matrix if provided
    if(is.matrix(ACOV)){
	require(MASS)
	#Matrix of values, need to be converted to a list
	dat <- mvrnorm(n=rep, mu=unlist(means), Sigma=ACOV)
	#Convert to a list,
	vecs<-as.list(as.data.frame(dat))
	}
	if(is.null(ACOV)){
  #Simulate values of all parameters for all reps
  for(i in 1:length(means)){vecs[[i]]<-rnorm(rep)*sds[[i]]+means[[i]]}
  }
  #Give names to it works with assign
  for(i in 1:length(vecs)){assign(names(means)[i],vecs[[i]])}
  
  #Apply the expression to compute the indirect effect
  indirect<-eval(parse(text=expression))
  #Get the CI
  low=(1-CI/100)/2
  upp=((1-CI/100)/2)+(CI/100)
  LL=round(quantile(indirect,low),digits=4)
  UL=round(quantile(indirect,upp),digits=4)
  interval<-rbind(LL,UL)
  dimnames(interval) <- list(c("LL", "UL"),paste(CI, "% Confidence Interval", sep=""))
  
  #Switch for outputting a plot
  if(plot) {
  hist(indirect,breaks='FD',col='skyblue',xlab=paste(CI,'% Confidence Interval ','LL',LL,'  UL',UL), main='Distribution of Indirect Effect')
  }
  
  #Switch to return simulated values
  if(outputValues) {
  interval <- list(interval, indirect)
  }
  
  return(interval)
}

#test simple mediation
 med <- 'a*b'
 ac <- 1
 ase <- .01
 bc<-2
 bse <- .02


 monteCarloMed(med, coef1=ac, se1=ase, coef2=bc, se2=bse, outputValues=F, plot=T)
 
AC <- matrix(c(.01,.00002,.00002,.02), nrow=2)
#If assymptotic covariance matrix is provided, SEs do not need to be provided seperately.
monteCarloMed(med, coef1=ac, coef2=bc, outputValues=F, plot=F, ACOV=AC)


#test complex mediation
# med <- 'a1*b2 + a1*b1'
# ac <- 1
# ase <- 1
# b1c<-2
# b1se <- 2
# b2c<-1
# b2se<-.005


# monteCarloMed(med, coef1=ac, se1=ase, coef2=b1c, se2=b1se, coef3=b2c, se3=b2se)
