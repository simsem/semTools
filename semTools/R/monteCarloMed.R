## Monte Carlo test of mediation for complex mediation cases
## Corbin Quick, Alex Schoemann, James Selig
## Function that takes an expression for an indirect effect, related parameter estimates and SEs and outputs a Monte Carlo SE
##Output: matrix of LL and UL, optional plot of indirect effect, or values of indirect effect.

monteCarloMed<-function(expression, ..., ACM=NULL, object = NULL, rep=20000, CI=95, plot=FALSE, outputValues=FALSE){
     	 
	 input<- c(...)
	 
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
	
  paramnames<-uniquepar(expression)

  #If input is a lavaan object pull out coefs and ACM
  if(class(object)=="lavaan"){
  input <- lavaan::coef(object)[paramnames]
  ACM <- lavaan::vcov(object)[paramnames,paramnames]  
  }

  vecs<-list()
	#Matrix of values, need to be converted to a list
	dat <- MASS::mvrnorm(n=rep, mu=input, Sigma=ACM)
	#Add parameters as the first row
	dat <-rbind(input, dat)
	#Convert to a list,
	vecs<-as.list(as.data.frame(dat))
  #Give names to it works with assign
  for(i in 1:length(vecs)){assign(paramnames[i],vecs[[i]])}
  
  #Apply the expression to compute the indirect effect
  indirect<-eval(parse(text=expression))
  #Get the CI
  low=(1-CI/100)/2
  upp=((1-CI/100)/2)+(CI/100)
  LL=round(quantile(indirect[-1],low),digits=4)
  UL=round(quantile(indirect[-1],upp),digits=4)
  interval<-list(indirect[1],rbind(LL,UL))
  dimnames(interval[[2]]) <- list(c("LL", "UL"),c(" "))
  names(interval) <- c("Point Estimate", paste(CI, "% Confidence Interval", sep=""))
  
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
