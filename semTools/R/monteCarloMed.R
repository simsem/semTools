## Monte Carlo test of mediation for complex mediation cases
## Corbin Quick, Alex Schoemann, Kris Preacher, James Selig
## Function that takes an expression for an indirect effect, related parameter estimates and SEs and outputs a Monte Carlo SE
##Output: matrix of LL and UL, optional plot of indirect effect, or values of indirect effect.

monteCarloMed<-function(expression, ..., rep=20000, CI=95, plot=FALSE, outputValues=FALSE){
  input<- c(...)
  
  #Get coefficents and SE from user input
  means<-as.list(input[seq(1,length(input),2)])
  sds<-as.list(input[seq(2,length(input),2)])
  
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
  names(means)<-uniquepar(expression)
  names(sds)<-paste(names(means),"sd",sep="")
  vecs<-list()
  
  #Simulate values of all parameters for all reps
  for(i in 1:length(means)){vecs[[i]]<-rnorm(rep)*sds[[i]]+means[[i]]}
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
  hist(indirect,breaks='FD',col='skyblue',xlab=paste(conf,'% Confidence Interval ','LL',LL,'  UL',UL), main='Distribution of Indirect Effect')
  }
  
  #Switch to return simulated values
  if(outputValues) {
  interval <- list(interval, indirect)
  }
  
  return(interval)
}

#test simple mediation
# med <- 'a*b'
# ac <- 1
# ase <- .01
# bc<-2
# bse <- .02

# monteCarloMed(med, coef1=ac, se1=ase, coef2=bc, se2=bse, outputValues=F, plot=T)

#test complex mediation
# med <- 'a1*b2 + a1*b1'
# ac <- 1
# ase <- 1
# b1c<-2
# b1se <- 2
# b2c<-1
# b2se<-.005


# monteCarloMed(med, coef1=ac, se1=ase, coef2=b1c, se2=b1se, coef3=b2c, se3=b2se)
