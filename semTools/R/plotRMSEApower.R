#Plot power of RMSEA over a range of possible sample sizes
#input: rmea of null and alternative model, degress of freedom, lower sampel size, upper sample sample, sample size steps, alpha
#Output: plot of power
#Alexander M. Schoemann, Kristopher J. Preacher, Donna Coffman
#5/30/2012


plotRMSEApower <- function(rmsea0, rmseaA, df, nlow, nhigh, steps, alpha=.05){ 

pow1<-0
nseq<-seq(nlow,nhigh, by=steps)
for(i in nseq){

ncp0 <- (i-1)*df*rmsea0^2
ncpa <- (i-1)*df*rmseaA^2
#Compute power
if(rmsea0<rmseaA) {
cval <- qchisq(alpha,df,ncp=ncp0,lower.tail=F)
pow <- pchisq(cval,df,ncp=ncpa,lower.tail=F)
}
if(rmsea0>rmseaA) {
cval <- qchisq(1-alpha,df,ncp=ncp0,lower.tail=F)
pow <- 1-pchisq(cval,df,ncp=ncpa,lower.tail=F)
}
pow1<-c(pow1,pow)
}
pow1<-pow1[-1]

plot(nseq,pow1,xlab="Sample Size",ylab="Power",main="Compute Power for RMSEA",type="l")
}

#Example Code
#plotRMSEApower(.025, .075, 23, 100, 500, 10)