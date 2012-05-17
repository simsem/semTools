## Title: Orthogonalize data for 2-way interaction in SEM
## Author: Alexander M. Schoemann <schoemann@ku.edu>
## Description: Orthogonalize data for 2-way interaction in SEM
##----------------------------------------------------------------------------##

## dat:  data frame or matrix with original variables
## xvars: vector of column numbers corresponding to indicators of the focal predictor (x).
## zvars: vector of column numbers corresponding to indicators of the moderator(z).


orthogonalize<-function(dat, xvars, zvars) {
vars <- c(xvars, zvars)
data<-dat[,vars]
#data2 <- dat[,-vars]
#gets total number of variables in data set#
nvar<-dim(data)[2]
#determine number of product terms#
numx <- length(xvars)
numz <- length(zvars)
nxz<-numx*numz

#Make data into a matrix#
mat<-as.matrix(data,ncol=nvar,nrow=dim(dat)[1], byrow=TRUE)

#sequences for variables to use in for loops#
xvar<-c(seq(1,numx,1))
zvar<-c(seq(1,numz,1))

#####create product terms for each variable#####

##Create matrix for product terms to end up in##
prodmat <- matrix(NA,ncol=nxz, nrow=dim(dat)[1])
i = 1
orthName <- NA
for(j in xvar){
   for(k in zvar){
	 prodmat[,i]<-mat[,j]*mat[,numx+k]
	 orthName <- c(orthName, paste(names(dat[j]),names(dat[numx+k]), sep=""))
     i <- i +1
}
}
orthName <- orthName[-1]


###Create Orthogonalized term ###

#determine the number of product terms#
nprod<-dim(prodmat)[2] 

#sequences for product terms to use in for loops#
xzseq<-c(seq(1,nxz,1))

form<-"prodmat[,i]~1"
for(i in xvar){
	ii<-as.character(i)
	formx<-paste("+mat[,",ii,"]", sep="")
	form<-paste(form,formx,sep="")
	}
for(i in zvar){
	ii<-as.character(i)
	formz<-paste("+mat[,numx+",ii,"]", sep="")
	form<-paste(form,formz,sep="")
}
formu<-as.formula(form)

#here's where we can use an apply statement.... or maybe not. Going back to a for loop for now

# residReg <- function(cols) {
	# forma <- paste(cols, form, sep="")
	# formu<-as.formula(forma)
	# reg<-lm(formu)
	# resid<-reg$residuals
	# return(resid)
	# }

# residDat <- apply(prodmat, 2, residReg)
residDat <- matrix(NA, ncol=nxz, nrow=dim(dat)[1])
for(i in 1:nxz){
	reg<-lm(formu)
	residDat[,i]<-residuals(reg)
}
colnames(residDat) <- orthName

mat <- cbind(dat, residDat)
return(mat)

}




