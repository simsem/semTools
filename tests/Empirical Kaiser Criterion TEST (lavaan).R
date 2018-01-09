#Empirical Kaiser criterion
#Load the EKC function

library(lavaan)
model <- 'f1=~.7*i1+.7*i2+.7*i3+.7*i4+.7*i5
f2=~.6*i6+.6*i7+.6*i8+.6*i9+.6*i10
f3=~.7*i11+.7*i12+.7*i13+.7*i14+.7*i15
f4=~.7*i16+.7*i17+.7*i18+.7*i19+.7*i20
f1~~.2*f2
f1~~.3*f3
f1~~.3*f4
f2~~.4*f3
f2~~.4*f4
f3~~.3*f4'

set.seed(123)
d <- simulateData(model=model, model.type = "cfa", sample.nobs=200)
dim(d)

#load efa.ekc function

efa.ekc(d, plot=FALSE)

x <- efa.ekc(d) # from data
x1 <- efa.ekc(sample.cov=cov(d), sample.nobs=200) #from covariance matrix
x2 <- efa.ekc(sample.cov=cor(d), sample.nobs=200) #from correlation matrix
x
x1
x2

#The original Keiser's criterion suggests 5 factors, but EKC correctly suggests 4

# check that results are consistent with data, covariance matrix and correlation matrix
round(x$Eigenvalues, 10) == round(x1$Eigenvalues, 10)
round(x$Eigenvalues, 10) == round(x2$Eigenvalues, 10)

#compare eigenvalues from ekc function with eigenvalues from eigen function
round(x$Eigenvalues[,1], 3) == round(eigen(cor(d))$values, 3)

#Compare reference eigenvalues from EKC to the ones produced by the shiny. Go to, 
# https://cemo.shinyapps.io/EKCapp/
#Insert sample PCA eigenvalues, change sample size to 200 and go to "Table":
x$Eigenvalues[,1]

#Compare the website EKC results to those from the EKC function
#Click on EKC unconstrained to show all reference eigenvalues
round(x$Eigenvalues[,2], 3)
