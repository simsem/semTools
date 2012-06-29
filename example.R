library(tools)
dirMan <- "C:/Users/student/Dropbox/semTools/semTools/man/runMI.Rd"
showNonASCIIfile(dirMan)

sourceDir <- function(path, trace = TRUE, ...) {
     for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
		if(nm != "AllClass.R" & nm != "AllGenerics.R") {
        if(trace) cat(nm,":") 
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
		}
     }
}

#get
#assign
#dir <- "C:/Users/Sunthud/Desktop/My Dropbox/simsem/simsem/R/"
#dir <- "C:/Users/Sunthud/simsem_backup/simsem/R/"
dir <- "C:/Users/student/Dropbox/semTools/semTools/R/"
sourceDir(dir)

######### Distribution

skew(1:5)
kurtosis(1:5)

######### measurementInvariance

library(lavaan)
HW.model <- ' visual =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed =~ x7 + x8 + x9 '

measurementInvariance(HW.model, data=HolzingerSwineford1939, group="school")

######### moreFitIndices

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
moreFitIndices(fit)

######### orthogonalize

library(MASS)

n <- 500
means <- c(0,0)
covmat <- matrix(c(1, 0.3, 0.3, 1),nrow=2)

data <- mvrnorm(n,means,covmat)

x<-as.vector(data[,1])
z<-as.vector(data[,2])

y<-rnorm(n,0,1)+.4*x+.4*z+.2*x*z

x1<-rnorm(n,0.2,.2)+.7*x
x2<-rnorm(n,0.2,.2)+.7*x
x3<-rnorm(n,0.2,.2)+.7*x
z1<-rnorm(n,0.2,.2)+.7*z
z2<-rnorm(n,0.2,.2)+.7*z
z3<-rnorm(n,0.2,.2)+.7*z
y1<-rnorm(n,0.2,.2)+.7*y
y2<-rnorm(n,0.2,.2)+.7*y
y3<-rnorm(n,0.2,.2)+.7*y

dat<-data.frame(cbind(x1,x2,x3,z1,z2,z3,y1,y2,y3))

datOrth <-orthogonalize(dat,(1:3), (4:6))

########### Better have the runMI example

library(lavaan)

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

HSMiss <- HolzingerSwineford1939[,paste("x", 1:9, sep="")]
randomMiss <- rbinom(prod(dim(HSMiss)), 1, 0.1)
randomMiss <- matrix(as.logical(randomMiss), nrow=nrow(HSMiss))
HSMiss[randomMiss] <- NA

out <- runMI(HSMiss, HS.model, m = 3)

HSMiss2 <- cbind(HSMiss, school = HolzingerSwineford1939[,"school"])
out2 <- runMI(HSMiss2, HS.model, m = 3, group="school", noms="school")

library(Amelia)

modsim <- '
f1 =~ 0.7*y1+0.7*y2+0.7*y3
f2 =~ 0.7*y4+0.7*y5+0.7*y6
f3 =~ 0.7*y7+0.7*y8+0.7*y9'

mod <- '
f1 =~ y1+y2+y3
f2 =~ y4+y5+y6
f3 =~ y7+y8+y9'

datsim <- simulateData(modsim,model.type="cfa",meanstructure=T,std.lv=T,
                     sample.nobs=c(200,200))
randomMiss2 <- rbinom(prod(dim(datsim)), 1, 0.1)
randomMiss2 <- matrix(as.logical(randomMiss2), nrow=nrow(datsim))
datsim[randomMiss2] <- NA
datsimMI <- amelia(datsim,m=3, noms="group")
datsim3 <- datsim3$imputations

out3 <- runMI(datsimMI$imputations, mod, group="group")

########### Saris, Satorra, and van der Veld (2009)

library(lavaan)

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939, group="sex", meanstructure=TRUE)
miPowerFit(fit)

model <- ' 
  # latent variable definitions
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # regressions
    dem60 ~ ind60
    dem65 ~ ind60 + dem60

  # residual correlations
    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
fit2 <- sem(model, data=PoliticalDemocracy, meanstructure=TRUE)
miPowerFit(fit2, stdLoad=0.3, cor=0.2, stdBeta=0.2, intcept=0.5)

############### Measurement Invariance ############################################

model <- ' f1t1 =~ y1t1 + y2t1 + y3t1
              f1t2 =~ y1t2 + y2t2 + y3t2
			  f1t3 =~ y1t3 + y2t3 + y3t3'

var1 <- c("y1t1", "y2t1", "y3t1")
var2 <- c("y1t2", "y2t2", "y3t2")
var3 <- c("y1t3", "y2t3", "y3t3")
constrainedVar <- list(var1, var2, var3)
longInvariance(model, auto=1, constrainAuto=TRUE, varList=constrainedVar, data=exLong)
longInvariance(model, auto=1, constrainAuto=TRUE, varList=constrainedVar, data=exLong, group="sex", group.equal=c("loadings", "intercepts"))
