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

sourceDirData <- function(path, trace = TRUE) {
     for (nm in list.files(path, pattern = "\\.[Rr]da$")) {
		if(nm != "AllClass.R" & nm != "AllGenerics.R") {
        if(trace) cat(nm,":") 
        load(file.path(path, nm))
        if(trace) cat("\n")
		}
     }
}

#get
#assign
#dir <- "C:/Users/Sunthud/Desktop/My Dropbox/simsem/simsem/R/"
#dir <- "C:/Users/Sunthud/simsem_backup/simsem/R/"
library(lavaan)

dir <- "C:/Users/student/Dropbox/semTools/semTools/R/"
sourceDir(dir)

dirData <- "C:/Users/student/Dropbox/semTools/semTools/data/"
sourceDirData(dirData)
######### Distribution

skew(1:5)
kurtosis(1:5)
mardiaSkew(HolzingerSwineford1939[,paste("x", 1:9, sep="")])
mardiaKurtosis(HolzingerSwineford1939[,paste("x", 1:9, sep="")])

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

################ Power Analysis ################################################

plotRMSEApower(rmsea0=.05, rmseaA=.075, df=23, nlow=100, nhigh=500, steps=10)

plotRMSEAdist(rmsea=c(.05, .08), n=200, df=20, ptile=0.95, rmseaScale = TRUE)
plotRMSEAdist(rmsea=c(.05, .01), n=200, df=20, ptile=0.05, rmseaScale = FALSE)

findRMSEApower(rmsea0=.05, rmseaA=.08, df=20, n=200)

findRMSEAsamplesize(rmsea0=.05, rmseaA=.08, df=20, power=0.80)

################ Probing interaction ###########################################

dat2wayRC <- orthogonalize(dat2way, 1:3, 4:6)
dat2wayMC <- indProd(dat2way, 1:3, 4:6)
dat2wayDMC <- indProd(dat2way, 1:3, 4:6, doubleMC=TRUE)

dat3wayRC <- orthogonalize(dat3way, 1:3, 4:6, 7:9)
dat3wayMC <- indProd(dat3way, 1:3, 4:6, 7:9)
dat3wayDMC <- indProd(dat3way, 1:3, 4:6, 7:9, doubleMC=TRUE)

library(lavaan) 
model1 <- "
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
f12 =~ x1.x4 + x2.x5 + x3.x6
f3 =~ x7 + x8 + x9
f3 ~ f1 + f2 + f12
f12 ~~0*f1
f12 ~~ 0*f2
x1 ~ 0*1
x4 ~ 0*1
x1.x4 ~ 0*1
x7 ~ 0*1
f1 ~ NA*1
f2 ~ NA*1
f12 ~ NA*1
f3 ~ NA*1
"
fitRC2way <- sem(model1, data=dat2wayRC, meanstructure=TRUE, std.lv=FALSE)
summary(fitRC2way)
result2wayRC <- probe2WayRC(fitRC2way, c("f1", "f2", "f12"), "f3", "f2", c(-1, 0, 1))
result2wayRC

fitMC2way <- sem(model1, data=dat2wayMC, meanstructure=TRUE, std.lv=FALSE)
summary(fitMC2way)
result2wayMC <- probe2WayMC(fitMC2way, c("f1", "f2", "f12"), "f3", "f2", c(-1, 0, 1))
result2wayMC

plotProbe(result2wayRC, xlim=c(-2, 2))
plotProbe(result2wayMC, xlim=c(-2, 2))

library(lavaan)
model3 <- "
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
f3 =~ x7 + x8 + x9
f12 =~ x1.x4 + x2.x5 + x3.x6
f13 =~ x1.x7 + x2.x8 + x3.x9
f23 =~ x4.x7 + x5.x8 + x6.x9
f123 =~ x1.x4.x7 + x2.x5.x8 + x3.x6.x9
f4 =~ x10 + x11 + x12
f4 ~ f1 + f2 + f3 + f12 + f13 + f23 + f123
f1 ~~ 0*f12
f1 ~~ 0*f13
f1 ~~ 0*f123
f2 ~~ 0*f12
f2 ~~ 0*f23
f2 ~~ 0*f123
f3 ~~ 0*f13
f3 ~~ 0*f23
f3 ~~ 0*f123
f12 ~~ 0*f123
f13 ~~ 0*f123
f23 ~~ 0*f123
x1 ~ 0*1
x4 ~ 0*1
x7 ~ 0*1
x10 ~ 0*1
x1.x4 ~ 0*1
x1.x7 ~ 0*1
x4.x7 ~ 0*1
x1.x4.x7 ~ 0*1
f1 ~ NA*1
f2 ~ NA*1
f3 ~ NA*1
f12 ~ NA*1
f13 ~ NA*1
f23 ~ NA*1
f123 ~ NA*1
f4 ~ NA*1
" 

fitRC3way <- sem(model3, data=dat3wayRC, meanstructure=TRUE, std.lv=FALSE)
summary(fitRC3way)
result3wayRC <- probe3WayRC(fitRC3way, c("f1", "f2", "f3", "f12", "f13", "f23", "f123"), "f4", c("f1", "f2"), c(-1, 0, 1), c(-1, 0, 1))
result3wayRC

fitMC3way <- sem(model3, data=dat3wayMC, meanstructure=TRUE, std.lv=FALSE)
summary(fitMC3way)
result3wayMC <- probe3WayMC(fitMC3way, c("f1", "f2", "f3", "f12", "f13", "f23", "f123"), "f4", c("f1", "f2"), c(-1, 0, 1), c(-1, 0, 1))
result3wayMC

plotProbe(result3wayRC, xlim=c(-2, 2))
plotProbe(result3wayMC, xlim=c(-2, 2))

################################# outClipboard #############################

library(lavaan)
HW.model <- ' visual  =~ c("c1", "c1")*x1 + NA*x1 + c("c2", "c2")*x2 + c("c3", "c3")*x3
              textual =~ c("c4", "c4")*x4 + NA*x4 + c("c5", "c5")*x5 + c("c6", "c6")*x6
               speed   =~ c("c7", "c7")*x7 + NA*x7 + c("c8", "c8")*x8 + c("c9", "c9")*x9 
			   visual ~~ c(1, NA)*visual
			   textual ~~ c(1, NA)*textual
			   speed ~~ c(1, NA)*speed
			   
			   '

fit <- cfa(HW.model, data=HolzingerSwineford1939, group="school", meanstructure=TRUE)

copy(fit)
copy(fit, "mifit")
copy(fit, "coef")
copy(fit, "se")
copy(fit, "samp")
copy(fit, "fit")

################################ rmseaNested ##############################

# plotPower
# plotDist
# findPower
# findDist
# equiPower


################################# group allocation power #################



################################### Transcribe code ###################

library(lavaan)
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

fit <- cfa(HS.model, data=HolzingerSwineford1939)
parTable(fit)



########### Auxiliary

library(lavaan)

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
			  
dat <- data.frame(HolzingerSwineford1939, z=rnorm(nrow(HolzingerSwineford1939), 0, 1))
			  
fit <- cfa(HS.model, data=dat) #, group="sex", meanstructure=TRUE)
fitaux <- auxiliary(fit, aux="z", data=dat)

fitgroup <- cfa(HS.model, data=dat, group="school")
fitgroupaux <- auxiliary(fitgroup, aux="z", data=dat, group="school")

mod <- ' x5 ~ x4
x4 ~ x3
x3 ~ x1 + x2'

fitpath <- sem(mod, data=dat)
fitpathaux <- auxiliary(fitpath, aux="z", data=dat)

dat2 <- data.frame(PoliticalDemocracy, z=rnorm(nrow(PoliticalDemocracy), 0, 1))
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
fitsem <- sem(model, data=dat2, meanstructure=TRUE)
fitsemaux <- auxiliary(fitsem, aux="z", data=dat2, meanstructure=TRUE)


# HS.model.cov <- ' visual  =~ x1 + x2 + x3
              # textual =~ x4 + x5 + x6
              # speed   =~ x7 + x8 + x9 
			  # visual ~ sex
			  # textual ~ sex
			  # speed ~ sex'
	  
# fitcov <- sem(HS.model.cov, data=dat) #, group="sex", meanstructure=TRUE)
# as.data.frame(fitcov@ParTable)
# fitcovaux <- auxiliary(fitcov, aux="z", data=dat)

# HS.model.covxx <- ' visual  =~ x1 + x2 + x3
              # textual =~ x4 + x5 + x6
              # speed   =~ x7 + x8 + x9 
			  # visual ~ sex
			  # textual ~ sex
			  # speed ~ sex
			  # z ~~ z
			  # z ~~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
			  # z ~ sex'
	  
# fitcovxx <- sem(HS.model.covxx, data=dat) #, group="sex", meanstructure=TRUE)
# as.data.frame(fitcovxx@ParTable)
# fitcovaux <- auxiliary(fitcov, aux="z", data=dat)


# HS.model.cov2 <- ' visual  =~ x1 + x2 + x3
              # textual =~ x4 + x5 + x6
              # x7 ~ visual + textual'
	  
# fitcov2 <- sem(HS.model.cov2, data=dat, fixed.x=FALSE) #, group="sex", meanstructure=TRUE)
# as.data.frame(fitcov2@ParTable)
# fitcov2aux <- auxiliary(fitcov2, aux="z", data=dat)



# HS.model.covyy <- ' visual  =~ x1 + x2 + x3
              # textual =~ x4 + x5 + x6
              # x7 ~ visual + textual
			  # z ~~ x1 + x2 + x3 + x4 + x5 + x6 + x7'
# fitcovyy <- sem(HS.model.covyy, data=dat) #, group="sex", meanstructure=TRUE)
		  
		  
#################################### runMI function ###########################################

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



test<-HolzingerSwineford1939[,-5]
HS.model <- ' x1 ~ x4
			  x4 ~ x5 + x6'
summary(cfa(HS.model,data=test), fit.measures=TRUE)

test[(test$x6>3 && test$x6<4)]
##Impose missing data to test
log.mat1 <- matrix(FALSE, nrow=dim(test)[1], ncol=dim(test)[2])
log.mat1[,9] <- test$x6>3
test[log.mat1] <- NA

runMI(test,HS.model,10, idvars='id')


test<-HolzingerSwineford1939[,-5]
HS.model <- ' visual  =~ NA*x1 + x2 + x3
               textual =~ NA*x4 + x5 + x6
               speed   =~ NA*x7 + x8 + x9 
			   visual ~~1*visual
			   textual ~~ 1*textual
			   speed ~~ 1*speed
			   visual ~ speed + textual'
summary(cfa(HS.model,data=test), fit.measures=TRUE)

test[(test$x6>3 && test$x6<4)]
##Impose missing data to test
log.mat1 <- matrix(FALSE, nrow=dim(test)[1], ncol=dim(test)[2])
log.mat1[,9] <- test$x6>3
test[log.mat1] <- NA

runMI(test,HS.model,3, idvars='id')






test <- HolzingerSwineford1939[-301,-5]
data.model <- ' visual  =~ NA*x1 + x2 + x3
		   textual =~ NA*x4 + x5 + x6
		   speed   =~ NA*x7 + x8 + x9 
		   visual ~~1*visual
		   textual ~~ 1*textual
		   speed ~~ 1*speed
		   visual ~ speed'

Impose missing data to test
log.mat1 <- matrix(FALSE, nrow=dim(test)[1], ncol=dim(test)[2])
log.mat1[,9] <- test$x6>3
test[log.mat1] <- NA
data.mat <- test
m <- 20
miArgs=list(); chi="all"; miPackage="Amelia"; digits=3; seed=12345; std.lv = FALSE; estimator = "ML"; group = "grade"; group.equal = ""
fun <- "sem"
	
y <- runMI(data=test,model=data.model, m=10, fun="sem", meanstructure=TRUE) 
x <- sem.mi(data=test,model=data.model, m=10, meanstructure=TRUE)
