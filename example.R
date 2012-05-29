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