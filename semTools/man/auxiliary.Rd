\name{auxiliary}
\alias{auxiliary}
\title{
	Analyzing data with full-information maximum likelihood with auxiliary variables
}
\description{
	Analyzing data with full-information maximum likelihood with auxiliary variables. The techniques used to account for auxiliary variables are both extra-dependent-variables and saturated-correlates approaches (Enders, 2008). The extra-dependent-variables approach is used for exogenous varibles in the model (see \code{fixed.x} attribute when fitting the \code{\link[lavaan]{lavaan}} model). For other variables, the saturated-correlates approach is used. WARNINGS: this function does not work for all models. See the examples below.
}
\usage{
auxiliary(object, aux, ...)
}
\arguments{
  \item{object}{
	The \code{lavaan} object or the parameter table
}
  \item{aux}{
	The list of auxiliary variable
  }
  \item{\dots}{
	The additional arguments in the \code{\link[lavaan]{lavaan}} function.
  }
}
\value{
	The \code{\linkS4class{lavaanStar}} object which contains the original \code{lavaan} object and the additional values of the null model, which need to be adjusted to account for auxiliary variables.
}
\references{
Enders, C. K. (2008). A note of the use of missing auxiliary variables in full information maximum likelihood-based structural equation models. \emph{Structural Equation Modeling, 15}, 434-448.
}
\seealso{
	\code{\linkS4class{lavaanStar}}
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
}
\examples{
# Example of confirmatory factor analysis

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
			  
dat <- data.frame(HolzingerSwineford1939, z=rnorm(nrow(HolzingerSwineford1939), 0, 1))
			  
fit <- cfa(HS.model, data=dat) #, group="sex", meanstructure=TRUE)
fitaux <- auxiliary(fit, aux="z", data=dat)

# Example of multiple groups confirmatory factor analysis

fitgroup <- cfa(HS.model, data=dat, group="school")
fitgroupaux <- auxiliary(fitgroup, aux="z", data=dat, group="school")

# Example of path analysis

mod <- ' x5 ~ x4
x4 ~ x3
x3 ~ x1 + x2'

fitpath <- sem(mod, data=dat)
fitpathaux <- auxiliary(fitpath, aux="z", data=dat)

# Example of full structural equation modeling

dat2 <- data.frame(PoliticalDemocracy, z=rnorm(nrow(PoliticalDemocracy), 0, 1))
model <- ' 
     ind60 =~ x1 + x2 + x3
     dem60 =~ y1 + a*y2 + b*y3 + c*y4
     dem65 =~ y5 + a*y6 + b*y7 + c*y8

    dem60 ~ ind60
    dem65 ~ ind60 + dem60

    y1 ~~ y5
    y2 ~~ y4 + y6
    y3 ~~ y7
    y4 ~~ y8
    y6 ~~ y8
'
fitsem <- sem(model, data=dat2, meanstructure=TRUE)
fitsemaux <- auxiliary(fitsem, aux="z", data=dat2, meanstructure=TRUE)

#########################################################
## These following codes show the models that do not work with the current function

##### 1. covariate at the factor level
## HS.model.cov <- ' visual  =~ x1 + x2 + x3
##              textual =~ x4 + x5 + x6
##             speed   =~ x7 + x8 + x9 
##			  visual ~ sex
##			  textual ~ sex
##			  speed ~ sex'
	  
## fitcov <- cfa(HS.model.cov, data=dat) 
## fitcovaux <- auxiliary(fitcov, aux="z", data=dat)

### The auxiliary code does not work when specifying manually.
## HS.model.covxx <- ' visual  =~ x1 + x2 + x3
##              textual =~ x4 + x5 + x6
##              speed   =~ x7 + x8 + x9 
##			  visual ~ sex
##			  textual ~ sex
##			  speed ~ sex
##			  z ~~ z
##			  z ~~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
##			  z ~ sex'
	  
## fitcovxx <- cfa(HS.model.covxx, data=dat) 
## fitcovaux <- auxiliary(fitcov, aux="z", data=dat)

##### 2. Endogenous variable with single indicator 
## HS.model.cov2 <- ' visual  =~ x1 + x2 + x3
##               textual =~ x4 + x5 + x6
##               x7 ~ visual + textual'
## 	  
## fitcov2 <- sem (HS.model.cov2, data=dat, fixed.x=FALSE) #, group="sex", meanstructure=TRUE)
## fitcov2aux <- auxiliary(fitcov2, aux="z", data=dat)


### The auxiliary code does not work when specifying manually.
## HS.model.covyy <- ' visual  =~ x1 + x2 + x3
##               textual =~ x4 + x5 + x6
##               x7 ~ visual + textual
## 			  z ~~ x1 + x2 + x3 + x4 + x5 + x6 + x7'
## fitcovyy <- sem(HS.model.covyy, data=dat) #, group="sex", meanstructure=TRUE)
		  
}
