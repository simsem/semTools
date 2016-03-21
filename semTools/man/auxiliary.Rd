\name{auxiliary}
\alias{auxiliary}
\alias{cfa.auxiliary}
\alias{sem.auxiliary}
\alias{growth.auxiliary}
\alias{lavaan.auxiliary}
\title{
	Analyzing data with full-information maximum likelihood with auxiliary variables
}
\description{
	Analyzing data with full-information maximum likelihood with auxiliary variables. The techniques used to account for auxiliary variables are both extra-dependent-variables and saturated-correlates approaches (Enders, 2008). The extra-dependent-variables approach is used for all single variables in the model (such as covariates or single-indicator dependent varaible) For variables that are belong to a multiple-indicator factor, the saturated-correlates approach is used. Note that all covariates are treated as endogenous varaibles in this model (fixed.x = FALSE) so multivariate normality is assumed for the covariates. CAUTION: (1) this function will automatically change the missing data handling method to full-information maximum likelihood and (2) this function is still not applicable for categorical variables (because the maximum likelhood method is not available in lavaan for estimating models with categorical variables currently).
}
\usage{
auxiliary(model, aux, fun, ...)
cfa.auxiliary(model, aux, ...)
sem.auxiliary(model, aux, ...)
growth.auxiliary(model, aux, ...)
lavaan.auxiliary(model, aux, ...)
}
\arguments{
  \item{model}{
	The \code{lavaan} object, the parameter table, or lavaan script. If the \code{lavaan} object is provided, the \code{lavaan} object must be evaluated with mean structure.
}
  \item{aux}{
	The list of auxiliary variable
  }
  \item{fun}{
	The character of the function name used in running lavaan model (\code{"cfa"}, \code{"sem"}, \code{"growth"}, \code{"lavaan"}). 
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
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
# Example of confirmatory factor analysis

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
			  
dat <- data.frame(HolzingerSwineford1939, z=rnorm(nrow(HolzingerSwineford1939), 0, 1))
			  
fit <- cfa(HS.model, data=dat, meanstructure=TRUE) 
fitaux <- auxiliary(HS.model, aux="z", data=dat, fun="cfa") # Use lavaan script
fitaux <- cfa.auxiliary(fit, aux="z", data=dat) # Use lavaan output

# Example of multiple groups confirmatory factor analysis

fitgroup <- cfa(HS.model, data=dat, group="school", meanstructure=TRUE)
fitgroupaux <- cfa.auxiliary(fitgroup, aux="z", data=dat, group="school")

\dontrun{
# Example of path analysis

mod <- ' x5 ~ x4
x4 ~ x3
x3 ~ x1 + x2'

fitpath <- sem(mod, data=dat, fixed.x=FALSE, meanstructure=TRUE) # fixed.x must be FALSE
fitpathaux <- sem.auxiliary(fitpath, aux="z", data=dat)

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
fitsemaux <- sem.auxiliary(fitsem, aux="z", data=dat2, meanstructure=TRUE)

# Example of covariate at the factor level

HS.model.cov <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
             speed   =~ x7 + x8 + x9 
			  visual ~ sex
			  textual ~ sex
			  speed ~ sex'
	  
fitcov <- cfa(HS.model.cov, data=dat, fixed.x=FALSE, meanstructure=TRUE) 
fitcovaux <- cfa.auxiliary(fitcov, aux="z", data=dat)

# Example of  Endogenous variable with single indicator 
HS.model.cov2 <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              x7 ~ visual + textual'
 	  
fitcov2 <- sem(HS.model.cov2, data=dat, fixed.x=FALSE, meanstructure=TRUE) 
fitcov2aux <- sem.auxiliary(fitcov2, aux="z", data=dat)

# Multiple auxiliary variables
HS.model2 <- ' visual  =~ x1 + x2 + x3
              speed   =~ x7 + x8 + x9'
fit <- cfa(HS.model2, data=HolzingerSwineford1939)
fitaux <- cfa.auxiliary(HS.model2, data=HolzingerSwineford1939, aux=c("x4", "x5")) 
}
}
