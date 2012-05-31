\name{orthogonalize}
\alias{orthogonalize}
\title{
Orthogonalize data for 2-way interaction in SEM
}
\description{
Orthogonalize indicators of a 2-way interaction between latent variables
}
\usage{
orthogonalize(dat, xvars, zvars)
}
\arguments{
  \item{dat}{Matrix or data frame of item level data.}
  \item{xvars}{A vector of column numbers corresponding to indicators of the focal predictor (x).}
  \item{zvars}{A vector of column numbers corresponding to indicators of the moderator (z).}
  }
\details{
This functions will take a data frame or matrix and create orthogonalized product terms to compute latent variable interactions based on the method proposed by Little, Bovaird, & Widaman. The orthogonalized product terms can be entered into a SEM as indicators of a latent interaction variable. This function will compute all possible orthogonalized product terms (e.g., x has 3 indicators and z has 4 indicators, the function will return 3*4=12 new orthognalized product terms)
}

\value{
\enumerate{
    \item{data} Original data with orthogonalized product terms appended.
}
}
\references{

Little, T. D., Bovaird, J. A., & Widaman, K. F. (2006). On the merits of orthogonalizing powered and product terms: Implications for modeling interactions among latent variables. \emph{Structural Equation Modeling, 13} 497-519.

}
\author{
    Alexander M. Schoemann (University of Kansas; \email{schoemann@ku.edu})
}
\examples{

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

#Fit model in Lavaan
library(lavaan)

syntax <- ' 
x =~ x1 + x2 +x3
z =~ z1 + z2 + z3
xz =~ x1z1 + x1z2 + x1z3 + x2z1 + x2z2 + x2z3 + x3z1 + x3z2 + x3z3
y =~ y1 + y2 + y3
x ~~ z
x ~~ 0*xz
z ~~ 0*xz
y ~ x + z +xz
'

fit <- sem(model = syntax, data=datOrth, std.lv=TRUE)
summary(fit, fit.measures=TRUE)
}