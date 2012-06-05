\name{monteCarloMed}
\alias{monteCarloMed}
Monte Carlo Confidence Intervals to Test Complex Indirect Effects
}
\description{
This function takes an expression for an indirect effect, the parameters and standard errors associated with the expression and returns a confidence interval based on a Monte Carlo test of mediation (MacKinnon, Lockwood, & Williams, 2004).
}
\usage{
monteCarloMed(expression, ..., rep=20000, CI=95, plot=FALSE, outputValues=FALSE, ACOV=NULL)
}
\arguments{
  \item{expression}{A character scalar representing the computation of an indirect effect. Different parameters in the expression should have different alphanumeric values. Expressions can use either addition (+) or multiplication (*) operators.}
  \item{\dots}{Parameter estimates and standard errors for all parameters named in the expression. The order of parameters should follow from the expression (the first parameter named in the expression should be the first parameter listed in \dots). Paramters and standard errors should be grouped together (e.g., parameter1, standard error1, parameter 2, standard error2, ...). If the asymptotic covariance matrix is provided in ACOV, the \dots should only contain parameter estimates.}
  \item{rep}{The number of replications to compute. Many thousand are reccomended.}	
  \item{CI}{Width of the confidence interval computed.}
  \item{plot}{Should the function output a plot of simulated values of the indirect effect?}
  \item{outputValues}{Should the function output all simulated values of the indirect effect?}
  \item{ACOV} {A matrix representing the asymptotic covariance matrix of the parameters described in the expression. If the ACOV is provided, then \dots should only contain parameter estimates, not standard errors)
  }
\details{
This function implements the Monte Carlo test of mediation first described in MacKinnon, Lockwood, & Williams (2004) and extends it to complex cases where the indirect effect is more than a function of two parameters. The function takes an expression for the indirect effect, randomly simulated values of the indirect effect based on the values of the parameters (and the associated standard errors) comprising the indirect effect, and outputs a confidence interval of the indirect effect based on the simulated values. For further information on the Monte Carlo test of mediation see MacKinnon, Lockwood, & Williams (2004), Preacher & Selig (in press), and Selig & Preacher (2008). For a Monte Carlo test of mediation with a random effects model see Selig & Preacher (2010).
}
\value{
A matrix with values for the upper and lower limits of the confidence interval generated from the Monte Carlo test of mediation. If outputValues=TRUE, output will be a list with a matrix with values for the upper and lower limits of the confidence interval as the first element and a vector of simulated values of the indirect effect as the second element.
}
\references{
Preacher, K. J., & Selig, J. P. (in press). Advantages of Monte Carlo confidence intervals for indirect effects. \emph{Communication Methods and Measures}.

Selig, J. P., & Preacher, K. J. (2008, June). Monte Carlo method for assessing mediation: An interactive tool for creating confidence intervals for indirect effects [Computer software]. Available from \url{http://quantpsy.org/}.

Preacher, K. J., & Selig, J. P. (2010, July). Monte Carlo method for assessing multilevel mediation: An interactive tool for creating confidence intervals for indirect effects in 1-1-1 multilevel models [Computer software]. Available from \url{http://quantpsy.org/}.
}
\authors{
    Corbin Quick (University of Kansas; \email{corbinq@ku.edu})
	Alexander M. Schoemann (University of Kansas; \email{schoemann@ku.edu})
	James P. Selig (University of New Mexico; \email{selig@unm.edu}}
	Kristopher J. Preacher {Vanderbilt University; \email{kris.preacher@vanderbilt.edu}}
}
\examples{
#Simple two path mediation
#Write expression of indirect effect
med <- 'a*b'
#Paramter values and standard errors from analyses
aparam <- 1
ase <- .01
bparam<-2
bse <- .02
#Compute CI include a plot
monteCarloMed(med, coef1=aparam, se1=ase, coef2=bparam, se2=bse, plot=T)


#complex mediation with two paths for the indirect effect
#Write expression of indirect effect
med <- 'a1*b1 + a1*b2'
#Paramter values and standard errors from analyses
aparam <- 1
ase <- .5
b1param<-2
b1se <- .2
b2param<-1
b2se<-.005
#Compute CI do not include a plot
monteCarloMed(med, coef1=aparam, se1=ase, coef2=b1c, se2=b1se, coef3=b2c, se3=b2se)
}
