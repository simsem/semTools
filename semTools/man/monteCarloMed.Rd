\name{monteCarloMed}
\alias{monteCarloMed}
\title{
Monte Carlo Confidence Intervals to Test Complex Indirect Effects
}
\description{
This function takes an expression for an indirect effect, the parameters and standard errors associated with the expression and returns a confidence interval based on a Monte Carlo test of mediation (MacKinnon, Lockwood, & Williams, 2004).
}
\usage{
monteCarloMed(expression, ..., ACM=NULL, object=NULL, rep=20000, CI=95, plot=FALSE, 
  outputValues=FALSE)
}
\arguments{
  \item{expression}{A character scalar representing the computation of an indirect effect. Different parameters in the expression should have different alphanumeric values. Expressions can use either addition (+) or multiplication (*) operators.}
  \item{\dots}{Parameter estimates for all parameters named in \code{expression}. The order of parameters should follow from \code{expression} (the first parameter named in \code{expression} should be the first parameter listed in \dots). Alternatively \dots can be a vector of parameter estimates.}
  \item{ACM}{A matrix representing the asymptotic covariance matrix of the parameters described in \code{expression}. This matrix should be a symetric matrix with dimensions equal to the number of parameters names in \code{expression}. Information on finding the ACOV is popular SEM software is described below.)}
  \item{object}{A lavaan model object fitted after running the running the \code{cfa}, \code{sem}, \code{growth}, or \code{lavaan} functions. The model must have parameters labelled with the same labels used in \code{expression}. When using this option do not specify values for \dots or \code{ACM}}
  \item{rep}{The number of replications to compute. Many thousand are reccomended.}	
  \item{CI}{Width of the confidence interval computed.}
  \item{plot}{Should the function output a plot of simulated values of the indirect effect?}
  \item{outputValues}{Should the function output all simulated values of the indirect effect?}
  }
\details{
This function implements the Monte Carlo test of mediation first described in MacKinnon, Lockwood, & Williams (2004) and extends it to complex cases where the indirect effect is more than a function of two parameters. The function takes an expression for the indirect effect, randomly simulated values of the indirect effect based on the values of the parameters (and the associated standard errors) comprising the indirect effect, and outputs a confidence interval of the indirect effect based on the simulated values. For further information on the Monte Carlo test of mediation see MacKinnon, Lockwood, & Williams (2004), Preacher & Selig (in press), and Selig & Preacher (2008). For a Monte Carlo test of mediation with a random effects model see Selig & Preacher (2010).

The asymptotic covariance matrix can be easily found in many popular SEM software applications.
	\itemize{
		\item{LISREL}{Including the EC option on the OU line will print the ACM to a seperate file. The file contains the lower triangular elements of the ACM in free format and scientific notation}
		\item{Mplus}{Include the command TECH3; in the OUTPUT section. The ACM will be printed in the output.}
		\item{lavaan} {Use the command \code{vcov} on the fitted lavaan object to print the ACM to the screen}
	}
}
\value{
A list with two elements. The first element is the point estimate for the indirect effect. The second element is a matrix with values for the upper and lower limits of the confidence interval generated from the Monte Carlo test of mediation. If \code{outputValues=TRUE}, output will be a list with a list with  the point estimate and values for the upper and lower limits of the confidence interval as the first element and a vector of simulated values of the indirect effect as the second element.
}
\references{
Preacher, K. J., & Selig, J. P. (2010, July). Monte Carlo method for assessing multilevel mediation: An interactive tool for creating confidence intervals for indirect effects in 1-1-1 multilevel models [Computer software]. Available from \url{http://quantpsy.org/}.

Preacher, K. J., & Selig, J. P. (2012). Advantages of Monte Carlo confidence intervals for indirect effects. \emph{Communication Methods and Measures, 6}, 77-98.

Selig, J. P., & Preacher, K. J. (2008, June). Monte Carlo method for assessing mediation: An interactive tool for creating confidence intervals for indirect effects [Computer software]. Available from \url{http://quantpsy.org/}.

}
\author{
    Corbin Quick (University of Michigan; \email{corbinq@umich.edu})
	Alexander M. Schoemann (East Carolina University; \email{schoemanna@ecu.edu})
	James P. Selig (University of New Mexico; \email{selig@unm.edu})
}
\examples{
#Simple two path mediation
#Write expression of indirect effect
med <- 'a*b'
#Paramter values from analyses
aparam <- 1
bparam<-2
#Asymptotic covariance matrix from analyses
AC <- matrix(c(.01,.00002,
               .00002,.02), nrow=2, byrow=TRUE)
#Compute CI, include a plot
monteCarloMed(med, coef1=aparam, coef2=bparam, outputValues=FALSE, plot=TRUE, ACM=AC)

#Use a vector of parameter estimates as input
aparam<-c(1,2)
monteCarloMed(med, coef1=aparam, outputValues=FALSE, plot=TRUE, ACM=AC)



#Complex mediation with two paths for the indirect effect
#Write expression of indirect effect
med <- 'a1*b1 + a1*b2'
#Paramter values and standard errors from analyses
aparam <- 1
b1param<-2
b2param<-1
#Asymptotic covariance matrix from analyses
AC <- matrix(c(1,.00002, .00003,
                    .00002,1, .00002,
					.00003, .00002, 1), nrow=3, byrow=TRUE)
#Compute CI do not include a plot
monteCarloMed(med, coef1=aparam, coef2=b1param, coef3=b2param, ACM=AC)
}
