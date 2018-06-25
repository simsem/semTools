### Corbin Quick, Alex Schoemann, James Selig, Terrence D. Jorgnensen
### Last updated: 25 June 2018

# FIXME: work out a path-analysis example like slide 25:
# http://www.da.ugent.be/cvs/pages/en/Presentations/Presentation%20Yves%20Rosseel.pdf
# add example to help page, to illustrate a complex function of parameters


#' Monte Carlo Confidence Intervals to Test Complex Indirect Effects
#'
#' This function takes an expression for an indirect effect, the parameters and
#' standard errors associated with the expression and returns a confidence
#' interval based on a Monte Carlo test of mediation (MacKinnon, Lockwood, &
#' Williams, 2004).
#'
#' This function implements the Monte Carlo test of mediation first described
#' in MacKinnon, Lockwood, & Williams (2004) and extends it to complex cases
#' where the indirect effect is more than a function of two parameters. The
#' function takes an expression for the indirect effect, randomly simulated
#' values of the indirect effect based on the values of the parameters (and the
#' associated standard errors) comprising the indirect effect, and outputs a
#' confidence interval of the indirect effect based on the simulated values.
#' For further information on the Monte Carlo test of mediation see MacKinnon,
#' Lockwood, & Williams (2004) and Preacher & Selig (2012).
#'
#' The asymptotic covariance matrix can be easily found in many popular SEM
#' software applications.
#' \itemize{
#'  \item LISREL: Including the EC option on the OU line will print the ACM
#'   to a seperate file. The file contains the lower triangular elements of
#'   the ACM in free format and scientific notation
#'  \item Mplus Include the command TECH3; in the OUTPUT section. The ACM will be
#'   printed in the output.
#'  \item lavaan: Use the command \code{vcov} on the fitted lavaan object to
#'   print the ACM to the screen
#' }
#'
#'
#' @importFrom stats quantile
#'
#' @param expression A character scalar representing the computation of an
#' indirect effect. Different parameters in the expression should have
#' different alphanumeric values. Expressions can use either addition (+) or
#' multiplication (*) operators.
#' @param \dots Parameter estimates for all parameters named in
#' \code{expression}. The order of parameters should follow from
#' \code{expression} (the first parameter named in \code{expression} should be
#' the first parameter listed in \dots{}). Alternatively \dots can be a
#' vector of parameter estimates.
#' @param ACM A matrix representing the asymptotic covariance matrix of the
#' parameters described in \code{expression}. This matrix should be a symetric
#' matrix with dimensions equal to the number of parameters names in
#' \code{expression}. Information on finding the ACOV is popular SEM software
#' is described below.)
#' @param object A lavaan model object fitted after running the running the
#' \code{cfa}, \code{sem}, \code{growth}, or \code{lavaan} functions. The model
#' must have parameters labelled with the same labels used in
#' \code{expression}. When using this option do not specify values for \dots
#' or \code{ACM}
#' @param rep The number of replications to compute. Many thousand are
#' reccomended.
#' @param CI Width of the confidence interval computed.
#' @param plot Should the function output a plot of simulated values of the
#' indirect effect?
#' @param outputValues Should the function output all simulated values of the
#' indirect effect?
#' @return A list with two elements. The first element is the point estimate
#' for the indirect effect. The second element is a matrix with values for the
#' upper and lower limits of the confidence interval generated from the Monte
#' Carlo test of mediation. If \code{outputValues = TRUE}, output will be a list
#' with a list with the point estimate and values for the upper and lower
#' limits of the confidence interval as the first element and a vector of
#' simulated values of the indirect effect as the second element.
#' @author
#' Corbin Quick (University of Michigan; \email{corbinq@@umich.edu})
#'
#' Alexander M. Schoemann (East Carolina University; \email{schoemanna@@ecu.edu})
#'
#' James P. Selig (University of New Mexico; \email{selig@@unm.edu})
#'
#' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
#' @references
#' MacKinnon, D. P., Lockwood, C. M., & Williams, J. (2004). Confidence limits
#' for the indirect effect: Distribution of the product and resampling methods.
#' \emph{Multivariate Behavioral Research, 39}(1) 99--128.
#' doi:10.1207/s15327906mbr3901_4
#'
#' Preacher, K. J., & Selig, J. P. (2010, July). Monte Carlo method
#' for assessing multilevel mediation: An interactive tool for creating
#' confidence intervals for indirect effects in 1-1-1 multilevel models
#' [Computer software]. Available from \url{http://quantpsy.org/}.
#'
#' Preacher, K. J., & Selig, J. P. (2012). Advantages of Monte Carlo confidence
#' intervals for indirect effects. \emph{Communication Methods and Measures,
#' 6}(2), 77--98. doi:10.1080/19312458.2012.679848
#'
#' Selig, J. P., & Preacher, K. J. (2008, June). Monte Carlo method for
#' assessing mediation: An interactive tool for creating confidence intervals
#' for indirect effects [Computer software]. Available from
#' \url{http://quantpsy.org/}.
#' @examples
#'
#' ## Simple two path mediation
#' ## Write expression of indirect effect
#' med <- 'a*b'
#' ## Paramter values from analyses
#' aparam <- 1
#' bparam <- 2
#' ## Asymptotic covariance matrix from analyses
#' AC <- matrix(c(.01,.00002,
#'                .00002,.02), nrow=2, byrow=TRUE)
#' ## Compute CI, include a plot
#' monteCarloMed(med, coef1 = aparam, coef2 = bparam, outputValues = FALSE,
#'               plot = TRUE, ACM = AC)
#'
#' ## Use a vector of parameter estimates as input
#' aparam <- c(1,2)
#' monteCarloMed(med, coef1 = aparam, outputValues = FALSE,
#'               plot = TRUE, ACM = AC)
#'
#'
#' ## Complex mediation with two paths for the indirect effect
#' ## Write expression of indirect effect
#' med <- 'a1*b1 + a1*b2'
#' ## Paramter values and standard errors from analyses
#' aparam <- 1
#' b1param <- 2
#' b2param <- 1
#' ## Asymptotic covariance matrix from analyses
#' AC <- matrix(c(1, .00002, .00003,
#'                .00002, 1, .00002,
#' 					      .00003, .00002, 1), nrow = 3, byrow = TRUE)
#' ## Compute CI do not include a plot
#' monteCarloMed(med, coef1 = aparam, coef2 = b1param,
#'               coef3 = b2param, ACM = AC)
#'
#' @export
monteCarloMed <- function(expression, ..., ACM = NULL, object = NULL,
                          rep = 20000, CI = 95, plot = FALSE,
                          outputValues = FALSE) {

  input <- c(...)

  ## Get names and the number of unique variables in the expression
  paramnames <- all.vars(stats::as.formula(paste("~", expression)))

  ## If input is a lavaan object pull out coefs and ACM
  if (class(object) == "lavaan"){
    input <- lavaan::coef(object)[paramnames]
    ACM <- lavaan::vcov(object)[paramnames, paramnames]
  }

  vecs <- list()
  ## Matrix of values, need to be converted to a list
  dat <- MASS::mvrnorm(n = rep, mu = input, Sigma = ACM)
  ## Add parameters as the first row
  dat <-rbind(input, dat)
  ## Convert to a list,
  vecs <- as.list(as.data.frame(dat))
  ## Give names to it works with assign
  for (i in 1:length(vecs)){assign(paramnames[i], vecs[[i]])}

  ## Apply the expression to compute the indirect effect
  indirect <- eval(parse(text = expression))
  ## Get the CI
  low  <- (1-CI/100)/2
  upp <- ((1-CI/100)/2) + (CI/100)
  LL <- round(quantile(indirect[-1], low), digits = 4)
  UL <- round(quantile(indirect[-1], upp), digits = 4)
  interval <- list(indirect[1], rbind(LL,UL))
  dimnames(interval[[2]]) <- list(c("LL", "UL"), c(" "))
  names(interval) <- c("Point Estimate",
                       paste(CI, "% Confidence Interval", sep = ""))

  ## Switch for outputting a plot
  if (plot) {
    hist(indirect, breaks = 'FD', col = 'skyblue',
         xlab = paste(CI, '% Confidence Interval ', 'LL', LL, '  UL', UL),
         main = 'Distribution of Indirect Effect')
  }

  ## Switch to return simulated values
  if (outputValues) {
    interval <- list(interval, indirect)
  }

  return(interval)
}


