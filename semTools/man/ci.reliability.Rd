\name{ci.reliability}
\alias{ci.reliability}
\title{ Confidence Interval for a Reliability Coefficient }
\description{
A function to calculate the confidence interval for a reliability coefficient: coefficient alpha or coefficient omega. 
}
\usage{
ci.reliability(data = NULL, S = NULL, N = NULL, aux = NULL, 
  type = NULL, inttype = NULL, B = 1000, conf.level = 0.95)
}
\arguments{
  \item{data}{The dataset that the reliability coefficient is obtained from. Real data set is required for categorical omega. Also, real data set is required for bootstrap confidence intervals or asymptotic distribution free confidence interval }
  \item{S}{Symmetric covariance matrix. Correlation matrix can be specified here but not recommended because, in the function, Confirmatory Factor Analysis (CFA) is analyzed based on covariance matrix. }
  \item{N}{The total sample size. Sample size is needed only that \code{S} is specified.}
  \item{aux}{The names of auxiliary variables. Auxiliary variables will not be used as a composite but they will use to handle missing observations. Note that full information maximum likelihood is used if auxiliary variables are specified. See \link[semTools]{auxiliary} for further details.}
  \item{type}{The type of reliability coefficient to be calculated: \code{"alpha"} or \code{1} for coefficient alpha analyzed by the formula proposed by Cronbach (1951), \code{"alpha-cfa"} or \code{2} for coefficient alpha analyzed by CFA with tau-equivalence (method of estimator depending on confidence interval method but none of them is unweighted least square so technically the result is not equal to the formula from Cronbach), \code{"omega"} for coefficient omega, \code{"hierarchical"} for hierarchical omega, \code{"categorical"} for categorical omega. The default is to use hierarchical omega for continuous items and categorical omega for categorical items.}
  \item{inttype}{There are 13 options for the methods. See \code{details} below. The default is to not provide any interval estimates. Based on our simulation studies (Kelley and Pornprasertmanit, in press), bias corrected and accelerated bootstrap, \code{"bca"}, is recommended for categorical omega. Any bootstrap approaches (e.g., \code{"bca"} or \code{"perc"}) are recommended for hierarchical omega, coefficient omega, and coefficient alpha.}
  \item{B}{the number of bootstrap replications}
  \item{conf.level}{ the confidence level (i.e., 1-Type I error rate)}
}
\details{
When the coefficient alpha is used, the measurement model is assumed to be true-score equivalent (or tau equivalent) model such that factor loadings are equal across items. When the coefficient omega, hierarchical omega, and categorical omega are used, the measurement model is assumed to be congeneric model (similar to one-factor confirmatory factor analysis model). Coefficient omega assumes that a model fits data perfectly so the variance of the composite scores is calculated from model-implied covariance matrix. However, hierarchical omega allows a model to not fit data perfectly (Kelley and Pornprasertmanit, in press). Categorical omega is a method to calculate coefficient omega for categorical items (Green and Yang, 2009). That is, categorical omega is estimated by the parameter estimates from CFA for categorical items. If coefficient omega or hierarchical omega is used, CFA for continuous items is used, which is not appropriate for categorical items.

If researchers wish to make the measurement model with all parallel items (equal factor loadings and equal error variances), users can specify it by setting \code{inttype = "parallel"} and \code{type = "alpha"} or \code{type = "alpha-cfa"}. See McDonald (1999) for the assumptions of each of these models.  

The list below shows all methods to find the confidence interval of reliability. 

\enumerate{
	\item \code{"none"} or \code{0} to not find any confidence interval
	\item \code{"parallel"} or \code{11} to assume that the items are parallel and analyze confidence interval based on wald confidence interval (see van Zyl, Neudecker, & Nel, 2000, Equation 22; also referred as the asymptotic method of Koning & Franses, 2003).
	\item \code{"feldt"} or \code{12} is based on that \eqn{\frac{1 - \alpha}{1 - \hat{alpha}}} is distributed as \eqn{F} distribution with the degree of freedoms of \eqn{N - 1} and \eqn{(N - 1) \times (p - 1)} (Feldt, 1965).
	\item \code{"siotani"} or \code{13} is the same as the \code{"feldt"} method but using the degree of freedoms of \eqn{N} and \eqn{N \times (p - 1)} (Siotani, Hayakawa, & Fujikoshi, 1985; van Zyl et al., 2000, Equations 7 and 8; also referred as the exact method of Koning & Franses, 2003).
	\item \code{"fisher"} or \code{21} for the Fisher's \eqn{z} transformation on the correlation coefficient approach, \eqn{z = 0.5 \times \log{\frac{1 + \alpha}{1 - \alpha}}}, directly on the coefficient alpha and find confidence interval of transformed scale (Fisher, 1950). The variance of the \eqn{z} is \eqn{\frac{1}{N - 3}} where \eqn{N} is the total sample size.
	\item \code{"bonett"} or \code{22} for the Fisher's \eqn{z} transformation on the intraclass correlation approach with the variance of \eqn{\frac{2p}{(N - 2)(p - 1)}} (Bonett, 2002, Equation 6).
	\item \code{"hakstian"} or \code{23} uses the cube root transformation and assumes normal distribution on the cube root transformation (Hakstian & Whalen, 1976). The variance of the transformed reliability is based on the degrees of freedom in the \code{"feldt"} method.
	\item \code{"hakstianbarchard"} or \code{24} uses a correction of the violation of compound symmetry of covariance matrix by adjusting the degrees of freedom in the \code{"hakstian"}. This correction is used for the inference in type 12 sampling (both persons and items are sampled from the population of persons and items) See Hakstian and Barchard (2000) for further details.
	\item \code{"icc"} or \code{25} for the Fisher's \eqn{z} transformation on the intraclass correlation approach, \eqn{z = \log{1 - \alpha}}. The variance of the \eqn{z} is \eqn{\frac{2p}{N(p - 1)}} where \eqn{p} is the number of items (Fisher, 1991, p. 221; van Zyl et al., 2000, p. 277).
	\item \code{"ml"} or \code{31} to analyze the confidence interval based on normal-theory approach (or multivariate delta method). See van Zyl, Neudecker, & Nel (2000, Equation 21) for the confidence interval of coefficient alpha (also be referred as Iacobucci & Duhachek's, 2003, method). See Raykov (2002) for details for coefficient omega. If users use \code{analytic.type="cfa"}, the \code{sem} package will be used to obtain parameter estimates and standard errors used for the formula proposed by Raykov (2002)
	\item \code{"mll"} or \code{32} to analyze the confidence interval based on normal-theory approach as above. However, the point estimate and standard error were used to build confidence interval using logistic transformation as the note below.
	\item \code{"mlr"} or \code{33} to analyze the confidence interval based on normal-theory approach (or multivariate delta method). However, the estimation method is to use robust standard error (Satorra and Bentler, 2000). 
	\item \code{"mlrl"} or \code{34} to analyze he confidence interval based on normal-theory approach using robust standard error and logistic transformation (see below). 
	\item \code{"adf"} or \code{35} for asymtotic distribution-free method (see Maydeu-Olivares, Coffman, & Hartman, 2007 for further details for coefficient omega; we use phantom variable approach, Cheung, 2009, and \code{"WLS"} estimator for coefficient omega, Browne, 1984, in the \code{lavaan} package, Rosseel, 2012).
	\item \code{"adfl"} or \code{36} to use asymptotic distribution-free method to derive standard error and parameter estimate. Then, logistic transformation is used to build confidence interval (see below).
	\item \code{"ll"} or \code{37} for profile likelihood-based confidence interval of both reliability coefficients (Cheung, 2009) analyzed by the \code{OpenMx} package (Boker et al., 2011)
	\item \code{"bsi"} or \code{41} for standard bootstrap confidence interval which find the standard deviation across the bootstrap estimates, multiply the standard deviation by critical value, and add and subtract from the reliability estimate.
	\item \code{"bsil"} or \code{42} to use standard bootstrap confidence interval. However, logistic transformation is used to build confidence interval.
	\item \code{"perc"} or \code{43} for percentile bootstrap confidence interval.
	\item \code{"bca"} or \code{44} for bias-corrected and accelerated bootstrap confidence interval.
}

The logistic transformation (Browne, 1982) is applicable for \code{"ml"}, \code{"mlr"}, \code{"adf"}, and \code{"bsi"} as \code{"mll"}, \code{"mlrl"}, \code{"adfl"}, and \code{"bsil"}. The logistic transformation does not assume that the sampling distribution of reliability is symmetric. It acknowledges the fact that reliability ranges from 0 and 1. Logistic transformation is applied to the reliability estimates. Confidence interval is established for the transformed value. The lower and upper bounds of the transformed value is translated back to the reliability estimates. See Browne (1982) or Kelley and Pornprasertmanit (in press) for further details.

Note that not all confidence interval methods are available for all types of reliability and all types of input. For example, bootstrap confidence intervals are not available for covariance matrix input. Parallel confidence intervals are not available for hierarchical omega. We provided appropriate error messages for all impossible combinations. 
}
\value{
  \item{est}{The estimated reliability coefficient}
  \item{se}{The standard error of the reliability coefficient. If the bootstrap methods are used, this value represents the standard deviation across bootstrap estimates.}
  \item{ci.lower }{The lower bound of the computed confidence interval}
  \item{ci.upper}{The upper bound of the computed confidence interval}
  \item{conf.Level}{The confidence level (i.e., 1 - Type I error rate)}
  \item{type}{The type of estimated reliability coefficient (alpha or omega)}
  \item{inttype}{The method used to find confidence interval}
}
\references{ 
Boker, S., M., N., Maes, H., Wilde, M., Spiegel, M., Brick, T., et al. (2011).
OpenMx: An open source extended structural equation modeling framework.
\emph{Psychometrika, 76,} 306-317.

Bonett, D. G. (2002). Sample size requirements for testing and estimating coefficient
alpha. \emph{Journal of Educational and Behavioral Statistics, 27,} 335-340.

Browne, M. W. (1982). Covariance structures. In D. M. Hawkins (Ed.), \emph{Topics 
in applied multivariate analysis} (pp. 72-141). Cambridge, UK: Cambridge University Press.

Browne, M. W. (1984). Asymptotic distribution free methods in the analysis of
covariance structures. \emph{British Journal of Mathematical and Statistical
Psychology, 24,} 445-455.

Cheung, M. W.-L. (2009). Constructing approximate confidence intervals for
parameters with structural constructing approximate confidence intervals for
parameters with structural equation models. \emph{Structural Equation Modeling,
16,} 267-294.

Feldt, L.S. (1965). The approximate sampling distribution of Kuder-Richardson reliability coefficient twenty. \emph{Psychometrika,
30}, 357-370.

Fisher, R. A. (1950). \emph{Statistical methods for research workers.} Edinburgh, UK:
Oliver & Boyd.

Fisher, R. A. (1991). Statistical methods for research workers. In J.H. Bennett (Ed.), \emph{Statistical methods, experimental
design, and scientific inference.} Oxford: Oxford University Press.

Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using structural equation modeling: 
An alternative to coefficient alpha. \emph{Psychometrika, 74,} 155-167. 

Hakstian, A. R., & Whalen, T. E. (1976). A k-sample significance test for independent alpha
coefficients. \emph{Psychometrika, 41,} 219-231.

Iacobucci, D., & Duhachek, A. (2003). Advancing alpha: measuring reliability with confidence.
\emph{Journal of Consumer Psychology, 13,} 478-487.

Kelley, K. & Pornprasertmanit, P. (in press). Confidence intervals for population reliability coefficients: Evaluation of methods, recommendations, and software for homogeneous composite measures. \emph{Psychological Methods.}

Koning, A. J., & Franses, P. H. (2003). \emph{Confidence intervals for Cronbach&grave;s coefficient
alpha values} (ERIM Report Series Ref. No. ERS-2003-041-MKT). Rotterdam, The Netherlands:
Erasmus Research Institute of Management.

Maydeu-Olivares, A., Coffman, D. L., & Hartmann, W. M. (2007). Asymptotically
distribution-free (ADF) interval estimation of coefficient alpha. \emph{Psychological
Methods, 12,} 157-176.

McDonald, R. P. (1999). \emph{Test theory: A unified approach}. Mahwah, New Jersey: Lawrence 
Erlbaum Associates, Publishers.

Raykov, T. (2002). Analytic estimation of standard error and confidence interval
for scale reliability. \emph{Multivariate Behavioral Research, 37,} 89-103.

Rosseel, Y. (2012). lavaan: An R package for structural equation modeling. \emph{Journal of Statistical Software, 48,} 1-36.

Satorra, A. & Bentler, P. M. (2001). A scaled difference chi-square test statistic for moment structure analysis. \emph{Psychometrika, 66,} 507-514.

Siotani, M., Hayakawa, T., & Fujikoshi, Y. (1985). \emph{Modem multivariate statistical analysis: A graduate course and
handbook}. Columbus, Ohio: American Sciences Press.

van Zyl, J. M., Neudecker, H., & Nel, D. G. (2000) On the distribution of the maximum likelihood estimator 
of Cronbach's alpha. \emph{Psychometrika, 65} (3), 271-280.

Yuan, K. & Bentler, P. M. (2002) On robustness of the normal-theory based asymptotic distributions of three
reliability coefficient estimates. \emph{Psychometrika, 67} (2), 251-259.
}
\author{Sunthud Pornprasertmanit (\email{psunthud@gmail.com}). The previous version was written by Keke Lai (University of California, Merced; \email{klai25@ucmerced.edu}), Leann J. Terry, and Ken Kelley (University of Notre Dame; \email{kkelley@nd.edu}); }
\examples{
# Use this function for the attitude dataset (ignoring the overall rating variable)
ci.reliability(data=attitude[,-1], type = "omega", inttype = "mlrl", B = 100)

## Forming a hypothetical population covariance matrix
Pop.Cov.Mat <- matrix(.3, 9, 9)
diag(Pop.Cov.Mat) <- 1
ci.reliability(S=Pop.Cov.Mat, N=50, type="alpha", inttype = "bonett")
}
