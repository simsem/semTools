\name{permuteMeasEq}
\alias{permuteMeasEq}
\title{
	Permutation Randomization Tests of Measurement Equivalence and Differential Item Functioning (DIF)
}
\description{
This function accepts a pair of nested lavaan objects, the less constrained of which freely estimates a set of measurement parameters (e.g., factor loadings, intercepts, or thresholds) in all groups, and the more constrained of which constrains those measurement parameters to equality across groups. Group assignment is repeatedly permuted and the model is fit to each permutation, in order to produce an empirical distribution of (a) changes in fit indices and (b) differences in measurement parameters, for which the null hypothesis of no group differences is true. This function is for testing measurement equivalence only across groups, not occasions.
}
\usage{
permuteMeasEq(nPermute, uncon, con, null = NULL, AFIs = NULL, moreAFIs = NULL,
              param = "loadings", maxSparse = 10, maxNonconv = 10)
}
\arguments{
  \item{nPermute}{
	An integer indicating the number of random permutations of group assignment used to form empirical distributions under the null hypothesis.
}
  \item{uncon}{
	The unconstrained \code{lavaan} object, in which a set of measurement parameters of interest (e.g., factor loadings for testing metric/weak invariance, item intercepts for testing scalar/strong invariance) are freely estimated in all groups.
}
  \item{con}{
  The constrained \code{lavaan} object, in which the measurement parameters of interest in \code{uncon} are constrained to equality across all groups.
}
  \item{null}{
  Optional.  A \code{lavaan} object, in which an alternative null model is fit (besides the default independence model specified by \code{lavaan}) for the calculation of incremental fit indices. See Widamin & Thompson (2003) for details. If \code{NULL}, \code{lavaan}'s default independence model is used.
}
  \item{AFIs}{
  A character vector indicating which alternative fit indices, returned by \code{lavaan::fitMeasures}, is to be used to test the global null hypothesis of no group differences in any parameters. If \code{NULL}, the default \code{AFIs} will be returned: \code{c("cfi","rni","tli","rmsea","srmr","mfi")}.
}
  \item{moreAFIs}{
  A character vector indicating which alternative fit indices, returned by \code{semTools::moreFitIndices}, is to be used to test the global null hypothesis of no group differences in any parameters. If \code{NULL}, the default \code{moreAFIs} will be returned: \code{c("gammaHat","adjGammaHat")}.
}
  \item{param}{
  A character vector indicating which parameters are to be tested for significant DIF. Parameter names must match those returned by \code{names(coef(uncon))}, but omitting any group-specific suffixes (e.g., \code{"f1~1"} rather than \code{"f1~1.g2"}). Also, do not specify user-defined labels from model syntax, but only lavaan parameter names. Alternatively, to test an entire set of a certain type of parameter, \code{param} may take any one of the following values: \code{"loadings"}, \code{"intercepts"}, \code{"thresholds"}, or \code{"residuals"}.
}
  \item{maxSparse}{
  An integer indicating the maximum number of consecutive times that randomly permuted group assignment can yield a sample in which at least one category (of an \code{ordered} indicator) is unobserved in at least one group, such that the same set of parameters cannot be estimated in each group. If such a sample occurs, group assignment is randomly permuted again, repeatedly until a sample is obtained with all categories observed in all groups. If \code{maxSparse} is exceeded, \code{NA} will be returned for that iteration of the permutation distribution.
}
  \item{maxNonconv}{
  An integer indicating the maximum number of consecutive times that randomly permuted group assignment can yield a sample for which the model does not converge on a solution. If such a sample occurs, group assignment is randomly permuted again, repeatedly until a sample is obtained for which the model does converge. If \code{maxNonconv} is exceeded, \code{NA} will be returned for that iteration of the permutation distribution, and a warning will be printed when using \code{show} or \code{summary}.
}
}
\details{
	The global null hypothesis of measurement equivalence/invariance is that there are no group differences in any measurement parameters. This can be tested using the \code{anova} method on nested lavaan objects, as seen in the output of \code{\link[semTools]{measurementInvariance}}, or by inspecting the change in alternative fit indices (AFIs) such as the CFI. See Cheung & Rensvold (2002) and Meade, Johnson, & Braddy (2008) for details.
  
  If the global null hypothesis is rejected, partial invariance can still be established by freeing parameters that differ across groups, while maintaining equality constraints for other indicators. DIF can be estimated using the less constrained model, in which the parameters are allowed to differ, but multiple testing leads to inflation of Type I error rates. The permutation randomization method creates a distribution of the maximum absolute-value of DIF if the null hypothesis is true, similar to Tukey's \emph{q} distribution for the Honestly Significant Difference (HSD) post hoc test, which allows the researcher to control the familywise Type I error rate.
}
\value{
	The \linkS4class{MeasEq.permute} object representing the results of testing measurement equivalence (global null hypothesis) and DIF (tests of individaul measurement parameters).
}
\author{
    Terrence D. Jorgensen (University of Kansas; \email{TJorgensen314@gmail.com})
}
\references{
Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit indexes for testing measurement invariance. \emph{Structural Equation Modeling, 9}(2), 233–255. doi:10.1207/S15328007SEM0902_5

Meade, A. W., Johnson, E. C., & Braddy, P. W. (2008). Power and sensitivity of alternative fit indices in tests of measurement invariance. \emph{Journal of Applied Psychology, 93}(3), 568–592. doi:10.1037/0021-9010.93.3.568

Widamin, K. F., & Thompson, J. S. (2003). On specifying the null model for incremental fit indices in structural equation modeling. \emph{Psychological Methods, 8}(1), 16–37. doi:10.1037/1082-989X.8.1.16
}
\seealso{
\code{\link[stats]{TukeyHSD}}, \code{\link[semTools]{measurementInvariance}}, \code{\link[semTools]{measurementInvarianceCat}}
}
\examples{
## fit indices of interest for global test of measurement equivalence
myAFIs <- c("cfi","rni","tli","rmsea","srmr","mfi","gfi","aic","bic")
moreAFIs <- c("gammaHat","adjGammaHat")

## run models to be compared
mod.null <- c(paste0("x", 1:9, " ~~ c(L", 1:9, ", L", 1:9, ")*x", 1:9),
              paste0("x", 1:9, " ~ c(T", 1:9, ", T", 1:9, ")*1"))
fit.null <- cfa(mod.null, data = HolzingerSwineford1939, group = "sex")

mod.config <- '
visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9
'

fit.config <- cfa(mod.config, data = HolzingerSwineford1939,
                  std.lv = TRUE, group = "sex")
AFI.config <- c(fitMeasures(fit.config, fit.measures = myAFIs,
                            baseline.model = fit.null),
                moreFitIndices(fit.config, fit.measures = moreAFIs))

mod.metric <- '
visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9
visual ~~ c(1, NA)*visual
textual ~~ c(1, NA)*textual
speed ~~ c(1, NA)*speed
'
fit.metric <- cfa(mod.metric, data = HolzingerSwineford1939,
                  std.lv = TRUE, group = "sex", group.equal = "loadings")
AFI.metric <- c(fitMeasures(fit.metric, fit.measures = myAFIs,
                            baseline.model = fit.null),
                moreFitIndices(fit.metric, fitmeasures = moreAFIs))

## calculate observed differences in fit indices
myDiffs <- AFI.config - AFI.metric
myDiffs
## compare to ANOVA result
anova(fit.config, fit.metric)


## Use only 20 permutations for a demo (in practice, use > 500)
set.seed(12345)
out <- permuteMeasEq(nPermute = 20, uncon = fit.config, con = fit.metric,
                     AFIs = myAFIs, moreAFIs = moreAFIs, param = "loadings",
                     null = fit.null) 

## Results object contains info about global and follow-up tests

## The "show" method gives results for global null hypothesis of equivalence
out

## The "summary" method gives details about "follow-up" DIF tests
summary(out, digits = 2)

## Not much to see with only one significant DIF. To illustrate, use individual
## maximum-DIF distributions for each parameter, and raise the alpha level.
summary(out, type = "each", alpha = .10, digits = 2)
}

