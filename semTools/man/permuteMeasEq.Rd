\name{permuteMeasEq}
\alias{permuteMeasEq}
\alias{calculateDIF}
\title{
	Permutation Randomization Tests of Measurement Equivalence and Differential Item Functioning (DIF)
}
\description{
The function \code{permuteMeasEq} accepts a pair of nested lavaan objects, the less constrained of which freely estimates a set of measurement parameters (e.g., factor loadings, intercepts, or thresholds) in all groups, and the more constrained of which constrains those measurement parameters to equality across groups. Group assignment is repeatedly permuted and the model is fit to each permutation, in order to produce an empirical distribution of (a) changes in fit indices and (b) differences in measurement parameters, for which the null hypothesis of no group differences is true. This function is for testing measurement equivalence only across groups, not occasions.

DIF estimates of specified parameters (\code{param}) are calculated from a single unconstrained model (\code{uncon}) using the function \code{calculateDIF}, which can be used to program more nuanced permutation methods if \code{permuteMeasEq} does not address the user's needs.
}
\usage{
permuteMeasEq(nPermute, uncon, con, null = NULL, AFIs = NULL, moreAFIs = NULL,
              param = "loadings", maxSparse = 10, maxNonconv = 10)
calculateDIF(uncon, param)
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
  A character vector indicating which alternative fit indices, returned by \code{lavaan::fitMeasures}, is to be used to test the multivariate omnibus null hypothesis of no group differences in any parameters. If \code{NULL}, the default \code{AFIs} will be returned: \code{c("cfi","rni","tli","rmsea","srmr","mfi")}. If \code{NA}, none of the AFIs from \code{lavaan::fitMeasures} will be returned.
}
  \item{moreAFIs}{
  A character vector indicating which alternative fit indices, returned by \code{semTools::moreFitIndices}, is to be used to test the multivariate omnibus null hypothesis of no group differences in any parameters. If \code{NULL}, the default \code{moreAFIs} will be returned: \code{c("gammaHat","adjGammaHat")}. If \code{NA}, none of the AFIs from \code{semTools::moreFitIndices} will be returned.
}
  \item{param}{
  A character vector indicating which parameters are to be tested for significant DIF. Parameter names must match those returned by \code{names(coef(uncon))}, but omitting any group-specific suffixes (e.g., \code{"f1~1"} rather than \code{"f1~1.g2"}). Alternatively, to test an entire set of a certain type of parameter, \code{param} may take any one of the following values: \code{"loadings"}, \code{"intercepts"}, \code{"thresholds"}, or \code{"residuals"}.
}
  \item{maxSparse}{
  An integer indicating the maximum number of consecutive times that randomly permuted group assignment can yield a sample in which at least one category (of an \code{ordered} indicator) is unobserved in at least one group, such that the same set of parameters cannot be estimated in each group. If such a sample occurs, group assignment is randomly permuted again, repeatedly until a sample is obtained with all categories observed in all groups. If \code{maxSparse} is exceeded, \code{NA} will be returned for that iteration of the permutation distribution.
}
  \item{maxNonconv}{
  An integer indicating the maximum number of consecutive times that randomly permuted group assignment can yield a sample for which the model does not converge on a solution. If such a sample occurs, group assignment is randomly permuted again, repeatedly until a sample is obtained for which the model does converge. If \code{maxNonconv} is exceeded, \code{NA} will be returned for that iteration of the permutation distribution, and a warning will be printed when using \code{show} or \code{summary}.
}
}
\details{
	The multivariate omnibus null hypothesis of measurement equivalence/invariance is that there are no group differences in any measurement parameters. This can be tested using the \code{anova} method on nested lavaan objects, as seen in the output of \code{\link[semTools]{measurementInvariance}}, or by inspecting the change in alternative fit indices (AFIs) such as the CFI. See Cheung & Rensvold (2002) and Meade, Johnson, & Braddy (2008) for details.

  If the multivariate omnibus null hypothesis is rejected using a global indicator of fit, partial invariance can still be established by freeing parameters that differ across groups, while maintaining equality constraints for other indicators. DIF can be estimated using the less constrained model, in which the parameters are allowed to differ, but multiple testing across items and across pairwise comparisons between groups) leads to inflation of Type I error rates. The permutation randomization method employed by \code{permuteMeasEq} creates a distribution of the maximum absolute-value of DIF if the null hypothesis is true, similar to Tukey's \emph{q} distribution for the Honestly Significant Difference (HSD) post hoc test, which allows the researcher to control the familywise Type I error rate. Two distributions are estimated: (1) the maximum absolute DIF across all pairwise comparisons within each item (using \code{type = "all"}), and (2) the maximum absolute DIF across all items and pairwise comparisons (using \code{type = "each"}).

  As an alternative to the multivariate omnibus test using global fit measures, \code{permuteMeasEq} also creates a distribution of the maximum sum-of-squared-DIF for each item. Thus, univariate omnibus tests can be performed for each item, with familywise Type I errors (across items) controlled analogously to Tukey's HSD. If the univariate omnibus null hypothesis for any item is rejected, then follow-up tests can be conducted using \code{type = "all"}.

  Because extreme observations become more likely to occur at least once as the number of tests grows, the familywise Type I error rate can only be controlled at the expense of power. Rather than trying to keep Type I errors close to a nominal alpha level, a linear step-up method (Maxwell & Delaney, 2004, pp. 230-234) can be used to keep the expected false discovery rate (FDR) close to a nominal alpha level. The expected FDR is the proportion of all rejected null hypotheses that are Type I errors. Controlling FDR guarantees controlling the familywise Type I error rate when all null hypotheses are true, but when at least one null hypothesis is false, controlling FDR is a compromise that keeps the number of incorrectly rejected null hypotheses (Type I errors) low relative to the total number of rejected null hypotheses, but with greater power to detect true group differences than methods that control the familywise Type I error rate. \code{permuteMeasEq} also employs the linear step-up procedure, and those results can be extracted by setting the argument \code{type = "step-up"}).
}
\value{
	The \linkS4class{permuteMeasEq} object representing the results of testing measurement equivalence (the multivariate omnibus null hypothesis) and DIF (tests of differences among individaul measurement parameters).
}
\author{
    Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\references{
Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit indexes for testing measurement invariance. \emph{Structural Equation Modeling, 9}(2), 233-255. doi:10.1207/S15328007SEM0902_5

Maxwell, S. E., & Delaney, H. D. (2004). \emph{Designing experiments and analyzing data: A model comparison perspective} (2nd ed.). Mahwah, NJ: Lawrence Erlbaum Associates, Inc.

Meade, A. W., Johnson, E. C., & Braddy, P. W. (2008). Power and sensitivity of alternative fit indices in tests of measurement invariance. \emph{Journal of Applied Psychology, 93}(3), 568-592. doi:10.1037/0021-9010.93.3.568

Widamin, K. F., & Thompson, J. S. (2003). On specifying the null model for incremental fit indices in structural equation modeling. \emph{Psychological Methods, 8}(1), 16-37. doi:10.1037/1082-989X.8.1.16
}
\seealso{
\code{\link[stats]{TukeyHSD}}, \code{\link[semTools]{measurementInvariance}}, \code{\link[semTools]{measurementInvarianceCat}}
}
\examples{
\dontrun{
## fit indices of interest for multivariate omnibus test of measurement equivalence
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
                moreFitIndices(fit.metric, fit.measures = moreAFIs))

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
## Results object contains info about multivariate omnibus and follow-up tests.
## It contains the observed DIF that can be calculated using calculateDIF:
out@observed@DIF
calculateDIF(uncon = fit.config, param = "loadings")

## The p values for each method's distribution can be inspected:
out@p.values

## The "show" method prints results for multivariate omnibus null hypothesis,
## as well as individual omnibus tests per item (controlling familywise errors)
out

## The "summary" method gives details about "follow-up" DIF tests.

## The default (type = "all") is a Tukey-type method appropriate if the
## multivariate omnibus null is reject AFIs
summary(out, digits = 2)
## Not much to see with only one significant DIF.

## Try using individual maximum-DIF distributions for each parameter, which is
## appropriate for a particular item when its univariate omnibus null is
## rejected using the maximum Sum-of-Squared DIF across items.
## Note that we can raise the alpha level and control the number of digits.
summary(out, type = "each", alpha = .10, digits = 2)

## To not control the Type I error rate at all (and use a liberal alpha):
outsum <- summary(out, type = "pairs", alpha = 0.30)
## notice that the returned object is a logical matrix: is.rejected?
outsum

## To control FDR (instead of familywise Type I errors) at level alpha:
summary(out, type = "step-up", digits = 2)
}
}

