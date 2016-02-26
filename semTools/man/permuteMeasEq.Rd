\name{permuteMeasEq}
\alias{permuteMeasEq}
\title{
	Permutation Randomization Tests of Measurement Equivalence and Differential Item Functioning (DIF)
}
\description{
The function \code{permuteMeasEq} accepts a pair of nested lavaan objects, the less constrained of which (\code{uncon}) freely estimates a set of measurement parameters (e.g., factor loadings, intercepts, or thresholds) in all groups, and the more constrained of which (\code{con}) constrains those measurement parameters to equality across groups. Group assignment is repeatedly permuted and the model is fit to each permutation, in order to produce an empirical distribution under the null hypothesis of no group differences, both for (a) changes in user-specified fit measures (see \code{AFIs} and \code{moreAFIs}) and for (b) the maximum modification index among the user-specified parameters (see \code{param}). This function is for testing measurement equivalence only across groups, not occasions. Configural invariance can also be tested by providing that fitted lavaan object to \code{con} and leaving \code{uncon = NULL}, in which case \code{param} must be \code{NULL} as well.

Modification indices for equality constraints on parameters specified in \code{param} are calculated from the constrained model (\code{con}) using the function \code{\link[lavaan]{lavTestScore}}, which can also be used to request expected parameter changes if the user has a need for them.
}
\usage{
permuteMeasEq(nPermute, con, uncon = NULL, null = NULL,
              param = NULL, AFIs = NULL, moreAFIs = NULL,
              maxSparse = 10, maxNonconv = 10, showProgress = TRUE)
}
\arguments{
  \item{nPermute}{
	An integer indicating the number of random permutations of group assignment used to form empirical distributions under the null hypothesis.
}
  \item{con}{
  The constrained \code{lavaan} object, in which the parameters specified in \code{param} are constrained to equality across all groups. In the case of testing \emph{configural} invariance, \code{con} is the configural model (implicitly, the unconstrained model is the saturated model, so use the default \code{uncon = NULL}).
}
  \item{uncon}{
	Optional.  The unconstrained \code{lavaan} object, in which the parameters specified in \code{param} are freely estimated in all groups. Only in the case of testing \emph{configural} invariance should this argument be \code{NULL}.
}
  \item{null}{
  Optional.  A \code{lavaan} object, in which an alternative null model is fit (besides the default independence model specified by \code{lavaan}) for the calculation of incremental fit indices. See Widamin & Thompson (2003) for details. If \code{NULL}, \code{lavaan}'s default independence model is used.
}
  \item{param}{
  A character vector indicating which parameters are constrained across groups in \code{con} and are unconstrained in \code{uncon}. Parameter names must match those returned by \code{names(coef(uncon))}, but omitting any group-specific suffixes (e.g., \code{"f1~1"} rather than \code{"f1~1.g2"}) or user-specified labels. Alternatively, to test an entire set of a certain type of parameter, \code{param} may take any one of the following values: \code{"loadings"}, \code{"intercepts"}, \code{"thresholds"}, or \code{"residuals"}.
}
  \item{AFIs}{
  A character vector indicating which alternative fit indices (or chi-squared itself) are to be used to test the multiparameter omnibus null hypothesis of no group differences in any parameters specified in \code{param}. Any fit measures returned by \code{\link[lavaan]{fitMeasures}} may be specified (including constants like \code{"df"}, which would be nonsensical). If both \code{AFIs} and \code{moreAFIs} are \code{NULL}, only \code{"chisq"} will be returned.

}
  \item{moreAFIs}{
  Optional. A character vector indicating which (if any) alternative fit indices returned by \code{\link[semTools]{moreFitIndices}} are to be used to test the multiparameter omnibus null hypothesis of no group differences in any parameters specified in \code{param}.
}
  \item{maxSparse}{
  An integer indicating the maximum number of consecutive times that randomly permuted group assignment can yield a sample in which at least one category (of an \code{ordered} indicator) is unobserved in at least one group, such that the same set of parameters cannot be estimated in each group. If such a sample occurs, group assignment is randomly permuted again, repeatedly until a sample is obtained with all categories observed in all groups. If \code{maxSparse} is exceeded, \code{NA} will be returned for that iteration of the permutation distribution.
}
  \item{maxNonconv}{
  An integer indicating the maximum number of consecutive times that randomly permuted group assignment can yield a sample for which the model does not converge on a solution. If such a sample occurs, group assignment is randomly permuted again, repeatedly until a sample is obtained for which the model does converge. If \code{maxNonconv} is exceeded, \code{NA} will be returned for that iteration of the permutation distribution, and a warning will be printed when using \code{show} or \code{summary}.
}
  \item{showProgress}{
	Logical. Indicating whether to display a progress bar while permuting.
}
}
\details{
	The multiparameter omnibus null hypothesis of measurement equivalence/invariance is that there are no group differences in any measurement parameters. This can be tested using the \code{anova} method on nested \code{lavaan} objects, as seen in the output of \code{\link[semTools]{measurementInvariance}}, or by inspecting the change in alternative fit indices (AFIs) such as the CFI. See Cheung & Rensvold (2002) or Meade, Johnson, & Braddy (2008) for details.

  If the multiparameter omnibus null hypothesis is rejected, partial invariance can still be established by freeing parameters that differ across groups, while maintaining equality constraints for at least two other indicators per factor. Modification indices can be calculated from the constrained model (\code{con}), but multiple testing leads to inflation of Type I error rates. The permutation randomization method employed by \code{permuteMeasEq} creates a distribution of the maximum modification index if the null hypothesis is true, which allows the user to control the familywise Type I error rate in a manner similar to Tukey's \emph{q} (studentized range) distribution for the Honestly Significant Difference (HSD) post hoc test.
}
\value{
	The \linkS4class{permuteMeasEq} object representing the results of testing measurement equivalence (the multiparameter omnibus test) and DIF (modification indices).
}
\author{
    Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\references{
Cheung, G. W., & Rensvold, R. B. (2002). Evaluating goodness-of-fit indexes for testing measurement invariance. \emph{Structural Equation Modeling, 9}(2), 233-255. doi:10.1207/S15328007SEM0902_5

Meade, A. W., Johnson, E. C., & Braddy, P. W. (2008). Power and sensitivity of alternative fit indices in tests of measurement invariance. \emph{Journal of Applied Psychology, 93}(3), 568-592. doi:10.1037/0021-9010.93.3.568

Widamin, K. F., & Thompson, J. S. (2003). On specifying the null model for incremental fit indices in structural equation modeling. \emph{Psychological Methods, 8}(1), 16-37. doi:10.1037/1082-989X.8.1.16
}
\seealso{
\code{\link[stats]{TukeyHSD}}, \code{\link[semTools]{measurementInvariance}}, \code{\link[semTools]{measurementInvarianceCat}}
}
\examples{
\dontrun{

########################
## Traditional Method ##
########################

## create 3-group data in lavaan example(cfa) data
HS <- lavaan::HolzingerSwineford1939
HS$ageGroup <- ifelse(HS$ageyr < 13, "preteen",
                      ifelse(HS$ageyr > 13, "teen", "thirteen"))

## specify and fit an appropriate null model for incremental fit indices
mod.null <- c(paste0("x", 1:9, " ~ c(T", 1:9, ", T", 1:9, ", T", 1:9, ")*1"),
              paste0("x", 1:9, " ~~ c(L", 1:9, ", L", 1:9, ", L", 1:9, ")*x", 1:9))
fit.null <- cfa(mod.null, data = HS, group = "ageGroup")

## fit target model with varying levels of measurement equivalence
mod.config <- '
visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9
'
miout <- measurementInvariance(mod.config, data = HS, std.lv = TRUE,
                               group = "ageGroup")

(fit.config <- miout[["fit.configural"]])
(fit.metric <- miout[["fit.loadings"]])
(fit.scalar <- miout[["fit.intercepts"]])

########################
## Permutation Method ##
########################

## fit indices of interest for multiparameter omnibus test
myAFIs <- c("chisq","cfi","rmsea","srmr","mfi","aic")
moreAFIs <- c("gammaHat","adjGammaHat")

## Use only 20 permutations for a demo.  In practice,
## use > 1000 to reduce sampling variability of estimated p values

## test configural invariance
set.seed(123)
out.config <- permuteMeasEq(nPermute = 20, con = fit.config)
out.config

## test metric equivalence
set.seed(456)
out.metric <- permuteMeasEq(nPermute = 20, uncon = fit.config, con = fit.metric,
                            param = "loadings", AFIs = myAFIs,
                            moreAFIs = moreAFIs, null = fit.null)
summary(out.metric, nd = 4)

## test scalar equivalence
set.seed(789)
out.scalar <- permuteMeasEq(nPermute = 20, uncon = fit.metric, con = fit.scalar,
                            param = "intercepts", AFIs = myAFIs,
                            moreAFIs = moreAFIs, null = fit.null)
summary(out.scalar)

## Not much to see without significant DIF.
## Try using an absurdly high alpha level for illustration.
outsum <- summary(out.scalar, alpha = .50)

## notice that the returned object is the table of DIF tests
outsum

## visualize permutation distribution
hist(out.config, AFI = "chisq")
hist(out.metric, AFI = "chisq", nd = 2, alpha = .01,
     legendArgs = list(x = "topright"))
hist(out.scalar, AFI = "cfi", printLegend = FALSE)

}
}

