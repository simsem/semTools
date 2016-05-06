\name{permuteMeasEq-class}
\docType{class}
\alias{permuteMeasEq-class}
\alias{show,permuteMeasEq-method}
\alias{summary,permuteMeasEq-method}
\alias{hist,permuteMeasEq-method}
\title{
  Class for the Results of Permutation Randomization Tests of Measurement Equivalence and DIF
}
\description{
	This class contains the results of tests of Measurement Equivalence and Differential Item Functioning (DIF).
}
\section{Objects from the Class}{
	Objects can be created via the \code{\link[semTools]{permuteMeasEq}} function.
}
\section{Slots}{
  \describe{
    \item{\code{PT}:}{A \code{data.frame} returned by a call to \code{\link[lavaan]{parTable}} on the constrained model}
    \item{\code{modelType}:}{A character indicating the specified \code{modelType} in the call to \code{permuteMeasEq}}
    \item{\code{ANOVA}:}{A vector indicating the results of the observed chi-squared (difference) test, based on the central chi-squared distribution}
    \item{\code{AFI.obs}:}{A vector of observed (changes in) user-selected fit measures}
    \item{\code{AFI.dist}:}{The permutation distribution(s) of user-selected fit measures. A \code{data.frame} with \code{n.Permutations} rows and one column for each \code{AFI.obs}.}
    \item{\code{AFI.pval}:}{A vector of \emph{p} values (one for each element in slot \code{AFI.obs}) calculated using slot \code{AFI.dist}, indicating the probability of observing a change at least as extreme as \code{AFI.obs} if the null hypothesis were true}
    \item{\code{MI.obs}:}{A \code{data.frame} of observed Lagrange Multipliers (modification indices) associated with the equality constraints or fixed parameters specified in the \code{param} argument. This is a subset of the output returned by a call to \code{\link[lavaan]{lavTestScore}} on the constrained model.}
    \item{\code{MI.dist}:}{The permutation distribution of the maximum modification index (among those seen in slot \code{MI.obs$X2}) at each permutation of group assignment or of \code{covariates}}
    \item{\code{extra.obs}:}{If \code{permuteMeasEq} was called with an \code{extra} function, the output when applied to the original data is concatenated into this vector}
    \item{\code{extra.dist}:}{A \code{data.frame}, each column of which contains the permutation distribution of the corresponding statistic in slot \code{extra.obs}}
    \item{\code{n.Permutations}:}{An integer indicating the number of permutations requested by the user}
    \item{\code{n.Converged}:}{An integer indicating the number of permuation iterations which yielded a converged solution}
    \item{\code{n.nonConverged}:}{A vector of length \code{n.Permutations} indicating how many times group assignment was randomly permuted (at each iteration) before converging on a solution}
    \item{\code{n.Sparse}:}{Only relevant with \code{ordered} indicators when \code{modelType == "mgcfa"}. A vector of length \code{n.Permutations} indicating how many times group assignment was randomly permuted (at each iteration) before obtaining a sample with all categories observed in all groups}
    \item{\code{oldSeed}:}{An integer vector storing the value of \code{.Random.seed} before running \code{permuteMeasEq}. Only relevant when using a parallel/multicore option and the original \code{RNGkind() != "L'Ecuyer-CMRG"}. This enables users to restore their previous \code{.Random.seed} state, if desired, by running:  \code{.Random.seed[-1] <- permutedResults@oldSeed[-1]}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "permuteMeasEq"):} The \code{show} function is used to summarize the results of the multiparameter omnibus test of measurement equivalence, using the user-specified AFIs. The parametric chi-squared (difference) test is also displayed.}
    \item{summary}{\code{signature(object = "permuteMeasEq", alpha = .05, nd = 3, extra = FALSE):} The summary function prints the same information from the \code{show} method, but when \code{extra = FALSE} (the default) it also provides a table summarizing any requested follow-up tests of DIF using modification indices in slot \code{MI.obs}. The user can also specify an \code{alpha} level for flagging modification indices as significant, as well as \code{nd} (the number of digits displayed). For each modification index, the \emph{p} value is displayed using a central chi-squared distribution with the \emph{df} shown in that column. Additionally, a \emph{p} value is displayed using the permutation distribution of the maximum index, which controls the familywise Type I error rate in a manner similar to Tukey's studentized range test. If any indices are flagged as significant using the \code{tukey.p.value}, then a message is displayed for each flagged index. The invisibly returned \code{data.frame} is the displayed table of modification indices, unless \code{\link[semTools]{permuteMeasEq}} was called with \code{param = NULL}, in which case the invisibly returned object is \code{object}. If \code{extra = TRUE}, the permutation-based \emph{p} values for each statistic returned by the \code{extra} function are displayed and returned in a \code{data.frame} instead of the modification indices requested in the \code{param} argument.}
    \item{hist}{\code{signature(x = "permuteMeasEq", ..., AFI, alpha = .05, nd = 3, printLegend = TRUE, legendArgs = list(x = "topleft")):} The \code{hist} function provides a histogram for the permutation distribution of the specified \code{AFI}, including observed and critical values from the specified \code{alpha} level. Distributions of modification indices and any extra output are not available with this method, but they can be created manually by accessing the distributions in slot \code{MI.dist} or \code{extra.dist}. The user can also specify additional graphical parameters to \code{\link[graphics]{hist}} via \code{...}, as well as pass a list of arguments to an optional \code{\link[graphics]{legend}} via \code{legendArgs}. If \code{AFI = "chisq"}, then the probability density and critical value from the central chi-squared distribution are also included in the plot. If the user wants more control over customization, \code{hist} returns a list of \code{length == 2}, containing the arguments for the call to \code{hist} and the arguments to the call for \code{legend}, respectively. This list may facilitate creating a customized histogram of \code{AFI.dist}, \code{MI.dist}, or \code{extra.dist}.}
  }
}
\author{
    Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\seealso{
\code{\link[semTools]{permuteMeasEq}}
}
\examples{
# See the example from the permuteMeasEq function
}

