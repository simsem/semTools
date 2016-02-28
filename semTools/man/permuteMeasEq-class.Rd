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
    \item{\code{ANOVA}:}{A vector indicating the results of the observed chi-squared (difference) test, based on the theoretical chi-squared distribution}
    \item{\code{AFI.obs}:}{A vector of observed changes in user-selected alternative fit indices (AFIs)}
    \item{\code{AFI.dist}:}{The permutation distribution(s) of AFI(s). A \code{data.frame} with \code{nrow == n.Permutations} and one column for each AFI.}
    \item{\code{AFI.pval}:}{A vector of \emph{p} values (one for each AFI in slot \code{AFI.obs}) calculated using slot \code{AFI.dist}, indicating the probability of observing a change at least as extreme as \code{AFI.obs} if the null hypothesis of no group differences were true}
    \item{\code{MI.obs}:}{A vector of observed Lagrange Multipliers (modification indices) associated with the equality constraints on the user-specified parameters. This is a subset of the output returned by a call to \code{\link[lavaan]{lavTestScore}} on the constrained model.}
    \item{\code{MI.dist}:}{The permutation distribution of the maximum modification index (among those seen in slot \code{MI.obs}) at each permutation of group assignment}
    \item{\code{n.Permutations}:}{An integer indicating the number of permutations requested by the user}
    \item{\code{n.Converged}:}{An integer indicating the number of permuation iterations which yielded a converged solution}
    \item{\code{n.nonConverged}:}{A vector of length \code{n.Permutations} indicating how many times group assignment was randomly permuted (at each iteration) before converging on a solution}
    \item{\code{n.Sparse}:}{Only relevant with \code{ordered} indicators. A vector of length \code{n.Permutations} indicating  how many times group assignment was randomly permuted (at each iteration) before obtaining a sample with all categories observed in all groups}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "permuteMeasEq"):} The \code{show} function is used to summarize the results of the multiparameter omnibus test of measurement equivalence, using the user-specified AFIs. The parametric chi-squared (difference) test is also displayed.}
    \item{summary}{\code{signature(object = "permuteMeasEq", alpha = .05, nd = 3):} The summary function prints the same information from the \code{show} method, but also provides a table summarizing follow-up tests of DIF using modification indices in slot \code{MI.obs}. The user can also specify an \code{alpha} level for flagging modification indices as significant, as well as \code{nd} (the number of digits displayed). For each modification index, the \emph{p} value is displayed using the parametric chi-squared distribution with \emph{df} = 1. Additionally, a \emph{p} value is displayed using the permutation distribution of the maximum index, which controls the familywise Type I error rate in a manner similar to Tukey's studentized range test. If any indices are flagged as significant (using the \code{tukey.p.value}), then a message is displayed for each flagged index. The invisibly returned \code{data.frame} is the displayed table of modification indices, unless \code{\link[semTools]{permuteMeasEq}} was called with \code{uncon = NULL} to test configural invariance, in which case the invisibly returned object is \code{object}.}
    \item{hist}{\code{signature(x = "permuteMeasEq", ..., AFI, alpha = .05, nd = 3, printLegend = TRUE, legendArgs = list(x = "topleft")):} The \code{hist} function provides a histogram for the permutation distribution of the specified \code{AFI}, including observed and critical values from the specified \code{alpha} level. The user can also specify additional graphical parameters to \code{\link[graphics]{hist}} via \code{...}, as well as pass a list of arguments to an optional \code{\link[graphics]{legend}} via \code{legendArgs}. If \code{AFI = "chisq"}, then the probability density and critical value from the theoretical chi-squared distribution are also included in the plot. If the user wants more control over customization, \code{hist} returns a list of \code{length == 2}, containing the arguments for the call to \code{hist} and the arguments to the call for \code{legend}, respectively.}
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

