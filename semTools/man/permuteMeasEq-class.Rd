\name{permuteMeasEq-class}
\docType{class}
\alias{permuteMeasEq-class}
\alias{show,permuteMeasEq-method}
\alias{summary,permuteMeasEq-method}
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
    \item{\code{observed}:}{An object of S4 class \code{MeasEq.observed}, which includes a vector of observed changes in user-selected alternative fit indices (AFIs) and a matrix of observed differences in user-specified measurement parameters (i.e., DIF estimates), which is returned by \code{\link[semTools]{calculateDIF}}.}
    \item{\code{p.values}:}{An object of S4 class \code{MeasEq.p.values}, which includes a vector of \emph{p} values associated with each observed change in user-selected AFIs under the null; a vectorof \emph{p} values associated with the maximum Sum-of-Squared-DIF across items; and three matrices of \emph{p} values associated with each observed DIF estimate under the null: one based on a distribution of the maximum DIF across all parameters and pairwise comparisons, and one based on separate distributions of the maximum DIF across all pairwise comparisons within each individual parameter, and one based on each pairwise comparison within each item (i.e., without controlling Type I errors).}
    \item{\code{ANOVA}:}{A vector indicating the results of the observed chi-squared difference test, based on the theoretical chi-squared distribution}
    \item{\code{n.Permutations}:}{An integer indicating the number of permutations requested by the user}
    \item{\code{n.Converged}:}{An integer indicating the number of permuation iterations which yielded a converged solution}
    \item{\code{n.nonConverged}:}{A vector of length \code{n.Permutations} indicating how many times group assignment was randomly permuted (at each iteration) before converging on a solution}
    \item{\code{n.Sparse}:}{Only relevant with \code{ordered} indicators. A vector of length \code{n.Permutations} indicating  how many times group assignment was randomly permuted (at each iteration) before obtaining a sample with all categories observed in all groups}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "permuteMeasEq"):} The show function is used to summarize the results of the multivariate omnibus test of measurement equivalence, using the user-specified AFIs.}
	  \item{summary}{\code{signature(object = "permuteMeasEq", alpha = .05, type = c("all","each","pairs","step-up"), digits = 3):} The summary function is used to summarize the results of the univariate omnibus tests (using the maximum Sum-of-Squared-DIF across items), as well as results of DIF tests for each user-specified parameter, using a specified \code{alpha} level and \code{type} of reference distribution for the maximum DIF under the null. The user can also specify the number of \code{digits} reported for DIF estimates. Results of the linear step-up procedure can also be summarized, which control the expected false discovery rate (FDR) to yield greater power than controlling the familywise Type I error rate.}
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

