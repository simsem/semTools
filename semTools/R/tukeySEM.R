### Alexander M. Schoemann
### Last updated: 9 March 2018


#' Tukey's WSD post-hoc test of means for unequal variance and sample size
#'
#' This function computes Tukey's WSD post hoc test of means when variances and
#' sample sizes are not equal across groups. It can be used as a post hoc test
#' when comparing latent means in multiple group SEM.
#'
#' After conducting an omnibus test of means across three of more groups,
#' researchers often wish to know which sets of means differ at a particular
#' Type I error rate. Tukey's WSD test holds the error rate stable across
#' multiple comparisons of means. This function implements an adaptation of
#' Tukey's WSD test from Maxwell & Delaney (2004), that allows variances and
#' sample sizes to differ across groups.
#'
#'
#' @importFrom stats ptukey
#'
#' @param m1 Mean of group 1.
#' @param m2 Mean of group 2.
#' @param var1 Variance of group 1.
#' @param var2 Variance of group 2.
#' @param n1 Sample size of group 1.
#' @param n2 Sample size of group 2.
#' @param ng Total number of groups to be compared (i.e., the number of groups
#' compared in the omnibus test).
#' @return A vector with three elements:
#' \enumerate{
#'  \item \code{q}: The \emph{q} statistic
#'  \item \code{df}: The degrees of freedom for the \emph{q} statistic
#'  \item \code{p}: A \emph{p} value based on the \emph{q} statistic, \emph{df},
#'   and the total number of groups to be compared
#' }
#' @author Alexander M. Schoemann (East Carolina University;
#' \email{schoemanna@@ecu.edu})
#' @references Maxwell, S. E., & Delaney, H. D. (2004). \emph{Designing
#' experiments and analyzing data: A model comparison perspective} (2nd ed.).
#' Mahwah, NJ: Lawrence Erlbaum Associates.
#' @examples
#'
#' ## For a case where three groups have been compared:
#' ## Group 1: mean = 3.91, var = 0.46, n = 246
#' ## Group 2: mean = 3.96, var = 0.62, n = 465
#' ## Group 3: mean = 2.94, var = 1.07, n = 64
#'
#' ## compare group 1 and group 2
#' tukeySEM(3.91, 3.96, 0.46, 0.62, 246, 425, 3)
#'
#' ## compare group 1 and group 3
#' tukeySEM(3.91, 2.94, 0.46, 1.07, 246, 64, 3)
#'
#' ## compare group 2 and group 3
#' tukeySEM(3.96, 2.94, 0.62, 1.07, 465, 64, 3)
#'
#' @export
tukeySEM <- function(m1, m2, var1, var2, n1, n2, ng) {
  qNum <- abs(m1 - m2)
  qDenom <- sqrt(((var1/n1) + (var2/n2))/2)
  Tukeyq <- qNum / qDenom
  Tukeydf <- ((var1/n1) + (var2/n2))^2 /
    (((var1/n1)^2 / (n1 - 1)) + ((var2/n2)^2 / (n2 - 2)))
  c(q = Tukeyq, df = Tukeydf, p = 1 - ptukey(Tukeyq, ng, Tukeydf))
}

##Example from Schoemann (2013)
##Bio vs. policial science on evo misconceptions
#tukeySEM(3.91, 3.96,.46, .62, 246, 425,3)

