\name{twostage}
\alias{twostage}
\alias{cfa.2stage}
\alias{sem.2stage}
\alias{growth.2stage}
\alias{lavaan.2stage}
\title{
	Fit a lavaan model using 2-Stage Maximum Likelihood (TSML) estimation for missing data.
}
\description{
	This function automates 2-Stage Maximum Likelihood (TSML) estimation, optionally with auxiliary variables.  Step 1 involves fitting a saturated model to the partially observed data set (to variables in the hypothesized model as well as auxiliary variables related to missingness).  Step 2 involves fitting the hypothesized model to the model-implied means and covariance matrix (also called the "EM" means and covariance matrix) as if they were complete data.  Step 3 involves correcting the Step-2 standard errors (\emph{SE}s) and chi-squared statistic to account for additional uncertainty due to missing data (using information from Step 1; see References section for sources with formulas).

	All variables (including auxiliary variables) are treated as endogenous varaibles in the Step-1 saturated model (\code{fixed.x = FALSE}), so data are assumed continuous, although not necessarily multivariate normal (dummy-coded auxiliary variables may be included in Step 1, but categorical endogenous variables in the Step-2 hypothesized model are not allowed).  To avoid assuming multivariate normality, request \code{se = "robust.huber.white"}.  CAUTION: In addition to setting \code{fixed.x = FALSE} and \code{conditional.x = FALSE} in \code{\link[lavaan]{lavaan}}, this function will automatically set \code{meanstructure = TRUE}, \code{estimator = "ML"}, \code{missing = "fiml"}, and \code{test = "standard"}.  \code{\link[lavaan]{lavaan}}'s \code{se} option can only be set to \code{"standard"} to assume multivariate normality or to \code{"robust.huber.white"} to relax that assumption.
}
\usage{
twostage(..., aux, fun, baseline.model = NULL)
cfa.2stage(..., aux = NULL, baseline.model = NULL)
sem.2stage(..., aux = NULL, baseline.model = NULL)
growth.2stage(..., aux = NULL, baseline.model = NULL)
lavaan.2stage(..., aux = NULL, baseline.model = NULL)
}
\arguments{
  \item{\dots}{
	Arguments passed to the \code{\link[lavaan]{lavaan}} function specified in the \code{fun} argument.  At a minimum, the user must supply the first two named arguments to \code{\link[lavaan]{lavaan}} (i.e., \code{model} and \code{data}).
  }
  \item{aux}{
	An optional character vector naming auxiliary variable(s) in \code{data}
  }
  \item{fun}{
	The character string naming the lavaan function used to fit the Step-2 hypothesized model (\code{"cfa"}, \code{"sem"}, \code{"growth"}, or \code{"lavaan"}).
  }
  \item{baseline.model}{
	An optional character string, specifying the lavaan \code{\link[lavaan]{model.syntax}} for a user-specified baseline model.  Interested users can use the fitted baseline model to calculate incremental fit indices (e.g., CFI and TLI) using the corrected chi-squared values (see the \code{anova} method in \code{\linkS4class{twostage}}).  If \code{NULL}, the default "independence model" (i.e., freely estimated means and variances, but all covariances constrained to zero) will be specified internally.
  }
}
\value{
	The \code{\linkS4class{twostage}} object contains 3 fitted lavaan models (saturated, target/hypothesized, and baseline) as well as the names of auxiliary variables.  None of the individual models provide the correct model results (except the point estimates in the target model are unbiased).  Use the methods in \code{\linkS4class{twostage}} to extract corrected \emph{SE}s and test statistics.
}
\references{
Savalei, V., \& Bentler, P. M. (2009). A two-stage approach to missing data: Theory and application to auxiliary variables. \emph{Structural Equation Modeling, 16}(3), 477-497. doi:10.1080/10705510903008238

Savalei, V., \& Falk, C. F. (2014). Robust two-stage approach outperforms robust full information maximum likelihood with incomplete nonnormal data. \emph{Structural Equation Modeling, 21}(2), 280-302. doi:10.1080/10705511.2014.882692

}
\seealso{
	\code{\linkS4class{twostage}}
}
\author{
    Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\examples{

## set some example data missing at random
dat1 <- HolzingerSwineford1939
dat1$x5 <- ifelse(dat1$x1 <= quantile(dat1$x1, .3), NA, dat1$x5)
dat1$age <- dat1$ageyr + dat1$agemo/12
dat1$x9 <- ifelse(dat1$age <= quantile(dat1$age, .3), NA, dat1$x9)

## fit CFA model from lavaan's ?cfa help page
model <- '
visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9
'
## use ageyr and agemo as auxiliary variables
out <- cfa.2stage(model = model, data = dat1, aux = c("ageyr","agemo"))

## two versions of a corrected chi-squared test results are shown
out
## see Savalei & Bentler (2009) and Savalei & Falk (2014) for details

## the summary additionally provides the parameter estimates with corrected
## standard errors, test statistics, and confidence intervals, along with
## any other options that can be passed to parameterEstimates()
summary(out, standardized = TRUE)



## use parameter labels to fit a more constrained model
modc <- '
visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + a*x8 + a*x9
'
outc <- cfa.2stage(model = modc, data = dat1, aux = c("ageyr","agemo"))


## use the anova() method to test this constraint
anova(out, outc)
## like for a single model, two corrected statistics are provided

}
