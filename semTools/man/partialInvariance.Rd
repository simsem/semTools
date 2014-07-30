\name{partialInvariance}
\alias{partialInvariance}
\title{
Partial Measurement Invariance Testing Across Groups
}
\description{
This test will provide partial invariance testing by (a) freeing a parameter one-by-one from nested model and compare with the original nested model or (b) fixing (or constraining) a parameter one-by-one from the parent model and compare with the original parent model. This function only works with continuous variables and congeneric models.
}
\usage{
partialInvariance(fit, type, free = NULL, fix = NULL, return.fit = FALSE)
}
\arguments{
  \item{fit}{A list of models for invariance testing. Each model should be assigned by appropriate names (see details). The result from \code{\link{measurementInvariance}} could be used in this argument directly.}
  \item{type}{The types of invariance testing: "metric", "scalar", "strict", or "means"}
  \item{free}{A vector of variable names that are free across groups in advance. If partial mean invariance is tested, this argument represents a vector of factor names that are free across groups.}
  \item{fix}{A vector of variable names that are constrained to be equal across groups in advance. If partial mean invariance is tested, this argument represents a vector of factor names that are fixed across groups.}
  \item{return.fit}{Return the submodels fitted by this function}
}
\details{
There are four types of partial invariance testing:

\itemize{
  \item{Partial weak invariance. The model named 'configural' from the list of models is compared with the model named 'metric'. Each loading will be freed or fixed from the metric and configural invariance models respectively. The modified models are compared with the original model. Note that the objects in the list of models must have the names of "configural" and "metric". Users may use "metric", "weak", "loading", or "loadings" in the \code{type} argument. Note that, for testing invariance on marker variables, other variables will be assigned as marker variables automatically.}
  \item{Partial strong invariance. The model named 'metric' from the list of models is compared with the model named 'scalar'. Each intercept will be freed or fixed from the scalar and metric invariance models respectively. The modified models are compared with the original model. Note that the objects in the list of models must have the names of "metric" and "scalar". Users may use "scalar", "strong", "intercept", or "intercepts" in the \code{type} argument. Note that, for testing invariance on marker variables, other variables will be assigned as marker variables automatically.}
  \item{Partial strict invariance. The model named 'scalar' from the list of models is compared with the model named 'strict'. Each residual variance will be freed or fixed from the strict and scalar invariance models respectively. The modified models are compared with the original model. Note that the objects in the list of models must have the names of "scalar" and "strict". Users may use "strict", "residual", "residuals", "error", or "errors" in the \code{type} argument.}
  \item{Partial mean invariance. The model named 'scalar' (or 'strict') from the list of models is compared with the model named 'means'. Each factor mean will be freed or fixed from the means and scalar (or strict) invariance models respectively. The modified models are compared with the original model. Note that the objects in the list of models must have the names of "means" and either "strict" or "scalar". Users may use "means" or "mean" in the \code{type} argument.}
}

Two types of comparisons are used in this function:
\enumerate{
  \item{\code{free}: The nested model is used as a template. Then, one parameter indicating the differences between two models is free. The new model is compared with the nested model. This process is repeated for all differences between two models.}
  \item{\code{fix}: The parent model is used as a template. Then, one parameter indicating the differences between two models is fixed or constrained to be equal to other parameters. The new model is then compared with the parent model. This process is repeated for all differences between two models.}
}

Note that this function does not adjust for the inflated Type I error rate from multiple tests. The degree of freedom of all tests would be the number of groups minus 1.
}
\value{
	If \code{return.fit = FALSE}, the result tables are provided. Chi-square and p-value are provided for all methods. 
	
	If \code{return.fit = TRUE}, a list with two elements are provided. The first element is the tabular result. The second element is the submodels used in the \code{free} and \code{fix} methods. 
}
\author{
	Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
}
\seealso{
	\code{\link{measurementInvariance}} for measurement invariance for continuous variables
}
\examples{
# Conduct weak invariance testing manually by using fixed-factor
# method of scale identification

conf <- "
f1 =~ NA*x1 + x2 + x3
f2 =~ NA*x4 + x5 + x6
f1 ~~ c(1, 1)*f1
f2 ~~ c(1, 1)*f2
"

weak <- "
f1 =~ NA*x1 + x2 + x3
f2 =~ NA*x4 + x5 + x6
f1 ~~ c(1, NA)*f1
f2 ~~ c(1, NA)*f2
"

configural <- cfa(conf, data = HolzingerSwineford1939, std.lv = TRUE, group="school")
weak <- cfa(weak, data = HolzingerSwineford1939, group="school", group.equal="loadings")
models <- list(configural = configural, metric = weak)
partialInvariance(models, "metric")

\dontrun{
partialInvariance(models, "metric", free = "x5") # "x5" is free across groups in advance
partialInvariance(models, "metric", fix = "x4") # "x4" is fixed across groups in advance

# Use the result from the measurementInvariance function
HW.model <- ' visual =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed =~ x7 + x8 + x9 '

models2 <- measurementInvariance(HW.model, data=HolzingerSwineford1939, group="school")
partialInvariance(models2, "scalar")
}
}
