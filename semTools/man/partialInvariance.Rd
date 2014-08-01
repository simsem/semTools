\name{partialInvariance}
\alias{partialInvariance}
\alias{partialInvarianceCat}
\title{
Partial Measurement Invariance Testing Across Groups
}
\description{
This test will provide partial invariance testing by (a) freeing a parameter one-by-one from nested model and compare with the original nested model or (b) fixing (or constraining) a parameter one-by-one from the parent model and compare with the original parent model. This function only works with congeneric models. The \code{partialInvariance} is used for continuous variable. The \code{partialInvarianceCat} is used for categorical variables.
}
\usage{
partialInvariance(fit, type, free = NULL, fix = NULL, return.fit = FALSE)
partialInvarianceCat(fit, type, free = NULL, fix = NULL, return.fit = FALSE)
}
\arguments{
  \item{fit}{A list of models for invariance testing. Each model should be assigned by appropriate names (see details). The result from \code{\link{measurementInvariance}} or \code{\link{measurementInvarianceCat}} could be used in this argument directly.}
  \item{type}{The types of invariance testing: "metric", "scalar", "strict", or "means"}
  \item{free}{A vector of variable names that are free across groups in advance. If partial mean invariance is tested, this argument represents a vector of factor names that are free across groups.}
  \item{fix}{A vector of variable names that are constrained to be equal across groups in advance. If partial mean invariance is tested, this argument represents a vector of factor names that are fixed across groups.}
  \item{return.fit}{Return the submodels fitted by this function}
}
\details{
There are four types of partial invariance testing:

\itemize{
  \item{Partial weak invariance. The model named 'configural' from the list of models is compared with the model named 'metric'. Each loading will be freed or fixed from the metric and configural invariance models respectively. The modified models are compared with the original model. Note that the objects in the list of models must have the names of "configural" and "metric". Users may use "metric", "weak", "loading", or "loadings" in the \code{type} argument. Note that, for testing invariance on marker variables, other variables will be assigned as marker variables automatically.}
  \item{Partial strong invariance. The model named 'metric' from the list of models is compared with the model named 'scalar'. Each intercept will be freed or fixed from the scalar and metric invariance models respectively. The modified models are compared with the original model. Note that the objects in the list of models must have the names of "metric" and "scalar". Users may use "scalar", "strong", "intercept", "intercepts", "threshold", or "thresholds" in the \code{type} argument. Note that, for testing invariance on marker variables, other variables will be assigned as marker variables automatically. Note that if all variables are dichotomous, scalar invariance testing is not available.}
  \item{Partial strict invariance. The model named 'scalar' (or 'metric') from the list of models is compared with the model named 'strict'. Each residual variance will be freed or fixed from the strict and scalar (or metric) invariance models respectively. The modified models are compared with the original model. Note that the objects in the list of models must have the names of "strict" and either "scalar" or "metric". Users may use "strict", "residual", "residuals", "error", or "errors" in the \code{type} argument.}
  \item{Partial mean invariance. The model named 'scalar' (or 'strict' or 'metric') from the list of models is compared with the model named 'means'. Each factor mean will be freed or fixed from the means and scalar (or strict or metric) invariance models respectively. The modified models are compared with the original model. Note that the objects in the list of models must have the names of "means" and either "strict", "scalar", or "metric". Users may use "means" or "mean" in the \code{type} argument.}
}

Two types of comparisons are used in this function:
\enumerate{
  \item{\code{free}: The nested model is used as a template. Then, one parameter indicating the differences between two models is free. The new model is compared with the nested model. This process is repeated for all differences between two models.}
  \item{\code{fix}: The parent model is used as a template. Then, one parameter indicating the differences between two models is fixed or constrained to be equal to other parameters. The new model is then compared with the parent model. This process is repeated for all differences between two models.}
  \item{\code{wald}: This method is similar to the \code{fix} method. However, instead of building a new model and compare them with likelihood-ratio test, multivariate wald test is used to compare equality between parameter estimates. See \code{\link{wald}} for further details. Note that if any rows of the contrast cannot be summed to 0, the Wald test is not provided, such as comparing two means where one of the means is fixed as 0. This test statistic is not as accurate as likelihood-ratio test provided in \code{fix}. I provide it here in case that likelihood-ratio test fails to converge.}
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
	\code{\link{measurementInvariance}} for measurement invariance for continuous variables; \code{\link{measurementInvarianceCat}} for measurement invariance for categorical variables; \code{\link{wald}} for multivariate Wald test
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

# Conduct weak invariance testing manually by using fixed-factor
# method of scale identification for dichotomous variables

f <- rnorm(1000, 0, 1)
u1 <- 0.9*f + rnorm(1000, 1, sqrt(0.19))
u2 <- 0.8*f + rnorm(1000, 1, sqrt(0.36))
u3 <- 0.6*f + rnorm(1000, 1, sqrt(0.64))
u4 <- 0.7*f + rnorm(1000, 1, sqrt(0.51))
u1 <- as.numeric(cut(u1, breaks = c(-Inf, 0, Inf)))
u2 <- as.numeric(cut(u2, breaks = c(-Inf, 0.5, Inf)))
u3 <- as.numeric(cut(u3, breaks = c(-Inf, 0, Inf)))
u4 <- as.numeric(cut(u4, breaks = c(-Inf, -0.5, Inf)))
g <- rep(c(1, 2), 500)
dat2 <- data.frame(u1, u2, u3, u4, g)

configural2 <- "
f1 =~ NA*u1 + u2 + u3 + u4
u1 | c(t11, t11)*t1 
u2 | c(t21, t21)*t1 
u3 | c(t31, t31)*t1 
u4 | c(t41, t41)*t1 
f1 ~~ c(1, 1)*f1
f1 ~ c(0, NA)*1
u1 ~~ c(1, 1)*u1
u2 ~~ c(1, NA)*u2
u3 ~~ c(1, NA)*u3
u4 ~~ c(1, NA)*u4
"

outConfigural2 <- cfa(configural2, data = dat2, group = "g", parameterization="theta", 
	estimator="wlsmv", ordered = c("u1", "u2", "u3", "u4"))

weak2 <- "
f1 =~ NA*u1 + c(f11, f11)*u1 + c(f21, f21)*u2 + c(f31, f31)*u3 + c(f41, f41)*u4
u1 | c(t11, t11)*t1 
u2 | c(t21, t21)*t1 
u3 | c(t31, t31)*t1 
u4 | c(t41, t41)*t1 
f1 ~~ c(1, NA)*f1
f1 ~ c(0, NA)*1
u1 ~~ c(1, 1)*u1
u2 ~~ c(1, NA)*u2
u3 ~~ c(1, NA)*u3
u4 ~~ c(1, NA)*u4
"

outWeak2 <- cfa(weak2, data = dat2, group = "g", parameterization="theta", estimator="wlsmv", 
	ordered = c("u1", "u2", "u3", "u4"))
modelsCat <- list(configural = outConfigural2, metric = outWeak2)

partialInvarianceCat(modelsCat, type = "metric") 

\dontrun{
partialInvarianceCat(modelsCat, type = "metric", free = "u2") 
partialInvarianceCat(modelsCat, type = "metric", fix = "u3") 

# Use the result from the measurementInvarianceCat function

model <- ' f1 =~ u1 + u2 + u3 + u4
           f2 =~ u5 + u6 + u7 + u8'

modelsCat2 <- measurementInvarianceCat(model, data = datCat, group = "g", 
	parameterization="theta", estimator="wlsmv", strict = TRUE)
	
partialInvarianceCat(modelsCat2, type = "scalar")	
}
}
