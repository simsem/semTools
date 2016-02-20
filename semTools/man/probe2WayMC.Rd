\name{probe2WayMC}
\alias{probe2WayMC}
\title{
Probing two-way interaction on the no-centered or mean-centered latent interaction
}
\description{
Probing interaction for simple intercept and simple slope for the no-centered or mean-centered latent two-way interaction
}
\usage{
probe2WayMC(fit, nameX, nameY, modVar, valProbe)
}
\arguments{
  \item{fit}{The lavaan model object used to evaluate model fit}
  \item{nameX}{The vector of the factor names used as the predictors. The first-order factor will be listed first. The last name must be the name representing the interaction term.}
  \item{nameY}{The name of factor that is used as the dependent variable.}
  \item{modVar}{The name of factor that is used as a moderator. The effect of the other independent factor on each moderator variable value will be probed.}
  \item{valProbe}{The values of the moderator that will be used to probe the effect of the other independent factor.}
}
\details{
Before using this function, researchers need to make the products of the indicators between the first-order factors using mean centering (Marsh, Wen, & Hau, 2004). Note that the double-mean centering may not be appropriate for probing interaction if researchers are interested in simple intercepts. The mean or double-mean centering can be done by the \code{\link{indProd}} function. The indicator products can be made for all possible combination or matched-pair approach (Marsh et al., 2004). Next, the hypothesized model with the regression with latent interaction will be used to fit all original indicators and the product terms. See the example for how to fit the product term below. Once the lavaan result is obtained, this function will be used to probe the interaction.

Let that the latent interaction model regressing the dependent variable (\eqn{Y}) on the independent varaible (\eqn{X}) and the moderator (\eqn{Z}) be
\deqn{
	Y = b_0 + b_1X + b_2Z + b_3XZ + r,
}
where \eqn{b_0} is the estimated intercept or the expected value of \eqn{Y} when both \eqn{X} and \eqn{Z} are 0, \eqn{b_1} is the effect of \eqn{X} when \eqn{Z} is 0, \eqn{b_2} is the effect of \eqn{Z} when \eqn{X} is 0, \eqn{b_3} is the interaction effect between \eqn{X} and \eqn{Z}, and \eqn{r} is the residual term.

For probing two-way interaction, the simple intercept of the independent variable at each value of the moderator (Aiken & West, 1991; Cohen, Cohen, West, & Aiken, 2003; Preacher, Curran, & Bauer, 2006) can be obtained by
\deqn{
	b_{0|X = 0, Z} = b_0 + b_2Z.
}

The simple slope of the independent varaible at each value of the moderator can be obtained by
\deqn{
	b_{X|Z} = b_1 + b_3Z.
}

The variance of the simple intercept formula is
\deqn{
	Var\left(b_{0|X = 0, Z}\right) = Var\left(b_0\right) + 2ZCov\left(b_0, b_2\right) + Z^2Var\left(b_2\right)
}
where \eqn{Var} denotes the variance of a parameter estimate and \eqn{Cov} denotes the covariance of two parameter estimates.

The variance of the simple slope formula is 
\deqn{
	Var\left(b_{X|Z}\right) = Var\left(b_1\right) + 2ZCov\left(b_1, b_3\right) + Z^2Var\left(b_3\right)
}

Wald statistic is used for test statistic.
}
\value{
A list with two elements:
\enumerate{
  \item{SimpleIntercept} The intercepts given each value of the moderator. This element will be shown only if the factor intercept is estimated (e.g., not fixed as 0).
  \item{SimpleSlope} The slopes given each value of the moderator. 
}
In each element, the first column represents the values of the moderators specified in the \code{valProbe} argument. The second column is the simple intercept or simple slope. The third column is the standard error of the simple intercept or simple slope. The fourth column is the Wald (\emph{z}) statistic. The fifth column is the \emph{p}-value testing whether the simple intercepts or slopes are different from 0.
}
\references{
Aiken, L. S., & West, S. G. (1991). Multiple regression: Testing and interpreting interactions. Newbury Park, CA: Sage.

Cohen, J., Cohen, P., West, S. G., & Aiken, L. S. (2003). Applied multiple regression/correlation analysis for the behavioral sciences (3rd ed.). New York: Routledge.

Marsh, H. W., Wen, Z., & Hau, K. T. (2004). Structural equation models of latent interactions: Evaluation of alternative estimation strategies and indicator construction. \emph{Psychological Methods, 9}, 275-300.

Preacher, K. J., Curran, P. J., & Bauer, D. J. (2006). Computational tools for probing interactions in multiple linear regression, multilevel modeling, and latent curve analysis. \emph{Journal of Educational and Behavioral Statistics, 31}, 437-448.
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\seealso{
	\itemize{
		\item \code{\link{indProd}} For creating the indicator products with no centering, mean centering, double-mean centering, or residual centering.
		\item \code{\link{probe3WayMC}} For probing the three-way latent interaction when the results are obtained from mean-centering, or double-mean centering.
		\item \code{\link{probe2WayRC}} For probing the two-way latent interaction when the results are obtained from residual-centering approach.
		\item \code{\link{probe3WayRC}} For probing the two-way latent interaction when the results are obtained from residual-centering approach.
		\item \code{\link{plotProbe}} Plot the simple intercepts and slopes of the latent interaction.
	}
}
\examples{
library(lavaan) 

dat2wayMC <- indProd(dat2way, 1:3, 4:6)

model1 <- "
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
f12 =~ x1.x4 + x2.x5 + x3.x6
f3 =~ x7 + x8 + x9
f3 ~ f1 + f2 + f12
f12 ~~0*f1
f12 ~~ 0*f2
x1 ~ 0*1
x4 ~ 0*1
x1.x4 ~ 0*1
x7 ~ 0*1
f1 ~ NA*1
f2 ~ NA*1
f12 ~ NA*1
f3 ~ NA*1
"

fitMC2way <- sem(model1, data=dat2wayMC, meanstructure=TRUE, std.lv=FALSE)
summary(fitMC2way)

result2wayMC <- probe2WayMC(fitMC2way, c("f1", "f2", "f12"), "f3", "f2", c(-1, 0, 1))
result2wayMC
}
