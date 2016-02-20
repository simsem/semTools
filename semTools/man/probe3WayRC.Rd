\name{probe3WayRC}
\alias{probe3WayRC}
\title{
Probing three-way interaction on the residual-centered latent interaction
}
\description{
Probing interaction for simple intercept and simple slope for the residual-centered latent three-way interaction (Pornprasertmanit, Schoemann, Geldhof, & Little, submitted)
}
\usage{
probe3WayRC(fit, nameX, nameY, modVar, valProbe1, valProbe2)
}
\arguments{
  \item{fit}{The lavaan model object used to evaluate model fit}
  \item{nameX}{The vector of the factor names used as the predictors. The three first-order factors will be listed first. Then the second-order factors will be listeed. The last element of the name will represent the three-way interaction. Note that the fourth element must be the interaction between the first and the second variables. The fifth element must be the interaction between the first and the third variables. The sixth element must be the interaction between the second and the third variables.}
  \item{nameY}{The name of factor that is used as the dependent variable.}
  \item{modVar}{The name of two factors that are used as the moderators. The effect of the independent factor on each combination of the moderator variable values will be probed.}
  \item{valProbe1}{The values of the first moderator that will be used to probe the effect of the independent factor.}
  \item{valProbe2}{The values of the second moderator that will be used to probe the effect of the independent factor.}
}
\details{
Before using this function, researchers need to make the products of the indicators between the first-order factors and residualize the products by the original indicators (Lance, 1988; Little, Bovaird, & Widaman, 2006). The process can be automated by the \code{\link{indProd}} function. Note that the indicator products can be made for all possible combination or matched-pair approach (Marsh et al., 2004). Next, the hypothesized model with the regression with latent interaction will be used to fit all original indicators and the product terms (Geldhof, Pornprasertmanit,  Schoemann, & Little, in press). To use this function the model must be fit with a mean structure. See the example for how to fit the product term below. Once the lavaan result is obtained, this function will be used to probe the interaction.

The probing process on residual-centered latent interaction is based on transforming the residual-centered result into the no-centered result. See Pornprasertmanit, Schoemann, Geldhof, and Little (submitted) for further details. Note that this approach based on a strong assumption that the first-order latent variables are normally distributed. The probing process is applied after the no-centered result (parameter estimates and their covariance matrix among parameter estimates) has been computed See the \code{\link{probe3WayMC}} for further details.
}
\value{
A list with two elements:
\enumerate{
  \item{SimpleIntercept} The intercepts given each value of the moderator. This element will be shown only if the factor intercept is estimated (e.g., not fixed as 0).
  \item{SimpleSlope} The slopes given each value of the moderator. 
}
In each element, the first column represents the values of the first moderator specified in the \code{valProbe1} argument. The second column represents the values of the second moderator specified in the \code{valProbe2} argument. The third column is the simple intercept or simple slope. The fourth column is the standard error of the simple intercept or simple slope. The fifth column is the Wald (\emph{z}) statistic. The sixth column is the \emph{p}-value testing whether the simple intercepts or slopes are different from 0.
}
\references{
Geldhof, G. J., Pornprasertmanit, S., Schoemann, A., & Little, T. D. (in press). Orthogonalizing through residual centering: Applications and caveats. \emph{Educational and Psychological Measurement.}

Lance, C. E. (1988). Residual centering, exploratory and confirmatory moderator analysis, and decomposition of effects in path models containing interactions. \emph{Applied Psychological Measurement, 12}, 163-175.

Little, T. D., Bovaird, J. A., & Widaman, K. F. (2006). On the merits of orthogonalizing powered and product terms: Implications for modeling interactions. \emph{Structural Equation Modeling, 13}, 497-519.

Marsh, H. W., Wen, Z., & Hau, K. T. (2004). Structural equation models of latent interactions: Evaluation of alternative estimation strategies and indicator construction. \emph{Psychological Methods, 9}, 275-300.

Pornprasertmanit, S., Schoemann, A. M., Geldhof, G. J., & Little, T. D. (submitted). \emph{Probing latent interaction estimated with a residual centering approach.} 

}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\seealso{
	\itemize{
		\item \code{\link{indProd}} For creating the indicator products with no centering, mean centering, double-mean centering, or residual centering.
		\item \code{\link{probe2WayMC}} For probing the two-way latent interaction when the results are obtained from mean-centering, or double-mean centering.
		\item \code{\link{probe3WayMC}} For probing the three-way latent interaction when the results are obtained from mean-centering, or double-mean centering.
		\item \code{\link{probe2WayRC}} For probing the two-way latent interaction when the results are obtained from residual-centering approach.
		\item \code{\link{plotProbe}} Plot the simple intercepts and slopes of the latent interaction.
	}
}
\examples{
library(lavaan)

dat3wayRC <- orthogonalize(dat3way, 1:3, 4:6, 7:9)

model3 <- "
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
f3 =~ x7 + x8 + x9
f12 =~ x1.x4 + x2.x5 + x3.x6
f13 =~ x1.x7 + x2.x8 + x3.x9
f23 =~ x4.x7 + x5.x8 + x6.x9
f123 =~ x1.x4.x7 + x2.x5.x8 + x3.x6.x9
f4 =~ x10 + x11 + x12
f4 ~ f1 + f2 + f3 + f12 + f13 + f23 + f123
f1 ~~ 0*f12
f1 ~~ 0*f13
f1 ~~ 0*f123
f2 ~~ 0*f12
f2 ~~ 0*f23
f2 ~~ 0*f123
f3 ~~ 0*f13
f3 ~~ 0*f23
f3 ~~ 0*f123
f12 ~~ 0*f123
f13 ~~ 0*f123
f23 ~~ 0*f123
x1 ~ 0*1
x4 ~ 0*1
x7 ~ 0*1
x10 ~ 0*1
x1.x4 ~ 0*1
x1.x7 ~ 0*1
x4.x7 ~ 0*1
x1.x4.x7 ~ 0*1
f1 ~ NA*1
f2 ~ NA*1
f3 ~ NA*1
f12 ~ NA*1
f13 ~ NA*1
f23 ~ NA*1
f123 ~ NA*1
f4 ~ NA*1
" 

fitRC3way <- sem(model3, data=dat3wayRC, meanstructure=TRUE, std.lv=FALSE)
summary(fitRC3way)

result3wayRC <- probe3WayRC(fitRC3way, c("f1", "f2", "f3", "f12", "f13", "f23", "f123"), 
	"f4", c("f1", "f2"), c(-1, 0, 1), c(-1, 0, 1))
result3wayRC
}
