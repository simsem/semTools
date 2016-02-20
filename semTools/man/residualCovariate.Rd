\name{residualCovariate}
\alias{residualCovariate}
\title{
	Residual centered all target indicators by covariates
}
\description{
This function will regress target variables on the covariate and replace the target variables by the residual of the regression analysis. This procedure is useful to control the covariate from the analysis model (Geldhof, Pornprasertmanit, Schoemann, & Little, in press).
}
\usage{
residualCovariate(data, targetVar, covVar)
}
\arguments{
  \item{data}{
	The desired data to be transformed.
}
  \item{targetVar}{
	Varible names or the position of indicators that users wish to be residual centered (as dependent variables)
}
  \item{covVar}{
	Covariate names or the position of the covariates using for residual centering (as independent variables) onto target variables
}
}
\value{
	The data that the target variables replaced by the residuals
}
\references{
Geldhof, G. J., Pornprasertmanit, S., Schoemann, A. M., & Little, T. D. (2013). Orthogonalizing through residual centering: Applications and caveats. \emph{Educational and Psychological Measurement, 73}, 27-46.
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\seealso{
	\code{\link{indProd}} For creating the indicator products with no centering, mean centering, double-mean centering, or residual centering.	
}
\examples{
dat <- residualCovariate(attitude, 2:7, 1)
}
