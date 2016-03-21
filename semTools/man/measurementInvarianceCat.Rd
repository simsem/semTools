\name{measurementInvarianceCat}
\alias{measurementInvarianceCat}
\title{
Measurement Invariance Tests for Categorical Items
}
\description{
Testing measurement invariance across groups using a typical sequence of
model comparison tests.
}
\usage{
measurementInvarianceCat(..., std.lv = FALSE, strict = FALSE, quiet = FALSE,
fit.measures = "default", method = "satorra.bentler.2001")
}
\arguments{
  \item{...}{The same arguments as for any lavaan model.
    See \code{\link{cfa}} for more information.}
  \item{std.lv}{If \code{TRUE}, the fixed-factor method of scale identification is used. If \code{FALSE}, the first variable for each factor is used as marker variable.}
  \item{strict}{If \code{TRUE}, the sequence requires `strict' invariance.
    See details for more information.}
  \item{quiet}{If \code{TRUE}, a summary is printed out containing an
    overview of the different models that are fitted, together with some
    model comparison tests.}
  \item{fit.measures}{Fit measures used to calculate the differences between nested models.}
  \item{method}{The method used to calculate likelihood ratio test. See \code{\link[lavaan]{lavTestLRT}} for available options}
}
\details{
Theta parameterization is used to represent SEM for categorical items. 
That is, residual variances are modeled instead of the total variance of underlying normal variate for each item. 
Five models can be tested based on different constraints across groups.
\enumerate{
  \item{Model 1: configural invariance. The same factor structure is imposed
        on all groups.}
  \item{Model 2: weak invariance. The factor loadings are constrained to be
        equal across groups.}
  \item{Model 3: strong invariance. The factor loadings and thresholds are
        constrained to be equal across groups.}
  \item{Model 4: strict invariance. The factor loadings, thresholds and
        residual variances are constrained to be equal across groups. For 
		categorical variables, all residual variances are fixed as 1.}
  \item{Model 5: The factor loadings, threshoulds, residual variances and means
        are constrained to be equal across groups.}
}

However, if all items have two items (dichotomous), scalar invariance and 
weak invariance cannot be separated because thresholds need to be equal across 
groups for scale identification. Users can specify \code{strict} option to 
include the strict invariance model for the invariance testing. See the further details 
of scale identification and different parameterization in Millsap and Yun-Tein (2004).
}
\value{
Invisibly, all model fits in the sequence are returned as a list.
}
\references{
Millsap, R. E., & Yun-Tein, J. (2004). Assessing factorial invariance in ordered-categorical measures. \emph{Multivariate Behavioral Research, 39}, 479-515.
}
\author{
	Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
    Yves Rosseel (Ghent University; \email{Yves.Rosseel@UGent.be}) 
}
\seealso{
	\code{\link{measurementInvariance}} for measurement invariance for continuous variables; 
	\code{\link{longInvariance}} For the measurement invariance test within person with continuous variables; 
	\code{partialInvariance} for the automated function for finding partial invariance models
}
\examples{
\dontrun{
model <- ' f1 =~ u1 + u2 + u3 + u4'

measurementInvarianceCat(model, data = datCat, group = "g", parameterization="theta", 
    estimator="wlsmv", ordered = c("u1", "u2", "u3", "u4"))
}
}
