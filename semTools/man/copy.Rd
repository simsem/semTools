\name{copy}
\alias{copy}
\title{
	Copy the result of \code{lavaan} object into the clipboard
}
\description{
Copy the result of \code{lavaan} object into the clipboard. The users may paste the result into the Microsoft Excel or spreadsheet application to create a table of the output.
}
\usage{
copy(object, what="summary", ...)
}
\arguments{
  \item{object}{
	The \code{lavaan} object
}
  \item{what}{
	The attributes of the \code{lavaan} object to be copied in the clipboard. \code{"summary"} is to copy the screen provided from the \code{summary} function. \code{"mifit"} is to copy the result from the \code{\link{miPowerFit}} function. Other attributes listed in the \code{inspect} method in the \link[lavaan]{lavaan-class} could also be used, such as \code{"coef"}, \code{"se"}, \code{"fit"}, \code{"samp"}, and so on.
}
  \item{\dots}{
	Additional argument listed in the \code{\link{miPowerFit}} function.
}
}
\value{
	The resulting output will be saved in the object. Users can paste it in the other application.
}
\author{
    Sunthud Pornprasertmanit (University of Kansas; \email{psunthud@ku.edu})
}
\examples{
library(lavaan)
HW.model <- ' visual  =~ x1 + c1*x2 + x3
              textual =~ x4 + c1*x5 + x6
               speed   =~ x7 + x8 + x9 '

fit <- cfa(HW.model, data=HolzingerSwineford1939, group="school", meanstructure=TRUE)

# Copy the summary of the lavaan object
copy(fit)

# Copy the modification indices and the model fit from the miPowerFit function
copy(fit, "mifit")

# Copy the parameter estimates
copy(fit, "coef")

# Copy the standard errors
copy(fit, "se")

# Copy the sample statistics
copy(fit, "samp")

# Copy the fit measures
copy(fit, "fit")
}

