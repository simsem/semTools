\name{clipboard_saveFile}
\alias{clipboard}
\alias{saveFile}
\title{
	Copy or save the result of \code{lavaan} or \code{FitDiff} objects into a clipboard or a file
}
\description{
Copy or save the result of \code{lavaan} or \code{\linkS4class{FitDiff}} object into a clipboard or a file. From the clipboard, users may paste the result into the Microsoft Excel or spreadsheet application to create a table of the output.
}
\usage{
clipboard(object, what="summary", ...)
saveFile(object, file, what="summary", tableFormat=FALSE, ...)
}
\arguments{
  \item{object}{
	The \code{lavaan} or \code{\linkS4class{FitDiff}} object
}
  \item{what}{
	The attributes of the \code{lavaan} object to be copied in the clipboard. \code{"summary"} is to copy the screen provided from the \code{summary} function. \code{"mifit"} is to copy the result from the \code{\link{miPowerFit}} function. Other attributes listed in the \code{inspect} method in the \link[lavaan]{lavaan-class} could also be used, such as \code{"coef"}, \code{"se"}, \code{"fit"}, \code{"samp"}, and so on. For the The \code{\linkS4class{FitDiff}} object, this argument is not active yet.
}
  \item{file}{
	A file name used for saving the result
}
  \item{tableFormat}{
	If \code{TRUE}, save the result in the table format using tabs for seperation. Otherwise, save the result as the output screen printed in the R console.
}
  \item{\dots}{
	Additional argument listed in the \code{\link{miPowerFit}} function (for \code{lavaan} object only).
}
}
\value{
	The resulting output will be saved into a clipboard or a file. If using the \code{clipboard} function, users may paste it in the other applications.
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
\dontrun{
library(lavaan)
HW.model <- ' visual  =~ x1 + c1*x2 + x3
              textual =~ x4 + c1*x5 + x6
               speed   =~ x7 + x8 + x9 '

fit <- cfa(HW.model, data=HolzingerSwineford1939, group="school", meanstructure=TRUE)

# Copy the summary of the lavaan object
clipboard(fit)

# Copy the modification indices and the model fit from the miPowerFit function
clipboard(fit, "mifit")

# Copy the parameter estimates
clipboard(fit, "coef")

# Copy the standard errors
clipboard(fit, "se")

# Copy the sample statistics
clipboard(fit, "samp")

# Copy the fit measures
clipboard(fit, "fit")

# Save the summary of the lavaan object
saveFile(fit, "out.txt")

# Save the modification indices and the model fit from the miPowerFit function
saveFile(fit, "out.txt", "mifit")

# Save the parameter estimates
saveFile(fit, "out.txt", "coef")

# Save the standard errors
saveFile(fit, "out.txt", "se")

# Save the sample statistics
saveFile(fit, "out.txt", "samp")

# Save the fit measures
saveFile(fit, "out.txt", "fit")
}
}

