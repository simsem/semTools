% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clipboard.R
\name{clipboard}
\alias{clipboard}
\alias{saveFile}
\title{Copy or save the result of \code{lavaan} or \code{FitDiff} objects into a
clipboard or a file}
\usage{
clipboard(object, what = "summary", ...)

saveFile(object, file, what = "summary", tableFormat = FALSE,
  fit.measures = "default", writeArgs = list(), ...)
}
\arguments{
\item{object}{An object of class \link[lavaan:lavaan-class]{lavaan::lavaan} or
\linkS4class{FitDiff}.}

\item{what}{The attributes of the \code{lavaan} object to be copied in the
clipboard. \code{"summary"} is to copy the screen provided from the
\code{summary} function. \code{"mifit"} is to copy the result from the
\code{\link[=miPowerFit]{miPowerFit()}} function. Other attributes listed in the
\code{inspect} method in the \link[lavaan:lavaan-class]{lavaan::lavaan} could also be
used, such as \code{"coef"}, \code{"se"}, \code{"fit"}, \code{"samp"}, and
so on.  Ignored for \linkS4class{FitDiff}-class objects.}

\item{\dots}{Additional arguments when passing a \code{lavaan} object to the
\code{summary} or \code{\link[=miPowerFit]{miPowerFit()}} function.}

\item{file}{A file name used for saving the result.}

\item{tableFormat}{If \code{TRUE}, save the result in the table format using
tabs for separation. Otherwise, save the result as the output screen
printed in the R console.}

\item{fit.measures}{\code{character} vector specifying names of fit measures
returned by \code{\link[lavaan:fitMeasures]{lavaan::fitMeasures()}} to be copied/saved.  Only
relevant if \code{object} is class \linkS4class{FitDiff}.}

\item{writeArgs}{\code{list} of additional arguments to be passed to
\code{\link[utils:write.table]{utils::write.table()}}}
}
\value{
The resulting output will be saved into a clipboard or a file. If
using the \code{clipboard} function, users may paste it in the other
applications.
}
\description{
Copy or save the result of \code{lavaan} or \linkS4class{FitDiff}
object into a clipboard or a file. From the clipboard, users may paste the
result into the Microsoft Excel or spreadsheet application to create a table
of the output.
}
\examples{

library(lavaan)
HW.model <- ' visual  =~ x1 + c1*x2 + x3
              textual =~ x4 + c1*x5 + x6
              speed   =~ x7 +    x8 + x9 '

fit <- cfa(HW.model, data = HolzingerSwineford1939, group = "school")

if(interactive()){
# Copy the summary of the lavaan object
clipboard(fit)

# pass additional arguments to summary() method for class?lavaan
clipboard(fit, rsquare = TRUE, standardized = TRUE, fit.measures = TRUE)

# Copy modification indices and fit stats from the miPowerFit() function
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

# Save modification indices and fit stats from the miPowerFit() function
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
\author{
Sunthud Pornprasertmanit (\email{psunthud@gmail.com})

Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
