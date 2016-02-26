\name{BootMiss-class}
\docType{class}
\alias{BootMiss-class}
\alias{show,BootMiss-method}
\alias{summary,BootMiss-method}
\alias{hist,BootMiss-method}
\title{
	Class For the Results of Bollen-Stine Bootstrap with Incomplete Data
}
\description{
	This class contains the results of Bollen-Stine bootstrap with missing data.
}
\section{Objects from the Class}{
	Objects can be created via the \code{\link{bsBootMiss}} function.
}
\section{Slots}{
  \describe{
    \item{\code{time}:}{A list containing 2 \code{difftime} objects (\code{transform} and \code{fit}), indicating the time elapsed for data transformation and for fitting the model to bootstrap data sets, respectively.}
    \item{\code{transData}:}{Transformed data}
    \item{\code{bootDist}:}{The vector of chi-square values from Bootstrap data sets fitted by the target model}
    \item{\code{origChi}:}{The chi-square value from the original data set}
    \item{\code{df}:}{The degree of freedom of the model}
    \item{\code{bootP}:}{The p-value comparing the original chi-square with the bootstrap distribution}
  }
}
\section{methods}{
  \describe{
    \item{show}{\code{signature(object = "BootMiss"):} The \code{show} function is used to display the results of the Bollen-Stine bootstrap.}
    \item{summary}{\code{signature(object = "BootMiss"):} The summary function prints the same information from the \code{show} method, but also provides information about the time elapsed, as well as the expected (theoretical) and observed (bootstrap) mean and variance of the chi-squared distribution.}
    \item{hist}{\code{signature(x = "BootMiss", ..., alpha = .05, nd = 2, printLegend = TRUE, legendArgs = list(x = "topleft")):} The \code{hist} function provides a histogram for the bootstrap distribution of chi-squared, including observed and critical values from the specified \code{alpha} level. The user can also specify additional graphical parameters to \code{\link[graphics]{hist}} via \code{...}, as well as pass a list of arguments to an optional \code{\link[graphics]{legend}} via \code{legendArgs}. If the user wants more control over customization, \code{hist} returns a list of \code{length == 2}, containing the arguments for the call to \code{hist} and the arguments to the call for \code{legend}, respectively.}
  }
}
\author{
    Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
}
\seealso{
\code{\link{bsBootMiss}}
}
\examples{
# See the example from the bsBootMiss function
}
