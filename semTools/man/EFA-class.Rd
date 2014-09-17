\name{EFA-class}
\docType{class}
\alias{EFA-class}
\alias{show,EFA-method}
\alias{summary,EFA-method}
\title{
	Class For Rotated Results from EFA
}
\description{
	This class contains the results of rotated exploratory factor analysis 
}
\section{Objects from the Class}{
	Objects can be created via the \code{\link{orthRotate}} or \code{\link{oblqRotate}} function.
}
\section{Slots}{
  \describe{
    \item{\code{loading}:}{Rotated standardized factor loading matrix}
    \item{\code{rotate}:}{Rotation matrix}
    \item{\code{gradRotate}:}{The gradient of the objective function at the rotated loadings}
    \item{\code{convergence}:}{Convergence status}
    \item{\code{phi}:}{Factor correlation. Will be an identity matrix if orthogonal rotation is used.}
    \item{\code{se}:}{Standard errors of the rotated standardized factor loading matrix}
    \item{\code{method}:}{Method of rotation}
    \item{\code{call}:}{The command used to generate this object}
  }
}
\section{methods}{
  \itemize{
	\item \code{summary} The \code{summary} function shows the detailed results of the rotated solution. This function has two arguments: \code{suppress} and \code{sort}. The \code{suppress} argument is used to not show the standardized loading values that less than the specified value. The default is 0.1. The \code{sort} is used to sort the factor loadings by the sizes of factor loadings in each factor. The default is \code{TRUE}.  
  }
}
\seealso{
\code{\link{efaUnrotate}}; \code{\link{orthRotate}}; \code{\link{oblqRotate}}
}
\author{
    Sunthud Pornprasertmanit (Texas Tech University; \email{sunthud.pornprasertmanit@ttu.edu})
}
\examples{
unrotated <- efaUnrotate(HolzingerSwineford1939, nf=3, varList=paste0("x", 1:9), estimator="mlr")
summary(unrotated, std=TRUE)
inspect(unrotated, "std")

# Rotated by Quartimin
rotated <- oblqRotate(unrotated, method="quartimin")
summary(rotated)
}
