\name{FitDiff-class}
\docType{class}
\alias{FitDiff-class}
\alias{show,FitDiff-method}
\alias{summary,FitDiff-method}
\title{
	Class For Representing A Template of Model Fit Comparisons
}
\description{
	This class contains model fit measures and model fit comparisons among multiple models
}
\section{Objects from the Class}{
	Objects can be created via the \code{\link{compareFit}} function.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{The name of each model}
    \item{\code{nested}:}{Model fit comparisons between adjacent nested models that are ordered based on their degrees of freedom}
    \item{\code{ordernested}:}{The order of nested models regarding to their degrees of freedom}
    \item{\code{fit}:}{Fit measures of all models specified in the \code{name} slot}
  }
}
\section{methods}{
  \itemize{
	\item \code{summary} The summary function is used to provide the nested model comparison results and the summary of the fit indices across models. This function has one argument: \code{fit.measures}. If \code{"default"} is specified, chi-square values, degree of freedom, \emph{p} value, CFI, TLI, RMSEA, SRMR, AIC, and BIC are provided. If \code{"all"} is specified, all information given in the \code{\link[lavaan]{fitMeasures}} function is provided. Users may specify a vector of the name of fit indices that they wish. 
  }
}
\author{
    Sunthud Pornprasertmanit (Texas Tech University; \email{sunthud.pornprasertmanit@ttu.edu})
}
\seealso{
\code{\link{compareFit}}; \code{\link{clipboard}}
}
\examples{
HW.model <- ' visual =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed =~ x7 + x8 + x9 '

out <- measurementInvariance(HW.model, data=HolzingerSwineford1939, group="school", quiet=TRUE)
modelDiff <- compareFit(out)
summary(modelDiff)
summary(modelDiff, fit.measures="all")
summary(modelDiff, fit.measures=c("aic", "bic"))

\dontrun{
# Save results to a file 
saveFile(modelDiff, file="modelDiff.txt")

# Copy to a clipboard
clipboard(modelDiff)
}
}
