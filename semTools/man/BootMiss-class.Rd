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
    \item{\code{timeTrans}:}{Time (in seconds) using for data transformation}
    \item{\code{timeFit}:}{Time (in seconds) using for fitting Bootstrap data sets}
    \item{\code{transData}:}{Transformed data}
    \item{\code{bootDist}:}{The vector of chi-square values from Bootstrap data sets fitted by the target model}
    \item{\code{origChi}:}{The chi-square value from the original data set}
    \item{\code{df}:}{The degree of freedom of the model}
    \item{\code{bootP}:}{The p-value comparing the original chi-square with the bootstrap distribution}
  }
}
\section{methods}{
  \itemize{
	\item \code{summary} The summary function is used to provide the detailed results of the Bollen-Stine bootstrap.  
	\item \code{hist} The hist function is used to provided the bootstrap distribution compared with the original chi-square value.
  }
}
\author{
    Terrence D. Jorgensen (University of Kansas; \email{TJorgensen314@gmail.com})
}
\seealso{
\code{\link{bsBootMiss}}
}
\examples{
# See the example from the bsBootMiss function
}
