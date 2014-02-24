\name{Net-class}
\docType{class}
\alias{Net-class}
\alias{show,Net-method}
\alias{summary,Net-method}
\title{
	Class For the Result of Nesting and Equivalence Testing 
}
\description{
	This class contains the results of nesting and equivalence testing among multiple models
}
\section{Objects from the Class}{
	Objects can be created via the \code{\link{net}} function.
}
\section{Slots}{
  \describe{
    \item{\code{test}:}{Logical matrix of results of nesting and equivalence testing across models}
    \item{\code{df}:}{The degrees of freedom of tested models}
  }
}
\section{methods}{
  \itemize{
	\item \code{summary} The summary function is used to provide the results in narrative.  
  }
}
\author{
    Terrence D. Jorgensen (University of Kansas; \email{TJorgensen314@gmail.com})
}
\seealso{
\code{\link{net}}
}
\examples{
# See the example in the net function.
}
