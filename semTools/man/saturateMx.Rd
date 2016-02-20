\name{saturateMx}
\alias{saturateMx}
\title{
	Analyzing data using a saturate model
}
\description{
	Analyzing data using a saturate model by full-information maximum likelihood. In the saturate model, all means and covariances are free if items are continuous. For ordinal variables, their means are fixed as 0 and their variances are fixed as 1--their covariances and thresholds are estimated. In multiple-group model, all means are variances are separately estimated. 
}
\usage{
saturateMx(data, groupLab = NULL)
}
\arguments{
  \item{data}{
	The target data frame
}
  \item{groupLab}{
	The name of grouping variable
  }
}
\value{
	The \code{MxModel} object which contains the analysis result of the saturate model.
}
\seealso{
	\code{\link{nullMx}}, \code{\link{fitMeasuresMx}}, \code{\link{standardizeMx}}
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
\dontrun{
library(OpenMx)
data(demoOneFactor)
satModel <- saturateMx(demoOneFactor)
}
}
