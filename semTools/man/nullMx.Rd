\name{nullMx}
\alias{nullMx}
\title{
	Analyzing data using a null model
}
\description{
	Analyzing data using a null model by full-information maximum likelihood. In the null model, all means and covariances are free if items are continuous. All covariances are fixed to 0. For ordinal variables, their means are fixed as 0 and their variances are fixed as 1 where their thresholds are estimated. In multiple-group model, all means are variances are separately estimated. 
}
\usage{
nullMx(data, groupLab = NULL)
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
	The \code{MxModel} object which contains the analysis result of the null model.
}
\seealso{
	\code{\link{saturateMx}}, \code{\link{fitMeasuresMx}}, \code{\link{standardizeMx}}
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
\dontrun{
library(OpenMx)
data(demoOneFactor)
nullModel <- nullMx(demoOneFactor)
}
}
