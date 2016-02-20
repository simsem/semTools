\name{standardizeMx}
\alias{standardizeMx}
\title{
	Find standardized estimates for OpenMx output
}
\description{
	Find standardized estimates for OpenMx output. This function is applicable for the \code{MxRAMObjective} only. 
}
\usage{
standardizeMx(object, free = TRUE)
}
\arguments{
  \item{object}{
	Target OpenMx output using \code{MxRAMObjective}
}
  \item{free}{
	If \code{TRUE}, the function will show only standardized values of free parameters. If \code{FALSE}, the function will show the results for fixed and free parameters.
  }
}
\value{
	A vector of standardized estimates
}
\seealso{
	\code{\link{saturateMx}}, \code{\link{nullMx}}, \code{\link{fitMeasuresMx}}
}
\author{
    Sunthud Pornprasertmanit (\email{psunthud@gmail.com})
}
\examples{
\dontrun{
library(OpenMx)
data(myFADataRaw)
myFADataRaw <- myFADataRaw[,c("x1","x2","x3","x4","x5","x6")]
oneFactorModel <- mxModel("Common Factor Model Path Specification", 
	type="RAM",
	mxData(
		observed=myFADataRaw, 
		type="raw"
	),
	manifestVars=c("x1","x2","x3","x4","x5","x6"),
	latentVars="F1",
	mxPath(from=c("x1","x2","x3","x4","x5","x6"),
		arrows=2,
		free=TRUE,
		values=c(1,1,1,1,1,1),
		labels=c("e1","e2","e3","e4","e5","e6")
	), 
	# residual variances
	# -------------------------------------
	mxPath(from="F1",
		arrows=2,
		free=TRUE,
		values=1,
		labels ="varF1"
	), 
	# latent variance
	# -------------------------------------
	mxPath(from="F1",
		to=c("x1","x2","x3","x4","x5","x6"),
		arrows=1,
		free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
		values=c(1,1,1,1,1,1),
		labels =c("l1","l2","l3","l4","l5","l6")
	), 
	# factor loadings
	# -------------------------------------
	mxPath(from="one",
		to=c("x1","x2","x3","x4","x5","x6","F1"),
		arrows=1,
		free=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE),
		values=c(1,1,1,1,1,1,0),
		labels =c("meanx1","meanx2","meanx3","meanx4","meanx5","meanx6",NA)
	) 
	# means
	# -------------------------------------
) # close model
# Create an MxModel object
# -----------------------------------------------------------------------------
oneFactorFit <- mxRun(oneFactorModel)      
standardizeMx(oneFactorFit)

# Compare with lavaan
library(lavaan)
script <- "f1 =~ x1 + x2 + x3 + x4 + x5 + x6"
fit <- cfa(script, data=myFADataRaw, meanstructure=TRUE)
standardizedSolution(fit)
}
}
