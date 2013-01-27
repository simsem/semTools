\name{lisrel2lavaan}
\alias{lisrel2lavaan}
\title{
	Latent variable modeling in \code{\linkS4class{lavaan}} using LISREL syntax
}
\description{
	Estimate a structural equation model in \code{\linkS4class{lavaan}} using LISREL syntax. Data will be imported from the same directory as the syntax file itself if data files names are provided within LISREL syntax, or from within syntax directly as per LISREL conventions. 
}
\usage{
lisrel2lavaan(filename = NULL, analyze = TRUE, ...)
}
\arguments{
  \item{filename}{
	Filename of the LISREL syntax file. If not specified, user will be prompted with a file browser with which LISREL syntax file can be selected. 
}
  \item{analyze}{
	Logical. If true (default) the data will be analyzed and output displayed. If false, data will not be analyzed; a \code{\linkS4class{lavaan}} parameter table containing the model specifications will be returned. 
}
  \item{\dots}{ 
	Additional arguments to be passed to \code{\link[lavaan]{lavaan}}.
}
}
\value{
	Output summary is printed to screen and \code{\linkS4class{lavaan}} fit object is returned.
}
\author{
    Corbin Quick (University of Kansas; \email{corbinq@ku.edu})
}
\examples{
\dontrun{
## calling lisrel2lavaan without specifying the filename argument will open
## a file browser window with which a LISREL syntax file can be selected. 

## any additional arguments to be passed to lavaan for data analysis can be
## specified normally. 

lisrel2lavaan(se="standard")
## lavaan output summary printed to screen
## lavaan fit object returned silently

## manual file specification 

lisrel2lavaan(filename="myFile.LS8", se="standard")
## lavaan output summary printed to screen
## lavaan fit object returned silently
}
}

