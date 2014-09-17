\name{lisrel2lavaan}
\alias{lisrel2lavaan}
\title{
	Latent variable modeling in \code{\linkS4class{lavaan}} using LISREL syntax
}
\description{
	This function can be used to estimate a structural equation model in \code{\linkS4class{lavaan}} using LISREL syntax. Data are automatically imported from the LISREL syntax file, or, if data files names are provided within LISREL syntax, from the same directory as the syntax itself, as per standard LISREL data importation. 
}
\usage{
lisrel2lavaan(filename = NULL, analyze = TRUE, silent = FALSE, ...)
}
\arguments{
  \item{filename}{
	Filename of the LISREL syntax file. If the \code{filename} arguement is not specified, the user will be prompted with a file browser with which LISREL syntax file can be selected (recommended). 
}
  \item{analyze}{
	Logical. If \code{analyze==TRUE} (default), data will be automatically imported and analyzed; \code{\linkS4class{lavaan}} summary output displayed and fit object will be returned silently. If \code{analyze==FALSE}, data will not be imported or analyzed; instead, a \code{\linkS4class{lavaan}} parameter table containing the model specifications will be returned. 
}
  \item{silent}{
	Logical. If false (default) the data will be analyzed and output displayed. If true, a fit object will be returned and summary output will not be displayed. 
}
  \item{\dots}{ 
	Additional arguments to be passed to \code{\link[lavaan]{lavaan}}.
}
}
\value{
	Output summary is printed to screen and \code{\linkS4class{lavaan}} fit object is returned.
}
\note{
	\code{lisrel2lavaan} is still in development, and not all LISREL commands are currently functional. A number of known limitations are outlined below. If an error is encountered that is not listed, please contact \email{corbinq@ku.edu}.
	\enumerate{
		\item{data importation}{
	\code{lisrel2lavaan} currently supports .csv, .dat, and most other delimited data formats. However, formats that are specific to LISREL or PRELIS (e.g., the .PSF file format) cannot be imported. \code{lisrel2lavaan} supports raw data, covariance matrices, and correlation matrices (accompanied by a variance vector). Symmetric matrices can either contain lower triangle or full matrix. For MACS structure models, either raw data or summary statistics (that include a mean vector) are supported.
}
	\item{variable labels}{
	Certain variable labels that are permitted in LISREL cannot be supported in \code{lisrel2lavaan}. 
		\item{duplicate labels}{
		Most importantly, no two variables of any kind (including phantom variables) should be given the same label when using \code{lisrel2lavaan}. If multiple variables are given the same label, \code{\link[lavaan]{lavaan}} will estimate an incorrect model. 
}
		\item{numeric character labels}{
		All variable labels are recommended to include non-numeric characters. In addition, the first character in each variable label is recommended to be non-numeric. 
}
		\item{labels not specified}{
		If variable labels are not provided by the user, names will be generated reflecting variable assignment (e.g. 'eta1', 'ksi1'); manifest variables will be in lower case and latent variables in upper case. 
}
}
	\item{OU paragraph}{
	Not all commands in the OU paragraph are presently supported in \code{lisrel2lavaan}. The ME command can be used to specify estimation method; however, not all estimations available in LISREL are currently supported by \code{\link[lavaan]{lavaan}}. If the specified ME is unsupported, \code{lisrel2lavaan} will revert to default estimation. The AD, EP, IT, ND and NP keywords will be ignored. Requests for text files containing starting values (e.g., \code{OU BE}) will also be ignored. 
}
	\item{starting values}{
	Certain functionalities related to starting values in LISREL are not yet operational in \code{lisrel2lavaan}. Note that due to differences in estimation, starting values are not as important in \code{\link[lavaan]{lavaan}} model estimation as in LISREL. 
	\item{text file output}{
	Requests for text files containing starting values for individual matrices in the in the \code{OU} command (e.g., \code{OU BE}) are not currently supported. These requests will be ignored.
}
	\item{MA paragraph}{
	Specification of matrix starting values using the MA command is permitted by providing starting values within syntax directly. However, \code{lisrel2lavaan} has sometimes encountered problems with importation when files are specified following the MA paragraph. 
}
}
}
}
\author{
    Corbin Quick (University of Michigan; \email{corbinq@umich.edu})
}
\examples{
\dontrun{
	## calling lisrel2lavaan without specifying the filename argument will  
	## open a file browser window with which LISREL syntax can be selected. 
	
	## any additional arguments to be passed to lavaan for data analysis can
	## be specified normally. 
	
	lisrel2lavaan(se="standard")
	## lavaan output summary printed to screen
	## lavaan fit object returned silently
	
	## manual file specification 
	
	lisrel2lavaan(filename="myFile.LS8", se="standard")
	## lavaan output summary printed to screen
	## lavaan fit object returned silently
}
}

