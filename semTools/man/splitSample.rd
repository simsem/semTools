\name{splitSample}
\alias{splitSample}
\title{
Randomly Split a Data Set into Halves
}
\description{
This function randomly splits a data set into two halves, and saves the resulting data sets to the same folder as the original. 
}
\usage{
splitSample(dataset,path="default", div=2, type="default", name="splitSample")
}
\arguments{
  \item{dataset}{The original data set to be divided. Can be a file path to a .csv or .dat file (headers will automatically be detected) or an R object (matrix or dataframe). (Windows users: file path must be specified using FORWARD SLASHES ONLY.)}
  \item{path}{File path to folder for output data sets. NOT REQUIRED if dataset is a filename. Specify ONLY if dataset is an R object, or desired output folder is not that of original data set. If path is specified as "object", output data sets will be returned as a list, and not saved to hard drive. }
  \item{div}{Number of output data sets. NOT REQUIRED if default, 2 halves.}	
  \item{type}{Output file format ("dat" or "csv"). NOT REQUIRED unless desired output formatting differs from that of input, or dataset is an R object and csv formatting is desired.}	
  \item{name}{Output file name. NOT REQUIRED unless desired output name differs from that of input, or input dataset is an R object. (If input is an R object and name is not specified, name will be "splitSample".)}	
}
\details{
This function randomly orders the rows of a data set, divides the data set into two halves, and saves the halves to the same folder as the original data set, preserving the original formatting. Data set type (.csv or .dat) and formatting (headers) are automatically detected, and output data sets will preserve input type and formatting unless specified otherwise. Input can be in the form of a file path (.dat or .csv), or an R object (matrix or dataframe). If input is an R object and path is default, output data sets will be returned as a list object.  
}
\value{
\item{dataL}{List of output data sets. ONLY IF dataset is an R object and path is default. Otherwise, output will saved to hard drive with the same formatting as input.}
}
\author{
    Corbin Quick (University of Michigan; \email{corbinq@umich.edu})
}
\examples{
#### Input is .dat file
#splitSample("C:/Users/Default/Desktop/MYDATA.dat")
#### Output saved to "C:/Users/Default/Desktop/" in .dat format
#### Names are "MYDATA_s1.dat" and "MYDATA_s2.dat"

#### Input is R object
##Split C02 dataset from the datasets package
library(datasets)
splitMyData <- splitSample(CO2, path="object")
summary(splitMyData[[1]])
summary(splitMyData[[2]])
#### Output object splitMyData becomes list of output data sets

#### Input is .dat file in "C:/" folder
#splitSample("C:/testdata.dat", path = "C:/Users/Default/Desktop/", type = "csv")
#### Output saved to "C:/Users/Default/Desktop/" in .csv format
#### Names are "testdata_s1.csv" and "testdata_s2.csv"

#### Input is R object
#splitSample(myData, path = "C:/Users/Default/Desktop/", name = "splitdata")
#### Output saved to "C:/Users/Default/Desktop/" in .dat format
#### Names are "splitdata_s1.dat" and "splitdata_s2.dat"
}
