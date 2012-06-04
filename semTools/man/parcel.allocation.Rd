\name{parcelAllocation}
\alias{parcelAllocation}
\title{
Random Allocation of Items to Parcels in a Structural Equation Model
}
\description{
This function generates a given number of randomly generated item-to-parcel allocations, fits a model to each allocation, and provides averaged results over all allocations. 
}
\usage{
parcelAllocation(Nperpar, Facplc, Nalloc=100, syntax=output='', Dataset, names='default', leaveout=0, ...)
}
\arguments{
  \item{Nperpar}{A list in which each element is a vector corresponding to each factor indicating sizes of parcels. If variables are left out of parceling, they need not be accounted for here (i.e., parcels of size "1").}
  \item{Facpl}{A list in which each element is a vector corresponding to each factor indicating manifest variables in respective factor. Either variable names or column numbers. Variables not included listed will not be included in output datasets. }
  \item{Nalloc}{The number of random allocations of items to parcels to generate.)	
  \item{syntax}{\link{lavaan} syntax. If substituted with a file name, parcelAllocation will print output data sets to a specified folder rather than analyzing using lavaan (note for Windows users: file path must be specified using forward slashes).}
  \item{Dataset}{Data set. Can be file path or R object (matrix or dataframe). If the data has missing values multiple imputation before parceling is recommended.}
  \item{names}{(Optional) A character vector containing the names of parceled variables.}
  \item{leaveout}{A vector of variables to be left out of randomized parceling. Either variable names or column numbers are allowed.}
   \item{...}{Additional arguments to be passed to \link{lavaan}}
  }
 
\details{
This function implements the random item to parcel allocation procedure described in Sterba (2011) and Sterba and MccCallum (2010). The function takes a single data set with item level data, randomly assigns items to parcels, fits a structural equation model to the parceled data (using \link{lavaan}), and repeats this process for a user specified number of random allocations. Results from all fitted models are summarized and output. For further details on the benefits of the random allocation of itesm to parcels see Sterba (2011) and Sterba and MccCallum (2010). 
}

\values{
A list with the following elements:
\enumerate{
  \item{Estimates}{A data frame containing results related to parameter estimates with columns corresponding to parameter names, average parameter estimates across allocations, the standard deviation of parameter estimates across allocations, the minimum parameter estimate across allocations, the maximum parameter estimate across allocations, the range of parameter estimates across allocations, and the proportions of allocations in which the parameter estimate is significant.}
  \item{SE}{A data frame containing results related to standard errors with columns corresponding to parameter names, average standard errors across allocations, the standard deviation of standard errors across allocations, the minimum standard error across allocations, the maximum standard error across allocations, and the range of standard errors across allocations.}
  \item{Fit}{A data frame containing results related to model fit with columns corresponding to fit index names, the average of each index across allocations, the standard deviation of each fit index across allocations, the minimum of each fit index across allocations, the maximum of each fit index across allocations, and the range of each fit index across allocations.}
  }
  }
 
\references{
Sterba, S.K. (2011). Implications of parcel-allocation variability for comparing fit of item-solutions and parcel-solutions. \emph{Structural Equation Modeling, 18,} 554-577.
Sterba, S.K. & MacCallum, R.C. (2010). Variability in parameter estimates and model fit across random allocations of items to parcels. \emph{Multivariate Behavioral Research, 45,} 322-358.
}

\authors{
    Corbin Quick (University of Kansas; \email{corbinq@ku.edu})
	Alexander M. Schoemann (University of Kansas; \email{schoemann@ku.edu})
}

\examples{

#Fit 3 factor CFA to simulated data.
#Each factor has 9 indicators that are randomly parceled into 3 parcels
#Lavaan syntax for the model to be fit to parceled data
syntax <- 'La =~ V1 + V2 + V3 
           Lb =~ V4 + V5 + V6
'

#Parcel and fit data 100 times 
parcelAllocation(list(c(3,3,3),c(3,3,3)), list(colnames[1:9],colnames[10:18]), Nalloc=100, syntax=syntax, Dataset=simParcel)
}