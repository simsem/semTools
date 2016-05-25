\name{PAVranking}
\alias{PAVranking}
\title{
Parcel-Allocation Variability in Model Ranking
}
\description{
	This function quantifies and assesses the consequences of parcel-allocation variability for model ranking of structural equation models (SEMs) that differ in their structural specification but share the same parcel-level measurement specification (see Sterba & Rights, 2016). This function is a modified version of \code{\link{parcelAllocation}} which can be used with only one SEM in isolation. The \code{PAVranking} function repeatedly generates a specified number of random item-to-parcel allocations, and then fits two models to each allocation. Output includes summary information about the distribution of model selection results (including plots) and the distribution of results for each model individually, across allocations within-sample. Note that this function can be used when selecting among more than two competing structural models as well (see instructions below involving \code{seed}).
}
\usage{
PAVranking(nPerPar, facPlc, nAlloc=100, parceloutput = 0,
   syntaxA, syntaxB, dataset, names = NULL, 
   leaveout=0, seed=NA, ...) 
}
\arguments{
  \item{nPerPar}{
	A list in which each element is a vector, corresponding to each factor, indicating sizes of parcels. If variables are left out of parceling, they should not be accounted for here (i.e., there should not be parcels of size "1").
}
  \item{facPlc}{
	A list of vectors, each corresponding to a factor, specifying the item indicators of that factor (whether included in parceling or not). Either variable names or column numbers. Variables not listed will not be modeled or included in output datasets.
  }
  \item{nAlloc}{
	The number of random allocations of items to parcels to generate. 
  }
  \item{syntaxA}{
	lavaan syntax for Model A. Note that, for likelihood ratio test (LRT) results to be interpreted, Model A should be nested within Model B (though the function will still provide results when Models A and B are nonnested).
  }
  \item{syntaxB}{
	lavaan syntax for Model B. Note that, for likelihood ratio test (LRT) results to be appropriate, Model A should be nested within Model B (though the function will still provide results when Models A and B are nonnested).
  }
  \item{dataset}{
	 Item-level dataset
  }
  \item{parceloutput}{
	folder where parceled data sets will be outputted (note for Windows users: file path must  specified using forward slashes).
  }
  \item{seed}{
	(Optional) Random seed used for parceling items. When the same random seed is specified and the program is re-run, the same allocations will be generated. The seed argument can be used to assess parcel-allocation variability in model ranking when considering more than two models. For each pair of models under comparison, the program should be rerun using the same random seed. Doing so ensures that multiple model comparisons will employ the same set of parcel datasets.
  }
  \item{names}{
	(Optional) A character vector containing the names of parceled variables.
  }
  \item{leaveout}{
	(Optional) A vector of variables to be left out of randomized parceling. Either variable names or column numbers are allowed.
  }
  \item{\dots}{
	Additional arguments to be passed to \code{\link[lavaan]{lavaan}} 
  }
}
\details{
This is a modified version of \code{\link{parcelAllocation}} which was, in turn, based on the SAS macro \code{ParcelAlloc} (Sterba & MacCallum, 2010). The \code{PAVranking} function produces results discussed in Sterba and Rights (2016) relevant to the assessment of parcel-allocation variability in model selection and model ranking. Specifically, the \code{PAVranking} function first uses a modified version of parcelAllocation to generate a given number (\code{nAlloc}) of item-to-parcel allocations. Then, \code{PAVranking} provides the following new developments: specifying more than one SEM and producing results for Model A and Model B separately that summarize parcel allocation variability in estimates, standard errors, and fit indices. \code{PAVranking} also newly produces results summarizing parcel allocation variability in model selection index values and model ranking between Models A and B. Additionally, \code{PAVranking} newly allows for nonconverged solutions and outputs the proportion of allocations that converged as well as the proportion of proper solutions (results are summarized for converged and proper allocations only).  

For further details on the benefits of the random allocation of items to parcels, see Sterba (2011) and Sterba and MacCallum (2010).

NOTE: This function requires the \code{lavaan} package. Missing data code needs to be \code{NA}. If function returns \code{"Error in plot.new() : figure margins too large,"} user may need to increase size of the plot window and rerun. 
}
\value{
\item{Estimates_A, Estimates_B}{A table containing results related to parameter estimates (in table Estimates_A for Model A and in table Estimates_B for Model B) with columns corresponding to parameter name, average parameter estimate across allocations, standard deviation of parameter estimate across allocations, the maximum parameter estimate across allocations, the minimum parameter estimate across allocations, the range of parameter estimates across allocations, and the percent of allocations in which the parameter estimate is significant.}
\item{SE_A, SE_B}{A table containing results related to standard errors (in table SE_A for Model A and in table SE_B for Model B) with columns corresponding to parameter name, average standard error across allocations, the standard deviation of standard errors across allocations, the maximum standard error across allocations, the minimum standard error across allocations, and the range of standard errors across allocations.}
\item{Fit_A, Fit_B}{A table containing results related to model fit (in table Fit_A for Model A and in table Fit_B for Model B) with columns corresponding to fit index name, the average of the fit index across allocations, the standard deviation of the fit index across allocations, the maximum of the fit index across allocations, the minimum of the fit index across allocations, the range of the fit index across allocations, and the percent of allocations where the chi-square test of absolute fit was significant.}
\item{LRT Summary, Model A vs. Model B}{A table with columns corresponding to: average likelihood ratio test (LRT) statistic for comparing Model A vs. Model B (null hypothesis is no difference in fit between Models A and B in the population), degrees of freedom (i.e. difference in the number of free parameters between Models A and B), as well as the standard deviation, maximum, and minimum of LRT statistics across allocations, and the percent of allocations where the LRT was significant (indicating preference for the more complex Model B). }
\item{LRT Summary, Model A vs. Model B}{A table with columns corresponding to: average likelihood ratio test (LRT) statistic for comparing Model A vs. Model B (null hypothesis is no difference in fit between Models A and B in the population), degrees of freedom (i.e. difference in the number of free parameters between Models A and B), as well as the standard deviation, maximum, and minimum of LRT statistics across allocations, and the percent of allocations where the LRT was significant (indicating preference for the more complex Model B). }
\item{Fit index differences}{A table containing percentage of allocations where Model A is preferred over Model B according to BIC, AIC, RMSEA, CFI, TLI and SRMR and where Model B is preferred over Model A according to the same indices. Also includes the average amount by which the given model is preferred (calculated only using allocations where it was preferred).}
\item{Fit index difference histograms}{Histograms are automatically outputted showing the distribution of the differences (Model A - Model B) for each fit index and for the p-value of the likelihood ratio difference test.}
\item{Percent of Allocations with | BIC Diff |  > 10}{A table containing the percentage of allocations with (BIC for Model A) - (BIC for Model B) < -10, indicating "very strong evidence" to prefer Model A over Model B and the percentage of allocations with (BIC for Model A) - (BIC for Model B) > 10, indicating "very strong evidence" to prefer Model B over Model A (Raftery, 1995).}
\item{Converged and proper}{A table containing the proportion of allocations that converged for Model A, Model B, and both models, and the proportion of allocations with converged and proper solutions for Model A, Model B, and both models.}
}
\references{
Raftery, A. E. (1995). Bayesian model selection in social research. \emph{Sociological Methodology, 25}, 111-163.

Sterba, S. K. (2011). Implications of parcel-allocation variability for comparing fit of item-solutions and parcel-solutions. \emph{Structural Equation Modeling: A Multidisciplinary Journal, 18}(4), 554-577.

Sterba, S. K., & MacCallum, R. C. (2010). Variability in parameter estimates and model fit across repeated allocations of items to parcels. \emph{Multivariate Behavioral Research, 45}(2), 322-358.

"Sterba, S. K., & Rights, J. D. (2016). Effects of parceling on model selection: Parcel-allocation variability in model ranking. \emph{Psychological Methods}. \url{http://dx.doi.org/10.1037/met0000067}
}
\seealso{
	\code{\link{parcelAllocation}}, \code{\link{poolMAlloc}}
}
\author{
    Jason D. Rights (Vanderbilt University; \email{jason.d.rights@vanderbilt.edu})
	
	The author would also like to credit Corbin Quick and Alexander Schoemann for providing the original parcelAllocation function on which this function is based. 
}
\examples{
\dontrun{
## Lavaan syntax for Model A: a 2 Uncorrelated 
## factor CFA model to be fit to parceled data

parmodelA <- '
   f1 =~ NA*p1f1 + p2f1 + p3f1
   f2 =~ NA*p1f2 + p2f2 + p3f2
   p1f1 ~ 1  
   p2f1 ~ 1
   p3f1 ~ 1
   p1f2 ~ 1
   p2f2 ~ 1
   p3f2 ~ 1
   p1f1 ~~ p1f1 
   p2f1 ~~ p2f1
   p3f1 ~~ p3f1
   p1f2 ~~ p1f2
   p2f2 ~~ p2f2
   p3f2 ~~ p3f2
   f1 ~~ 1*f1
   f2 ~~ 1*f2
   f1 ~~ 0*f2
'

## Lavaan syntax for Model B: a 2 Correlated 
## factor CFA model to be fit to parceled data

parmodelB <- '
   f1 =~ NA*p1f1 + p2f1 + p3f1
   f2 =~ NA*p1f2 + p2f2 + p3f2
   p1f1 ~ 1
   p2f1 ~ 1
   p3f1 ~ 1
   p1f2 ~ 1
   p2f2 ~ 1
   p3f2 ~ 1 
   p1f1 ~~ p1f1 
   p2f1 ~~ p2f1
   p3f1 ~~ p3f1
   p1f2 ~~ p1f2
   p2f2 ~~ p2f2
   p3f2 ~~ p3f2
   f1 ~~ 1*f1
   f2 ~~ 1*f2
   f1 ~~ f2
'

##specify items for each factor
f1name <- colnames(simParcel)[1:9]  
f2name <- colnames(simParcel)[10:18]

##run function
PAVranking(nPerPar=list(c(3,3,3),c(3,3,3)), 
  facPlc=list(f1name,f2name), nAlloc=100, 
  parceloutput=0, syntaxA=parmodelA, 
  syntaxB=parmodelB, dataset = simParcel, 
  names=list("p1f1","p2f1","p3f1","p1f2","p2f2","p3f2"),
  leaveout=0)
}
}
