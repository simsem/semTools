### Corbin Quick & Alex Schoemann
### Last updated: 4 April 2017


#' Random Allocation of Items to Parcels in a Structural Equation Model
#' 
#' This function generates a given number of randomly generated item-to-parcel
#' allocations, fits a model to each allocation, and provides averaged results
#' over all allocations.
#' 
#' This function implements the random item to parcel allocation procedure
#' described in Sterba (2011) and Sterba and MacCallum (2010). The function
#' takes a single data set with item level data, randomly assigns items to
#' parcels, fits a structural equation model to the parceled data (using
#' \link{lavaan}), and repeats this process for a user specified number of
#' random allocations. Results from all fitted models are summarized and
#' output. For further details on the benefits of the random allocation of
#' itesm to parcels see Sterba (2011) and Sterba and MccCallum (2010).
#' 
#' @importFrom stats sd runif
#' 
#' @param nPerPar A list in which each element is a vector corresponding to
#' each factor indicating sizes of parcels. If variables are left out of
#' parceling, they should not be accounted for here (there should NOT be
#' parcels of size "1").
#' @param facPlc A list of vectors, each corresponding to a factor, specifying
#' the variables in that factor (whether included in parceling or not). Either
#' variable names or column numbers. Variables not listed will not be modeled
#' or included in output datasets.
#' @param nAlloc The number of random allocations of items to parcels to
#' generate.
#' @param syntax \link{lavaan} syntax. If substituted with a file name,
#' parcelAllocation will print output data sets to a specified folder rather
#' than analyzing using lavaan (note for Windows users: file path must be
#' specified using forward slashes).
#' @param dataset Data set. Can be file path or R object (matrix or dataframe).
#' If the data has missing values multiple imputation before parceling is
#' recommended.
#' @param names (Optional) A character vector containing the names of parceled
#' variables.
#' @param leaveout A vector of variables to be left out of randomized
#' parceling. Either variable names or column numbers are allowed.
#' @param \dots Additional arguments to be passed to \link{lavaan}
#' @return \item{Estimates}{A data frame containing results related to
#' parameter estimates with columns corresponding to parameter names, average
#' parameter estimates across allocations, the standard deviation of parameter
#' estimates across allocations, the minimum parameter estimate across
#' allocations, the maximum parameter estimate across allocations, the range of
#' parameter estimates across allocations, and the proportions of allocations
#' in which the parameter estimate is significant.} \item{SE}{A data frame
#' containing results related to standard errors with columns corresponding to
#' parameter names, average standard errors across allocations, the standard
#' deviation of standard errors across allocations, the minimum standard error
#' across allocations, the maximum standard error across allocations, and the
#' range of standard errors across allocations.} \item{Fit}{A data frame
#' containing results related to model fit with columns corresponding to fit
#' index names, the average of each index across allocations, the standard
#' deviation of each fit index across allocations, the minimum of each fit
#' index across allocations, the maximum of each fit index across allocations,
#' and the range of each fit index across allocations.}
#' @author
#'  Corbin Quick (University of Michigan; \email{corbinq@@umich.edu})
#'  
#' Alexander M. Schoemann (East Carolina University; \email{schoemanna@@ecu.edu})
#' 
#' @seealso \code{\link{PAVranking}}, \code{\link{poolMAlloc}}
#' @references Sterba, S. K. (2011). Implications of parcel-allocation
#' variability for comparing fit of item-solutions and parcel-solutions.
#' \emph{Structural Equation Modeling, 18}(4), 554-577.
#' doi:10.1080/10705511.2011.607073
#' 
#' Sterba, S. K.m & MacCallum, R. C. (2010). Variability in parameter estimates
#' and model fit across random allocations of items to parcels.
#' \emph{Multivariate Behavioral Research, 45}(2), 322-358.
#' doi:10.1080/00273171003680302
#' @examples
#' 
#' ## Fit 3 factor CFA to simulated data.
#' ## Each factor has 9 indicators that are randomly parceled into 3 parcels
#' ## lavaan syntax for the model to be fit to parceled data
#' library(lavaan)
#' 
#' syntax <- 'La =~ V1 + V2 + V3 
#'            Lb =~ V4 + V5 + V6
#' '
#' ## Parcel and fit data 20 times.
#' ## The actual parcel number should be higher than 20 times.
#' name1 <- colnames(simParcel)[1:9]
#' name2 <- colnames(simParcel)[10:18]
#' parcelAllocation(list(c(3,3,3),c(3,3,3)), list(name1, name2), nAlloc = 20,
#'                  syntax = syntax, dataset = simParcel)
#' 
#' @export
parcelAllocation <- function(nPerPar, facPlc, nAlloc = 100, syntax, dataset,
                             names = "default", leaveout = 0, ...) {
  if (is.character(dataset)) dataset <- utils::read.csv(dataset)
  dataset <- as.matrix(dataset)
  
  if (nAlloc < 2) stop("Minimum of two allocations required.")
  
  if (is.list(facPlc)) {
    
    if (is.numeric(facPlc[[1]][1]) == FALSE) {
      facPlcb <- facPlc
      Namesv <- colnames(dataset)
      
      for (i in 1:length(facPlc)) {
        for (j in 1:length(facPlc[[i]])) {
          facPlcb[[i]][j] <- match(facPlc[[i]][j], Namesv)
        }
        facPlcb[[i]] <- as.numeric(facPlcb[[i]])
      }
      facPlc <- facPlcb
    }
    
    # facPlc2 <- rep(0, sum(sapply(facPlc, length)))
    facPlc2 <- rep(0, ncol(dataset))
    
    for (i in 1:length(facPlc)) {
      for (j in 1:length(facPlc[[i]])) {
        facPlc2[facPlc[[i]][j]] <- i
      }
    }
    facPlc <- facPlc2
  }
  
  if (leaveout != 0) {
    
    if (is.numeric(leaveout) == FALSE) {
      leaveoutb <- rep(0, length(leaveout))
      Namesv <- colnames(dataset)
      
      for (i in 1:length(leaveout)) {
        leaveoutb[i] <- match(leaveout[i], Namesv)
      }
      leaveout <- as.numeric(leaveoutb)
      
    }
    
    k1 <- 0.001
    for (i in 1:length(leaveout)) {
      facPlc[leaveout[i]] <- facPlc[leaveout[i]] + k1
      k1 <- k1 + 0.001
    }
  }
  
  if (0 %in% facPlc == TRUE) {
    Zfreq <- sum(facPlc == 0)
    for (i in 1:Zfreq) {
      Zplc <- match(0, facPlc)
      dataset <- dataset[, -Zplc]
      facPlc <- facPlc[-Zplc]
    }
    ## this allows for unused variables in dataset, which are specified by
    ## zeros, and deleted
  }
  
  if (is.list(nPerPar)) {
    nPerPar2 <- c()
    for (i in 1:length(nPerPar)) {
      Onesp <- sum(facPlc > i & facPlc < i + 1)
      nPerPar2 <- c(nPerPar2, nPerPar[i], rep(1, Onesp), recursive = TRUE)
    }
    nPerPar <- nPerPar2
  }
  
  Npp <- c()
  for (i in 1:length(nPerPar)) Npp <- c(Npp, rep(i, nPerPar[i]))
  
  Locate <- sort(round(facPlc))
  Maxv <- max(Locate) - 1
  
  if (length(Locate) != length(Npp)) stop("Parcels incorrectly specified.",
                                          " Check input!")
  
  if (Maxv > 0) {
    ## Bug was here. With 1 factor Maxv=0. Skip this with a single factor
    for (i in 1:Maxv) {
      Mat <- match(i + 1, Locate)
      if (Npp[Mat] == Npp[Mat - 1]) stop("Parcels incorrectly specified.",
                                         " Check input!")
    }
  }
  ## warning message if parcel crosses into multiple factors vector, parcel to
  ## which each variable belongs vector, factor to which each variables belongs
  ## if variables are in the same parcel, but different factors error message
  ## given in output
  
  Onevec <- facPlc - round(facPlc)
  NleaveA <- length(Onevec) - sum(Onevec == 0)
  NleaveP <- sum(nPerPar == 1)
  
  if (NleaveA < NleaveP) warning("Single-variable parcels have",
                                 " been requested. Check input!")
  
  if (NleaveA > NleaveP) warning("More non-parceled variables",
                                 " have been requested than provided for in",
                                 " parcel vector. Check input!")
  
  if (length(names) > 1) {
    if (length(names) != length(nPerPar))
      warning("Number of parcel names provided not equal to",
              " number of parcels requested. Check input!")
  }
  
  if (NA %in% dataset == TRUE) warning("Missing data detected. Prior",
                                       " multiple imputation recommended.")
  
  Data <- c(1:ncol(dataset))
  ## creates a vector of the number of indicators e.g. for three indicators, c(1, 2,
  ## 3)
  Nfactors <- max(facPlc)
  ## scalar, number of factors
  Nindicators <- length(Data)
  ## scalar, number of indicators
  Npar <- length(nPerPar)
  ## scalar, number of parcels
  Rmize <- runif(Nindicators, 1, Nindicators)
  ## create vector of randomly ordered numbers, length of number of indicators
  
  Data <- rbind(facPlc, Rmize, Data)
  ## 'Data' becomes object of three rows, consisting of 1) factor to which each
  ## indicator belongs (in order to preserve indicator/factor assignment during
  ## randomization) 2) randomly order numbers 3) indicator number
  
  Results <- matrix(numeric(0), nAlloc, Nindicators)
  ## create empty matrix for parcel allocation matrix
  
  Pin <- nPerPar[1]
  for (i in 2:length(nPerPar)) {
    
    Pin <- c(Pin, nPerPar[i] + Pin[i - 1])
    ## creates vector which indicates the range of columns (endpoints) in each parcel
  }
  
  for (i in 1:nAlloc) {
    Data[2, ] <- runif(Nindicators, 1, Nindicators)
    ## Replace second row with newly randomly ordered numbers
    
    Data <- Data[, order(Data[2, ])]
    ## Order the columns according to the values of the second row
    
    Data <- Data[, order(Data[1, ])]
    ## Order the columns according to the values of the first row in order to preserve
    ## factor assignment
    
    Results[i, ] <- Data[3, ]
    ## assign result to allocation matrix
  }
  
  Alpha <- rbind(Results[1, ], dataset)
  ## bind first random allocation to dataset 'Alpha'
  
  Allocations <- list()
  ## create empty list for allocation data matrices
  
  for (i in 1:nAlloc) {
    
    Ineff <- rep(NA, ncol(Results))
    Ineff2 <- c(1:ncol(Results))
    for (inefficient in 1:ncol(Results)) {
      Ineff[Results[i, inefficient]] <- Ineff2[inefficient]
    }
    
    Alpha[1, ] <- Ineff
    ## replace first row of dataset matrix with row 'i' from allocation matrix
    
    Beta <- Alpha[, order(Alpha[1, ])]
    ## arrangle dataset columns by values of first row assign to temporary matrix
    ## 'Beta'
    
    Temp <- matrix(NA, nrow(dataset), Npar)
    ## create empty matrix for averaged parcel variables
    
    TempAA <- if (length(1:Pin[1]) > 1) 
      Beta[2:nrow(Beta), 1:Pin[1]] else cbind(Beta[2:nrow(Beta), 1:Pin[1]],
                                              Beta[2:nrow(Beta), 1:Pin[1]])
    Temp[, 1] <- rowMeans(TempAA)
    ## fill first column with averages from assigned indicators
    for (al in 2:Npar) {
      Plc <- Pin[al - 1] + 1
      ## placeholder variable for determining parcel width
      TempBB <- if (length(Plc:Pin[al]) > 1) 
        Beta[2:nrow(Beta), Plc:Pin[al]] else cbind(Beta[2:nrow(Beta), Plc:Pin[al]],
                                                   Beta[2:nrow(Beta), Plc:Pin[al]])
      Temp[, al] <- rowMeans(TempBB)
      ## fill remaining columns with averages from assigned indicators
    }
    
    if (length(names) > 1) {
      colnames(Temp) <- names
    }
    
    Allocations[[i]] <- Temp
    ## assign result to list of parcel datasets
  }
  
  if (as.vector(regexpr("/", syntax)) != -1) {
    replist <- matrix(NA, nAlloc, 1)
    for (i in 1:nAlloc) {
      if (names != "default") {
        colnames(Allocations[[i]]) <- names
      } else {
        colnames(Allocations[[i]]) <- NULL
      }
      utils::write.table(Allocations[[i]],
                         paste(syntax, "parcelruns", i, ".dat", sep = ""),
                         row.names = FALSE, col.names = TRUE)
      replist[i, 1] <- paste("parcelrun", i, ".dat", sep = "")
    }
    utils::write.table(replist, paste(syntax, "parcelrunsreplist.dat", sep = ""),
                       quote = FALSE, row.names = FALSE, col.names = FALSE)
  } else {
    
    Param <- list()
    ## list for parameter estimated for each imputation
    Fitind <- list()
    ## list for fit indices estimated for each imputation
    
    for (i in 1:nAlloc) {
      data <- as.data.frame(Allocations[[i]], row.names = NULL, optional = FALSE)
      ## convert allocation matrix to dataframe for model estimation
      fit <- lavaan::sem(syntax, data = data, ...)
      ## estimate model in lavaan
      Param[[i]] <- lavaan::parameterEstimates(fit)
      ## assign allocation parameter estimates to list
      Fitind[[i]] <- lavaan::fitMeasures(fit, c("chisq", "df", "cfi", "tli", 
                                                "rmsea", "srmr"))
      ## assign allocation parameter estimates to list
    }
    
    Parmn <- Param[[1]]
    ## assign first parameter estimates to mean dataframe
    
    ParSE <- matrix(NA, nrow(Parmn), nAlloc)
    ParSEmn <- Parmn[, 5]
    
    Parsd <- matrix(NA, nrow(Parmn), nAlloc)
    ## assign parameter estimates for S.D. calculation
    
    Fitmn <- Fitind[[1]]
    ## assign first fit indices to mean dataframe
    
    Fitsd <- matrix(NA, length(Fitmn), nAlloc)
    ## assign fit indices for S.D. calculation
    
    Sigp <- matrix(NA, nrow(Parmn), nAlloc)
    ## assign p-values to calculate percentage significant
    
    for (i in 1:nAlloc) {
      Parsd[, i] <- Param[[i]][, 4]
      ## assign parameter estimates for S.D. estimation
      ParSE[, i] <- Param[[i]][, 5]
      if (i > 1) {
        ParSEmn <- ParSEmn + Param[[i]][, 5]
      }
      Sigp[, ncol(Sigp) - i + 1] <- Param[[i]][, 7]
      ## assign p values to calculate percentage significant
      Fitsd[, i] <- Fitind[[i]]
      ## assign fit indices for S.D. estimation
      if (i > 1) {
        Parmn[, 4:ncol(Parmn)] <- Parmn[, 4:ncol(Parmn)] + Param[[i]][, 4:ncol(Parmn)]
      }
      ## add together all parameter estimates
      if (i > 1) Fitmn <- Fitmn + Fitind[[i]]
      ## add together all fit indices
    }
    
    
    Sigp <- Sigp + 0.45
    Sigp <- apply(Sigp, c(1, 2), round)
    Sigp <- 1 - as.vector(rowMeans(Sigp))
    ## calculate percentage significant parameters
    
    Parsum <- cbind(apply(Parsd, 1, sd), apply(Parsd, 1, max),
                    apply(Parsd, 1, min),
                    apply(Parsd, 1, max) - apply(Parsd, 1, min), Sigp)
    colnames(Parsum) <- c("S.D.", "MAX", "MIN", "Range", "% Sig")
    ## calculate parameter SD, min, max, range, bind to percentage significant
    
    ParSEmn <- cbind(Parmn[, 1:3], ParSEmn / nAlloc)
    ParSEfn <- cbind(ParSEmn, apply(ParSE, 1, sd), apply(ParSE, 1, max),
                     apply(ParSE, 1, min),
                     apply(ParSE, 1, max) - apply(ParSE, 1, min))
    colnames(ParSEfn) <- c("lhs","op","rhs","Avg SE","S.D.","MAX","MIN","Range")
    
    Fitsum <- cbind(apply(Fitsd, 1, sd), apply(Fitsd, 1, max),
                    apply(Fitsd, 1, min),
                    apply(Fitsd, 1, max) - apply(Fitsd, 1, min))
    rownames(Fitsum) <- c("chisq", "df", "cfi", "tli", "rmsea", "srmr")
    ## calculate fit S.D., minimum, maximum, range
    
    Parmn[, 4:ncol(Parmn)] <- Parmn[, 4:ncol(Parmn)] / nAlloc
    ## divide totalled parameter estimates by number allocations
    Parmn <- Parmn[, 1:4]
    ## remove confidence intervals from output
    Parmn <- cbind(Parmn, Parsum)
    ## bind parameter average estimates to cross-allocation information
    Fitmn <- Fitmn/nAlloc
    ## divide totalled fit indices by number allocations
    
    Fitsum <- cbind(Fitmn, Fitsum)
    colnames(Fitsum) <- c("Avg Ind", "S.D.", "MAX", "MIN", "Range")
    ## bind to fit averages
    
    ParSEfn[, 4:8] <- apply(ParSEfn[, 4:8], 2, round, digits = 3)
    Parmn[, 4:9] <- apply(Parmn[, 4:9], 2, round, digits = 3)
    Fitsum <- apply(Fitsum, 2, round, digits = 3)
    ## round output to three digits
    
    Output <- list(Parmn, ParSEfn, Fitsum)
    names(Output) <- c("Estimates", "SE", "Fit")
    
    return(Output)
  }
}

#parcelAllocation(list(c(3,3,3)), list(name1), nAlloc=20, syntax=syntax, dataset=simParcel)
