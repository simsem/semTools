##PoolMAlloc

poolMAlloc <- function(nPerPar, facPlc, nAllocStart, 
   nAllocAdd=0, parceloutput=0, syntax,
   dataset, stopProp, stopValue,
   selectParam = NULL,
   double = FALSE, checkConv=FALSE, 
   names='default', leaveout=0, useTotalAlloc=FALSE, ...) 
{
 
  StartTimeFull <- proc.time()
  #### start clock for calculating loop runtime 
  
  if(is.character(dataset)){
    dataset <- read.csv(dataset)
  }  
  
  
  
  {
  nloop <- 0
  ### start loop counter
  
  nAllocStarttemp <- nAllocStart
  ### save initial nAllocStart for final calculation
  
  options(max.print=1000000) 
  ### allow many tables to be outputted 
  
  BreakCounter <- NA
  ### start break counter for double and checkConv options
  
repeat
{

 
  StartTime <- proc.time()
  #### start clock for calculating loop runtime 
  
  nloop <- nloop + 1
  ## add to loop counter
  
  if (double==TRUE & is.na(BreakCounter)==FALSE) BreakCounter <- BreakCounter + 1
  ### add to break counter after stopping criteria reached
  
  if (checkConv==TRUE & is.na(BreakCounter)==FALSE) BreakCounter <- BreakCounter + 1
  ### add to break counter after stopping criteria reached
  
  
  if (nloop > 1) {
  
##Final output##

 if (is.na(BreakCounter)==TRUE){
  
  Parmn_revFinal <- Parmn_rev[[nloop-1]]
  ## save parameter estimates and pooled se table from previous loop for final output
  
  nConvergedOutput <- nConverged
  ## save # allocations converged from previous loop for final output
  
  nConvergedProperOutput <- nConvergedProper
  ## save # allocations converged and proper from previous loop for final output
  
  PooledSEwithinvarFinal <- PooledSEwithinvar
  ## save pooled se within variance for final output
  
  PooledSEbetweenvarFinal <- PooledSEbetweenvar
  ## save pooled se between variance for final output
  
  PooledSEFinal <- PooledSE
  ## save pooled se between variance for final output
  
  FitsumOutput <- Fitsum
  ## save Fit table from previous loop for final output
  
  nAllocOutput <- nAllocStart - nAllocAdd
  #### save nAlloc for output
  
  AllocationsOutput <- Allocations
  ## save datasets from previous loop for final output
  
  ParamFinal <- Param

  
  }
  

  
  ParamPooledSE_temp <- ParamPooledSE
  ### make current "pre-loop" parameter estimates a temporary vector for comparison with "post-loop" estimates
  
  ParamTest_temp <- ParamTest
  #### make current "pre-loop" parameter estimates a temporary vector for comparison with "post-loop" estimates (parameter estimates only)
  
  PooledSE_temp <- PooledSE
  #### make current "pre-loop" parameter estimates a temporary vector for comparison with "post-loop" estimates (pooled SE only)
  
  ParamPoolSEdiffmin <- abs(ParamPooledSE_temp*stopProp)
  ### create vector of minimum differences to continue looping
  
  ParamPoolSEdiffmin[ParamPoolSEdiffmin<stopValue] <- stopValue
  ### change small values in minimum difference vector to a set minimum value
  
  ParamDiffMin <- abs(ParamTest*stopProp)
  ParamDiffMin[ParamDiffMin<stopValue] <- stopValue
  #### create vector of minimum differences to continue looping for parameter estimates
  
  PooledSEmin <- abs(PooledSE*stopProp)
  PooledSEmin[PooledSEmin<stopValue] <- stopValue
  #### create vector of minimum differences to continue looping for pooled se
  } 
  
## Write parceled datasets
  
    dataset <- as.matrix(dataset)
  
  if(nAllocStart<2) stop("Minimum of two allocations required.")
  
  if(is.list(facPlc)){ 
    
    if(is.numeric(facPlc[[1]][1])==FALSE){
      facPlcb <- facPlc
      Namesv <- colnames(dataset)
      
      for(i in 1:length(facPlc)){
        for(j in 1:length(facPlc[[i]])){
          facPlcb[[i]][j] <- match(facPlc[[i]][j],Namesv)
        }
        facPlcb[[i]] <- as.numeric(facPlcb[[i]])
      }
      facPlc <- facPlcb
      
    }
    
    # facPlc2 <- rep(0, sum(sapply(facPlc, length)))
    facPlc2 <- rep(0,ncol(dataset))
    
    for(i in 1:length(facPlc)){
      for(j in 1:length(facPlc[[i]])){
        facPlc2[facPlc[[i]][j]] <- i
      }
    }
    facPlc <- facPlc2
  }
  
  if(leaveout!=0){
    
    if(is.numeric(leaveout)==FALSE){
      leaveoutb <- rep(0,length(leaveout))
      Namesv <- colnames(dataset)
      
      for(i in 1:length(leaveout)){
        leaveoutb[i] <- match(leaveout[i],Namesv)
      }
      leaveout <- as.numeric(leaveoutb)
      
    }
    
    k1 <- .001
    for(i in 1:length(leaveout)){
      facPlc[leaveout[i]] <- facPlc[leaveout[i]] + k1
      k1 <- k1 +.001
    }
  }
  
  if(0 %in% facPlc == TRUE){
    Zfreq <- sum(facPlc==0)
    for (i in 1:Zfreq){
      Zplc <- match(0,facPlc)
      dataset <- dataset[ , -Zplc]
      facPlc <- facPlc[-Zplc]
    }
    ## this allows for unused variables in dataset, 
    ## which are specified by zeros, and deleted
  }

if(is.list(nPerPar)){
    
  nPerPar2 <- c()
  for (i in 1:length(nPerPar)){
    Onesp <- sum(facPlc>i & facPlc<i+1) 
    nPerPar2 <- c(nPerPar2, nPerPar[i], rep(1, Onesp), recursive = TRUE)
  }
    
  nPerPar <- nPerPar2
}
  
  Npp <- c()
  for (i in 1:length(nPerPar)){
    Npp <- c(Npp, rep(i, nPerPar[i]))
  }
  
  Locate <- sort(round(facPlc))
  Maxv <- max(Locate)-1
  
  if(length(Locate)!=length(Npp)){
    stop('** WARNING! ** Parcels incorrectly specified. Check input!')}
  
if(Maxv > 0){
  ##Bug was here. With 1 factor Maxv=0. Skip this with a single factor
  for (i in 1:Maxv){
    Mat <- match(i+1, Locate)
    if(Npp[Mat] == Npp[Mat-1]){ 
      stop('** WARNING! ** Parcels incorrectly specified. Check input!')} 
  }
  }
    ## warning message if parcel crosses into multiple factors
          ## vector, parcel to which each variable belongs
          ## vector, factor to which each variables belongs
      ## if variables are in the same parcel, but different factors
      ## error message given in output
    
  Onevec <- facPlc - round(facPlc)
  NleaveA <- length(Onevec) - sum(Onevec==0)
  NleaveP <- sum(nPerPar==1)
  
  if(NleaveA < NleaveP){
    print('** WARNING! ** Single-variable parcels have been requested. Check input!')}

  if(NleaveA > NleaveP)
    print('** WARNING! ** More non-parceled variables have been requested than provided for in parcel vector. Check input!')
  
  if(length(names)>1){
    if(length(names) != length(nPerPar)){
      print('** WARNING! ** Number of parcel names provided not equal to number of parcels requested. Check input!')}}
  

  Data <- c(1:ncol(dataset))
    ## creates a vector of the number of indicators
    ## e.g. for three indicators, c(1, 2, 3)
  Nfactors <- max(facPlc)
    ## scalar, number of factors
  Nindicators <- length(Data)
    ## scalar, number of indicators
  Npar <- length(nPerPar)
    ## scalar, number of parcels
  Rmize <- runif(Nindicators, 1, Nindicators)
    ## create vector of randomly ordered numbers, 
    ## length of number of indicators
  
  Data <- rbind(facPlc, Rmize, Data)
    ## "Data" becomes object of three rows, consisting of
    ##    1) factor to which each indicator belongs
    ##          (in order to preserve indicator/factor 
    ##           assignment during randomization)
    ##    2) randomly order numbers
    ##    3) indicator number
  
  Results <- matrix(numeric(0), nAllocStart, Nindicators)
    ##create empty matrix for parcel allocation matrix
  
  Pin <- nPerPar[1]
  for (i in 2:length(nPerPar)){

    Pin <- c(Pin, nPerPar[i]+Pin[i-1])
      ## creates vector which indicates the range
      ## of columns (endpoints) in each parcel
  }
  
  for (i in 1:nAllocStart) {
    Data[2,]<-runif(Nindicators, 1, Nindicators)
      ## Replace second row with newly randomly ordered numbers

    Data <- Data[, order(Data[2,])]
      ## Order the columns according 
      ## to the values of the second row

    Data <- Data[, order(Data[1,])]
      ## Order the columns according
      ## to the values of the first row
      ## in order to preserve factor assignment 

    Results[i,] <- Data[3,]
      ## assign result to allocation matrix
  }
  
  Alpha <- rbind(Results[1,], dataset)
      ## bind first random allocation to dataset "Alpha"
  
  Allocations <- list()
      ## create empty list for allocation data matrices

  for (i in 1:nAllocStart){
    
    Ineff <- rep(NA, ncol(Results))
    Ineff2 <- c(1:ncol(Results))
    for (inefficient in 1:ncol(Results)){
      Ineff[Results[i,inefficient]] <- Ineff2[inefficient]
    }
    
    Alpha[1,] <- Ineff
      ## replace first row of dataset matrix 
      ## with row "i" from allocation matrix
    
    Beta <- Alpha[, order(Alpha[1,])]
      ## arrangle dataset columns by values of first row 
      ## assign to temporary matrix "Beta"
    
    Temp <- matrix(NA, nrow(dataset), Npar)
      ## create empty matrix for averaged parcel variables
    
    TempAA <- if(length(1:Pin[1])>1) Beta[2:nrow(Beta) , 1:Pin[1]] else cbind(Beta[2:nrow(Beta) , 1:Pin[1]],Beta[2:nrow(Beta) , 1:Pin[1]])
    Temp[, 1] <- rowMeans(TempAA,na.rm = TRUE)
      ## fill first column with averages from assigned indicators 
    for (al in 2:Npar) {    
      Plc <- Pin[al-1]+1
        ## placeholder variable for determining parcel width
      TempBB <- if(length(Plc:Pin[al])>1) Beta[2:nrow(Beta) , Plc:Pin[al]] else cbind(Beta[2:nrow(Beta) , Plc:Pin[al]],Beta[2:nrow(Beta) , Plc:Pin[al]])
      Temp[, al] <- rowMeans(TempBB,na.rm = TRUE)
        ## fill remaining columns with averages from assigned indicators
    }
    
    if(length(names)>1){
      colnames(Temp) <- names
    }
    
    Allocations[[i]] <- Temp
      ## assign result to list of parcel datasets
	  
	 
  } 
  
  Param <- list()
  ## list for parameter estimated for each imputation
  Fitind <- list()
  ## list for fit indices estimated for each imputation
  Converged <- list()
  ## list for whether or not each allocation converged
  ProperSolution <- list()
  ## list for whether or not each allocation has proper solutions
  ConvergedProper <- list()
  ## list for whether or not each allocation is converged and proper


  
 for (i in 1:(nAllocStart)){
    data_parcel <- as.data.frame(Allocations[[i]], row.names = NULL, optional = FALSE)
    ## convert allocation matrix to dataframe for model estimation 
    fit <- lavaan::sem(syntax, control=list(iter.max=100), data=data_parcel, ...)
    ## estimate model in lavaan
	if (lavaan::lavInspect(fit, "converged")==TRUE){
	Converged[[i]] <- 1
	} else Converged[[i]] <- 0
	## determine whether or not each allocation converged
    Param[[i]] <- lavaan::parameterEstimates(fit)[,c("lhs","op","rhs","est","se","z","pvalue","ci.lower","ci.upper")]
    ## assign allocation parameter estimates to list
    if (lavaan::lavInspect(fit, "post.check") & Converged[[i]] == 1) {
	ProperSolution[[i]] <- 1
	} else ProperSolution[[i]] <- 0
	## determine whether or not each allocation has proper solutions
    if (any(is.na(Param[[i]][,5]==TRUE))) ProperSolution[[i]] <- 0
	## make sure each allocation has existing SE
	if (Converged[[i]]==1 & ProperSolution[[i]]==1) {
	ConvergedProper[[i]] <- 1
	} else ConvergedProper[[i]] <- 0
	## determine whether or not each allocation converged and has proper solutions

	if (ConvergedProper[[i]]==0) Param[[i]][,4:9] <- matrix(data=NA,nrow(Param[[i]]),6)
	## make parameter estimates null for nonconverged, improper solutions
	
    if (ConvergedProper[[i]]==1) {
	Fitind[[i]] <- lavaan::fitMeasures(fit,  c("chisq", "df", "cfi", "tli", "rmsea"))
    } else Fitind[[i]] <- c(NA,NA,NA,NA,NA)
	### assign allocation parameter estimates to list 
	
	}
  

  nConverged <- Reduce("+",Converged)
  ## count number of converged allocations
  
  nProperSolution <- Reduce("+",ProperSolution)
  ## count number of allocations with proper solutions
  
  nConvergedProper <- Reduce("+",ConvergedProper)
  ## count number of allocations with proper solutions
  
  if (nConvergedProper==0) stop("All allocations failed to converge and/or yielded improper solutions for a given loop.")
  ## stop program if no allocations converge
  
  Parmn <- Param[[1]]
  ## assign first parameter estimates to mean dataframe
  if(is.null(selectParam)) selectParam <- 1:nrow(Parmn)
  
  ParSE <- matrix(NA, nrow(Parmn), nAllocStart)
  ParSEmn <- Parmn[,5]
  
  Parsd <- matrix(NA, nrow(Parmn), nAllocStart)
  ## assign parameter estimates for S.D. calculation
  
  Fitmn <- Fitind[[1]]
  ## assign first fit indices to mean dataframe
  
  Fitsd <- matrix(NA, length(Fitmn), nAllocStart)
  ## assign fit indices for S.D. calculation
  
  Sigp <- matrix(NA, nrow(Parmn), nAllocStart)
  ## assign p-values to calculate percentage significant
  
  Fitind <- data.frame(Fitind)
  ## convert fit index table to dataframe

  ParamSEsquared <- list()
  #### create empty list for squared SE
  
  
  for (i in 1:nAllocStart){
    
	ParamSEsquared[[i]] <- cbind(Param[[i]][,5],Param[[i]][,5])
	if (any(is.na(ParamSEsquared[[i]])==TRUE)) ParamSEsquared[[i]] <- 0
	ParamSEsquared[[i]] <- apply(as.data.frame(ParamSEsquared[[i]]),1,prod) 
	### square SE for each allocation
	
    Parsd[,i] <- Param[[i]][,4]
    ## assign parameter estimates for S.D. estimation
    
    ParSE[,i] <- Param[[i]][,5]
    
    Sigp[,ncol(Sigp)-i+1] <- Param[[i]][,7]
    ## assign p-values to calculate percentage significant
    
    Fitsd[,i] <- Fitind[[i]]
    ## assign fit indices for S.D. estimation
  }
 
  
  Sigp <- Sigp + .45
  Sigp <- apply(Sigp, c(1,2), round)
  Sigp <- 1 - as.vector(rowMeans(Sigp, na.rm = TRUE))
  ## calculate percentage significant parameters
    
   
  Parsum <- cbind(apply(Parsd,1,mean,na.rm=TRUE),apply(Parsd,1,sd,na.rm=TRUE),apply(Parsd,1,max,na.rm=TRUE),apply(Parsd,1,min,na.rm=TRUE),apply(Parsd,1,max,na.rm=TRUE)-apply(Parsd,1,min,na.rm=TRUE), Sigp)
  colnames(Parsum) <- c("Avg Est.","S.D.","MAX","MIN","Range", "% Sig")
  ## calculate parameter S.D., minimum, maximum, range, bind to percentage significant
  
  ParSEmn <- Parmn[,1:3]
  ParSEfn <- cbind(ParSEmn,apply(ParSE,1,mean,na.rm=TRUE),apply(ParSE,1,sd,na.rm=TRUE),apply(ParSE,1,max,na.rm=TRUE),apply(ParSE,1,min,na.rm=TRUE),apply(ParSE,1,max,na.rm=TRUE)-apply(ParSE,1,min,na.rm=TRUE))
  colnames(ParSEfn) <- c("lhs", "op", "rhs", "Avg SE","S.D.","MAX","MIN","Range")
  
  Fitsum <- cbind(apply(Fitsd,1,mean,na.rm=TRUE),apply(Fitsd,1,sd,na.rm=TRUE),apply(Fitsd,1,max,na.rm=TRUE),apply(Fitsd,1,min,na.rm=TRUE),apply(Fitsd,1,max,na.rm=TRUE)-apply(Fitsd,1,min,na.rm=TRUE))
  rownames(Fitsum) <- c("chisq", "df", "cfi", "tli", "rmsea")
  ## calculate fit S.D., minimum, maximum, range
  
  Parmn[,4:ncol(Parmn)] <- Parmn[,4:ncol(Parmn)] / nConvergedProper
  ## divide totalled parameter estimates by number converged allocations
  Parmn <- Parmn[,1:3]
  ## remove confidence intervals from output
  Parmn <- cbind(Parmn, Parsum)
  ## bind parameter average estimates to cross-allocation information
  Fitmn <- Fitmn / nConvergedProper
  ## divide totalled fit indices by number converged allocations
  
  pChisq <- list()
  ## create empty list for Chi-square p-values
  sigChisq <- list()
  ## create empty list for Chi-square significance 
 
  for (i in 1:nAllocStart){
  
  pChisq[[i]] <- (1-pchisq(Fitsd[1,i],Fitsd[2,i]))
  ## calculate p-value for each Chi-square
  
  if (is.na(pChisq[[i]])==FALSE & pChisq[[i]]<.05) {
  sigChisq[[i]] <- 1
  } else sigChisq[[i]] <- 0 
  }
  ## count number of allocations with significant chi-square
    
  PerSigChisq <- (Reduce("+",sigChisq))/nConvergedProper*100
  PerSigChisq <- round(PerSigChisq,4)
  ## calculate percent of allocations with significant chi-square
   
  PerSigChisqCol <- c(PerSigChisq,"n/a","n/a","n/a","n/a")
  ## create list of Chi-square Percent Significant and "n/a" 
  
  options(stringsAsFactors=FALSE)
  ## set default option to allow strings into dataframe without converting to factors
  
  Fitsum <- data.frame(Fitsum,PerSigChisqCol)
  colnames(Fitsum) <- c("Avg Ind","S.D.","MAX","MIN","Range","% Sig")
  ### bind to fit averages (changed to dataframe)
  
  options(stringsAsFactors=TRUE)
  ## unset option to allow strings into dataframe without converting to factors;
  
  PooledSEwithinvar <- Reduce("+",ParamSEsquared)/nConvergedProper
  #### calculate within variance for pooled SE
  
  PooledSEbetweenvar <- Parmn[,5]^2
  ## calculate between variance for pooled SE
 
  PooledSE <- sqrt(PooledSEwithinvar + PooledSEbetweenvar + PooledSEbetweenvar/nConvergedProper)
  ### calculate pooled SE
  
  ParamPooledSE <- c(Parmn[,4],PooledSE)
  ### create vector of "post-loop" paramater estimates and pooled SE
  

  
  
  ParamTest <- Parmn[,4]
  #### create vector of parameter estimates
  
  if (nloop>1){
  
  ParamPoolSEdiff <- abs(ParamPooledSE_temp - ParamPooledSE)
  ### create vector of absolute differences between "pre-loop" and "post-loop" vectors
  
  Paramdiff <- abs(ParamTest_temp - ParamTest)
  #### create vector of absolute differences between "pre-loop" and "post-loop" vectors (parameter estimates only)
  
  PooledSEdiff <- abs(PooledSE - PooledSE_temp)
  #### create vector of absolute differences between "pre-loop" and "post-loop" vectors (pooled SE only)
   
  ParamPoolSEdifftest <- ParamPoolSEdiff - ParamPoolSEdiffmin
  ParamPoolSEdifftest[ParamPoolSEdifftest<=0] <- 0
  ParamPoolSEdifftest[ParamPoolSEdifftest>0] <- 1
  ##create vector of difference between (absolute differences between "pre-loop" and "post-loop" vectors) 
  ##and (minimum differences required to continue looping) and set all negative values to 0
     
  Paramdifftest <- Paramdiff - ParamDiffMin
  Paramdifftest[Paramdifftest<=0] <- 0
  Paramdifftest[Paramdifftest>0] <- 1
  PooledSEdifftest <- PooledSEdiff - PooledSEmin
  PooledSEdifftest[PooledSEdifftest<=0] <- 0
  PooledSEdifftest[PooledSEdifftest>0] <- 1
  ##create vector of difference between (absolute differences between "pre-loop" and "post-loop" vectors) 
  ##and (minimum differences required to continue looping) and set all negative values to 0 (parameter estimates and pooled SE separately) 
  
  if (nloop==2){
  ParamPoolSEdifftesttable <- cbind(ParamPoolSEdifftest)
  Paramdifftesttable <- cbind(Paramdifftest)
  PooledSEdifftesttable <- cbind(PooledSEdifftest)
  ### create table of whether or not parameter estimates/ pooled se met stopping criteria for each parameter 
  }
  
  if (nloop>2){
  ParamPoolSEdifftesttable <- cbind(ParamPoolSEdifftesttable,ParamPoolSEdifftest)
  Paramdifftesttable <- cbind(Paramdifftesttable,Paramdifftest)
  PooledSEdifftesttable <- cbind(PooledSEdifftesttable,PooledSEdifftest)
  ##create table indicating whether or not parameter estimates/ pooled se met stopping criteria for each parameter 
  }
  
  PropStopParam <- 1-(Reduce("+",Paramdifftesttable[selectParam,nloop-1])/length(selectParam))
  PropStopPooled <- 1-(Reduce("+",PooledSEdifftesttable[selectParam,nloop-1])/length(selectParam))
  PropStopParamPooled <- 1-(Reduce("+",ParamPoolSEdifftesttable[c(selectParam,selectParam+nrow(Parmn)),nloop-1])/(2*length(selectParam)))
  ##calculate proportion of values meeting stopping criteria
 	 
	 
  if (checkConv==TRUE & is.na(BreakCounter)==TRUE) {
  print(nAllocStart)
  print("Proportion of pooled estimates meeting stop criteria:")
  print(PropStopParam)
  print("Proportion of pooled SE meeting stop criteria:")
  print(PropStopPooled)
  #### print number of allocations, proportion of parameters meeting stop criteria, and proportion of pooled SE meeting stop criteria
  }
  
  if (checkConv==FALSE){
  print(nAllocStart)
  print("Proportion of pooled estimates meeting stop criteria:")
  print(PropStopParam)
  print("Proportion of pooled SE meeting stop criteria:")
  print(PropStopPooled)
  #### print number of allocations, proportion of parameters meeting stop criteria, and proportion of pooled SE meeting stop criteria
  }
  
  
  }
		 
  nAllocStart <- nAllocStart + nAllocAdd
  ### update # allocations for potential next loop
  
  StopTime <- proc.time() - StartTime
  #### calculate time taken to run loop
  
  print("Runtime:")
  print(StopTime)
  #### print time needed for loop
  
  Parmn_rev <- list()
  Parmn_rev[[nloop]] <- cbind(Parmn[,1:4],PooledSE)
  Parmn_rev[[nloop]][,4:5] <- sapply(Parmn_rev[[nloop]][,4:5],as.numeric)
  colnames(Parmn_rev[[nloop]]) <- c("lhs","op","rhs","Estimate","Pooled SE")
  #### calc estimates + pooled SE table
  

  if (nloop==1){
  Param_revTemp <- cbind(Parmn[,1:3],Parmn_rev[[nloop]][,4])
  Param_revTemp[,4] <- as.numeric(Param_revTemp[,4])
  Param_revTotal <- cbind(Param_revTemp)
  PooledSE_revTemp <- cbind(Parmn[,1:3],Parmn_rev[[nloop]][,5])
  PooledSE_revTemp[,4] <- as.numeric(PooledSE_revTemp[,4])
  PooledSE_revTotal <- cbind(PooledSE_revTemp)
  }
  
  if (nloop>1){
  Param_revTemp <- cbind(Parmn_rev[[nloop]][,4])
  Param_revTemp <- as.numeric(Param_revTemp)
  Param_revTotal <- cbind(Param_revTotal,Param_revTemp)
  PooledSE_revTemp <- cbind(Parmn_rev[[nloop]][,5])
  PooledSE_revTemp <- as.numeric(PooledSE_revTemp)
  PooledSE_revTotal <- cbind(PooledSE_revTotal,PooledSE_revTemp)  
  }
  ## create table of parameter estimates and pooled se for each loop
  
if (nloop==1){
  
  ParamTotal <- Param
  
  FitindTotal <- Fitind
  
  AllocationsTotal <- Allocations
  
    nAllocTotal <- nAllocStart - nAllocAdd
  
  nConvergedTotal <- nConverged
  
  nProperSolutionTotal <- nProperSolution
  
  nConvergedProperTotal <- nConvergedProper
  }
  
  if (nloop>1){
  
  ParamTotal <- c(ParamTotal, Param)
   
  FitindTotal <- c(FitindTotal, Fitind)
  
  AllocationsTotal <- c(AllocationsTotal, Allocations)
  
  nAllocTotal <- nAllocTotal + nAllocStart - nAllocAdd
  
  nConvergedTotal <- nConverged + nConvergedTotal
  
  nProperSolution <- nProperSolution + nProperSolutionTotal
  
  nConvergedProperTotal <- nConvergedProper + nConvergedProperTotal  

  }

  
  #print(Parmn_rev[[nloop]])
  #print(ParSEfn)
  #### print all relevant tables
  

  if (nloop>1 & double==TRUE & is.na(BreakCounter)==FALSE & BreakCounter==2){
  if (Reduce("+",ParamPoolSEdifftesttable[c(selectParam,selectParam+nrow(Parmn_rev[[nloop]])),nloop-1])==0)
  break;
  ### with double option selected, break loop after two consecutive hits
  }
  
  if (nloop>1 & double==TRUE){
  if (Reduce("+",ParamPoolSEdifftesttable[c(selectParam,selectParam+nrow(Parmn_rev[[nloop]])),nloop-1])==0){
  BreakCounter <- 1
  } else BreakCounter <- NA
  ### with double option selected, start break counter if stopping criteria are met, otherwise reset BreakCounter to NA
  }
  
  
  if (nloop>1 & checkConv==TRUE & is.na(BreakCounter)==TRUE){
  if (Reduce("+",ParamPoolSEdifftesttable[c(selectParam,selectParam+nrow(Parmn_rev[[nloop]])),nloop-1])==0)
  BreakCounter <- 0
  ### with checkConv option, start break counter if stopping criteria are met
  }
  
  if (nloop>1 & double==FALSE & checkConv==FALSE){
  if (Reduce("+",ParamPoolSEdifftesttable[c(selectParam,selectParam+nrow(Parmn_rev[[nloop]])),nloop-1])==0)
  break;
  }
  ### break loop if differences between "pre-loop" and "post-loop" estimates are sufficiently small
  
  if (nAllocAdd==0)
  break;
  
  ### break loop if nAllocAdd is 0
  

  
  if (checkConv==TRUE & is.na(BreakCounter)==FALSE & BreakCounter==9)
  break;
  ### for checkConv option, break loop after 9 loops after stopping criteria met

  
  }
  

##Write objects for Output when nAllocAdd is set to 0

if (nAllocAdd==0){

  Parmn_revFinal <- Parmn_rev[[nloop]]
  ## save parameter estimates and pooled se table from previous loop for final output
  
  nConvergedOutput <- nConverged
  ## save # allocations converged from previous loop for final output
  
  nConvergedProperOutput <- nConvergedProper
  ## save # allocations converged and proper from previous loop for final output
  
  PooledSEwithinvarFinal <- PooledSEwithinvar
  ## save pooled se within variance for final output
  
  PooledSEbetweenvarFinal <- PooledSEbetweenvar
  ## save pooled se between variance for final output
  
  PooledSEFinal <- PooledSE
  ## save pooled se between variance for final output
  
  FitsumOutput <- Fitsum
  ## save Fit table from previous loop for final output
  
  nAllocOutput <- nAllocStart - nAllocAdd
  #### save nAlloc for output
  
  AllocationsOutput <- Allocations
  ## save datasets from previous loop for final output
}
  
  
##Write parceled datasets

if(as.vector(regexpr("/",parceloutput))!=-1){
    replist<-matrix(NA,nAllocOutput,1)
	for (i in 1:(nAllocOutput)){
	  colnames(AllocationsOutput[[i]])<-names
      write.table(AllocationsOutput[[i]],paste(parceloutput,'/parcelruns',i,'.dat',sep=''),row.names=FALSE,col.names=TRUE)
      replist[i,1]<-paste('parcelruns',i,'.dat',sep='')
    }
    write.table(replist,paste(parceloutput,"/parcelrunsreplist.dat",sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE) }
	
  
  
##Results for using all Allocations  


  if (useTotalAlloc==TRUE) 
{
 
   ParmnTotal <- ParamTotal[[1]]
  ## assign first parameter estimates to mean dataframe
  
  ParSETotal <- matrix(NA, nrow(ParmnTotal), nAllocTotal)
  ParSEmnTotal <- ParmnTotal[,5]
  
  ParsdTotal <- matrix(NA, nrow(ParmnTotal), nAllocTotal)
  ## assign parameter estimates for S.D. calculation
  
  FitmnTotal <- FitindTotal[[1]]
  ## assign first fit indices to mean dataframe
  
  FitsdTotal <- matrix(NA, length(FitmnTotal), nAllocTotal)
  ## assign fit indices for S.D. calculation
  
  SigpTotal <- matrix(NA, nrow(ParmnTotal), nAllocTotal)
  ## assign p-values to calculate percentage significant
  
  FitindTotal <- data.frame(FitindTotal)
  ## convert fit index table to dataframe

  ParamSEsquaredTotal <- list()
  #### create empty list for squared SE
  
  

    for (i in 1:nAllocTotal){
    
	ParamSEsquaredTotal[[i]] <- cbind(ParamTotal[[i]][,5],ParamTotal[[i]][,5])
	if (any(is.na(ParamSEsquaredTotal[[i]])==TRUE)) ParamSEsquaredTotal[[i]] <- 0
	ParamSEsquaredTotal[[i]] <- apply(as.data.frame(ParamSEsquaredTotal[[i]]),1,prod) 
	### square SE for each allocation
	
    ParsdTotal[,i] <- ParamTotal[[i]][,4]
    ## assign parameter estimates for S.D. estimation
    
    ParSETotal[,i] <- ParamTotal[[i]][,5]
    
    SigpTotal[,ncol(Sigp)-i+1] <- ParamTotal[[i]][,7]
    ## assign p-values to calculate percentage significant
    
    FitsdTotal[,i] <- FitindTotal[[i]]
    ## assign fit indices for S.D. estimation
    
  }
 
  
  SigpTotal <- SigpTotal + .45
  SigpTotal <- apply(SigpTotal, c(1,2), round)
  SigpTotal <- 1 - as.vector(rowMeans(SigpTotal, na.rm = TRUE))
  ## calculate percentage significant parameters
  

  
  ParsumTotal <- cbind(apply(ParsdTotal,1,mean,na.rm=TRUE),apply(ParsdTotal,1,sd,na.rm=TRUE),apply(ParsdTotal,1,max,na.rm=TRUE),apply(ParsdTotal,1,min,na.rm=TRUE),apply(ParsdTotal,1,max,na.rm=TRUE)-apply(ParsdTotal,1,min,na.rm=TRUE), SigpTotal)
  colnames(ParsumTotal) <- c("Avg Est.","S.D.","MAX","MIN","Range", "% Sig")
  ## calculate parameter S.D., minimum, maximum, range, bind to percentage significant
  
  ParSEmnTotal <- ParmnTotal[,1:3]
  ParSEfnTotal <- cbind(ParSEmnTotal,apply(ParSETotal,1,mean,na.rm=TRUE),apply(ParSETotal,1,sd,na.rm=TRUE),apply(ParSETotal,1,max,na.rm=TRUE),apply(ParSETotal,1,min,na.rm=TRUE),apply(ParSETotal,1,max,na.rm=TRUE)-apply(ParSETotal,1,min,na.rm=TRUE))
  colnames(ParSEfnTotal) <- c("lhs", "op", "rhs", "Avg SE","S.D.","MAX","MIN","Range")
  
  FitsumTotal <- cbind(apply(FitsdTotal,1,mean,na.rm=TRUE),apply(FitsdTotal,1,sd,na.rm=TRUE),apply(FitsdTotal,1,max,na.rm=TRUE),apply(FitsdTotal,1,min,na.rm=TRUE),apply(FitsdTotal,1,max,na.rm=TRUE)-apply(FitsdTotal,1,min,na.rm=TRUE))
  rownames(FitsumTotal) <- c("chisq", "df", "cfi", "tli", "rmsea")
  ## calculate fit S.D., minimum, maximum, range
  
  ParmnTotal[,4:ncol(ParmnTotal)] <- ParmnTotal[,4:ncol(Parmn)] / nConvergedProperTotal
  ## divide totalled parameter estimates by number converged allocations
  ParmnTotal <- ParmnTotal[,1:3]
  ## remove confidence intervals from output
  ParmnTotal <- cbind(ParmnTotal, ParsumTotal)
  ## bind parameter average estimates to cross-allocation information
  FitmnTotal <- FitmnTotal / nConvergedProperTotal
  ## divide totalled fit indices by number converged allocations
  
  pChisqTotal <- list()
  ## create empty list for Chi-square p-values
  sigChisqTotal <- list()
  ## create empty list for Chi-square significance 
 
  for (i in 1:nAllocTotal){
  
  pChisqTotal[[i]] <- (1-pchisq(FitsdTotal[1,i],FitsdTotal[2,i]))
  ## calculate p-value for each Chi-square
  
  if (is.na(pChisqTotal[[i]])==FALSE & pChisqTotal[[i]]<.05) {
  sigChisqTotal[[i]] <- 1
  } else sigChisqTotal[[i]] <- 0 
  }
  ## count number of allocations with significant chi-square
    
  PerSigChisqTotal <- (Reduce("+",sigChisqTotal))/nConvergedProperTotal*100
  PerSigChisqTotal <- round(PerSigChisqTotal,4)
  ## calculate percent of allocations with significant chi-square
   
  PerSigChisqColTotal <- c(PerSigChisqTotal,"n/a","n/a","n/a","n/a")
  ## create list of Chi-square Percent Significant and "n/a" (used for fit summary table)
  
  options(stringsAsFactors=FALSE)
  ## set default option to allow strings into dataframe without converting to factors
  
  FitsumTotal <- data.frame(FitsumTotal,PerSigChisqColTotal)
  colnames(FitsumTotal) <- c("Avg Ind","S.D.","MAX","MIN","Range","% Sig")
  ### bind to fit averages (changed to dataframe)
  
  options(stringsAsFactors=TRUE)
  ## unset option to allow strings into dataframe without converting to factors;
  
  PooledSEwithinvarTotal <- Reduce("+",ParamSEsquaredTotal)/nConvergedProperTotal
  #### calculate within variance for pooled SE
  
  PooledSEbetweenvarTotal <- ParmnTotal[,5]^2
  ## calculate between variance for pooled SE
 
  PooledSETotal <- sqrt(PooledSEwithinvarTotal + PooledSEbetweenvarTotal + PooledSEbetweenvarTotal/nConvergedProperTotal)
  ### calculate pooled SE
  
  ParamPooledSETotal <- c(ParmnTotal[,4],PooledSETotal)
  ### create vector of "post-loop" paramater estimates and pooled SE
  

  
  
  ParamTestTotal <- ParmnTotal[,4]
  #### create vector of parameter estimates
  
  
  
  #Parmn_revTotal <- list()
  Parmn_revTotal <- cbind(ParmnTotal[,1:4],PooledSETotal)
  Parmn_revTotal[,4:5] <- sapply(Parmn_revTotal[,4:5],as.numeric)
  colnames(Parmn_revTotal) <- c("lhs","op","rhs","Estimate","Pooled SE")
  #### calc estimates + pooled SE table
  

  
  df_tTotal <- (nConvergedProperTotal-1)*(1 + (nConvergedProperTotal*PooledSEwithinvarTotal)/(nConvergedProperTotal*PooledSEbetweenvarTotal + PooledSEbetweenvarTotal))^2
  crit_tTotal <- abs(qt(0.05/2, df_tTotal)) 
  ### compute degrees of freedom and critical value for t 
  
  pval_zTotal <- 2*(1-pnorm(abs(Parmn_revTotal[,4]/PooledSETotal)))
  pval_tTotal <- 2*(1-pt(abs(Parmn_revTotal[,4]/PooledSETotal),df=df_tTotal))
  ### calc p-value for z and t distribution 
  
  
  CI95_Lower_zTotal <- Parmn_revTotal[,4]-1.959963985*PooledSETotal
  CI95_Upper_zTotal <- Parmn_revTotal[,4]+1.959963985*PooledSETotal
  ## compute confidence interval for z-tests
  
  CI95_Lower_tTotal <- Parmn_revTotal[,4]-crit_tTotal*PooledSETotal
  CI95_Upper_tTotal <- Parmn_revTotal[,4]+crit_tTotal*PooledSETotal
  ## compute confidence interval for t-tests
  
  Parmn_revTotal <- cbind(Parmn_revTotal,pval_zTotal,CI95_Lower_zTotal,CI95_Upper_zTotal,pval_tTotal,CI95_Lower_tTotal,CI95_Upper_tTotal)
  colnames(Parmn_revTotal) <- c("lhs","op","rhs","Pooled Est","Pooled SE","pval_z","CI95_LB_z","CI95_UB_z","pval_t","CI95_LB_t","CI95_UB_t")
  ## add confidence intervals to final output table
  
  for (i in 1:nrow(Parmn_revTotal)){
  if (Parmn_revTotal[i,5]==0) Parmn_revTotal[i,6:11] <- NA
  }
  ## make all z/t p-values and CI's NA for fixed parameters (or when pooled se = 0) 
  
  RPAVTotal <- (PooledSEbetweenvarTotal+(PooledSEbetweenvarTotal/(nConvergedProperTotal)))/PooledSEwithinvarTotal
  PPAVTotal <- (((nConvergedProperTotal+1)/(nConvergedProperTotal))*PooledSEbetweenvarTotal)/(PooledSEwithinvarTotal+(((nConvergedProperTotal+1)/(nConvergedProperTotal))*PooledSEbetweenvarTotal))
  PAVtableTotal <- cbind(ParmnTotal[1:3],RPAVTotal,PPAVTotal)
  ### create table for RPAV and PPAV
  

  
  Parmn_revTotal[,4:11] <- apply(Parmn_revTotal[,4:11], 2, round, digits = 4)
  FitsumTotal[,1:5] <- apply(FitsumTotal[,1:5], 2, round, digits = 4)
 
 
  PAVtableTotal[,4:5] <- apply(PAVtableTotal[,4:5], 2, round, digits = 4)
  ### round output to three digits
  
  FitsumTotal[2,2:5] <- c("n/a","n/a","n/a","n/a")
  ## Change df row to "n/a" for sd, max, min, and range
  
 
  ConvergedProperSumTotal <- rbind((nConvergedTotal)/(nAllocTotal),(nConvergedProperTotal)/(nAllocTotal))
  rownames(ConvergedProperSumTotal) <- c("Converged","Converged and Proper")
  colnames(ConvergedProperSumTotal) <- "Proportion of Allocations"
  ### create table summarizing proportions of converged allocations and allocations with proper solutions
  
  }
  
 
  
  
##Output results     

  if (nAllocAdd!=0){
  if (nloop==2) PropParamMet <- matrix(data=1,nrow(Parmn),1)
  if (nloop==2) PropPooledSEMet <- matrix(data=1,nrow(Parmn),1)
  if (nloop !=2) PropParamMet <- (1-apply(Paramdifftesttable[,1:nloop-1],1,mean))*100
  if (nloop !=2) PropPooledSEMet <- (1-apply(PooledSEdifftesttable[,1:nloop-1],1,mean))*100
  #### calc percent of loops where stopping criteria were met for parameters and pooledse
  
  FirstParamMet <- apply(Paramdifftesttable==0,1,which.max)
  FirstPooledSEMet <- apply(PooledSEdifftesttable==0,1,which.max)
  #### determine first loop in which stopping criteria were met for parameters and pooledse
  }
  
  if (nAllocAdd==0){
  PropParamMet <- matrix(data=NA,nrow(Parmn),1)
  PropPooledSEMet <- matrix(data=NA,nrow(Parmn),1)
  FirstParamMet <- matrix(data=NA,nrow(Parmn),1)
  FirstPooledSEMet <- matrix(data=NA,nrow(Parmn),1)
  }
  ### if only running one loop, change columns regarding stopping criteria to NA
	 

  PerLoops <- cbind(Parmn[,1:3],PropParamMet,PropPooledSEMet)
  colnames(PerLoops) <- c("lhs","op","rhs","Param Criteria Met","PooledSE Criteria Met")  
  FirstLoops <- cbind(Parmn[,1:3],FirstParamMet,FirstPooledSEMet)
  colnames(FirstLoops) <- c("lhs","op","rhs","Param Criteria Met","PooledSE Criteria Met")
  NumbAllocations <- cbind(Parmn[,1:3],(FirstParamMet-1)*nAllocAdd+nAllocStarttemp,(FirstPooledSEMet-1)*nAllocAdd+nAllocStarttemp)
  colnames(NumbAllocations) <- c("lhs","op","rhs","Param Criteria Met","PooledSE Criteria Met")
  ### create tables with parameter estimates, pooled SE, and critical value
  
  if (nAllocAdd!=0){
  for (i in 1:nrow(Parmn)){
  if ((i %in% selectParam)==FALSE) PerLoops[i,4:5] <- NA
  if ((i %in% selectParam)==FALSE) FirstLoops[i,4:5] <- NA
  if ((i %in% selectParam)==FALSE) NumbAllocations[i,4:5] <- NA
  ### if  parameter is not used for stopping criteria, change "percent of loops when met" and "loop when first met" to NA
  }
  } 
  

  df_t <- (nConvergedProperOutput-1)*(1 + (nConvergedProperOutput*PooledSEwithinvarFinal)/(nConvergedProperOutput*PooledSEbetweenvarFinal + PooledSEbetweenvarFinal))^2
  crit_t <- abs(qt(0.05/2, df_t)) 
  ### compute degrees of freedom and critical value for t 
  
  pval_z <- 2*(1-pnorm(abs(Parmn_revFinal[,4]/PooledSEFinal)))
  pval_t <- 2*(1-pt(abs(Parmn_revFinal[,4]/PooledSEFinal),df=df_t))
  ### calc p-value for z and t distribution 
  
  
  CI95_Lower_z <- Parmn_revFinal[,4]-1.959963985*PooledSEFinal
  CI95_Upper_z <- Parmn_revFinal[,4]+1.959963985*PooledSEFinal
  ## compute confidence interval for z-tests
  
  CI95_Lower_t <- Parmn_revFinal[,4]-crit_t*PooledSEFinal
  CI95_Upper_t <- Parmn_revFinal[,4]+crit_t*PooledSEFinal
  ## compute confidence interval for t-tests
  
  Parmn_revFinal <- cbind(Parmn_revFinal,pval_z,CI95_Lower_z,CI95_Upper_z,pval_t,CI95_Lower_t,CI95_Upper_t)
  colnames(Parmn_revFinal) <- c("lhs","op","rhs","Pooled Est","Pooled SE","pval_z","CI95_LB_z","CI95_UB_z","pval_t","CI95_LB_t","CI95_UB_t")
  ## add confidence intervals to final output table
  
  for (i in 1:nrow(Parmn_revFinal)){
  if (Parmn_revFinal[i,5]==0) Parmn_revFinal[i,6:11] <- NA
  }
  ## make all z/t p-values and CI's NA for fixed parameters (or when pooled se = 0) 
  
  RPAV <- (PooledSEbetweenvarFinal+(PooledSEbetweenvarFinal/(nConvergedProperOutput)))/PooledSEwithinvarFinal
  PPAV <- (((nConvergedProperOutput+1)/(nConvergedProperOutput))*PooledSEbetweenvarFinal)/(PooledSEwithinvarFinal+(((nConvergedProperOutput+1)/(nConvergedProperOutput))*PooledSEbetweenvarFinal))
  PAVtable <- cbind(Parmn[1:3],RPAV,PPAV)
  ### create table for RPAV and PPAV
  
  colnames(Param_revTotal) <- c("lhs","op","rhs",c(1:nloop))
  colnames(PooledSE_revTotal) <- c("lhs","op","rhs",c(1:nloop))
  ### create column names for tables with parameters estimates and pooled se for each loop

  Param_revTotal[,4:(nloop+3)] <- sapply(Param_revTotal[,4:(nloop+3)], as.numeric)
  PooledSE_revTotal[,4:(nloop+3)] <- sapply(PooledSE_revTotal[,4:(nloop+3)], as.numeric)
  
  Parmn_revFinal[,4:11] <- apply(Parmn_revFinal[,4:11], 2, round, digits = 4)
  FitsumOutput[,1:5] <- apply(FitsumOutput[,1:5], 2, round, digits = 4)
  if (nAllocAdd!=0) Param_revTotal[,4:(nloop+3)] <- apply(Param_revTotal[,4:(nloop+3)], 2, round, digits = 8)
  if (nAllocAdd==0) Param_revTotal[,4] <- round(Param_revTotal[,4],8)
  if (nAllocAdd!=0) PooledSE_revTotal[,4:(nloop+3)] <- apply(PooledSE_revTotal[,4:(nloop+3)], 2, round, digits = 8)
  if (nAllocAdd==0) PooledSE_revTotal[,4] <- round(PooledSE_revTotal[,4],8)
  PAVtable[,4:5] <- apply(PAVtable[,4:5], 2, round, digits = 4)
  ### round output to three digits
  
  FitsumOutput[2,2:5] <- c("n/a","n/a","n/a","n/a")
  ## Change df row to "n/a" for sd, max, min, and range
  
 
  ConvergedProperSum <- rbind((nConvergedOutput)/(nAllocOutput),(nConvergedProperOutput)/(nAllocOutput))
  rownames(ConvergedProperSum) <- c("Converged","Converged and Proper")
  colnames(ConvergedProperSum) <- "Proportion of Allocations"
  ### create table summarizing proportions of converged allocations and allocations with proper solutions
  
  

  #Output_mod <-  list(Parmn_revFinal,FitsumOutput,ConvergedProperSum,nAllocOutput,PAVtable,Param_revTotal,PooledSE_revTotal)
  #names(Output_mod) <- c("Estimates","Fit","Proportion of Converged and Proper Allocations", "Allocations needed for stability (M)",
  #"Indices to quantify uncertainty in estimates due to sampling vs. allocation variability","Pooled Estimates by Loop","Pooled SE by Loop")
  ### output summary for model estimation when checkConv is true (includes results by loop)

    
  StopTimeFull <- proc.time() - StartTimeFull
  #### calculate time taken to run loop
  
  
if (useTotalAlloc==FALSE){
  
 Output_mod <-  list(Parmn_revFinal,FitsumOutput,ConvergedProperSum,nAllocOutput,PAVtable,StopTimeFull[[3]]/60)
  names(Output_mod) <- c("Estimates","Fit","Proportion of Converged and Proper Allocations", "Allocations needed for stability (M)","Indices to quantify uncertainty in estimates due to sampling vs. allocation variability","Total runtime (minutes)")
  ### output summary for model estimation
 }
  
if (useTotalAlloc==TRUE){

 Output_mod <-  list(Parmn_revFinal,FitsumOutput,ConvergedProperSum,nAllocOutput,PAVtable,Parmn_revTotal,FitsumTotal,ConvergedProperSumTotal,nAllocTotal,PAVtableTotal,StopTimeFull[[3]]/60)
  names(Output_mod) <- c("Estimates (using M allocations)","Fit (using M allocations)","Proportion of Converged and Proper Allocations (using M allocations)", "Allocations needed for stability (M)","Indices to quantify uncertainty in estimates due to sampling vs. allocation variability (using M allocations)",
  "Estimates (using all allocations)","Fit (using all allocations)","Proportion of Converged and Proper Allocations (using all allocations)", "Total Allocations used by algorithm","Indices to quantify uncertainty in estimates due to sampling vs. allocation variability (using all allocations)","Total runtime (minutes)")
  ### output summary for model estimation
}

  
  }
 
  return(Output_mod)
  ### returns output for model

}

