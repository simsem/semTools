##Parcel Allocation
##Corbin Quick
##02/29/12
##Vector of numbers of indicators in each parcel, vector assigning each indicator to its factor, Number allocations, lavaan syntax, Data set, parcel names, variables left out of parceling, additional arguments to be passed to lavaan

parcelAllocation <- function(Nperpar,Facplc,Nalloc=100,syntax,Dataset,names='default',leaveout=0, ...) {
  
  require(lavaan)
  
  if(is.character(Dataset)){
    Dataset <- read.csv(Dataset)
  }  

  Dataset <- as.matrix(Dataset)
  
  if(Nalloc<2) stop("Minimum of two allocations required.")
  
  if(is.list(Facplc)){ 
    
    if(is.numeric(Facplc[[1]][1])==FALSE){
      Facplcb <- Facplc
      Namesv <- colnames(Dataset)
      
      for(i in 1:length(Facplc)){
        for(j in 1:length(Facplc[[i]])){
          Facplcb[[i]][j] <- match(Facplc[[i]][j],Namesv)
        }
        Facplcb[[i]] <- as.numeric(Facplcb[[i]])
      }
      Facplc <- Facplcb
      
    }
    
    # Facplc2 <- rep(0, sum(sapply(Facplc, length)))
    Facplc2 <- rep(0,ncol(Dataset))
    
    for(i in 1:length(Facplc)){
      for(j in 1:length(Facplc[[i]])){
        Facplc2[Facplc[[i]][j]] <- i
      }
    }
    Facplc <- Facplc2
  }
  
  if(leaveout!=0){
    
    if(is.numeric(leaveout)==FALSE){
      leaveoutb <- rep(0,length(leaveout))
      Namesv <- colnames(Dataset)
      
      for(i in 1:length(leaveout)){
        leaveoutb[i] <- match(leaveout[i],Namesv)
      }
      leaveout <- as.numeric(leaveoutb)
      
    }
    
    k1 <- .001
    for(i in 1:length(leaveout)){
      Facplc[leaveout[i]] <- Facplc[leaveout[i]] + k1
      k1 <- k1 +.001
    }
  }
  
  if(0 %in% Facplc == TRUE){
    Zfreq <- sum(Facplc==0)
    for (i in 1:Zfreq){
      Zplc <- match(0,Facplc)
      Dataset <- Dataset[ , -Zplc]
      Facplc <- Facplc[-Zplc]
    }
    ## this allows for unused variables in dataset, 
    ## which are specified by zeros, and deleted
  }

if(is.list(Nperpar)){
    
  Nperpar2 <- c()
  for (i in 1:length(Nperpar)){
    Onesp <- sum(Facplc>i & Facplc<i+1) 
    Nperpar2 <- c(Nperpar2, Nperpar[i], rep(1, Onesp), recursive = TRUE)
  }
    
  Nperpar <- Nperpar2
}
  
  Npp <- c()
  for (i in 1:length(Nperpar)){
    Npp <- c(Npp, rep(i, Nperpar[i]))
  }
  
  Locate <- sort(round(Facplc))
  Maxv <- max(Locate)-1
  
  if(length(Locate)!=length(Npp)){
    stop('** WARNING! ** Parcels incorrectly specified. Check input!')}
  
  for (i in 1:Maxv){
    Mat <- match(i+1, Locate)
    if(Npp[Mat] == Npp[Mat-1]){ 
      stop('** WARNING! ** Parcels incorrectly specified. Check input!')} 
  }
    ## warning message if parcel crosses into multiple factors
          ## vector, parcel to which each variable belongs
          ## vector, factor to which each variables belongs
      ## if variables are in the same parcel, but different factors
      ## error message given in output
    
  Onevec <- Facplc - round(Facplc)
  NleaveA <- length(Onevec) - sum(Onevec==0)
  NleaveP <- sum(Nperpar==1)
  
  if(NleaveA < NleaveP){
    print('** WARNING! ** Single-variable parcels have been requested. Check input!')}

  if(NleaveA > NleaveP)
    print('** WARNING! ** More non-parceled variables have been requested than provided for in parcel vector. Check input!')
  
  if(length(names)>1){
    if(length(names) != length(Nperpar)){
      print('** WARNING! ** Number of parcel names provided not equal to number of parcels requested. Check input!')}}
  
  if(NA %in% Dataset == TRUE){ 
    print('** WARNING! ** Missing data detected. Prior multiple imputation recommended.')}

  Data <- c(1:ncol(Dataset))
    ## creates a vector of the number of indicators
    ## e.g. for three indicators, c(1, 2, 3)
  Nfactors <- max(Facplc)
    ## scalar, number of factors
  Nindicators <- length(Data)
    ## scalar, number of indicators
  Npar <- length(Nperpar)
    ## scalar, number of parcels
  Rmize <- runif(Nindicators, 1, Nindicators)
    ## create vector of randomly ordered numbers, 
    ## length of number of indicators
  
  Data <- rbind(Facplc, Rmize, Data)
    ## "Data" becomes object of three rows, consisting of
    ##    1) factor to which each indicator belongs
    ##          (in order to preserve indicator/factor 
    ##           assignment during randomization)
    ##    2) randomly order numbers
    ##    3) indicator number
  
  Results <- matrix(numeric(0), Nalloc, Nindicators)
    ##create empty matrix for parcel allocation matrix
  
  Pin <- Nperpar[1]
  for (i in 2:length(Nperpar)){

    Pin <- c(Pin, Nperpar[i]+Pin[i-1])
      ## creates vector which indicates the range
      ## of columns (endpoints) in each parcel
  }
  
  for (i in 1:Nalloc) {
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
  
  Alpha <- rbind(Results[1,], Dataset)
      ## bind first random allocation to dataset "Alpha"
  
  Allocations <- list()
      ## create empty list for allocation data matrices

  for (i in 1:Nalloc){
    
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
    
    Temp <- matrix(NA, nrow(Dataset), Npar)
      ## create empty matrix for averaged parcel variables
    
    TempAA <- if(length(1:Pin[1])>1) Beta[2:nrow(Beta) , 1:Pin[1]] else cbind(Beta[2:nrow(Beta) , 1:Pin[1]],Beta[2:nrow(Beta) , 1:Pin[1]])
    Temp[, 1] <- rowMeans(TempAA)
      ## fill first column with averages from assigned indicators 
    for (al in 2:Npar) {    
      Plc <- Pin[al-1]+1
        ## placeholder variable for determining parcel width
      TempBB <- if(length(Plc:Pin[al])>1) Beta[2:nrow(Beta) , Plc:Pin[al]] else cbind(Beta[2:nrow(Beta) , Plc:Pin[al]],Beta[2:nrow(Beta) , Plc:Pin[al]])
      Temp[, al] <- rowMeans(TempBB)
        ## fill remaining columns with averages from assigned indicators
    }
    
    if(length(names)>1){
      colnames(Temp) <- names
    }
    
    Allocations[[i]] <- Temp
      ## assign result to list of parcel datasets
  }
  
  if(as.vector(regexpr("/",syntax))!=-1){
    replist<-matrix(NA,Nalloc,1)
    for (i in 1:Nalloc){
      if(names!='default'){colnames(Allocations[[i]])<-names}else{colnames(Allocations[[i]])<-NULL}
      write.table(Allocations[[i]],paste(syntax,'parcelruns',i,'.dat',sep=''),row.names=FALSE,col.names=TRUE)
      replist[i,1]<-paste('parcelrun',i,'.dat',sep='')
    }
    write.table(replist,paste(syntax,"parcelrunsreplist.dat",sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE)
  }
  else{
    
  Param <- list()
  ## list for parameter estimated for each imputation
  Fitind <- list()
  ## list for fit indices estimated for each imputation
  
  for (i in 1:Nalloc){
    data <- as.data.frame(Allocations[[i]], row.names = NULL, optional = FALSE)
    ## convert allocation matrix to dataframe for model estimation 
    fit <- sem(syntax, data=data, ...)
    ## estimate model in lavaan
    Param[[i]] <- parameterEstimates(fit)
    ## assign allocation parameter estimates to list
    Fitind[[i]] <- fitMeasures(fit)
    ## assign allocation parameter estimates to list
  }
  
  Parmn <- Param[[1]]
  ## assign first parameter estimates to mean dataframe
  
  ParSE <- matrix(NA, nrow(Parmn), Nalloc)
  ParSEmn <- Parmn[,5]
  
  Parsd <- matrix(NA, nrow(Parmn), Nalloc)
  ## assign parameter estimates for S.D. calculation
  
  Fitmn <- Fitind[[1]]
  ## assign first fit indices to mean dataframe
  
  Fitsd <- matrix(NA, 5, Nalloc)
  ## assign fit indices for S.D. calculation
  
  Sigp <- matrix(NA, nrow(Parmn), Nalloc)
  ## assign p-values to calculate percentage significant
  
  for (i in 1:Nalloc){
    
    Parsd[,i] <- Param[[i]][,4]
    ## assign parameter estimates for S.D. estimation
    
    ParSE[,i] <- Param[[i]][,5]
    
    if(i>1){ParSEmn <- ParSEmn + Param[[i]][,5]}
    
    Sigp[,ncol(Sigp)-i+1] <- Param[[i]][,7]
    ## assign p-values to calculate percentage significant
    
    
    Fitsd[1:5,i] <- c(rbind(Fitind[[i]][1],Fitind[[i]][7],Fitind[[i]][8],Fitind[[i]][16],Fitind[[i]][20]))
    ## assign fit indices for S.D. estimation
    
    if(i>1){Parmn[,4:ncol(Parmn)] <- Parmn[,4:ncol(Parmn)] + Param[[i]][,4:ncol(Parmn)]}
    ## add together all parameter estimates
    
    if(i>1){Fitmn <- Fitmn + Fitind[[i]]}
    ## add together all fit indices
  }
  
  
  Sigp <- Sigp + .45
  Sigp <- apply(Sigp, c(1,2), round)
  Sigp <- 1 - as.vector(rowMeans(Sigp))
  ## calculate percentage significant parameters
  
  Parsum <- cbind(apply(Parsd,1,sd),apply(Parsd,1,max),apply(Parsd,1,min),apply(Parsd,1,max)-apply(Parsd,1,min), Sigp)
  colnames(Parsum) <- c("S.D.","MAX","MIN","Range", "% Sig")
  ## calculate parameter S.D., minimum, maximum, range, bind to percentage significant
  
  ParSEmn <- cbind(Parmn[,1:3], ParSEmn/Nalloc)
  ParSEfn <- cbind(ParSEmn,apply(ParSE,1,sd),apply(ParSE,1,max),apply(ParSE,1,min),apply(ParSE,1,max)-apply(ParSE,1,min))
  colnames(ParSEfn) <- c("lhs", "op", "rhs", "Avg SE","S.D.","MAX","MIN","Range")
  
  Fitsum <- cbind(apply(Fitsd,1,sd),apply(Fitsd,1,max),apply(Fitsd,1,min),apply(Fitsd,1,max)-apply(Fitsd,1,min))
  rownames(Fitsum) <- c("chisq","cfi","tli","rmsea","srmr")
  ## calculate fit S.D., minimum, maximum, range
  
  Parmn[,4:ncol(Parmn)] <- Parmn[,4:ncol(Parmn)] / Nalloc
  ## divide totalled parameter estimates by number allocations
  Parmn <- Parmn[,1:4]
  ## remove confidence intervals from output
  Parmn <- cbind(Parmn, Parsum)
  ## bind parameter average estimates to cross-allocation information
  Fitmn <- Fitmn / Nalloc
  ## divide totalled fit indices by number allocations
  
  Fitsum <- cbind(rbind(Fitmn[1],Fitmn[7],Fitmn[8],Fitmn[16],Fitmn[20]),Fitsum)
  colnames(Fitsum) <- c("Avg Ind","S.D.","MAX","MIN","Range")
  ## bind to fit averages
  
  ParSEfn[,4:8] <- apply(ParSEfn[,4:8], 2, round, digits = 3)
  Parmn[,4:9] <- apply(Parmn[,4:9], 2, round, digits = 3)
  Fitsum <- apply(Fitsum, 2, round, digits = 3)
  ## round output to three digits
  
  Output <-  list(Parmn,ParSEfn,Fitsum)
  names(Output) <- c('Estimates', 'SE', 'Fit')
  
  return(Output)
  
}}