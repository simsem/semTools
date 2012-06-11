splitSample<-function(dataset,path="default",div=2,type="default",name="splitSample"){

  type1<-type
  hea=FALSE
  file<-dataset
  
  if(is.character(file)){
    temp <- strsplit(file,'/',fixed=TRUE)
    if(path=="default"){
      path<-paste(temp[[1]][1:(length(temp[[1]])-1)],"/",sep='',collapse="")
    }
    fileN <- temp[[1]][length(temp[[1]])]
    temp <- strsplit(fileN,'.',fixed=TRUE)
    type <- temp[[1]][2]
    name <- temp[[1]][1]
    if(type=='dat'){
      if(is.numeric(as.matrix(read.table(file, nrows=1)))==FALSE){
        data <- as.matrix(read.table(file,header=TRUE))
        hea=TRUE
      }
      else{data <- as.matrix(read.table(file))}
    }
    if(type=='csv'){
      if(is.numeric(as.matrix(read.table(file, nrows=1)))==FALSE){
        data <- as.matrix(read.csv(file,header=TRUE))
        hea=TRUE
      }else{data <- as.matrix(read.csv(file))}
    }
  }else{
    if(is.matrix(file) | is.data.frame(file)){
      data <- as.matrix(file)
    }else{stop("PROVIDE DATA IN .DAT OR .CSV FORMAT")}
  }
  
  if(type1!="default"){
    type<-type1
  }
  
  if(is.character(colnames(data))){
    hea=TRUE
  }
  
  random <- runif(nrow(data),1,nrow(data))
  data <- cbind(random, data)
  data <- data[order(random),]
  data <- data[,2:ncol(data)]
  
  size<-split((1:nrow(data)),cut((1:nrow(data)),div,labels=FALSE))
  size<-as.matrix(as.data.frame(lapply(size,length)))
  dataL <- list()
  dataL[[1]] <- data[1:size[1,1],]
  for(i in 2:div){
    size[1,i]<-size[1,(i-1)]+size[1,i]
    dataL[[i]] <- data[(size[1,(i-1)]+1):size[1,i],]
  }
  
  if(path=='default'){
    return(dataL)} 
  else{
    if(path=="object"){
      return(dataL)}
    else{
      for(i in 1:div){
        if(type=="dat"){
          write.table(dataL[[i]],paste(path,name,"_s",i,".dat",sep=''),sep='  ',row.names=FALSE,col.names=hea)}
        if(type=="csv"){
          write.table(dataL[[i]],paste(path,name,"_s",i,".csv",sep=''),sep=",",row.names=FALSE,col.names=hea)}
        if(type=="default"){
          write.table(dataL[[i]],paste(path,name,"_s",i,".dat",sep=''),sep='  ',row.names=FALSE,col.names=hea)}
      }
    }
  }
}