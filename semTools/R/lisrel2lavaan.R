##lisrel2lavaan
##Corbin Quick
##02/12/13
##file path/name of LS8 LISREL syntax file

lisrel2lavaan <- function(filename=NULL, analyze=TRUE, silent=FALSE, ...){

## if filename == null, prompt user with file browser
  
  if(is.null(filename)){
    reverseSlash <- function (x, pat = "\\", rep = "/") {
      x <- gsub(pat, rep, x, fixed = T)
      x <- gsub("'", "", x, fixed = T)
      x <- gsub('"', "", x, fixed = T)
      paste(x, collapse = " ")
    }
    filename <- reverseSlash(file.choose())
  }
  
## if a file path is included in 'filename', set the working directory
## to that path so that data files will be searched for in the same 
## directory as the syntax file regardless of the current directory.
## working directory is restored at the end of the function. 
  temp <- unlist(strsplit(filename,'/',fixed=T))
  restore.wd <- getwd()
  if(length(temp)>1){
    path <- paste(temp[1:(length(temp)-1)],"/",sep='',collapse="")
    filename <- temp[length(temp)]
    setwd(path)
  }
  
lisrel<-function(filename, analyze, ...){
  
## "find" function for manipulating syntax  
  
  find <- function(pat = 0, sou, n = 1) {
  	flag <- function(vec, pat){
      vec <- unlist(vec)
      if(is.null(vec[1])){
          FALSE
        }else{
        if(is.na(vec[1])){
          FALSE
        } else {
  	  	  if(vec[1]==pat){TRUE}else{FALSE}
        }
      }
  	}
  	if (is.data.frame(sou) | is.matrix(sou)) {
      out <- 1:nrow(sou)
  		out <- out[unlist(apply(sou, 1, flag, pat = pat))]
  		out <- out[length(out)]
  	} else if (is.list(sou) | is.vector(sou)) {
      if(is.vector(sou)){
        sou <- as.list(sou)
      }
      out <- 1:length(sou)
  		out <- out[unlist(lapply(sou, flag, pat = pat))]
  		if(n!=0) {
  			out <- out[n]
  		} else {
  			out <- out[length(out)]
  		}
  	} else {
  		out <- NULL
  	}
  	if(!is.null(out)){
  	  out <- out[!is.na(out)]
  	  if(length(out)<1){
  	    out <- NULL
  	  } else if(is.na(out)){
  	    out <- NULL
  	  }
  	}
  	out
  }

as.numeric.s <- function(x){
  suppressWarnings(as.numeric(x))
}

## function to evaluate MO matrix commands; creates pseudo-class for matrices

modMat <- function(name, line) {
  ## obtain row/col numbers using ref table (external)
		row <- eval(parse(text=paste(ref[find(name,ref),2])))
		col <- eval(parse(text=paste(ref[find(name,ref),3])))
	## constraint and misc are blank by default
		constraint <- matrix(0, row, col)
		misc <- matrix("", row, col)
	## if mode specified then obtain mode, else mode='de' (default)
		if(length(unlist(strsplit(line,",")))>1){
			form <- unlist(strsplit(line,","))[1]
			mode <- unlist(strsplit(line,","))[2]
		} else {
      if(any(line==c("fi","fr"))){
			  mode <- line
        if(any(name==c("lx","ly","ga"))){
          form <- "fu"
        }else if(any(name==c("ps","te","td"))){
          form <- "di"
        }else if(any(name==c("be","th"))){
          form <- "ze"
        }else if(any(name==c("ph"))){
          form <- "sy"
        }else {
          form <- "fu"
        }
      }else{
        form <- line
			  mode <- "de"
      }
		}
	## determine matrix type (properties differ)
	if(any(name== c("lx","ly") )){
		if(any(form==c("fu","ze"))){
			if(mode=="fr"){
				start <- matrix(NA, row, col)
				free <- matrix(1, row, col)
			}else{
				start <- matrix(0, row, col)
				free <- matrix(0, row, col)
			}
		}else if(any(form==c("sd","sy","st","iz","zi"))){
			if(mode=="fr"){
				start <- matrix(NA, row, col)
				free <- matrix(1, row, col)
			}else{
				start <- matrix(0, row, col)
				free <- matrix(0, row, col)
			}
		}else if(form=="di"){
			if(name=="ly"){
				if(ny==ne){
					if(mode=="fr"){
						start <- as.matrix(diag(NA, row))
						free <- as.matrix(diag(1, row))
					}else{
						start <- as.matrix(diag(1, row))
						free <- as.matrix(diag(0, row))
					}
				}else {
					stop("syntax error: LY matrix cannot be form DI when NY is not equal to NE")
				}
			}
			if(name=="lx"){
				if(nx==nk){
					if(mode=="fr"){
						start <- as.matrix(diag(NA, row))
						free <- as.matrix(diag(1, row))
					}else {
						start <- as.matrix(diag(1, row))
						free <- as.matrix(diag(0, row))
					}
				}else {
					stop("syntax error: LX matrix cannot be form DI when NX is not equal to NK")
				}
			}
		}else if(any(form==c("id"))){
			start <- matrix(0, row, col)
      diag(start) <- 1
			free <- matrix(0, row, col)
		}
	}else if(name=="ga") {
		if(form=="fu"){
			if(mode=="fr" | mode=="de"){
				start <- matrix(NA, row, col)
				free <- matrix(1, row, col)
			} else {
				start <- matrix(0, row, col)
				free <- matrix(0, row, col)
			}
		}else if(form=="ze"){
			if(mode=="fr"){
				start <- matrix(NA, row, col)
				free <- matrix(1, row, col)
			} else {
				start <- matrix(0, row, col)
				free <- matrix(0, row, col)
			}
		}else if(any(form==c("sd","sy","st","iz","zi"))){
			if(mode=="fr"){
				start <- matrix(NA, row, col)
				free <- matrix(1, row, col)
			} else {
				start <- matrix(0, row, col)
				free <- matrix(0, row, col)
			}
		}
		if(form=="di"){
			if(ny==nx){
				if(mode=="fr"){
					start <- as.matrix(diag(NA, row))
					free <- as.matrix(diag(1, row))
				} else {
					start <- as.matrix(diag(1, row))
					free <- as.matrix(diag(0, row))
				}
			} else {
				stop("syntax error: GA matrix cannot be form DI when NY is not equal to NX")
			}
		}
		if(form=="id"){
			start <- matrix(0, row, col)
			free <- matrix(0, row, col)
		}
	}else if(name=="be") {
		if(any(form==c("fu","ze"))){
			if(mode=="fr"){
				start <- matrix(NA, row, col)
				free <- matrix(1, row, col)
			}else {
				start <- matrix(0, row, col)
				free <- matrix(0, row, col)
			}
		} else if(form=="sy"){
			if(mode=="fr"){
				start <- matrix(NA, row, col)
				free <- matrix(1, row, col)
			}else if(mode=="de") {
				start <- as.matrix(diag(NA, row))
				free <- as.matrix(diag(1, row))
			}else {
				start <- matrix(0, row, col)
				free <- matrix(0, row, col)
			}
		}else if(any(form==c("sd","st","iz","zi","id", "di"))){
			if(mode=="fi"){
				start <- as.matrix(diag(1, row))
				free <- as.matrix(diag(0, row))
			}else {
				start <- as.matrix(diag(NA, row))
				free <- as.matrix(diag(1, row))
			}
		}
	}else if(any(name==c("td", "te", "th", "ph", "ps"))) {
		if(any(form==c("fu","ze"))){
			if(mode=="fr"){
				start <- matrix(NA, row, col)
				free <- matrix(1, row, col)
			}else {
				start <- matrix(0, row, col)
				free <- matrix(0, row, col)
			}
		}else if(form=="sy"){
			if(mode=="fr"){
				start <- matrix(NA, row, col)
				free <- matrix(1, row, col)
			}else if(mode=="de") {
				start <- as.matrix(diag(NA, row))
				free <- as.matrix(diag(1, row))
			}else {
				start <- matrix(0, row, col)
				free <- matrix(0, row, col)
			}
		}else if(any(form==c("sd","st","iz","zi","id","di"))){
			if(mode=="fi"){
				start <- as.matrix(diag(1, row))
				free <- as.matrix(diag(0, row))
			}else {
				start <- as.matrix(diag(NA, row))
				free <- as.matrix(diag(1, row))
			}
		}
	}else if(any(name==c("ty","tx","ka","kl","al"))){
		if(any(mode==c("fi","ze")) | any(form==c("fi","ze"))){
			start <- matrix(0, row, col)
			free <- matrix(0, row, col)
		}else{
			start <- matrix(NA, row, col)
			free <- matrix(1, row, col)
		}
	}
	list(start=start,free=free,constraint=constraint,misc=misc)
}


## function to format LISREL syntax
  
  doc <- scan(filename, "", sep="\n")

  format <- function(doc) {
      doc <- gsub("\t"," ",doc)
      doc <- gsub("(^ +)|( +$)", "", doc)
      doc <- gsub("\\(","[",doc)
      doc <- gsub("\\)","]",doc)
      doc <- gsub("]","] ",doc)
      doc <- gsub(" \\[","\\[",doc)
      doc <- gsub("/","",doc)
      doc <- lapply(doc, function(x){gsub("!","`!",x)})
      if(length(grep("!",doc))>0) {
        doc <- lapply(lapply(doc,strsplit,split="`",fixed=TRUE),unlist)
        doc <- lapply(doc, function(x){if(length(grep("!",x))>0){x[1:((grep("!",x))[1]-1)]}else{x}})
      }
      del <- lapply(doc,function(x){if(is.null(find("",x))){NULL}else{find("",x)+1}})
      for(i in seq_along(del)){
        if(!is.null(del[[i]])){
          doc[[i]] <- doc[[i]][doc[[i]]!=doc[[i]][del[[i]]]]
        }
      }
      doc<-unlist(doc)
      doc<-lapply(doc,gsub,pattern="\t",replacement="")
      doc<-lapply(doc,gsub,pattern="(^ +)|( +$)",replacement="")
      doc<-lapply(lapply(doc,strsplit,split=" ",fixed=TRUE),unlist)
      doc<-doc[doc!=""]
      doc<-lapply(doc,function(x){x[x!=""]})
      doc<-doc[!unlist(lapply(doc,is.null))]
      doc<-doc[unlist(lapply(doc,function(x){if(length(x)==0){FALSE}else{TRUE}}))]
      
    doc
  }
  
  doc0 <- format(doc)
  doc <- format(tolower(doc))
  
## OU output commands ...
  if(!is.null(find("ou",doc))){
    ou <- unlist(doc[[find("ou",doc)]])
    ou <- ou[ou!="ou"]
  }else{
    ou <- NULL
  }
  if(length(grep("me",ou))>0){
    estimator <- unlist(strsplit(ou[grep("me",ou)],"="))[2]
    if(estimator=="gl"){
      estimator <- "GLS"
    }else if(estimator=="wl"){
      estimator <- "WLS"
    }else if(estimator=="ul"){
      estimator <- "ULS"
    }else if(estimator=="dw"){
      estimator <- "DWLS"
    }
  }else{
    estimator <- "default"
  }
#   if(length(grep("se",ou))>0){
#     me <- 
#   }else{
#     me <- "default"
#   }
  
## Multiple-Group Models

  groupN <- 1
  da <- doc[[find("da",doc,1)]]
  da <- t(as.data.frame(strsplit(da[2:length(da)],"=")))
  if(!is.null(find("ng",da))){
    ng <- as.numeric.s(da[find("ng",(da)),2])
    if(ng>1){
      for(i in 2:ng){
        if(i==ng){
          tx <- ")):length(doc)]"
        }else{
          tx <- paste(")):(find('da',doc,",(i+1),")-1)]",sep="")
        }
        eval(parse(text=paste("doc",i,"<-doc[(find('da',doc,",i,tx,sep="")))
        eval(parse(text=paste("doc0",i,"<-doc0[(find('da',doc,",i,tx,sep="")))
      }
      doc0 <- doc0[1:(find("da",doc,2)-1)]
      doc <- doc[1:(find("da",doc,2)-1)]
    }
  }else{
    ng <- 1
  }
  
## FUNCTION TO EXTRACT DATA

  ## get # variables
  ## must be global environment
  ni <- doc[[find("da",doc)]][[grep("ni",doc[[find("da",doc)]])]]
  ni <- as.numeric.s(gsub("ni=","",ni))

  ## the 'makeSym' function is primarily used in 'getData';
  ## however, it has uses elsewhere (e.g. PA commands), and
  ## therefore must be left out of 'getData' itself.
  makeSym <- function(dat, ni){
    dat <- unlist(dat)
    lapply(1:ni,function(x,dat){
      if(x==1){
        return(dat[1])
      }else{
        dat[(sum(1:(x-1))+1):sum(1:x)]
      }
    },dat=dat)
  }
  
getData <- function(doc, doc0, ngroup = 1) {
  
## below is an unfortunate work-around ..
  if(length(grep("cm=", doc))>0){
    doc[[grep("cm=", doc)]] <- unlist(strsplit(unlist(doc[grep("cm=", doc)]),"="))
    doc0[[grep("cm=", doc0, ignore.case=T)]] <- unlist(strsplit(unlist(doc0[grep("cm=", doc0, ignore.case=T)]),"="))
  }
  if(length(grep("km=", doc))>0){
    doc[[grep("km=", doc)]] <- unlist(strsplit(unlist(doc[grep("km=", doc)]),"="))
    doc0[[grep("km=", doc0, ignore.case=T)]] <- unlist(strsplit(unlist(doc0[grep("km=", doc0, ignore.case=T)]),"="))
  }
  if(length(grep("me=", doc))>0){
    doc[[grep("me=", doc)]] <- unlist(strsplit(unlist(doc[grep("me=", doc)]),"="))
    doc0[[grep("me=", doc0, ignore.case=T)]] <- unlist(strsplit(unlist(doc0[grep("me=", doc0, ignore.case=T)]),"="))
  }
  if(length(grep("pm=", doc))>0){
    doc[[grep("pm=", doc)]] <- unlist(strsplit(unlist(doc[grep("pm=", doc)]),"="))
    doc0[[grep("pm=", doc0, ignore.case=T)]] <- unlist(strsplit(unlist(doc0[grep("pm=", doc0, ignore.case=T)]),"="))
  }
  if(length(grep("sd=", doc))>0){
    doc[[grep("sd=", doc)]] <- unlist(strsplit(unlist(doc[grep("sd=", doc)]),"="))
    doc0[[grep("sd=", doc0, ignore.case=T)]] <- unlist(strsplit(unlist(doc0[grep("sd=", doc0, ignore.case=T)]),"="))
  }
  if(length(grep("ra=", doc))>0){
    doc[[grep("ra=", doc)]] <- unlist(strsplit(unlist(doc[grep("ra=", doc)]),"="))
    doc0[[grep("ra=", doc0, ignore.case=T)]] <- unlist(strsplit(unlist(doc0[grep("ra=", doc0, ignore.case=T)]),"="))
  }
  
##paragraphs of interest ...
  paragraphs <- c("cm","km","me","pm","sd","ra")
  pValues <- 1:length(paragraphs)
  
  cm <- NULL
  km <- NULL
  me <- NULL
  pm <- NULL
  sd <- NULL
  ra <- NULL
  
  fLength <- c(((ni^2+ni)/2), ((ni^2+ni)/2), ni, ni, ni, NA)
  fExist <- rep(FALSE, length(paragraphs))
  fLocate <- rep(0, length(paragraphs))
  fNames <- rep(NA, length(paragraphs))
  fData <- list()
  
  dataList <- list()
  dN <- 1
  existList <- list()
  orderList <- list()
  
  charTest <-function(line){
    line <- unlist(line)
    if(is.na(as.numeric.s(line[1])))
      {TRUE} else {FALSE}
  }
  
  ## read in data for individual paragraphs
  for(i in seq_along(paragraphs)){
    p <- paragraphs[[i]]
    if(!is.null(find(p,doc))){
      line1a <- length(unlist(doc[find(p,doc)]))>1
      line1b <- if(line1a){charTest(doc[[find(p,doc)]][2])}else{FALSE}
      line2 <- charTest(unlist(doc[find(p,doc)+1]))
      if(line1a && line1b && line2){
        fExist[i] <- TRUE
        if(length(grep("=",doc0[[find(p,doc)]][2]))>0){
          if(doc0[[find(p,doc)]][2]=="="){
            fname <- doc0[[find(p,doc)]][3]
          }else{
            fname <- unlist(strsplit(doc0[[find(p,doc)]][2],"="))[2]
          }
        }else{
          fname <- doc0[[find(p,doc)]][2]
        }
        if(length(find(fname, fNames))==0){
          if(paragraphs[i]=="ra"){
              type <- tolower(substr(fname, (nchar(fname)-2), nchar(fname)))
              if(type=='dat'){
                if(is.numeric(as.matrix(read.table(fname, nrows=1)))==FALSE){
                  dataList[[dN]] <- as.matrix(read.table(fname,header=TRUE))
                }
                else{dataList[[dN]] <- as.matrix(read.table(fname))}
              } else if(type=='csv'){
                if(is.numeric(as.matrix(read.table(fname, nrows=1)))==FALSE){
                  dataList[[dN]] <- as.matrix(read.csv(fname,header=TRUE))
                }else{dataList[[dN]] <- as.matrix(read.csv(fname))}
              } else if (type=='psf'){
                stop("Please use a different data format: .PSF files are compatible only with PRELIS.")
              } else {
                if(is.numeric(as.matrix(read.table(fname, nrows=1)))==FALSE){
                  dataList[[dN]] <- as.matrix(read.table(fname,header=TRUE))
                }
                else{dataList[[dN]] <- as.matrix(read.table(fname))}
              }
              
          }else{
            dataList[[dN]] <- unlist(format(scan(fname,"",sep="\n")))
          }
          fLocate[i] <- dN
          existList[[dN]] <- rep(FALSE, length(paragraphs))
          existList[[dN]][i] <- TRUE
          dN <- dN + 1
        } else {
          existList[[(dN-1)]][i] <- TRUE
          fLocate[i] <- fLocate[find(fname, fNames)]
        }
        fNames[[i]] <- fname
      }
    }
  }
  ## determine order: which paragraphs are found in data files first?
  if(length(dataList)>0){
    for(x in 1:length(dataList)){
      tempFrame <- matrix(NA,2,length(pValues[existList[[x]]]))
      tempFrame[1,] <- pValues[existList[[x]]]
      for(i in pValues[existList[[x]]]){
        tempFrame[2,tempFrame[1,]==i] <- find(paragraphs[i], doc)
      }
      orderList[[x]] <- tempFrame[1,order(tempFrame[2,])]
    }
  ## assign appropriate data to paragraph list 
    for(x in 1:length(dataList)){
      for(i in orderList[[x]]){
            ## TEST: IS FULL (NON-SYMMETRIC) MATRIX??
          if(paragraphs[[i]]=="ra"){
            fData[[i]] <- dataList[[x]]
          }else{
            if(fLength[i]==((ni^2+ni)/2)){
              if(dataList[[x]][2]==dataList[[x]][(ni+1)] && dataList[[x]][3]==dataList[[x]][(2*ni+1)]){
                fLength[i] <- ni^2
                fData[[i]] <- dataList[[x]][1:fLength[i]]
              }else{
                fData[[i]] <- makeSym(dataList[[x]][1:fLength[i]], ni=ni)
              }
            } else {
              fData[[i]] <- dataList[[x]][1:fLength[i]]
            }
            
            dataList[[x]] <- dataList[[x]][(fLength[i]+1):length(dataList[[x]])]
          }
      }
    }
  }
    
  excerpt <- function(para, doc){
    ## determine whether or not paragraph is specified
    if(!is.null(find(para, doc))){
      if(fExist[[grep(para, paragraphs)]]){
        fData[[grep(para, paragraphs)]]
      } else {
        out <- find(para, doc):length(doc)
        out <- out[unlist(lapply(doc0[find(para, doc):length(doc)], charTest))]
        doc[(find(para, doc)+1):(out[2]-1)]
      }
    } else {
      return(NULL)
    }
  }
  makeMatrix <- function(x) {
    if(is.null(x)){
      NULL
    } else {
      if(length(x)>1 && length(x[[1]])!=length(x[[2]])){
        for(i in 1:(length(x)-1)){
          d <- unlist(lapply(x[(i+1):length(x)],function(z,i){z[i]}, i=i))
          x[[i]] <- c(x[[i]],d)
        }
        do.call(rbind,lapply(x, as.numeric.s))
      } else { 
        if(!is.matrix(x) && !is.data.frame(x)){
          sapply(x, as.numeric.s,simplify="vector")
        } else{
          apply(x, 2, as.numeric.s)
        }  
      }
    }
  }
  for(i in paragraphs){
    if(i=="me"|i=="sd"){
      assign(i,makeMatrix(unlist(excerpt(i, doc))))
    }else{
      assign(i,makeMatrix(excerpt(i, doc)))
    }
  }
  if(!is.null(cm)){
    if(length(var)>ncol(cm)){
      var <- var[1:ncol(cm)]
    }
  }
  if(!is.null(km)){
    if(length(var)>ncol(km)){
      var <- var[1:ncol(km)]
    }
  }
  output <- list(cm=cm,km=km,me=me,pm=pm,sd=sd,ra=ra)
  rows <- list(var,var,NULL,var,NULL,NULL)
  for(i in paragraphs[!sapply(output,is.null)]){
    if(i=="ra"){
      if(length(var)<ncol(output$ra)){
        output$ra <- output$ra[,1:length(var)]
      }
    }
    if(is.null(rows[paragraphs==i][[1]]) && i!="ra"){
      eval(parse(text=paste("output$",i,"<-matrix(output$",i,",nrow=1)")))
    }
    eval(parse(text=paste("colnames(output$",i,")<-var",sep="")))
    eval(parse(text=paste("rownames(output$",i,")<-unlist(rows[paragraphs==i])")))
  }
  output
}

## Model Defaults
  
  ne <- NULL
  nk <- NULL
  nx <- NULL
  ny <- NULL
    
  al <- NULL
  be <- NULL
  ga <- NULL
  ka <- NULL
  lx <- NULL
  ly <- NULL
  ph <- NULL
  ps <- NULL
  td <- NULL
  te <- NULL
  th <- NULL
  tx <- NULL
  ty <- NULL

  key.ref <- c("al","be","ga","ka","lx","ly","td","te","ph","ps","th","tx","ty")
  row.ref <- c("01","ne","nk","01","nx","ny","nx","ny","nk","ne","ne","01","01")
  col.ref <- c("ne","ne","ne","nk","nk","ne","nx","ny","nk","ne","nk","nx","ny")
  rel.ref <- c("~1","=~","~","~1","=~","=~","~~","~~","~~","~~","~~","~1","~1")

  ref <- data.frame(key.ref,row.ref,col.ref,rel.ref)
  
## establish basic model components: matrices, X and Y variables

  mo <- doc[[find("mo",doc,1)]]
  mo <- t(as.data.frame(strsplit(mo[2:length(mo)],"=")))
  ## determine number latent & manifest variables, remove from mo
  if(length(find("ny",mo))!=0){
    ny <- as.numeric.s(mo[find("ny",mo),2])
    mo <- mo[-find("ny",mo),]
  }
  if(length(find("ne",mo))!=0){
    ne <- as.numeric.s(mo[find("ne",mo),2])
    mo <- mo[-find("ne",mo),]
  }
  if(length(find("nx",mo))!=0){
    nx <- as.numeric.s(mo[find("nx",mo),2])
    mo <- mo[-find("nx",mo),]
  }
  if(length(find("nk",mo))!=0){
    nk <- as.numeric.s(mo[find("nk",mo),2])
    mo <- mo[-find("nk",mo),]
  }
  if(length(find("me",doc))!=0){
    macs<-TRUE
    if(is.null(ny)==FALSE){
      if(length(find("ty",mo))==0){
        mo <- rbind(mo,c("ty","fr"))
      }
      if(length(find("al",mo))==0){
        mo <- rbind(mo,c("al","fi"))
      }
    }
  } else {
    macs <- FALSE
  }
  
## determine variable numbers, names, x / y assignment 
  
  ##extrapNames function for sequences such as 'variable3 - variable33'
    extrapNames <- function(x){
      if(length(grep("-",x))>0){
        w <- grep("-",x)
        if(nchar(x[w])>1){
          x <- gsub("-","`-`",x)
          x <- unlist(strsplit(x, "`"))
          w <- grep("-",x)
        }
        start <- x[(w-1)]
        end <- x[(w+1)]
        startV <- as.numeric(gsub("[^0-9]","",start))
        endV <- as.numeric(gsub("[^0-9]","",end))
        name <- gsub("[0-9]","",start)
        out <- paste(name, startV:endV, sep="")
        if(length(x)>(w+2)){
          if(w==2){
            c(out,x[(w+2):length(x)])
          }else{
            c(x[1:(w-2)],out,x[(w+2):length(x)])
          }
        }else{
          if(w==2){
            out
          }else{
            c(x[1:(w-2)],out)
          }
        }
      }else{
        x
      }
    }
  
  ##pullNames function to simplify obtaining names
    pullNames <- function(x){
      y <- find(x, doc)
      if(is.null(y)){
        NULL
      } else {
        coms <- c("mo","km","cm","se","la","lk","le","ou","pd","ra","fr","va","fi","eq","co")
        coms <- coms[coms!=x]
        is.com <- function(line){any(unlist(line)[1]==coms)}
        com.l <- c(1:length(doc))[sapply(doc,is.com)]
        names <- unlist(doc0[(y+1):(c(com.l[com.l>y])[1]-1)])
        extrapNames(names[names!=""])
      }
    }
    
    use <- pullNames("se")
    var <- pullNames("la")
    
    if(is.null(use)){use <- var}
    
    use.t <- use
    
    if(!is.null(var)){
      name.def <- FALSE
    } else {
      name.def <- TRUE
    }
    if(!is.null(ny)){
      NY <- use.t[1:ny]
      use.t <- use.t[(ny+1):length(use.t)]
    }else{
      NY <- NULL
    }
    if(!is.null(nx)){
      NX <- use.t[1:nx]
    }else{
      NX <- NULL
    }
    
    if(name.def){
      ## names not specified
      if(is.numeric(nx)){
       NX <- paste("ksi",1:nx,sep="")
        if(is.null(ny)){
          NY <- NX
          var <- NX
        }
      }
      if(is.numeric(ny)){
       NY <- paste("eta",1:ny,sep="")
        if(is.null(nx)){
          NX <- NY
          var <- NY
        }
      }  
      use <- var
    }
      
    NK <- pullNames("lk")
    NE <- pullNames("le")
  
  ## for path analysis models...
  
  if(!is.null(nx)){
	if(nx>length(NX)){
		NX <- paste("ksi", 1:nx, sep="")
	}
  }
  if(!is.null(ny)){
	  if(ny>length(NY)){
		NY <- paste("eta", 1:ny, sep="")
	  }
  }
  if(is.null(NK)){
    NK<-NX
  }
  if(!is.null(nk)){
	  if(nk>length(NK)){
		NK <- paste("KSI", 1:nk, sep="")
	  }
  }
  if(is.null(NE)){
    NE<-NY
  }
  if(!is.null(ne)){
	  if(ne>length(NE)){
		NE <- paste("ETA", 1:ne, sep="")
	  }
  }

  if(is.null(nk)){nk<-nx}
  if(is.null(ne)){ne<-ny}
  
## generate model matrices
  
  for(i in 1:nrow(mo)){
    assign(mo[i,1],modMat(mo[i,1], mo[i,2]))
  }
  
  if((!is.null(ph) && !is.null(td))|(length(grep("lx",doc))>0)){
    if(is.null(find("lx",mo))){
      mo <- rbind(mo,c("lx","fu,fi"))
      lx <- modMat("lx", "fu,fi")
    }
  }
  if((!is.null(ps) && !is.null(te))|(length(grep("ly",doc))>0)){
    if(is.null(find("ly",mo))){
      mo <- rbind(mo,c("ly","fu,fi"))
      lx <- modMat("ly", "fu,fi")
    }
  }
  
## PA paragraph commands
  
while(!is.null(find("pa",doc))){
  if(!is.null(find("pa",doc))){
    loc.n <- find("pa",doc)
    nam.n <- unlist(doc[[loc.n]])[2]
    row.n <- (eval(parse(text=paste(ref[find(nam.n,ref),2]))))
    lis.n <- doc[(loc.n+1):(loc.n+row.n)]
    if(length(lis.n[[1]])!=length(lis.n[[length(lis.n)]])){
      lis.n <- lavaan::lower2full(lavaan::char2num(paste(lapply(lis.n,paste,collapse=", "),collapse="\n")))
    }else{
      lis.n <- do.call(rbind,lapply(lis.n,as.numeric.s))
    }
    eval(parse(text=paste(nam.n,"$free<-lis.n")))
    tex.n<-paste(nam.n,"$start[(is.na(",nam.n,"$start)|",nam.n,"$start=='NA')&(",nam.n,"$free==0)]<-0")
    eval(parse(text=tex.n))
    doc[(loc.n):(loc.n+row.n)] <- NULL
    doc <- doc[!is.null(doc)]
  }
}

## MA paragraph commands
  
while(!is.null(find("ma",doc))){
  if(!is.null(find("ma",doc))){
    if(length(grep("fi",doc[[find("ma",doc)]]))>0){
      loc.n <- find("ma",doc)
      nam.n <- unlist(doc[[loc.n]])[2]
      file <- unlist(strsplit(doc[[find("ma",doc)]],'='))
      file <- file[length(file)]
      lis.n <- read.table(file)
      if(length(grep('D',unlist(lis.n)))>0){
        s2n <- function(x){
          s2n.i <- function(x){
            x <- as.numeric(unlist(strsplit(x,'D')))
            x[1]*10^(x[2])
          }
          sapply(x, s2n.i)
        }
        lis.n <- apply(lis.n, 2, s2n)
      }
      if(is.list(lis.n)){
        if(length(lis.n[[1]])!=length(lis.n[[length(lis.n)]])){
          lis.n <- lavaan::lower2full(lavaan::char2num(paste(lapply(lis.n,paste,collapse=", "),collapse="\n")))
        }else{
          lis.n <- do.call(rbind,lapply(lis.n,as.numeric.s))
        }
      }
      eval(parse(text=paste(nam.n,"$start<-lis.n")))
      doc[(loc.n)] <- NULL
      doc <- doc[!is.null(doc)]
    }else{
      loc.n <- find("ma",doc)
      nam.n <- unlist(doc[[loc.n]])[2]
      row.n <- (eval(parse(text=paste(ref[find(nam.n,ref),2]))))
      lis.n <- doc[(loc.n+1):(loc.n+row.n)]
      if(length(lis.n[[1]])!=length(lis.n[[length(lis.n)]])){
        lis.n <- lavaan::lower2full(lavaan::char2num(paste(lapply(lis.n,paste,collapse=", "),collapse="\n")))
      }else{
        lis.n <- do.call(rbind,lapply(lis.n,as.numeric.s))
      }
      eval(parse(text=paste(nam.n,"$start<-lis.n")))
      doc[(loc.n):(loc.n+row.n)] <- NULL
      doc <- doc[!is.null(doc)]
    }
  }
}
  
## ensure that command lines have brackets
  
fixLazilyWrittenSyntax <- function(line){
  commands <- c("fr", "fi", "eq", "co", "va", "st", "pa")
  line <- unlist(line)
  if(any(line[1]==commands) && length(grep("\\[",line))==0){
    if(!is.na(as.numeric.s(line[2]))){
      l <- 2 } else {l <- 1}
    line[(l+1):length(line)]<-sapply(line[(l+1):length(line)],
      function(x){if(!is.na(as.numeric.s(x))){paste('[',x,']',sep='')}else{x}})
    temp <- paste(line[(l+1):length(line)],collapse="")
    temp <- strsplit(gsub("\\]","\\]:",gsub("\\]\\[",",",temp)),":")
    unlist(c(line[1:l],gsub(",,",",",unlist(temp))))
  } else { line }
}

doc<-lapply(doc,fixLazilyWrittenSyntax)


## function: process model commands
  
  eqN <- 1
  
  processCommands <- function(doc){
    
    fr1 <- doc[(find("mo",doc)+1):length(doc)]
    
    is.pertinent <- function(doc.l){
      commands <- c("fr", "fi", "eq", "co", "va", "st", "pa")
      if(any(doc.l[[1]][1]==commands)){
        return(TRUE)
      } else {
        return(FALSE)
      }
    }
    
    fr1 <- fr1[unlist(lapply(fr1,is.pertinent))]
    
    if(length(fr1)>0){
      
      comm <- c()
      commN <- 1
      eq <- list()
      co <- c()
      coN <- 1
      
      for(i in 1:length(fr1)){
        fr1[[i]] <- fr1[[i]][fr1[[i]] != ""]
        if(length(fr1[[i]])>1){
          for(z in 2:length(fr1[[i]])){
            if(length(unlist(strsplit(fr1[[i]][z],",")))>2){
              temp <- unlist(strsplit(fr1[[i]][z],","))
              fr1[[i]][z] <- paste(substr(temp[1],1,3),temp[2],",",temp[3],sep="")
            }
          }
        }
        if(fr1[[i]][1]=="fr"){
          comm[[commN]] <- c(1,NA,fr1[[i]][2:length(fr1[[i]])])
          commN <- commN+1
        }
        if(any(fr1[[i]][1]==c("va","st"))){
          comm[[commN]] <- c("X",fr1[[i]][2:length(fr1[[i]])])
          commN <- commN+1
        }
        if(fr1[[i]][1]=="fi"){
          comm[[commN]] <- c(0,0,fr1[[i]][2:length(fr1[[i]])])
          commN <- commN+1
        }
        if(fr1[[i]][1]=="eq"){
          eq[[eqN]] <- fr1[[i]][2:length(fr1[[i]])]
          eqN <- eqN+1
        }
        if(fr1[[i]][1]=="co"){
          tempc <- fr1[[i]][2:length(fr1[[i]])]
          if(length(tempc>1)){
            tempc <- paste(tempc,sep="",collapse="")
            tempc <- strsplit(tempc,"=")
          }
          co[[coN]] <- tempc[[1]]
          coN <- coN+1
        }
        if(fr1[[i]][1]=="pa"){
          a <- (i+1)
          b <- eval(parse(text=paste(ref[find(fr1[[i]][2],ref),2]))) + a - 1
          c <- 1
          matr <- list()
          for(d in a:b){
            matr[[c]] <- fr1[[d]]
            c <- c + 1
          }
          matr <- lapply(matr, as.numeric.s)
          matr <- do.call(rbind,matr)
          eval(parse(text=paste(fr1[[i]][2],"<-matr")))
        }
        
      }
      simpl <- function(x) {
        x <- lapply(x,unlist)
        x <- x[unlist(lapply(x, length) != 0)] 
        x <- x[x != ""]
      }
      comm <- simpl(comm)
      eq   <- simpl(eq)
      co   <- simpl(co)
      
      outp <- list(comm,eq,co)
      
      return(list(outp,eqN))
      
    } else {
      return(list(list(NULL,NULL,NULL,NULL,NULL),eqN))
    }
  }
  
  ## function: apply model commands to matrices
  
  eqID <-  1
  
  applyCommands <- function(commList){
    
    comm <- commList[[1]]
    eq   <- commList[[2]]
    co   <- commList[[3]]
    
    outList <- list()
    t0 <- 1
    
    ## apply fixed parameter values
    
    if(length(comm)>0){
      for(i in 1:length(comm)){
        for(z in 3:length(comm[[i]])){
          if(comm[[i]][1]=="X"){
            outList[[t0]] <- paste(gsub("\\[","$start[",comm[[i]][z]),"<-",comm[[i]][2],sep="")
            t0 <- t0 + 1
          }else{
            outList[[t0]] <- paste(gsub("\\[","$free[",comm[[i]][z]),"<-",comm[[i]][1],sep="")
            outList[[t0+1]] <- paste(gsub("\\[","$start[",comm[[i]][z]),"<-",comm[[i]][2],sep="")
            t0 <- t0 + 2
          }
        }
      }
    }
    
    ## apply equality constraints to matrices 
    
    if(length(eq)>0){
      for(i in 1:length(eq)){
        for(z in 1:length(eq[[i]])){
          qtest <- nrow(eval(parse(text=gsub("\\[","$start[",eq[[i]][z]))))
          if(is.null(qtest)){
            qtest <- 1
          }
          if(qtest > 1){
            subN <- gsub("\\D", "", eq[[i]][z])
            eq[[i]][z] <- gsub(subN, paste(subN,",",subN), eq[[i]][z])
          }
          if(z==1){
            outList[[t0]] <- paste(gsub("\\[","$constraint[",eq[[i]][z]),"<-",eqID,sep="")
            outList[[(t0+1)]] <- paste(gsub("\\[","$misc[",eq[[i]][z]),"<-",1,sep="")
            t0 <- t0 + 2
          }else{
            outList[[t0]] <- paste(gsub("\\[","$constraint[",eq[[i]][z]),"<-",eqID,sep="")
            t0 <- t0 + 1
          }
        }
        eqID <- eqID + 1
      }
    }
    
    ## apply CO parameter constraints to matrices 
    
    if(length(co)>0){
      for(i in 1:length(co)){
        temp <- gsub('\\+','?',co[[i]][2])
        temp <- gsub('\\-','?',temp)
        temp <- gsub('\\*','?',temp)
        temp <- gsub('/','?',temp)
        temp <- strsplit(temp,'\\?')
        temp <- unlist(c(co[[i]][1],temp))
        is.mat <- function(input){
          if(nchar(input)<4){
            return(FALSE)
          } else{
            return(TRUE)
          }
        }
        temp <- temp[unlist(lapply(temp,is.mat))]
      }
    }
    
    outList
    
  }

## execute processCommands, applyCommands
  
  commList1 <- processCommands(doc)
	eqN<-commList1[[2]]
	commList1<-commList1[[1]]
  
  commands1 <- applyCommands(commList1)
  
  if(length(commands1)>0){
    for(i in 1:length(commands1)){
      eval(parse(text=commands1[[i]]))
    }
  }
  
## matrix-to-parameter-table function

  toPara <- function(name)
  {
    ob <- eval(parse(text=name))

    if(nrow(ob$free)!=0 && ncol(ob$free)!=0){
      ID <- find(name, ref)

      ROW <- eval(parse(text=toupper(ref[ID,2])))
      COL <- eval(parse(text=toupper(ref[ID,3])))
      if(length(ROW)==1){
        if(ROW==1){
          ROW <- rep("",length(COL))
        }
      }

      if(paste(ref[ID,2])==paste(ref[ID,3]) && name!="be"){
        ob$start[upper.tri(ob$start)] <- 0
      }
      
      if(name=="ps"){
        if(is.null(be) && all(ob$free==diag(1,nrow(ob$free),ncol(ob$free)))){
          diag(ob$start) <- diag(ob$start) + 99
          ob$start[is.na(ob$start)] <- "NA"
          ob$start[ob$start[lower.tri(ob$start)]!="0"]<-(as.numeric.s(ob$start[ob$start[lower.tri(ob$start)]!="0"])+99)
        } else { 
          ob$start[lower.tri(ob$start,diag=T)] <- ob$start[lower.tri(ob$start,diag=T)]+99
          if(!is.null(be)){
            test <- (be$free+be$start+t(be$free)+t(be$start))[lower.tri(ob$start,diag=T)]
            test[is.na(test)] <- 100
            tmpPS <- ob$start[lower.tri(ob$start,diag=T)] 
            tmpPS[test!=0] <- tmpPS[test!=0] - 99
            ob$start[lower.tri(ob$start,diag=T)] <- tmpPS
          }
        }
      }
    
      if(name=="al" | name=="ka" | name=="ty"){
        ob$start <- ob$start + 99
        ob$start[is.na(ob$start)]<-"NA"
        ob$start <-  t(as.matrix(as.character(ob$start)))
      }else {
        ob$start <- apply(ob$start,2,function(x){
          x[is.na(x)]<-"NA"
          as.matrix(as.character(x))
        })
      }

      OP <- paste(ref[ID,4])

      lhs <- c()
      op  <- c()
      rhs <- c()
      user <- c()
      group <- c()
      free <- c()
      ustart <- c()
      exo <- c()
      label <- c()
      eq.id <- c()
      unco <- c()
    
      ob <- lapply(ob, as.matrix)
      
      correctPosition <- function(x){
        if(length(ROW)!=nrow(x) && length(ROW)==ncol(x)){
          t(x)
        }else{
          x
        }
      }
      
      

      ob <- lapply(ob, correctPosition)
      
      if(any(name==c("al","ka","tx","ty"))){
	ROW <- COL
	COL <- rep("", length(ROW))
      }

      for(i in 1:ncol(ob$start)){
        non <- 1:nrow(ob$start)
        non <- non[unlist(ob$start[,i]!=0)]
        if(length(non)!=0){
          for(z in non){
            lhs <- c(lhs, COL[i])
            op  <- c(op, OP)
            rhs <- c(rhs, paste(ROW[z]))
            user <- c(user, 1)
            group <- c(group, groupN)
            free <- c(free, ob$free[z,i])
            ustart <- c(ustart, ob$start[z,i])
            exo <- c(exo, 0)
            label <- c(label, paste(name,"_",z,"_",i,sep=""))
            eq.id <- c(eq.id, ob$constraint[z,i])
            unco <- c(unco, ob$free[z,i])
          }
        }
      }
      if(name=="al" | name=="ka" | name=="ty"){
        ustart<-as.character(as.numeric.s(ustart)-99)
        ustart[is.na(ustart)]<-"NA"
      }
      if(name=="ps"){
        ustart<-as.character(as.numeric.s(ustart)-99)
        ustart[is.na(ustart)]<-"NA"
      }
      if(any(name==c("al","ka","tx","ty"))){
	lhs <- rhs
	rhs <- rep("", length(lhs))
      }
      data.frame(lhs,op,rhs,user,group,free,ustart,exo,label,eq.id,unco)
      
    }else{
      NULL
    }
      
  }
  
## CHECK FOR EQUALITY CONSTRAINTS BEFORE GROUP 1 PARAMETER TABLE 
  
if(ng>1){
  
	moEQ <- mo
	moEQ[,2] <- 0
	
    for(gN in 2:ng){

      eval(parse(text=paste("docN<-doc",gN,sep="")))
      eval(parse(text=paste("docN0<-doc0",gN,sep="")))

      moN <- docN[[find("mo",docN,1)]]
      moN <- t(as.data.frame(strsplit(moN[2:length(moN)],"=")))
      	  
  	  ## check each matrix for constraints
  		for(i in 1:nrow(moEQ)){
  			if(!is.null(find(moEQ[i,1],moN))){
  				if(moN[find(moEQ[i,1],moN),2]=="in"){
  					moEQ[i,2]<-1
  				}
  			}
  		}
    }
     
  multi.grp.eq <- function(x){
    if(x[1]=="eq"){
      if(any(sapply(x,function(y)
        {if(length(unlist(strsplit(y,",")))>1){TRUE}else{FALSE}}
        ))){
          TRUE
        }else{
          FALSE
        }
    } else {
      FALSE
    }
  }
  
     eq <- docN[sapply(docN,multi.grp.eq)]
    
  if(length(eq)>0){
	t0G1 <- 1
	t0GN <- 1
	listG1 <- list()
	listGN <- list()
	if(is.list(eq)){
	  eq <- lapply(eq, function(x){x[2:length(x)]})
	}else{
	  eq <- list(eq[2:length(eq)])
	}
	for(i in seq_along(eq)){
	    for(z in 1:length(eq[[i]])){
		if(length(unlist(strsplit(eq[[i]][[z]],",")))>1){
			tmp <- unlist(strsplit(eq[[i]][[z]],","))
			wh.gr <- as.numeric(strsplit(tmp[1],'\\[')[[1]][2])
			eq[[i]][[z]] <- paste(strsplit(tmp,"\\[")[[1]][1],"[",tmp[2],",",tmp[3])
		}else{
			wh.gr <- 2
		}
         	qtest <- nrow(eval(parse(text=gsub("\\[","$start[",eq[[i]][z]))))
                if(is.null(qtest)){
            		qtest <- 1
                }
                if(qtest > 1){
                    subN <- gsub("\\D", "", eq[[i]][z])
                    eq[[i]][z] <- gsub(subN, paste(subN,",",subN), eq[[i]][z])
                }
             if(z==1){
		if(wh.gr==1){
	                listG1[[t0G1]] <- paste(gsub("\\[","$constraint[",eq[[i]][z]),"<-",eqID,sep="")
                  listG1[[(t0G1+1)]] <- paste(gsub("\\[","$misc[",eq[[i]][z]),"<-",1,sep="")
                  t0G1 <- t0G1 + 2
		} else {
		  listGN[[t0GN]] <- paste(gsub("\\[","$constraint[",eq[[i]][z]),"<-",eqID,sep="")
                  listG1[[(t0GN+1)]] <- paste(gsub("\\[","$misc[",eq[[i]][z]),"<-",1,sep="")
                  t0GN <- t0GN + 2
		}
             }else{
               if(wh.gr==1){
		  listG1[[t0G1]] <- paste(gsub("\\[","$constraint[",eq[[i]][z]),"<-",eqID,sep="")
		  t0G1 <- t0G1 + 1
		} else {
		  listGN[[t0GN]] <- paste(gsub("\\[","$constraint[",eq[[i]][z]),"<-",eqID,sep="")
		  t0GN <- t0GN + 1
		}
             }
        }
        eqID <- eqID + 1
      }
    listG1<-unlist(listG1)
    listGN<-unlist(listGN)
  ## apply these constraints to group 1
    if(length(listG1)>0){
      for(i in listG1){eval(parse(text=i))}
    }
  }else{
    listGN <- NULL
  }

  ## enter constraint requests into group 1 matrices
  
	if(sum(as.numeric.s(moEQ[,2]))>0){
	
		g1.eq <- list()
		gn.eq <- list()
	
	  moEQ <- unlist(moEQ[moEQ[,2]!=0,1])
	  liEQ <- list()
	  t0 <- 1
	  
	  for(i in seq_along(moEQ)){
		# matN is the current matrix for which multiple group
		# invariance constraints are being processed
	  
  		matN <- eval(parse(text=moEQ[i]))$start
  		matV <- 1:nrow(matN)
  		for(z in 1:ncol(matN)){
  			for(y in matV[is.na(matN[,z]) | matN[,z]!="0"]){
  				g1.eq[[t0]] <- list(paste(moEQ[i],"$constraint[",y,",",z,"]<-",eqID,sep=""))
				g1.eq[[t0]][[2]] <- paste(moEQ[i],"$misc[", y, ",", z, "]<-", 1, sep="")
  				gn.eq[[t0]] <- paste(moEQ[i],"$constraint[",y,",",z,"]<-",eqID,sep="")
  				t0 <- t0 + 1
  				eqID <- eqID + 1
  			}
		}
	  }
	  g1.eq <- unlist(g1.eq)
  } else {
  
	g1.eq <- NULL
	gn.eq <- NULL
	
  }
  
	## APPLY CONSTRAINTS TO GROUP 1 MATRICES
		## create global matrices so that constraints 
		## are not carried to other groups
	for(i in 1:nrow(mo)){
	  eval(parse(text=paste(mo[i,1],"G<-",mo[i,1],sep="")))
	}
	if(!is.null(g1.eq) && length(g1.eq)>1){
	  for(i in 1:length(g1.eq)){
	    eval(parse(text=g1.eq[i]))
	  }
	}
}

## PROCESS MATRICES TO PARAMETER TABLE, GROUP 1

  endL <- c()
  allL <- c("lx","ly","td","te","al","ka","tx","ty","be","ps","ph","th","ga")
  for(i in 1:13){
    if(is.null(find(allL[i],mo))==FALSE){
      endL <- c(endL, allL[i])
    }
  }
  
  tableList <- lapply(endL, toPara)
  parTable <- do.call(rbind,tableList)
  
## multiple group models, parameter table

  if(ng>1){
    
    for(groupN in 2:ng){
      
      eval(parse(text=paste("docN<-doc",groupN,sep="")))
      eval(parse(text=paste("docN0<-doc0",groupN,sep="")))
      
      docN <- docN[!sapply(docN,multi.grp.eq)]
      
      mo <- docN[[find("mo",docN,1)]]
      mo <- t(as.data.frame(strsplit(mo[2:length(mo)],"=")))
      
	  if(macs==TRUE){
  			if(length(find("ty",mo))==0){
  				mo <- rbind(mo,c("ty","fr"))
  			}
  			if(length(find("al",mo))==0){
  				mo <- rbind(mo,c("al","fi"))
  			}
		 }
	  
      for(i in 1:nrow(mo)){
        m.typ <- unlist(strsplit(mo[i,2],","))
        if(length(m.typ)>1){
          m.form <- m.typ[2]
          m.typ <- m.typ[1]
        }else{
          m.form <- "de"
        }
        if(m.typ=="fi"){
          ## fixed
          eval(parse(text=(paste(mo[i,1],"$start[]","<- 0"))))
          eval(parse(text=(paste(mo[i,1],"$free[]","<- 0"))))
          eval(parse(text=(paste(mo[i,1],"$misc[]","<- '' "))))
        }
        if(m.typ=="fr"){
          ## free
          eval(parse(text=(paste(mo[i,1],"$start[]","<-NA"))))
          eval(parse(text=(paste(mo[i,1],"$free[]","<-1"))))
          eval(parse(text=(paste(mo[i,1],"$misc[]","<- '' "))))
        }
        if(m.typ=="ps"){
          ## same pattern & starting values
          eval(parse(text=(paste(mo[i,1],"$misc[]","<- '' "))))
        }
        if(m.typ=="sp"){
          ## same pattern
          eval(parse(text=(paste(mo[i,1],"$start[",mo[i,1],"$start!=0]","<-NA"))))
          eval(parse(text=(paste(mo[i,1],"$misc[]","<- '' "))))
        }
        if(m.typ=="ss"){
          ## same starting values ... ?
          eval(parse(text=(paste(mo[i,1],"$free[]","<- 0"))))
          eval(parse(text=(paste(mo[i,1],"$misc[]","<- '' "))))
        }
        if(m.typ=="in"){
			    eval(parse(text=(paste(mo[i,1],"<-",mo[i,1],"G",sep=""))))
        }
        if(any(m.typ==c("fu","sy","ze"))){
          if(any(m.form==c("fi","de"))){
            eval(parse(text=(paste(mo[i,1],"$start[]","<- 0"))))
            eval(parse(text=(paste(mo[i,1],"$free[]","<- 0"))))
            eval(parse(text=(paste(mo[i,1],"$misc[]","<- '' "))))
          }else{
            eval(parse(text=(paste(mo[i,1],"$start[]","<-NA"))))
            eval(parse(text=(paste(mo[i,1],"$free[]","<-1"))))
            eval(parse(text=(paste(mo[i,1],"$misc[]","<- '' "))))
          }
        }
		
      }
      
	  if(!is.null(gn.eq) && length(gn.eq)>1){
	  	for(i in unlist(gn.eq)){
			  eval(parse(text=i))
		  }
	  }
	  
    if(length(listGN)>0){
      for(i in listGN){eval(parse(text=i))}
    }
      
      commListN <- processCommands(docN)
      
      commandsN <- applyCommands(commListN[[1]])
      
		if(length(commandsN>0)){
			for(i in 1:length(commandsN)){
				eval(parse(text=commandsN[[i]]))
			}
		}
 
      ## PROCESS MATRICES TO PARAMETER TABLE, MULTIPLE GROUPS
	  
		tableListN <- lapply(as.vector(mo[,1]), toPara)
		parTableN <- do.call(rbind,tableListN)
      
		colnames(parTableN) <- colnames(parTable)

		row.names(parTableN) <- NULL
      
    parTable <- rbind(parTable, parTableN)
		
    }
	
  }

   ## CO command constraints

  if(length(commList1[[3]])>1){
    for(i in 1:length(commList1[[3]])){
        exp <- gsub("\\[","_",commList1[[3]][[i]])
        exp <- gsub(",","_",exp)
        exp <- gsub("\\]","",exp)
        exp <- data.frame(exp[1], "==",exp[2],1,0,0,"NA",0,"",0,0)
        colnames(exp) <- colnames(parTable)
        parTable <- rbind(parTable,exp)
    }
  }

## format final parameter table
  
  if(sum(parTable[,6])!=0){
    if(sum(parTable$eq.id!=0)>0){
      j <- 1
      for(i in unique(parTable$eq.id[parTable$eq.id!=0])){
        test <- parTable$free[((parTable$eq.id==i)+(parTable$free!=0))==2]
        if(length(test)>0){
          parTable$free[((parTable$eq.id==i)+(parTable$free!=0))==2] <- j
          j <- j + 1
        }
      }
    }else{
      j <- 1
    }
    parTable[(((parTable$free==1)+(parTable$eq.id==0))==2),6] <- 
      j:(j-1+sum(((parTable$free==1)+(parTable$eq.id==0))==2))
  }

  parTable$unco[parTable$unco!=0] <- 1:length(parTable$unco[parTable$unco!=0])

  id <- 1:nrow(parTable)
  parTable <- data.frame(id,parTable)

  row.names(parTable) <- NULL
  
  parTable$id <- as.integer(as.numeric.s(parTable$id))
  parTable$lhs <- as.character(parTable$lhs)
  parTable$op <- as.character(parTable$op)
  parTable$rhs <- as.character(parTable$rhs)
  parTable$user <- as.integer(as.numeric.s(parTable$user))
  parTable$group <- as.integer(as.numeric.s(parTable$group))
  parTable$free <- as.integer(as.numeric.s(parTable$free))
  parTable$ustart <- as.numeric.s(as.character(parTable$ustart))
  parTable$exo <- as.integer(as.numeric.s(parTable$exo))
  parTable$label <- as.character(parTable$label)
  parTable$eq.id <- as.integer(as.numeric.s(as.character(parTable$eq.id)))
  parTable$unco <- as.integer(as.numeric.s(parTable$unco))
  
  if(analyze){
    for(i in 1:ng){
      if(i==1){
        data <- getData(doc, doc0)
      }else{
        data <- mapply("list",data,getData(
                         doc=eval(parse(text=paste("doc",i,sep=""))),
                         doc0=eval(parse(text=paste("doc0",i,sep="")))
                        ),SIMPLIFY=FALSE)
      }
    }
    if(ng>1){
      for(i in names(data)){
        if(is.null(data[[i]][[1]])){
          data[[i]] <- NULL
        }
      }
    }
    if(is.null(data$ra) | is.null(data$ra[[1]])){
      for(i in 1:ng){
        if(i==1){
          if(ng==1){
            n <- doc[[find("da",doc)]][[grep("no",doc[[find("da",doc)]])]]
            n <- as.numeric.s(gsub("no=","",n))
          } else {
            n2 <- doc[[find("da",doc)]][[grep("no",doc[[find("da",doc)]])]]
            n <- list()
            n[[i]] <- as.numeric.s(gsub("no=","",n2))
          }
        }else{
          n2 <- eval(parse(text=paste("doc",i,"[[find('da',doc",i,")]][[grep('no',doc",i,"[[find('da',doc",i,")]])]]",sep="")))
          n[[i]] <- as.numeric.s(gsub("no=","",n2))
        }
      }
    } else{
      n <- NULL
    }
    if(!is.null(data$sd) && is.null(data$ra)){
      cr2cv <- function(x,sd){
        na <- colnames(x)
        sd <- diag(as.vector(sd))
        ou <- sd%*%x%*%t(sd)
        colnames(ou) <- na
        rownames(ou) <- na
        ou
      }
      for(i in 1:ng){
        if(i==1){
          if(ng==1){
            data$cm <- cr2cv(data$km, data$sd)
          }else{
            data$cm <- list()
            data$cm[[i]] <- cr2cv(data$km[[i]], data$sd[[i]])
          }
        }else{
          data$cm[[i]] <- cr2cv(data$km[[i]], data$sd[[i]])
        }
      }
    }
    if(!is.null(data$ra)){
      if(ng==1){
        data$ra <- as.data.frame(data$ra)
      }else{
        data$ra <- lapply(data$ra, as.data.frame)
      }
    }
    if(is.null(data$cm) && is.null(data$ra)){
      invisible(suppressWarnings(lavaan::lavaan(model=parTable)))
      stop("lisrel2lavaan requires either 1) raw data (specified in the RA paragraph in LISREL syntax), 2) a variance-covariance matrix (the CM paragraph in LISREL syntax), or 3) a correlation matrix AND standard deviation vector (the KM and SD paragraphs respectively) in order to fit models.")
    } else {
      if(!is.null(data$me)){
        macs <- T
      } else {
        macs <- F
      }
#      return(parTable)
	fit <- lavaan::lavaan(model=parTable,data=data$ra,sample.cov=data$cm,sample.mean=data$me,estimator=estimator,sample.nobs=n,...)
	  if(silent==F){
		lavaan::summary(fit, standardized=TRUE, fit.measures=TRUE)
		invisible(fit)
	  }else{
		return(fit)
	  }
    }
  }else{
    invisible(parTable)
  }
 
}

return(suppressWarnings(lisrel(filename=filename, analyze=analyze, ...=...)))

setwd(restore.wd)

}
