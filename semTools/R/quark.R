quark <- function(data, id, order = 1, silent = FALSE){
  if(!is.data.frame(data) && !is.matrix(data)) {
    stop("Inappropriate data file provided.")
  }
  if(!silent) cat("Data Check Passed.\n")
  
  if(is.character(id)) id <- match(id, colnames(data))
  for(i in 1:length(id)){
    if(id[i] > ncol(data) || id[i] < 1){
      stop("At least one of the IDs is out of bounds.")
    }
  }
  if(!silent) cat("ID Check Passed.\n")
  if(!(order %in% 1:2)) stop("Currently, the order argument can take either 1 or 2.")
  
  final.collect <- list()
  final.collect$ID_Columns <- id
  final.collect$ID_Vars <- data[,id]
  final.collect$Used_Data <- data[,-c(id)]
  final.collect$Imputed_Data <- imputequark(data = final.collect$Used_Data, order = order, silent = silent)
  final.collect$Big_Data_Matrix <- bigquark(data = final.collect$Imputed_Data, silent = silent)
  cmp <- compquark(data = final.collect$Big_Data_Matrix, silent = silent)
  final.collect$Prin_Components <- cmp[[1]]
  final.collect$Prin_Components_Prcnt <- cmp[[2]]
  
  return(final.collect)
}

imputequark <- function(data, order, silent = FALSE){
  if(order==1){
    data <- aImp(data=data, silent = silent)
    data <- gImp(data=data, silent = silent)
  } else if(order==2) {
    data <- gImp(data=data, silent = silent)
    if(length(which(is.na(data > 0)))){
      data <- aImp(data=data, silent = silent)
    }
  }
  return(data)
}

gImp <- function(data, silent = FALSE){
  imputed_data <- data
  num_adds <- vector(length=ncol(data)) #number of columns combined into one for averaging.
  data.cor <- cor(data,use="pairwise",method="pearson")
  if(!silent) printCor(data.cor)
  #populate multiple matrices that can then be utilized to determine if one column should enhance another based upon
  #the correlations they share...
  if(!silent) cat("Imputing Column... \n")
  
  for(a in 1:ncol(data)){
    temp_mat <- matrix(ncol=ncol(data),nrow=nrow(data))
    list <- unique(sort(data[,a]))
    if(length(list)>1 && length(list)<=10){
      for(b in 1:nrow(data)){
        for(c in 1:length(list)){
          if(data[b,a]==list[c] && !is.na(data[b,a])){
            temp_mat[b,] <- round(colMeans(subset(data,data[,a]==list[c]),na.rm=T),digits=1)
          }
          else if(is.na(data[b,a])){
            for(p in 1:ncol(data)){
              temp_mat[b,p] <- data[b,p]
            }
          }
        }
      }
      
      #Here I need to determine if the other columns are correlated enough with the reference to ensure accuracy
      #of predictions
      temp_cor <- data.cor[,a]
      #if(countNA(temp_cor)==0){
      for(i in 1:length(temp_cor)){
        if(i!=a){
          if(abs(temp_cor[i])>=.5&&!is.na(temp_cor[i])){#Using a moderate effect size, column a, will inform other columns.
            for(x in 1:nrow(imputed_data)){
              imputed_data[x,i] <- sum(imputed_data[x,i],temp_mat[x,a],na.rm=T)
            }
            num_adds[i] <- num_adds[i] + 1
          }
        }
      }
      #}
      if(!silent) cat("\t", colnames(data)[a])
    }   
  }
  if(!silent) cat("\n")
  imputed_data <- cleanMat(m1=data,m2=imputed_data,impact=num_adds)
  imputed_data <- fixData(imputed_data)
  return(imputed_data)
}

cleanMat <- function(m1,m2,impact){
  #Impact is the number of influences on each column...
  #We need to clean up and then try to determine what final values should be...
  #Go through each of the cells...
  new_mat <- m2
  for(a in 1:ncol(m1)){
    for(b in 1:nrow(m1)){
      if(!is.na(m1[b,a])){
        new_mat[b,a] <- m1[b,a]
      }
      else if(is.na(m1[b,a])){
        new_mat[b,a] <- new_mat[b,a]/impact[a]
      }
    }
  }
  return(new_mat)
}

fixData <- function(data){
  for(a in 1:ncol(data)){
    for(b in 1:nrow(data)){
      data[b,a] <- round(data[b,a],digits=1)
    }
  }
  
  return(data)
}

aImp <- function(data, silent = FALSE){
  requireNamespace("mice")
  if(!("package:mice" %in% search())) attachNamespace("mice")
  if(!silent) cat("Starting Algorithm Imputation...\n")
  data <- mice::mice(data,maxit=1,m=1, printFlag = !silent)
  data <- mice::complete(data)
  if(!silent) cat("Ending Algorithm Imputation...\n")
  return(data)
}

bigquark <- function(data, silent = FALSE){
  if(!silent) cat("Calculating Polynomial Effects.\n")
  poly <- ((data^2)+(data^3))/2
  if(!silent) cat("Creating Matrix for Interaction Effects.\n")
  prod <- matrix(ncol=(ncol(data)-1),nrow=nrow(data))
  if(!silent) cat("Calculating Interaction Effects...0%..")
  for(i in 1:nrow(data)){
    if(!silent) printpct(percent=i/nrow(data))
    for(j in 1:(ncol(data)-1)){
      prod[i,j] <- mean(as.numeric(data[i,j])*as.numeric(data[i,(j+1):ncol(data)]))
    }
  }
  cat("\n")
  data <- cbind(data,poly,prod)
  return(data)
}

compquark <- function(data, silent = FALSE){
  if(!silent) cat("Calculating values for the PCA\n")
  pcam <- pcaquark(data, ncp=ncol(data))
  cmp <- list()
  cmp$pca <- pcam$ind$coord
  cmp$var <- pcam$eig[,3]
  colnames(cmp$pca) <- c(paste0("AuxVar",1:ncol(cmp$pca)))
  return(cmp)
}

printpct <- function(percent){
  if(round(percent,digits=10)==0)
    cat("0%..")
  if(round(percent,digits=10)==.10)
    cat("10%..")
  if(round(percent,digits=10)==.20)
    cat("20%..")
  if(round(percent,digits=10)==.30)
    cat("30%..")
  if(round(percent,digits=10)==.40)
    cat("40%..")
  if(round(percent,digits=10)==.50)
    cat("50%..")
  if(round(percent,digits=10)==.60)
    cat("60%..")
  if(round(percent,digits=10)==.70)
    cat("70%..")
  if(round(percent,digits=10)==.80)
    cat("80%..")
  if(round(percent,digits=10)==.90)
    cat("90%..")
  if(round(percent,digits=10)==1)
    cat("100%..")
}

combinequark <- function(quark,percent){
  data <- cbind(quark$ID_Vars,quark$Used_Data)
  pct <- quark$Prin_Components_Prcnt
  comp <- quark$Prin_Components
  
  for(i in 1:length(pct)){
    if(pct[i]>=percent){
      num <- i
      break
    }
  }
  return(cbind(data,comp[,1:num]))
}

# This function is modified from the FactoMinoR package.
pcaquark <- function (X, ncp = 5) {
    moy.p <- function(V, poids) {
        res <- sum(V * poids)/sum(poids)
    }
    ec <- function(V, poids) {
        res <- sqrt(sum(V^2 * poids)/sum(poids))
    }
    X <- as.data.frame(X)
    if (any(is.na(X))) {
        warnings("Missing values are imputed by the mean of the variable: you should use the imputePCA function of the missMDA package") 
        X[is.na(X)] <- matrix(apply(X,2,mean,na.rm=TRUE),ncol=ncol(X),nrow=nrow(X),byrow=TRUE)[is.na(X)]
    }
    if (is.null(rownames(X))) rownames(X) <- 1:nrow(X)
    if (is.null(colnames(X))) colnames(X) <- paste("V", 1:ncol(X), sep = "")
    colnames(X)[colnames(X) == ""] <- paste("V", 1:sum(colnames(X)==""),sep="")
    rownames(X)[is.null(rownames(X))] <- paste("row",1:sum(rownames(X)==""),sep="")
    Xtot <- X
    if (any(!sapply(X, is.numeric))) {
        auxi <- NULL
        for (j in 1:ncol(X)) if (!is.numeric(X[, j])) auxi <- c(auxi, colnames(X)[j])
        stop(paste("\nThe following variables are not quantitative: ", auxi))
    }
    ncp <- min(ncp, nrow(X) - 1, ncol(X))
    row.w <- rep(1, nrow(X))
    row.w.init <- row.w
    row.w <- row.w/sum(row.w)
    col.w <- rep(1, ncol(X))
    centre <- apply(X, 2, moy.p, row.w)
    X <- as.matrix(sweep(as.matrix(X), 2, centre, FUN = "-"))
	ecart.type <- apply(X, 2, ec, row.w)
	ecart.type[ecart.type <= 1e-16] <- 1
	X <- sweep(as.matrix(X), 2, ecart.type, FUN = "/")
	dist2.ind <- apply(sweep(X,2,sqrt(col.w),FUN="*")^2,1,sum)
	dist2.var <- apply(sweep(X,1,sqrt(row.w),FUN="*")^2,2,sum)
    tmp <- svd.triplet.quark(X, row.w = row.w, col.w = col.w, ncp = ncp)
    eig <- tmp$vs^2
    vp <- as.data.frame(matrix(NA, length(eig), 3))
    rownames(vp) <- paste("comp", 1:length(eig))
    colnames(vp) <- c("eigenvalue", "percentage of variance", 
        "cumulative percentage of variance")
    vp[, "eigenvalue"] <- eig
    vp[, "percentage of variance"] <- (eig/sum(eig)) * 100
    vp[, "cumulative percentage of variance"] <- cumsum(vp[, "percentage of variance"])
    V <- tmp$V
    U <- tmp$U
	eig <- eig[1:ncp]
    coord.ind <- sweep(as.matrix(U), 2, sqrt(eig), FUN = "*")
    coord.var <- sweep(as.matrix(V), 2, sqrt(eig), FUN = "*")
    contrib.var <- sweep(as.matrix(coord.var^2), 2, eig, "/")
    contrib.var <- sweep(as.matrix(contrib.var), 1, col.w, "*")
	dist2 <- dist2.var
    cor.var <- sweep(as.matrix(coord.var), 1, sqrt(dist2), FUN = "/")
    cos2.var <- cor.var^2
    rownames(coord.var) <- rownames(cos2.var) <- rownames(cor.var) <- rownames(contrib.var) <- colnames(X)
    colnames(coord.var) <- colnames(cos2.var) <- colnames(cor.var) <- colnames(contrib.var) <- paste("Dim", 
        c(1:ncol(V)), sep = ".")
    res.var <- list(coord = coord.var[, 1:ncp], cor = cor.var[, 
        1:ncp], cos2 = cos2.var[, 1:ncp], contrib = contrib.var[, 
        1:ncp] * 100)
	dist2 <- dist2.ind
    cos2.ind <- sweep(as.matrix(coord.ind^2), 1, dist2, FUN = "/")
    contrib.ind <- sweep(as.matrix(coord.ind^2), 1, row.w/sum(row.w), FUN = "*")
    contrib.ind <- sweep(as.matrix(contrib.ind), 2, eig, FUN = "/")
    rownames(coord.ind) <- rownames(cos2.ind) <- rownames(contrib.ind) <- names(dist2) <- rownames(X)
    colnames(coord.ind) <- colnames(cos2.ind) <- colnames(contrib.ind) <- paste("Dim", 
        c(1:ncol(U)), sep = ".")
    res.ind <- list(coord = coord.ind[, 1:ncp], cos2 = cos2.ind[, 
        1:ncp], contrib = contrib.ind[, 1:ncp] * 100, dist = sqrt(dist2))
    res <- list(eig = vp, var = res.var, ind = res.ind, svd = tmp)
    class(res) <- c("PCA", "list")
    return(res)
}

# This function is modified from the FactoMinoR package.
svd.triplet.quark <- function (X, row.w = NULL, col.w = NULL,ncp=Inf) {
	tryCatch.W.E <- function(expr){  ## function proposed by Maechlmr
		W <- NULL
		w.handler <- function(w){ # warning handler
			W <<- w
			invokeRestart("muffleWarning")
		}
		list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
										 warning = w.handler),
			 warning = W)
	}
   ncp <- min(ncp,nrow(X)-1,ncol(X))
   row.w <- row.w / sum(row.w)
    X <- sweep(X, 2, sqrt(col.w), FUN = "*")
    X <- sweep(X, 1, sqrt(row.w), FUN = "*")
	if (ncol(X) < nrow(X)){
		svd.usuelle <- tryCatch.W.E(svd(X,nu=ncp,nv=ncp))$val
		if (names(svd.usuelle)[[1]]=="message") {
		  svd.usuelle<- tryCatch.W.E(svd(t(X),nu=ncp,nv=ncp))$val
		  if (names(svd.usuelle)[[1]]=="d"){
			aux=svd.usuelle$u
			svd.usuelle$u=svd.usuelle$v
			svd.usuelle$v=aux
		  } else{
			  bb=eigen(t(X)%*%X,symmetric=TRUE)
			  svd.usuelle <- vector(mode = "list", length = 3)
			  svd.usuelle$d[svd.usuelle$d<0]=0
			  svd.usuelle$d=sqrt(svd.usuelle$d)
			  svd.usuelle$v=bb$vec[,1:ncp]
			  svd.usuelle$u=sweep(X%*%svd.usuelle$v,2,svd.usuelle$d[1:ncp],FUN="/")
		  }
		}
		U <- svd.usuelle$u
		V <- svd.usuelle$v
		if (ncp >1){
		  mult <- sign(apply(V,2,sum))
		  mult[mult==0] <- 1
		  U <- sweep(U,2,mult,FUN="*")
		  V <- sweep(V,2,mult,FUN="*")
		}
		U <- sweep(as.matrix(U), 1, sqrt(row.w), FUN = "/")
		V <- sweep(as.matrix(V), 1, sqrt(col.w), FUN = "/")
	} else {
		svd.usuelle=tryCatch.W.E(svd(t(X),nu=ncp,nv=ncp))$val
		if (names(svd.usuelle)[[1]]=="message"){
		  svd.usuelle=tryCatch.W.E(svd(X,nu=ncp,nv=ncp))$val
		  if (names(svd.usuelle)[[1]]=="d"){
			aux=svd.usuelle$u
			svd.usuelle$u=svd.usuelle$v
			svd.usuelle$v=aux
		  } else{
			  bb=eigen(X%*%t(X),symmetric=TRUE)
			  svd.usuelle <- vector(mode = "list", length = 3)
			  svd.usuelle$d[svd.usuelle$d<0]=0
			  svd.usuelle$d=sqrt(svd.usuelle$d)
			  svd.usuelle$v=bb$vec[,1:ncp]
			  svd.usuelle$u=sweep(t(X)%*%svd.usuelle$v,2,svd.usuelle$d[1:ncp],FUN="/")
		  }
		}
		U <-  svd.usuelle$v
		V <- svd.usuelle$u
		mult <- sign(apply(V,2,sum))
		mult[mult==0] <- 1
		V <- sweep(V,2,mult,FUN="*")
		U <- sweep(U,2,mult,FUN="*")
		U <- sweep(U, 1, sqrt(row.w), FUN = "/")
		V <- sweep(V, 1, sqrt(col.w), FUN = "/")
	}
    vs <- svd.usuelle$d[1:min(ncol(X),nrow(X)-1)]
	num <- which(vs[1:ncp]<1e-15)
    if (length(num)==1){
	  U[,num] <- U[,num]*vs[num]
      V[,num] <- V[,num]*vs[num]
	} 
    if (length(num)>1){
	  U[,num] <- sweep(U[,num],2,vs[num],FUN="*")
      V[,num] <- sweep(V[,num],2,vs[num],FUN="*")
	}
    res <- list(vs = vs, U = U, V = V)
    return(res)
}

# This function is copied from the psych package: lowerMat
printCor <- function (R, digits = 2) {
    lowleft <- lower.tri(R, diag = TRUE)
    nvar <- ncol(R)
    nc <- digits + 3
    width <- getOption("width")
    k1 <- width/(nc + 2)
    if (is.null(colnames(R))) {
        colnames(R) <- paste("C", 1:nvar, sep = "")
    }
    if (is.null(rownames(R))) {
        rownames(R) <- paste("R", 1:nvar, sep = "")
    }
    colnames(R) <- abbreviate(colnames(R), minlength = digits + 
        3)
    nvar <- ncol(R)
    nc <- digits + 3
    if (k1 * nvar < width) {
        k1 <- nvar
    }
    k1 <- floor(k1)
    fx <- format(round(R, digits = digits))
    if (nrow(R) == ncol(R)) {
        fx[!lowleft] <- ""
    }
    for (k in seq(0, nvar, k1)) {
        if (k < nvar) {
            print(fx[(k + 1):nvar, (k + 1):min((k1 + k), nvar)], 
                quote = FALSE)
        }
    }
}
