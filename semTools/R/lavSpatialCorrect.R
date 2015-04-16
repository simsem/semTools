spatialCorrect <- function(obj, xvar, yvar, alpha=0.05){
  require(lavaan)
  require(ape)
  
  #first, get the residuals from the model
  resids <- as.data.frame(residuals(obj, "casewise"))
  
  #get only endogenous variables
  resids <- resids[,which(apply(resids, 2, function(x) length(unique(x))) !=1)]
  

  #make a distance matrix
  distMat <- as.matrix(dist(cbind(xvar, yvar)))
  
  #invert this matrix for weights
  distsInv <- 1/distMat
  diag(distsInv) <- 0
  
  morans_i <- lapply(resids,  function(x){
    mi <- Moran.I(x, distsInv, na.rm = TRUE)
    if(mi$p.value>alpha){
      mi$n.eff <- nrow(resids) #don't correct sample size
    }else{
      #large sample size approximation
      mi$n.eff <- nrow(resids)*(1-mi$observed)/(1+mi$observed)
    }
    
    #return the list
    mi
    
    })
  
  
  #get the vcov matrix
  v <- diag(vcov(obj))
  n <- nrow(resids)
  #using new sample sizes, for each variable, calculate new Z-scores
  params <- lapply(names(morans_i), function(acol){
    idx <- grep(paste0(acol, "~"),names(v))  #regression or covariances
    idx <- c(idx, grep(paste0("=~",acol),names(v)))  #latent variable definitions
    v_idx <- v[idx]*n/morans_i[[acol]]$n.eff
    ret <-  data.frame(Parameter = names(v)[idx], Estimate=coef(obj)[idx], 
                       n.eff = morans_i[[acol]]$n.eff, Std.err = sqrt(v_idx))
    ret[["Z-value"]] <- ret$Estimate/ret$Std.err
    ret[["P(>|z|)"]] <- 2*pnorm(abs(ret[["Z-value"]]), lower.tail=F)
  
    ret
    
  })
  names(params) <- names(morans_i)

  mi <- lapply(morans_i, function(m) {
         data.frame(observed=m$observed, expected=m$expected, 
                    sd=m$sd, p.value = m$p.value, n.eff = m$n.eff)
  })

 return(list(Morans_I = mi, parameters=params)) 
  
}

#lavSpatialCorrect(obj, xvar, yvar)