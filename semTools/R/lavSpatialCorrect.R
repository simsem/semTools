spatialCorrect <- function(obj, xvar, yvar, alpha=0.05){
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
    mi <- moran(x, distsInv, na.rm = TRUE)
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

# Copy from the Moran.I function from the ape package
moran <- function (x, weight, scaled = FALSE, na.rm = FALSE, alternative = "two.sided") 
{
    if (dim(weight)[1] != dim(weight)[2]) 
        stop("'weight' must be a square matrix")
    n <- length(x)
    if (dim(weight)[1] != n) 
        stop("'weight' must have as many rows as observations in 'x'")
    ei <- -1/(n - 1)
    nas <- is.na(x)
    if (any(nas)) {
        if (na.rm) {
            x <- x[!nas]
            n <- length(x)
            weight <- weight[!nas, !nas]
        }
        else {
            warning("'x' has missing values: maybe you wanted to set na.rm=TRUE?")
            return(list(observed = NA, expected = ei, sd = NA, 
                p.value = NA))
        }
    }
    ROWSUM <- rowSums(weight)
    ROWSUM[ROWSUM == 0] <- 1
    weight <- weight/ROWSUM
    s <- sum(weight)
    m <- mean(x)
    y <- x - m
    cv <- sum(weight * y %o% y)
    v <- sum(y^2)
    obs <- (n/s) * (cv/v)
    if (scaled) {
        i.max <- (n/s) * (sd(rowSums(weight) * y)/sqrt(v/(n - 
            1)))
        obs <- obs/i.max
    }
    S1 <- 0.5 * sum((weight + t(weight))^2)
    S2 <- sum((apply(weight, 1, sum) + apply(weight, 2, sum))^2)
    s.sq <- s^2
    k <- (sum(y^4)/n)/(v/n)^2
    sdi <- sqrt((n * ((n^2 - 3 * n + 3) * S1 - n * S2 + 3 * s.sq) - 
        k * (n * (n - 1) * S1 - 2 * n * S2 + 6 * s.sq))/((n - 
        1) * (n - 2) * (n - 3) * s.sq) - 1/((n - 1)^2))
    alternative <- match.arg(alternative, c("two.sided", "less", 
        "greater"))
    pv <- pnorm(obs, mean = ei, sd = sdi)
    if (alternative == "two.sided") 
        pv <- if (obs <= ei) 
            2 * pv
        else 2 * (1 - pv)
    if (alternative == "greater") 
        pv <- 1 - pv
    list(observed = obs, expected = ei, sd = sdi, p.value = pv)
}