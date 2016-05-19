### HTMT function
#Written by Ylenio Longo

htmt <- function(data, model, ...){
	R <- lavaan::lavCor(object = data, ...)
	R <- abs(R) #this helps avoid errors
    diag(R) <- NA
    m <- lavaan::lavaanify(model)
    m <- m[m$op%in% "=~",]
    
    ##variable names for each scale / factor
    factors <- unique(m$lhs)
    var <- list()
    for(i in 1:length(factors)){
        var[[i]] <- m$rhs[which(m$lhs %in% factors[i])]
    }
    var
    
    ##mean correlation within scales
    m.cor.w <- list()
    for(i in 1:length(factors)){
        m.cor.w[[i]] <-  mean(R[var[[i]],var[[i]]], na.rm=TRUE)
    }
    m.cor.w <- as.numeric(m.cor.w)
    m.cor.w
    
    ##geometric mean correlations within scale pairs
    #all possible correlation combinations
    comb <- expand.grid(1:length(factors), 1:length(factors))
    g <- list() 
    for(i in 1:nrow(comb)){
        g[[i]] <- sqrt(m.cor.w[comb[i,2]]*m.cor.w[comb[i,1]])
    }
    g <- as.numeric(g)
    g #geometric mean results
    
    paste(comb[,2], comb[,1])
    
    ##mean correlations among items across scales
    m.cor.a <- list()
    for(i in 1:nrow(comb)){
        m.cor.a[[i]] <-  mean(R[var[[comb[i,2]]],  var[[comb[i,1]]]], na.rm=TRUE)
    }
    m.cor.a <- as.numeric(m.cor.a)
    m.cor.a
    
    ##htmt values
    outhtmt <- m.cor.a / g 
    
    ##results
    res <- matrix(outhtmt, nrow=length(factors), ncol=length(factors), dimnames=list(factors))
    colnames(res) <- factors
	class(res) <- c("lavaan.matrix.symmetric", "matrix")
    res
}
