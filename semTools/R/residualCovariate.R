# residualCovariate: Residual centered all target indicators by covariates

residualCovariate <- function(data, targetVar, covVar) {
    x <- as.list(match.call())
    cov <- eval(x$covVar)
    target <- eval(x$targetVar)
    if (all(is.numeric(cov))) 
        cov <- colnames(data)[cov]
    if (all(is.numeric(target))) 
        target <- colnames(data)[target]
    express <- paste("cbind(", paste(target, collapse = ", "), ") ~ ", paste(cov, collapse = " + "), sep = "")
    data[, target] <- lm(express, data = data)$residuals
    return(data)
} 
