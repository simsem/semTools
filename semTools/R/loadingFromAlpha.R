# loadingFromAlpha: Find a standardized factor loading that provide a specified
# alpha value

loadingFromAlpha <- function(alpha, ni) {
    denominator <- ni - ((ni - 1) * alpha)
    result <- sqrt(alpha/denominator)
    return(result)
} 
