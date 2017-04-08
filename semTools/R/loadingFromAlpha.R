### Sunthud Pornprasertmanit
### Last updated: 3 April 2017


#' Find standardized factor loading from coefficient alpha
#'
#' Find standardized factor loading from coefficient alpha assuming that all
#' items have equal loadings.
#'
#' @param alpha A desired coefficient alpha value.
#' @param ni A desired number of items.
#' @return \item{result}{The standardized factor loadings that make desired
#' coefficient alpha with specified number of items.}
#' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#' @examples
#'
#' loadingFromAlpha(0.8, 4)
#'
#' @export
loadingFromAlpha <- function(alpha, ni) {
    denominator <- ni - ((ni - 1) * alpha)
    result <- sqrt(alpha/denominator)
    return(result)
}


