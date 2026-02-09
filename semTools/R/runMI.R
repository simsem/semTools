### Terrence D. Jorgensen
### Last updated: 12 March 2025
### DEPRECATED: 16 June 2024, supplanted by lavaan.mi package
### DEFUNCT: 9 February 2026

## -------------
## Main function
## -------------


##' Fit a lavaan Model to Multiple Imputed Data Sets
##'
##' This functionality has been moved to the `lavaan.mi` package.
##'
##'
##' @aliases runMI-deprecated lavaan.mi-deprecated cfa.mi-deprecated sem.mi-deprecated growth.mi-deprecated
##'
##' @param \dots Functionality has ceased to be supported, replaced by `lavaan.mi` package
##'
##' @return An error message
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##'
##' @examples
##'
##' ## See the new lavaan.mi package
##'
##' @name runMI-deprecated
##' @usage
##' runMI(...)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL


##' @rdname semTools-deprecated
##' @section `runMI()` and the `lavaan.mi` class functionality:
##' The `runMI()` function and support for `lavaan.mi-class` objects became
##' such a large part of semTools that it made sense to move that functionality
##' to its own package.  The \pkg{lavaan.mi} package is now available for users
##' to fit `lavaan` models to their multiply imputed data.  The new package
##' already fixes many bugs and provides many new features that make the old
##' semTools `lavaan.mi-class` obsolete, so the latter are now defunct.
##'
##' @export
runMI <- function(...) {
  .Defunct(msg = c("\nThe runMI() function and associated functions have been ",
                   "replaced by the new lavaan.mi package."))
}


##' @name runMI-deprecated
##' @usage
##' lavaan.mi(...)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL

##' @rdname semTools-deprecated
##' @export
lavaan.mi <- function(...) {
  .Defunct(msg = c("\nThe lavaan.mi() function and associated functions have been ",
                   "replaced by the new lavaan.mi package."))
}

##' @name runMI-deprecated
##' @usage
##' cfa.mi(...)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL

##' @rdname semTools-deprecated
##' @export
cfa.mi <- function(...) {
  .Defunct(msg = c("\nThe cfa.mi() function and associated functions have been ",
                   "replaced by the new lavaan.mi package."))
}


##' @name runMI-deprecated
##' @usage
##' sem.mi...)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL

##' @rdname semTools-deprecated
##' @export
sem.mi <- function(...) {
  .Defunct(msg = c("\nThe sem.mi() function and associated functions have been ",
                   "replaced by the new lavaan.mi package."))
}

##' @name runMI-deprecated
##' @usage
##' growth.mi(...)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL

##' @rdname semTools-deprecated
##' @export
growth.mi <- function(...) {
  .Defunct(msg = c("\nThe growth.mi() function and associated functions have ",
                   "been replaced by the new lavaan.mi package."))
}



