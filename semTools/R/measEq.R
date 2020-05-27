### Terrence D. Jorgensen
### Last updated: 27 May 2020
### lavaan model syntax-writing engine for new measEq() to replace
### measurementInvariance(), measurementInvarianceCat(), and longInvariance()



## ----------------------
## Model-Fitting Function
## ----------------------


measEq <- function(configural.model,
                   ID.fac = c("std.lv","auto.fix.first","effects.coding"),
                   ID.cat = c("Wu.Estabrook.2016","Mplus","Millsap.Tein.2004","LISREL"),
                   ID.thr = c(1L, 2L), # only for ID.cat == "Millsap.Tein.2004"
                   group = NULL, longFacNames = list(), longIndNames = list(),
                   #group.equal = "", long.equal = "",
                   group.partial = "", long.partial = "",
                   auto = "all", extra = NULL,
                   test.seq = c("thresholds","loadings","intercepts","means",
                                "lv.variances","residuals"), # optional resid/lv.autocov, residual/lv.covariances
                   #fixed.x = TRUE, strict = FALSE, quiet = FALSE,
                   warn = TRUE, debug = FALSE,
                   alpha = .05, fit.measures = "default", argsLRT = list(), ...) {

  #TODO: check GLIST structure for multilevel (with single and multiple groups)
  #TODO: compatibility with auxiliary(), runMI(), parcelAllocation(), permuteMeasEq()
  #TODO: develop automated anchorSelection() and measEq.partial()
  #TODO: if (inherits(configural.model, "lavaan.measEq")) {continue sequence}?

  ## This example might help:  https://groups.google.com/d/msg/lavaan/LvALeUpJBDg/2zD1CoikAQAJ

  #TODO: add argument to accept measEq.partial output, to continue sequence (or make and update() method?)
  if (is.character(group.partial)) {
    if (group.partial == "" && length(group.partial) == 1L) {
      group.partial <- data.frame(stringsAsFactors = FALSE, lhs = character(0),
                                  op = character(0), rhs = character(0))
    } else {
      group.partial <- lavaan::lavParseModelString(group.partial,
                                                   as.data.frame. = TRUE,
                                                   warn = warn, debug = debug)
    }
  } #TODO: else {extract information from a measEq.partial object}
  if (is.character(long.partial)) {
    if (long.partial == "" && length(long.partial) == 1L) {
      long.partial <- data.frame(stringsAsFactors = FALSE, lhs = character(0),
                                 op = character(0), rhs = character(0))
    } else {
      long.partial <- lavaan::lavParseModelString(long.partial,
                                                  as.data.frame. = TRUE,
                                                  warn = warn, debug = debug)
    }
  } #TODO: else {extract information from a measEq.partial object}

  ## pass arguments to measEq.syntax(), which performs checks


}



## -----------------
## Class and Methods
## -----------------

##' Class for Representing a Measurement-Equivalence Model
##'
##' This class of object stores information used to automatically generate
##' lavaan model syntax to represent user-specified levels of measurement
##' equivalence/invariance across groups and/or repeated measures. See
##' \code{\link{measEq.syntax}} for details.
##'
##'
##' @name measEq.syntax-class
##' @aliases measEq.syntax-class show,measEq.syntax-method
##'   summary,measEq.syntax-method as.character,measEq.syntax-method
##'   update,measEq.syntax-method
##' @docType class
##'
##' @slot package \code{character} indicating the software package used to
##'   represent the model. Currently, only \code{"lavaan"} is available, which
##'   uses the LISREL representation (see \code{\link[lavaan]{lavOptions}}).
##'   In the future, \code{"OpenMx"} may become available, using RAM
##'   representation.
##' @slot model.type \code{character}. Currently, only "cfa" is available.
##'   Future versions may allow for MIMIC / RFA models, where invariance can be
##'   tested across levels of exogenous variables explicitly included as
##'   predictors of indicators, controlling for their effects on (or correlation
##'   with) the common factors.
##' @slot call The function call as returned by \code{match.call()}, with
##'   some arguments updated if necessary for logical consistency.
##' @slot meanstructure \code{logical} indicating whether a mean structure is
##'   included in the model.
##' @slot numeric \code{character} vector naming \code{numeric} manifest indicators.
##' @slot ordered \code{character} vector naming \code{ordered} indicators.
##' @slot parameterization \code{character}. See \code{\link[lavaan]{lavOptions}}.
##' @slot specify \code{list} of parameter matrices, similar in form to the
##'   output of \code{\link[lavaan]{lavInspect}(fit, "free")}. These matrices
##'   are \code{logical}, indicating whether each parameter should be specified
##'   in the model syntax.
##' @slot values \code{list} of parameter matrices, similar in form to the
##'   output of \code{\link[lavaan]{lavInspect}(fit, "free")}. These matrices
##'   are \code{numeric}, indicating whether each parameter should be freely
##'   estimated (indicated by \code{NA}) or fixed to a particular value.
##' @slot labels \code{list} of parameter matrices, similar in form to the
##'   output of \code{\link[lavaan]{lavInspect}(fit, "free")}. These matrices
##'   contain \code{character} labels used to constrain parameters to equality.
##' @slot constraints \code{character} vector containing additional equality
##'   constraints used to identify the model when \code{ID.fac = "fx"}.
##' @slot ngroups \code{integer} indicating the number of groups.
##'
##' @param x,object an object of class \code{measEq.syntax}
##' @param package \code{character} indicating the package for which the model
##'   syntax should be generated.  Currently, only \code{"lavaan"} and
##'   \code{"mplus"} are supported.
##' @param params \code{character} vector indicating which type(s) of parameter
##'   to print syntax for. Must match a type that can be passed to
##'   \code{group.equal} or \code{long.equal}, but \code{"residual.covariances"}
##'   and \code{"lv.covariances"} will be silently ignored. Instead, requesting
##'   \code{"residuals"} or \code{"lv.variances"} will return covariances along
##'   with variances. By default (\code{NULL}), all types are printed.
##' @param single \code{logical} indicating whether to concatenate lavaan
##'   \code{\link[lavaan]{model.syntax}} into a single \code{character} string.
##'   Setting \code{FALSE} will return a vector of strings, which may be
##'   convenient (or even necessary to prevent an error) in
##'   models with long variable names, many variables, or many groups.
##' @param groups.as.blocks \code{logical} indicating whether to write lavaan
##'   \code{\link[lavaan]{model.syntax}} using vectors of labels and values
##'   for multiple groups (the default: \code{FALSE}), or whether to write
##'   a separate "block" of syntax per group. The block structure could allow
##'   users to apply the generated multigroup syntax (after some editing) to
##'   test invariance across levels in a multilevel SEM (see final example on
##'   \code{\link{measEq.syntax}} help page).
##' @param verbose \code{logical} indicating whether to print a summary to the
##'   screen (default). If \code{FALSE}, only a pattern matrix is returned.
##' @param ... Additional arguments to the \code{call}, or arguments with
##'   changed values.
##' @param evaluate If \code{TRUE}, evaluate the new \code{call}; otherwise,
##'   return the new \code{call}.
##' @param change.syntax \code{lavaan \link[lavaan]{model.syntax}} specifying
##'   labels or fixed/free values of parameters in \code{object}.
##'   These provide some flexibility to customize
##'   existing parameters without having to copy/paste the output of
##'   \code{as.character(object)} into an R script. For example,
##'   \code{group.partial} will free a parameter across all groups, but
##'   \code{update} allows users to free the parameter in just one group
##'   while maintaining equality constraints among other groups.
##'
##' @return
##'   \item{summary}{\code{signature(object = "measEq.syntax", verbose = TRUE)}:
##'     A \code{character} matrix indicating the pattern of \code{numeric},
##'     \code{ordered}, or latent indicators loading on common factors.
##'     By default (\code{verbose = TRUE}), \code{summary} also prints descriptive
##'     details about the model, including the numbers of indicators and factors,
##'     and which parameters are constrained to equality.}
##'   \item{show}{\code{signature(object = "measEq.syntax")}: Prints a message
##'     about how to use the \code{object} for model fitting. Invisibly
##'     returns the \code{object}.}
##'   \item{update}{\code{signature(object = "measEq.syntax"), ...,
##'     evaluate = TRUE, change.syntax = NULL}: Creates a new
##'     \code{object} with updated arguments in \code{...}, or updated
##'     parameter labels or fixed/free specifications in \code{object}.}
##'   \item{as.character}{\code{signature(x = "measEq.syntax", package = "lavaan")}:
##'     Converts the \code{measEq.syntax} object to model syntax that can be
##'     copy/pasted or written to a syntax file to be edited before analysis,
##'     or simply passed to \code{\link[lavaan]{lavaan}} to fit the model to
##'     data. Generated M\emph{plus} syntax could also be utilized using the
##'     \pkg{MplusAuthomation} package.}
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##' \email{TJorgensen314@@gmail.com})
##'
##' @examples
##' ## See ?measEq.syntax help page for examples using lavaan
##'
## ## Here, we illustrate how measEq.syntax() objects can be used in
## ## tandem with MplusAutomation.
##
## \dontrun{
## ## borrow example data from Mplus user guide
## myData <- read.table("http://www.statmodel.com/usersguide/chap5/ex5.16.dat")
## names(myData) <- c("u1","u2","u3","u4","u5","u6","x1","x2","x3","g")
## bin.mod <- '
##   FU1 =~ u1 + u2 + u3
##   FU2 =~ u4 + u5 + u6
## '
## ## pretend the 2 factors actually measure the same factor (FU) twice
## longFacNames <- list(FU = c("FU1","FU2"))
## syntax.scalar <- measEq.syntax(configural.model = bin.mod,
##               data = myData, ordered = paste0("u", 1:6),
##               parameterization = "theta",
##               ID.fac = "std.lv", ID.cat = "Wu.Estabrook.2016",
##               group = "g", longFacNames = longFacNames,
##               group.equal = c("thresholds","loadings","intercepts"),
##               long.equal = c("thresholds","loadings","intercepts"))
## library(MplusAutomation)
## mpInp <- mplusObject(rdata = myData, TITLE = "Scalar Invariance",
##                      VARIABLE = "GROUPING = g (1=g1 2=g2);",
##                      usevariables = c(paste0("u", 1:6), "g"),
##                      ANALYSIS = "ESTIMATOR = WLSMV;",
##                 ## model specification from measEq.syntax():
##                      MODEL = as.character(syntax.scalar, package = "mplus"))
## ## add details for Mplus script:
## mpInp <- update(mpInp, ANALYSIS = ~ . + "PARAMETERIZATION = THETA;",
##                 VARIABLE = ~ . + "CATEGORICAL = u1 u2 u3 u4 u5 u6;")
## ## fit model
## mpOut <- mplusModeler(mpInp, modelout = "scalar.inp", run = 1L)
## }
#TODO: add configural and DIFFTEST example
##'
setClass("measEq.syntax", slots = c(package = "character", # lavaan, OpenMx in the future?
                                    model.type = "character", # cfa, extend to mimic/rfa?
                                    call = "call",
                                    meanstructure = "logical",
                                    numeric = "character",
                                    ordered = "character",
                                    parameterization = "character",
                                    specify = "list",
                                    values = "list",
                                    labels = "list",
                                    constraints = "character",
                                    updates = "list", # 2 data.frames: labels and values
                                    ngroups = "integer"))

##' @rdname measEq.syntax-class
##' @aliases as.character,measEq.syntax-method
##' @export
setMethod("as.character", "measEq.syntax", function(x, package = "lavaan",
                                                    params = NULL, single = TRUE,
                                                    groups.as.blocks = FALSE) {

  package <- tolower(package)[1]

  if (package == "mplus") {
    LL <- x@specify[[1]]$lambda
    nn <- c(rownames(LL), colnames(LL))
    over8 <- nchar(nn) > 8L
    if (any(over8)) warning('Mplus only allows variable names to have 8 ',
                            'characters. The following variable names in ',
                            'your model exceed 8 characters:\n',
                            paste(nn[over8], collapse = ", "), '\n',
                            'Consider shortening variable names before ',
                            'printing an Mplus MODEL statement.')
    ## print everything leading up to the MODEL command
    script <- c("MODEL:\n")
    if (length(x@ordered)) {
      script <- c(paste0("!Make sure your VARIABLE command indicates the ",
                         "following variables as CATEGORICAL:\n!",
                         paste(x@ordered, collapse = ", "), '\n'), script)
    }
    if (x@ngroups > 1L) {
      script <- c(script, "!This is the first group's MODEL.\n!Group 2's MODEL",
                  "!will be labeled as 'g2', and so on for any other groups.\n")
    }
    script <- c(script, write.mplus.syntax(object = x, group = 1L,
                                           params = params))
    if (x@ngroups > 1L) for (g in 2:x@ngroups) {
      script <- c(script, paste0("\nMODEL g", g, ":\n"),
                  write.mplus.syntax(object = x, group = g, params = params))
    }

    return(paste(script, collapse = "\n")) # always return a single string


  } else if (package == "lavaan") {
    script <- character(0)

    pmatList <- c("lambda","tau","nu","delta","theta","alpha","psi")
    names(pmatList) <- c("loadings","thresholds","intercepts","scales",
                         "residuals","means","lv.variances")
    ## selected parameter types?
    if (!is.null(params)) {
      requested <- intersect(names(pmatList), params)
      if (!length(requested)) stop('invalid choice: params = c("',
                                   paste(params, collapse = '", "'), '")\n',
                                   'Valid choices include: ',
                                   paste(names(pmatList), collapse = ", "))
      pmatList <- pmatList[requested]
    }

    if (groups.as.blocks) {
      ## loop over groups
      for (gg in 1:x@ngroups) {
        script <- c(script, paste("group:", gg, "\n", collapse = ""))

        ## loop over pmats
        for (pm in pmatList) {
          if (!pm %in% names(x@specify[[gg]])) next

          if (pm == "lambda" && "beta" %in% names(x@specify[[gg]])) {
            ## add higher-order loadings to lambda matrix
            specify <- list(rbind(x@specify[[gg]]$lambda, x@specify[[gg]]$beta))
            value   <- list(rbind(x@values[[gg]]$lambda, x@values[[gg]]$beta))
            label   <- list(rbind(x@labels[[gg]]$lambda, x@labels[[gg]]$beta))
          } else {
            specify <- list(x@specify[[gg]][[pm]])
            value   <- list(x@values[[gg]][[pm]])
            label   <- list(x@labels[[gg]][[pm]])
          }

          script <- c(script, write.lavaan.syntax(pmat = pm, specify = specify,
                                                  value = value, label = label))
        } # end pm
      } # end gg

    } else {
      ## the usual multigroup lavaan syntax:
      ## loop over pmats, send all groups together
      for (pm in pmatList) {
        if (!pm %in% names(x@specify[[1]])) next

        if (pm == "lambda" && "beta" %in% names(x@specify[[1]])) {
          ## add higher-order loadings to lambda matrix
          specify.l <- lapply(x@specify, "[[", i = "lambda")
          value.l   <- lapply(x@values , "[[", i = "lambda")
          label.l   <- lapply(x@labels , "[[", i = "lambda")
          specify.b <- lapply(x@specify, "[[", i = "beta")
          value.b   <- lapply(x@values , "[[", i = "beta")
          label.b   <- lapply(x@labels , "[[", i = "beta")
          specify   <- mapply(rbind, specify.l, specify.b, SIMPLIFY = FALSE)
          value     <- mapply(rbind, value.l  , value.b  , SIMPLIFY = FALSE)
          label     <- mapply(rbind, label.l  , label.b  , SIMPLIFY = FALSE)
        } else {
          specify   <- lapply(x@specify, "[[", i = pm)
          value     <- lapply(x@values, "[[", i = pm)
          label     <- lapply(x@labels, "[[", i = pm)
        }

        script <- c(script, write.lavaan.syntax(pmat = pm, specify = specify,
                                                value = value, label = label))
      }
    }

    if (length(x@constraints)) script <- c(script,
                                           "## MODEL CONSTRAINTS:\n",
                                           x@constraints, "")
  }
  #TODO: else if (package == "openmx") # concatenate matrices for RAM specification

  ## convert GLIST objects to a character string
  if (single) return(paste(script, collapse = "\n"))
  script
})

##' @rdname measEq.syntax-class
##' @aliases show,measEq.syntax-method
##' @export
setMethod("show", "measEq.syntax", function(object) {
  cat('This object contains information for specifying a CFA using lavaan',
      'model syntax.\nTo print the syntax (to copy/paste it into an R script),',
      'use the as.character() method:\n\n\tcat(as.character(object))\n\nTo fit',
      'this model to data, save the syntax to an object and pass it to lavaan:',
      '\n\n\tmodel <- as.character(object)\n\tfit <- lavaan(model, ...)',
      '\n\nTo view some key features of the model use: \n\n\tsummary(object)')
  invisible(object)
})

##' @rdname measEq.syntax-class
##' @aliases summary,measEq.syntax-method
##' @export
setMethod("summary", "measEq.syntax", function(object, verbose = TRUE) {

  nG <- object@ngroups
  nOrd <- length(object@ordered)
  higher <- !is.null(object@specify[[1]]$beta)

  ## create pattern matrix
  lambda <- object@specify[[1]]$lambda
  lambda[!lambda] <- ""
  for (RR in 1:nrow(lambda)) {
    if (rownames(lambda)[RR] %in% object@ordered) {
      lambda[RR, object@specify[[1]]$lambda[RR, ] ] <- "ord"
    } else lambda[RR, object@specify[[1]]$lambda[RR, ] ] <- "num"
  }
  if (higher) {
    beta <- object@specify[[1]]$beta
    beta[!beta] <- ""
    for (RR in 1:nrow(beta)) {
      beta[RR, object@specify[[1]]$beta[RR, ] ] <- "lat"
    }
    rownames(beta) <- paste("**", rownames(beta), "**")
    lambda <- rbind(lambda, beta)
  }
  if (!verbose) return(lambda)


  ## Basics: number of groups, factors, and indicators (higher order?); ID.fac
  nHigher <- if (higher) sum(apply(object@specify[[1]]$beta, 2, any)) else 0L
  if (object@call$ID.fac == "ul" && !object@meanstructure) {
    ID.fac.text <- 'first indicator`s factor loading was fixed to 1.'
  } else if (object@call$ID.fac == "ul" && object@meanstructure) {
    ID.fac.text <- paste('first indicator`s intercept and factor loading were',
                         'fixed to 0 and 1, respectively.')
  } else if (object@call$ID.fac == "uv" && !object@meanstructure) {
    ID.fac.text <- paste('factor variances were fixed to 1, unless equality',
                         'constraints on factor loadings allow them to be freed.')
  } else if (object@call$ID.fac == "uv" && object@meanstructure) {
    ID.fac.text <- paste('factor means and variances were fixed to 0 and 1,',
                         'respectively, unless equality constraints on',
                         'measurement parameters allow them to be freed.')
  } else if (object@call$ID.fac == "fx") {
    ID.fac.text <- paste('factor loadings were constrained to average 1',
                         if (object@meanstructure) 'and intercepts were constrained to average 0',
                         'within each factor. In models with partial',
                         'invariance, only the factor loadings',
                         if (object@meanstructure) 'and intercepts',
                         'that were constrained to equality across ALL groups',
                         'and repeated measures (when applicable) are used to',
                         'identify the common-factor distribution.')
  }
  cat('This lavaan model syntax specifies a CFA with ',
      nrow(object@specify[[1]]$lambda), ' manifest indicators ',
      if (nOrd == 1L) {
        paste0('(', nOrd, ' of which is ordinal) ')
      } else if (nOrd > 1L) {
        paste0('(', nOrd, ' of which are ordinal) ')
      }, 'of ', ncol(object@specify[[1]]$lambda), ' common factor(s)',
      if (nHigher == 1L) {
        paste(',', nHigher, 'of which is a higher-order factor. ')
      } else if (nHigher > 1L) {
        paste(',', nHigher, 'of which are higher-order factors. ')
      } else '.\n\n', 'To identify the ',
      if (object@meanstructure) 'location and ',
      'scale of each common factor, the ', ID.fac.text, "\n\n", sep = '')

  ## if (ordered) ID.cat and parameterization
  if (nOrd) {
    if (object@call$ID.cat == "wu") {
      ID.cat.author <- 'recommended by Wu & Estabrook (2016). '
      ID.cat.DOI <- 'https://doi.org/10.1007/s11336-016-9506-0 \n\n'
    } else if (object@call$ID.cat == "millsap") {
      ID.cat.author <- 'recommended by Millsap & Tein (2004). '
    } else if (object@call$ID.cat == "mplus") {
      ID.cat.author <- 'used by default in the Mplus (and lavaan) software. '
    } else if (object@call$ID.cat == "lisrel") {
      ID.cat.author <- 'used by default in the LISREL software. '
    }
    if (object@call$ID.cat != "wu") ID.cat.DOI <- 'http://dx.doi.org/10.1207/S15327906MBR3903_4 \n\n'
    cat('The location and scale of each latent item-response underlying ', nOrd,
        ' ordinal indicators were identified using the "', object@parameterization,
        '" parameterization, and the identification constraints ',
        ID.cat.author, 'For details, read:\n\n\t', ID.cat.DOI, sep = '')
  }

  ## number of occassions per longitudinal construct
  if (length(object@call$longFacNames)) {
    cat('The following factors were measured on multiple occasions:\n')
    for (f in names(object@call$longFacNames)) {
      cat('\t"', f, '" was measured on ', length(object@call$longFacNames[[f]]),
          ' occasions\n', sep = '')
    }
    cat('\n')
  }

  ## print pattern matrix
  cat('Pattern matrix indicating num(eric), ord(ered), and lat(ent)',
      'indicators per factor:\n\n')
  print(lambda, quote = FALSE)
  cat('\n')

  ## without any constraints, call it the configural model
  no.group.equal <- length(object@call$group.equal) == 1L && object@call$group.equal == ""
  no.long.equal <- length(object@call$long.equal) == 1L && object@call$long.equal == ""
  if (no.group.equal && no.long.equal) {
    cat('\nThis model hypothesizes only configural invariance.\n\n')
    ## return pattern matrix
    return(invisible(lambda))
  }
  ## otherwise, print the constraints & exceptions

  ## constrained parameters across groups (+ partial exceptions)
  if (nG > 1L) {
    if (no.group.equal) {
      cat('No parameters were constrained to equality across groups.\n')
    } else {
      cat('The following types of parameter were constrained to',
          'equality across groups:\n\n')
      for (i in object@call$group.equal) {
        group.partial <- object@call$group.partial
        ## first, check for exceptions
        if (i == "loadings") {
          man.ind <- group.partial$rhs %in% rownames(object@specify[[1]]$lambda)
          group.partial <- group.partial[group.partial$op == "=~" & man.ind, ]

        } else if (i == "regressions") {
          lat.ind <- group.partial$rhs %in% colnames(object@specify[[1]]$lambda)
          group.partial <- group.partial[group.partial$op == "=~" & lat.ind, ]

        } else if (i == "thresholds") {
          man.ind <- group.partial$lhs %in% rownames(object@specify[[1]]$lambda)
          group.partial <- group.partial[group.partial$op == "|" & man.ind, ]

        } else if (i == "residuals") {
          man.ind <- group.partial$rhs %in% rownames(object@specify[[1]]$lambda)
          same.ind <- group.partial$rhs == group.partial$lhs
          group.partial <- group.partial[group.partial$op == "~~" & man.ind & same.ind, ]

        } else if (i == "residual.covariances") {
          man.ind <- group.partial$rhs %in% rownames(object@specify[[1]]$lambda)
          same.ind <- group.partial$rhs == group.partial$lhs
          group.partial <- group.partial[group.partial$op == "~~" & man.ind & !same.ind, ]

        } else if (i == "lv.variances") {
          lat <- group.partial$rhs %in% colnames(object@specify[[1]]$lambda)
          same <- group.partial$rhs == group.partial$lhs
          group.partial <- group.partial[group.partial$op == "~~" & lat & same, ]

        } else if (i == "lv.covariances") {
          lat <- group.partial$rhs %in% colnames(object@specify[[1]]$lambda)
          same <- group.partial$rhs == group.partial$lhs
          group.partial <- group.partial[group.partial$op == "~~" & lat & !same, ]

        } else if (i == "intercepts") {
          man.ind <- group.partial$lhs %in% rownames(object@specify[[1]]$lambda)
          group.partial <- group.partial[group.partial$op == "~1" & man.ind, ]

        } else if (i == "means") {
          lat <- group.partial$lhs %in% colnames(object@specify[[1]]$lambda)
          group.partial <- group.partial[group.partial$op == "~1" & lat, ]
        }

        ## then print a message
        cat('\t', i,
            if (nrow(group.partial)) ', with the exception of:\n',
            '\n', sep = '')
        if (nrow(group.partial)) {
          rownames(group.partial) <- paste("            row-",
                                           rownames(group.partial), ":   ",
                                           sep = "")
          print(group.partial)
          cat('\n')
        }
      }

    }
    cat('\n')
  }
  ## constrained parameters across repeated measures (+ partial exceptions)
  if (length(object@call$longFacNames)) {
    if (no.long.equal) {
      cat('No parameters were constrained to equality across repeated measures:\n')
    } else {
      cat('The following types of parameter were constrained to equality',
          'across repeated measures:\n\n')
      for (i in object@call$long.equal) {
        long.partial <- object@call$long.partial
        ## first, check for exceptions
        if (i == "loadings") {
          man.ind <- long.partial$rhs %in% names(object@call$longIndNames)
          long.partial <- long.partial[long.partial$op == "=~" & man.ind, ]

        } else if (i == "regressions") {
          lat.ind <- long.partial$rhs %in% names(object@call$longFacNames)
          long.partial <- long.partial[long.partial$op == "=~" & lat.ind, ]

        } else if (i == "thresholds") {
          man.ind <- long.partial$lhs %in% names(object@call$longIndNames)
          long.partial <- long.partial[long.partial$op == "|" & man.ind, ]

        } else if (i == "residuals") {
          man.ind <- long.partial$rhs %in% names(object@call$longIndNames)
          same.ind <- long.partial$rhs == long.partial$lhs
          long.partial <- long.partial[long.partial$op == "~~" & man.ind & same.ind, ]

        } else if (i == "lv.variances") {
          lat <- long.partial$rhs %in% names(object@call$longFacNames)
          same <- long.partial$rhs == long.partial$lhs
          long.partial <- long.partial[long.partial$op == "~~" & lat & same, ]

        } else if (i == "intercepts") {
          man.ind <- long.partial$lhs %in% names(object@call$longIndNames)
          long.partial <- long.partial[long.partial$op == "~1" & man.ind, ]

        } else if (i == "means") {
          lat <- long.partial$lhs %in% names(object@call$longFacNames)
          long.partial <- long.partial[long.partial$op == "~1" & lat, ]
        }

        ## then print a message
        cat('\t', i,
            if (nrow(long.partial)) ', with the exception of:\n',
            '\n', sep = '')
        if (nrow(long.partial)) {
          rownames(long.partial) <- paste("            row-",
                                          rownames(long.partial), ":   ",
                                          sep = "")
          print(long.partial)
          cat('\n')
        }
      }

    }
    cat('\n')
  }

  ## return pattern matrix
  invisible(lambda)
})

updateMeasEqSyntax <- function(object, ..., evaluate = TRUE,
                               change.syntax = NULL) {
  # data.frame(stringsAsFactors = FALSE, extras = c(TRUE, FALSE),
  #            override = c(TRUE, TRUE, FALSE, FALSE),
  #            eval = c(TRUE, TRUE, TRUE, TRUE,
  #                     FALSE, FALSE, FALSE, FALSE),
  #            TODO = c("extras; eval; transfer, augment, and apply @updates",
  #                     "apply and augment @updates",
  #                     "extras; eval; transfer and apply @updates", "return object",
  #                     "nothing, can't add to call", "nothing, can't add to call",
  #                     "extras, return call", "return call")) -> foo
  # foo[order(foo$extras), ]

  #  extras override  eval                                                TODO
  # 1 FALSE     TRUE  TRUE                          apply and augment @updates
  # 2 FALSE    FALSE  TRUE                                       return object *
  # 3 FALSE     TRUE FALSE                          nothing, can't add to call *
  # 4 FALSE    FALSE FALSE                                         return call *
  # 5  TRUE     TRUE  TRUE extras; eval; transfer, augment, and apply @updates
  # 6  TRUE    FALSE  TRUE                     extras, eval, transfer @updates
  # 7  TRUE     TRUE FALSE                          nothing, can't add to call *
  # 8  TRUE    FALSE FALSE                                 extras, return call *

  #extras <- match.call(expand.dots = FALSE)$...
  extras <- list(...)
  custom <- !is.null(change.syntax)

  ## regardless of customization/evaluation, extras can be added to call first
  if (length(extras)) {
    ## prep 5:8
    call <- object@call
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
    if (!evaluate) {
      if (custom) warning('cannot apply "change.syntax" ',
                          'argument when evaluate=FALSE.')
      ## finish 7:8
      return(call)
    }
  } else if (!evaluate) {
    if (custom) warning('cannot apply "change.syntax" ',
                        'argument when evaluate=FALSE.')
    ## finish 3:4
    return(object@call)
  } else if (!custom) {
    ## finish 2
    return(object)
  }

  #  extras override  eval                                        TODO
  # 1 FALSE     TRUE  TRUE                  apply and augment @updates
  # 5  TRUE     TRUE  TRUE eval; transfer, augment, and apply @updates
  # 6  TRUE    FALSE  TRUE           eval; transfer and apply @updates

  if (length(extras)) {
    ## prep 5:6
    out <- eval(call, parent.frame())
    if (nrow(object@updates$values)) out@updates$values <- object@updates$values
    if (nrow(object@updates$labels)) out@updates$labels <- object@updates$labels
  } else out <- object # "prep" 1

  #  extras override  eval                        TODO
  # 1 FALSE     TRUE  TRUE  apply and augment @updates
  # 5  TRUE     TRUE  TRUE augment, and apply @updates
  # 6  TRUE    FALSE  TRUE              apply @updates

  ## check before augmenting to prep 1 and 5
  if (!is.null(change.syntax)) {
    stopifnot(is.character(change.syntax))
    ## convert syntax to data.frame of updates to make
    UPDATES <- char2update(object, change.syntax, return.object = FALSE)
    out@updates$values <- rbind(out@updates$values, UPDATES$values)
    out@updates$labels <- rbind(out@updates$labels, UPDATES$labels)
  }

  ## nothing left to do but apply @updates

  ## loop over any values/labels (now stored) to finish 1 and 5:6
  if (nrow(out@updates$values)) for (RR in 1:nrow(out@updates$values)) {
    valueArgs <- c(list(object = out, slotName = "values"),
                   as.list(out@updates$values[RR, ]))
    BB <- out@updates$values$group[RR]
    matName <- out@updates$values$matName[RR]
    out@values[[BB]][[matName]] <- do.call(override, valueArgs)
  }
  if (nrow(out@updates$labels)) for (RR in 1:nrow(out@updates$labels)) {
    labelArgs <- c(list(object = out, slotName = "labels"),
                   as.list(out@updates$labels[RR, ]))
    BB <- out@updates$labels$group[RR]
    matName <- out@updates$labels$matName[RR]
    out@labels[[BB]][[matName]] <- do.call(override, labelArgs)
  }
  ## user-specified parameters to override MUST include:
  ## - group (eventually block), defaults to 1 (convenient for longitudinal CFA)
  ## - matName (e.g., lambda)
  ## - row and col (integer or character indices)
  ## - replacement (NA or numeric for values, character for labels)

  out
}
##' @rdname measEq.syntax-class
##' @aliases update,measEq.syntax-method
##' @importFrom stats update
##' @export
setMethod("update", "measEq.syntax", updateMeasEqSyntax)



## -----------------------
## Syntax-Writing Function
## -----------------------


##' Syntax for measurement equivalence
##'
##' Automatically generates \code{lavaan} model syntax to specify a confirmatory
##' factor analysis (CFA) model with equality constraints imposed on
##' user-specified measurement (or structural) parameters. Optionally returns
##' the fitted model (if data are provided) representing some chosen level of
##' measurement equivalence/invariance across groups and/or repeated measures.
##'
##' This function is a pedagogical and analytical tool to generate model syntax
##' representing some level of measurement equivalence/invariance across any
##' combination of multiple groups and/or repeated measures. Support is provided
##' for confirmatory factor analysis (CFA) models with simple or complex
##' structure (i.e., cross-loadings and correlated residuals are allowed).
##' For any complexities that exceed the limits of automation, this function is
##' intended to still be useful by providing a means to generate syntax that
##' users can easily edit to accommodate their unique situations.
##'
##' Limited support is provided for bifactor models and higher-order constructs.
##' Because bifactor models have cross-loadings by definition, the option
##' \code{ID.fac = "effects.code"} is unavailable. \code{ID.fac = "UV"} is
##' recommended for bifactor models, but \code{ID.fac = "UL"} is available on
##' the condition that each factor has a unique first indicator in the
##' \code{configural.model}. In order to maintain generality, higher-order
##' factors may include a mix of manifest and latent indicators, but they must
##' therefore require \code{ID.fac = "UL"} to avoid complications with
##' differentiating lower-order vs. higher-order (or mixed-level) factors.
##' The keyword \code{"loadings"} in \code{group.equal} or \code{long.equal}
##' constrains factor loadings of all manifest indicators (including loadings on
##' higher-order factors that also have latent indicators), whereas the keyword
##' \code{"regressions"} constrains factor loadings of latent indicators. Users
##' can edit the model syntax manually to adjust constraints as necessary, or
##' clever use of the \code{group.partial} or \code{long.partial} arguments
##' could make it possible for users to still automated their model syntax.
##' The keyword \code{"intercepts"} constrains the intercepts of all manifest
##' indicators, and the keyword \code{"means"} constrains intercepts and means
##' of all latent common factors, regardless of whether they are latent
##' indicators of higher-order factors.  To test equivalence of lower-order and
##' higher-order intercepts/means in separate steps, the user can either
##' manually edit their generated syntax or conscientiously exploit the
##' \code{group.partial} or \code{long.partial} arguments as necessary.
##'
##' \strong{\code{ID.fac}:} If the \code{configural.model} fixes any (e.g.,
##' the first) factor loadings, the generated syntax object will retain those
##' fixed values. This allows the user to retain additional constraints that
##' might be necessary (e.g., if there are only 1 or 2 indicators). Some methods
##' must be used in conjunction with other settings:
##' \itemize{
##'   \item \code{ID.cat = "Millsap"} requires \code{ID.fac = "UL"} and
##'         \code{parameterization = "theta"}.
##'   \item \code{ID.cat = "LISREL"} requires \code{parameterization = "theta"}.
##'   \item \code{ID.fac = "effects.code"} is unavailable when there are any
##'         cross-loadings.
##' }
##'
##' \strong{\code{ID.cat}:} Wu & Estabrook (2016) recommended constraining
##' thresholds to equality first, and doing so should allow releasing any
##' identification constraints no longer needed. For each \code{ordered}
##' indicator, constraining one threshold to equality will allow the item's
##' intercepts to be estimated in all but the first group or repeated measure.
##' Constraining a second threshold (if applicable) will allow the item's
##' (residual) variance to be estimated in all but the first group or repeated
##' measure. For binary data, there is no independent test of threshold,
##' intercept, or residual-variance equality. Equivalence of thresholds must
##' also be assumed for three-category indicators. These guidelines provide the
##' least restrictive assumptions and tests, and are therefore the default.
##'
##' The default setting in M\emph{plus} is similar to Wu & Estabrook (2016),
##' except that intercepts are always constrained to zero (so they are assumed
##' to be invariant without testing them). Millsap & Tein (2004) recommended
##' \code{parameterization = "theta"} and identified an item's residual variance
##' in all but the first group (or occasion; Liu et al., 2017) by constraining
##' its intercept to zero and one of its thresholds to equality. A second
##' threshold for the reference indicator (so \code{ID.fac = "UL"}) is used to
##' identify the common-factor means in all but the first group/occasion. The
##' LISREL software fixes the first threshold to zero and (if applicable) the
##' second threshold to 1, and assumes any remaining thresholds to be equal
##' across groups / repeated measures; thus, the intercepts are always
##' identified, and residual variances (\code{parameterization = "theta"}) are
##' identified except for binary data, when they are all fixed to one.
##'
##' \strong{Repeated Measures:} If each repeatedly measured factor is measured
##' by the same indicators (specified in the same order in the
##' \code{configural.model}) on each occasion, without any cross-loadings, the
##' user can let \code{longIndNames} be automatically generated. Generic names
##' for the repeatedly measured indicators are created using the name of the
##' repeatedly measured factors (i.e., \code{names(longFacNames)}) and the
##' number of indicators. So the repeatedly measured first indicator
##' (\code{"ind"}) of a longitudinal construct called "factor" would be
##' generated as \code{"._factor_ind.1"}.
##'
##' The same types of parameter can be specified for \code{long.equal} as for
##' \code{group.equal} (see \code{\link[lavaan]{lavOptions}} for a list), except
##' for \code{"residual.covariances"} or \code{"lv.covariances"}. Instead, users
##' can constrain \emph{auto}covariances using keywords \code{"resid.autocov"}
##' or \code{"lv.autocov"}. Note that \code{group.equal = "lv.covariances"} or
##' \code{group.equal = "residual.covariances"} will constrain any
##' autocovariances across groups, along with any other covariances the user
##' specified in the \code{configural.model}. Note also that autocovariances
##' cannot be specified as exceptions in \code{long.partial}, so anything more
##' complex than the \code{auto} argument automatically provides should instead
##' be manually specified in the \code{configural.model}.
##'
##' When users set \code{orthogonal=TRUE} in the \code{configural.model} (e.g.,
##' in bifactor models of repeatedly measured constructs), autocovariances of
##' each repeatedly measured factor will still be freely estimated in the
##' generated syntax.
##'
##' \strong{Missing Data:} If users wish to utilize the \code{\link{auxiliary}}
##' function to automatically include auxiliary variables in conjunction with
##' \code{missing = "FIML"}, they should first generate the hypothesized-model
##' syntax, then submit that syntax as the model to \code{auxiliary()}.
##' If users utilized \code{\link{runMI}} to fit their \code{configural.model}
##' to multiply imputed data, that model can also be passed to the
##' \code{configural.model} argument, and if \code{return.fit = TRUE}, the
##' generated model will be fitted to the multiple imputations.
##'
##' @importFrom lavaan lavInspect lavNames cfa
##'
##' @param configural.model A model with no measurement-invariance constraints
##'   (i.e., representing only configural invariance), unless required for model
##'   identification. \code{configural.model} can be either:
##'   \itemize{
##'     \item lavaan \code{\link[lavaan]{model.syntax}} or a parameter table
##'           (as returned by \code{\link[lavaan]{parTable}}) specifying the
##'           configural model. Using this option, the user can also provide
##'           either raw \code{data} or summary statistics via \code{sample.cov}
##'           and (optionally) \code{sample.mean}. See argument descriptions in
##'           \code{\link[lavaan]{lavaan}}. In order to include thresholds in
##'           the generated syntax, either users must provide raw \code{data},
##'           or the \code{configural.model} syntax must specify all thresholds
##'           (see first example). If raw \code{data} are not provided, the
##'           number of blocks (groups, levels, or combination) must be
##'           indicated using an arbitrary \code{sample.nobs} argument (e.g.,
##'           3 groups could be specified using \code{sample.nobs=rep(1, 3)}).
##'     \item a fitted \code{\linkS4class{lavaan}} model (e.g., as returned by
##'           \code{\link[lavaan]{cfa}}) estimating the configural model
##'   }
##'   Note that the specified or fitted model must not contain any latent
##'   structural parameters (i.e., it must be a CFA model), unless they are
##'   higher-order constructs with latent indicators (i.e., a second-order CFA).
##'
##' @param ... Additional arguments (e.g., \code{data}, \code{ordered}, or
##'   \code{parameterization}) passed to the \code{\link[lavaan]{lavaan}}
##'   function. See also \code{\link[lavaan]{lavOptions}}.
##'
##' @param ID.fac \code{character}. The method for identifying common-factor
##'   variances and (if \code{meanstructure = TRUE}) means. Three methods are
##'   available, which go by different names in the literature:
##'   \itemize{
##'     \item Standardize the common factor (mean = 0, \emph{SD} = 1) by
##'           specifying any of: \code{"std.lv"}, \code{"unit.variance"},
##'           \code{"UV"}, \code{"fixed.factor"},
##'           \code{"fixed-factor"}
##'     \item Choose a reference indicator by specifying any of:
##'           \code{"auto.fix.first"}, \code{"unit.loading"}, \code{"UL"},
##'           \code{"marker"}, \code{"ref"},  \code{"ref.indicator"},
##'           \code{"reference.indicator"}, \code{"reference-indicator"},
##'           \code{"marker.variable"}, \code{"marker-variable"}
##'     \item Apply effects-code constraints to loadings and intercepts by
##'           specifying any of: \code{"FX"}, \code{"EC"}, \code{"effects"},
##'           \code{"effects.coding"}, \code{"effects-coding"},
##'           \code{"effects.code"}, \code{"effects-code"}
##'   }
##'   See Kloessner & Klopp (2019) for details about all three methods.
##'
##' @param ID.cat \code{character}. The method for identifying (residual)
##'   variances and intercepts of latent item-responses underlying any
##'   \code{ordered} indicators. Four methods are available:
##'   \itemize{
##'     \item To follow Wu & Estabrook's (2016) guidelines (default), specify
##'           any of: \code{"Wu.Estabrook.2016"}, \code{"Wu.2016"},
##'           \code{"Wu.Estabrook"}, \code{"Wu"}, \code{"Wu2016"}.
##'     \item To use the default settings of M\emph{plus} and \code{lavaan},
##'           specify any of: \code{"default"}, \code{"Mplus"}, \code{"Muthen"}.
##'           Details provided in Millsap & Tein (2004).
##'     \item To use the constraints recommended by Millsap & Tein (2004; see
##'           also Liu et al., 2017, for the longitudinal case)
##'           specify any of: \code{"millsap"}, \code{"millsap.2004"},
##'           \code{"millsap.tein.2004"}
##'     \item To use the default settings of LISREL, specify \code{"LISREL"}
##'           or \code{"Joreskog"}. Details provided in Millsap & Tein (2004).
##'   }
##'   See \strong{Details} and \strong{References} for more information.
##'
##' @param ID.thr \code{integer}. Only relevant when
##'   \code{ID.cat = "Millsap.Tein.2004"}. Used to indicate which thresholds
##'   should be constrained for identification. The first integer indicates the
##'   threshold used for all indicators, the second integer indicates the
##'   additional threshold constrained for a reference indicator (ignored if
##'   binary).
##'
##' @param group optional \code{character} indicating the name of a grouping
##'   variable. See \code{\link[lavaan]{cfa}}.
##'
##' @param group.equal optional \code{character} vector indicating type(s) of
##'   parameter to equate across groups. Ignored if \code{is.null(group)}.
##'   See \code{\link[lavaan]{lavOptions}}.
##'
##' @param group.partial optional \code{character} vector or a parameter table
##'   indicating exceptions to \code{group.equal} (see
##'   \code{\link[lavaan]{lavOptions}}). Any variables not appearing in the
##'   \code{configural.model} will be ignored, and any parameter constraints
##'   needed for identification (e.g., two thresholds per indicator when
##'   \code{ID.cat = "Millsap"}) will be removed.
##'
##' @param longFacNames optional named \code{list} of \code{character} vectors,
##'   each indicating multiple factors in the model that are actually the same
##'   construct measured repeatedly. See \strong{Details} and \strong{Examples}.
##'
##' @param longIndNames optional named \code{list} of \code{character} vectors,
##'   each indicating multiple indicators in the model that are actually the
##'   same indicator measured repeatedly. See \strong{Details} and
##'   \strong{Examples}.
##'
##' @param long.equal optional \code{character} vector indicating type(s) of
##'   parameter to equate across repeated measures. Ignored if no factors are
##'   indicated as repeatedly measured in \code{longFacNames}.
##'
##' @param long.partial optional \code{character} vector or a parameter table
##'   indicating exceptions to \code{long.equal}. Any longitudinal variable
##'   names not  appearing in \code{names(longFacNames)} or
##'   \code{names(longIndNames)} will be ignored, and any parameter constraints
##'   needed for identification will be removed.
##'
##' @param auto Used to automatically included autocorrelated measurement errors
##'   among repeatedly measured indicators in \code{longIndNames}. Specify a
##'   single \code{integer} to set the maximum order (e.g., \code{auto = 1L}
##'   indicates that an indicator's unique factors should only be correlated
##'   between adjacently measured occasions). \code{auto = TRUE} or \code{"all"}
##'   will specify residual covariances among all possible lags per repeatedly
##'   measured indicator in \code{longIndNames}.
##'
##' @param warn,debug \code{logical}. Passed to \code{\link[lavaan]{lavaan}}
##'   and \code{\link[lavaan]{lavParseModelString}}.
##'   See \code{\link[lavaan]{lavOptions}}.
##'
##' @param return.fit \code{logical} indicating whether the generated syntax
##'   should be fitted to the provided \code{data} (or summary statistics, if
##'   provided via \code{sample.cov}). If \code{configural.model} is a fitted
##'   lavaan model, the generated syntax will be fitted using the \code{update}
##'   method (see \code{\linkS4class{lavaan}}), and \dots will be passed to
##'   \code{\link[lavaan]{lavaan}}. If neither data nor a fitted lavaan model
##'   were provided, this must be \code{FALSE}. If \code{TRUE}, the generated
##'   \code{measEq.syntax} object will be included in the \code{lavaan} object's
##'   \code{@@external} slot, accessible by \code{fit@@external$measEq.syntax}.
##'
##' @return By default, an object of class \code{\linkS4class{measEq.syntax}}.
##'   If \code{return.fit = TRUE}, a fitted \code{\link[lavaan]{lavaan}}
##'   model, with the \code{measEq.syntax} object stored in the
##'   \code{@@external} slot, accessible by \code{fit@@external$measEq.syntax}.
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' @seealso  \code{\link{compareFit}}
##'
##' @references
##'   Kloessner, S., & Klopp, E. (2019). Explaining constraint interaction: How
##'   to interpret estimated model parameters under alternative scaling methods.
##'   \emph{Structural Equation Modeling, 26}(1), 143--155.
##'   doi:10.1080/10705511.2018.1517356
##'
##'   Liu, Y., Millsap, R. E., West, S. G., Tein, J.-Y., Tanaka, R., & Grimm,
##'   K. J. (2017). Testing measurement invariance in longitudinal data with
##'   ordered-categorical measures. \emph{Psychological Methods, 22}(3),
##'   486--506. doi:10.1037/met0000075
##'
##'   Millsap, R. E., & Tein, J.-Y. (2004). Assessing factorial invariance in
##'   ordered-categorical measures. \emph{Multivariate Behavioral Research, 39}(3),
##'   479--515. doi:10.1207/S15327906MBR3903_4
##'
##'   Wu, H., & Estabrook, R. (2016). Identification of confirmatory factor
##'   analysis models of different levels of invariance for ordered categorical
##'   outcomes. \emph{Psychometrika, 81}(4), 1014--1045.
##'   doi:10.1007/s11336-016-9506-0
##'
##' @examples
##' mod.cat <- ' FU1 =~ u1 + u2 + u3 + u4
##'              FU2 =~ u5 + u6 + u7 + u8 '
##' ## the 2 factors are actually the same factor (FU) measured twice
##' longFacNames <- list(FU = c("FU1","FU2"))
##'
##' ## CONFIGURAL model: no constraints across groups or repeated measures
##' syntax.config <- measEq.syntax(configural.model = mod.cat,
##'                                # NOTE: data provides info about numbers of
##'                                #       groups and thresholds
##'                                data = datCat,
##'                                ordered = paste0("u", 1:8),
##'                                parameterization = "theta",
##'                                ID.fac = "std.lv", ID.cat = "Wu.Estabrook.2016",
##'                                group = "g", longFacNames = longFacNames)
##' ## print lavaan syntax to the Console
##' cat(as.character(syntax.config))
##' ## print a summary of model features
##' summary(syntax.config)
##'
##' ## THRESHOLD invariance:
##' ## only necessary to specify thresholds if you have no data
##' mod.th <- '
##'   u1 | t1 + t2 + t3 + t4
##'   u2 | t1 + t2 + t3 + t4
##'   u3 | t1 + t2 + t3 + t4
##'   u4 | t1 + t2 + t3 + t4
##'   u5 | t1 + t2 + t3 + t4
##'   u6 | t1 + t2 + t3 + t4
##'   u7 | t1 + t2 + t3 + t4
##'   u8 | t1 + t2 + t3 + t4
##' '
##' syntax.thresh <- measEq.syntax(configural.model = c(mod.cat, mod.th),
##'                                # NOTE: data not provided, so syntax must
##'                                #       include thresholds, and number of
##'                                #       groups == 2 is indicated by:
##'                                sample.nobs = c(1, 1),
##'                                parameterization = "theta",
##'                                ID.fac = "std.lv", ID.cat = "Wu.Estabrook.2016",
##'                                group = "g", group.equal = "thresholds",
##'                                longFacNames = longFacNames,
##'                                long.equal = "thresholds")
##' ## notice that constraining 4 thresholds allows intercepts and residual
##' ## variances to be freely estimated in all but the first group & occasion
##' cat(as.character(syntax.thresh))
##' ## print a summary of model features
##' summary(syntax.thresh)
##'
##'
##' ## Fit a model to the data either in a subsequent step (recommended):
##' mod.config <- as.character(syntax.config)
##' fit.config <- cfa(mod.config, data = datCat, group = "g",
##'                   ordered = paste0("u", 1:8), parameterization = "theta")
##' ## or in a single step (not generally recommended):
##' fit.thresh <- measEq.syntax(configural.model = mod.cat, data = datCat,
##'                             ordered = paste0("u", 1:8),
##'                             parameterization = "theta",
##'                             ID.fac = "std.lv", ID.cat = "Wu.Estabrook.2016",
##'                             group = "g", group.equal = "thresholds",
##'                             longFacNames = longFacNames,
##'                             long.equal = "thresholds", return.fit = TRUE)
##' ## compare their fit to test threshold invariance
##' anova(fit.config, fit.thresh)
##'
##'
##' ## --------------------------------------------------------
##' ## RECOMMENDED PRACTICE: fit one invariance model at a time
##' ## --------------------------------------------------------
##'
##' ## - A downside of setting return.fit=TRUE is that if the model has trouble
##' ##   converging, you don't have the opportunity to investigate the syntax,
##' ##   or even to know whether an error resulted from the syntax-generator or
##' ##   from lavaan itself.
##' ## - A downside of automatically fitting an entire set of invariance models
##' ##   (like the old measurementInvariance() function did) is that you might
##' ##   end up testing models that shouldn't even be fitted because less
##' ##   restrictive models already fail (e.g., don't test full scalar
##' ##   invariance if metric invariance fails! Establish partial metric
##' ##   invariance first, then test equivalent of intercepts ONLY among the
##' ##   indicators that have invariate loadings.)
##'
##' ## The recommended sequence is to (1) generate and save each syntax object,
##' ## (2) print it to the screen to verify you are fitting the model you expect
##' ## to (and potentially learn which identification constraints should be
##' ## released when equality constraints are imposed), and (3) fit that model
##' ## to the data, as you would if you had written the syntax yourself.
##'
##' ## Continuing from the examples above, after establishing invariance of
##' ## thresholds, we proceed to test equivalence of loadings and intercepts
##' ##   (metric and scalar invariance, respectively)
##' ## simultaneously across groups and repeated measures.
##'
##' \dontrun{
##'
##' ## metric invariance
##' syntax.metric <- measEq.syntax(configural.model = mod.cat, data = datCat,
##'                                ordered = paste0("u", 1:8),
##'                                parameterization = "theta",
##'                                ID.fac = "std.lv", ID.cat = "Wu.Estabrook.2016",
##'                                group = "g", longFacNames = longFacNames,
##'                                group.equal = c("thresholds","loadings"),
##'                                long.equal  = c("thresholds","loadings"))
##' summary(syntax.metric)                    # summarize model features
##' mod.metric <- as.character(syntax.metric) # save as text
##' cat(mod.metric)                           # print/view lavaan syntax
##' ## fit model to data
##' fit.metric <- cfa(mod.metric, data = datCat, group = "g",
##'                   ordered = paste0("u", 1:8), parameterization = "theta")
##' ## test equivalence of loadings, given equivalence of thresholds
##' anova(fit.thresh, fit.metric)
##'
##' ## scalar invariance
##' syntax.scalar <- measEq.syntax(configural.model = mod.cat, data = datCat,
##'                                ordered = paste0("u", 1:8),
##'                                parameterization = "theta",
##'                                ID.fac = "std.lv", ID.cat = "Wu.Estabrook.2016",
##'                                group = "g", longFacNames = longFacNames,
##'                                group.equal = c("thresholds","loadings",
##'                                                "intercepts"),
##'                                long.equal  = c("thresholds","loadings",
##'                                                "intercepts"))
##' summary(syntax.scalar)                    # summarize model features
##' mod.scalar <- as.character(syntax.scalar) # save as text
##' cat(mod.scalar)                           # print/view lavaan syntax
##' ## fit model to data
##' fit.scalar <- cfa(mod.scalar, data = datCat, group = "g",
##'                   ordered = paste0("u", 1:8), parameterization = "theta")
##' ## test equivalence of intercepts, given equal thresholds & loadings
##' anova(fit.metric, fit.scalar)
##'
##'
##' ## For a single table with all results, you can pass the models to
##' ## summarize to the compareFit() function
##' compareFit(fit.config, fit.thresh, fit.metric, fit.scalar)
##'
##'
##'
##' ## ------------------------------------------------------
##' ## NOT RECOMMENDED: fit several invariance models at once
##' ## ------------------------------------------------------
##' test.seq <- c("thresholds","loadings","intercepts","means","residuals")
##' meq.list <- list()
##' for (i in 0:length(test.seq)) {
##'   if (i == 0L) {
##'     meq.label <- "configural"
##'     group.equal <- ""
##'     long.equal <- ""
##'   } else {
##'     meq.label <- test.seq[i]
##'     group.equal <- test.seq[1:i]
##'     long.equal <- test.seq[1:i]
##'   }
##'   meq.list[[meq.label]] <- measEq.syntax(configural.model = mod.cat,
##'                                          data = datCat,
##'                                          ordered = paste0("u", 1:8),
##'                                          parameterization = "theta",
##'                                          ID.fac = "std.lv",
##'                                          ID.cat = "Wu.Estabrook.2016",
##'                                          group = "g",
##'                                          group.equal = group.equal,
##'                                          longFacNames = longFacNames,
##'                                          long.equal = long.equal,
##'                                          return.fit = TRUE)
##' }
##'
##' compareFit(meq.list)
##'
##'
##' ## -----------------
##' ## Binary indicators
##' ## -----------------
##'
##' ## borrow example data from Mplus user guide
##' myData <- read.table("http://www.statmodel.com/usersguide/chap5/ex5.16.dat")
##' names(myData) <- c("u1","u2","u3","u4","u5","u6","x1","x2","x3","g")
##' bin.mod <- '
##'   FU1 =~ u1 + u2 + u3
##'   FU2 =~ u4 + u5 + u6
##' '
##' ## Must SIMULTANEOUSLY constrain thresholds, loadings, and intercepts
##' test.seq <- list(strong = c("thresholds","loadings","intercepts"),
##'                  means = "means",
##'                  strict = "residuals")
##' meq.list <- list()
##' for (i in 0:length(test.seq)) {
##'   if (i == 0L) {
##'     meq.label <- "configural"
##'     group.equal <- ""
##'     long.equal <- ""
##'   } else {
##'     meq.label <- names(test.seq)[i]
##'     group.equal <- unlist(test.seq[1:i])
##'     # long.equal <- unlist(test.seq[1:i])
##'   }
##'   meq.list[[meq.label]] <- measEq.syntax(configural.model = bin.mod,
##'                                          data = myData,
##'                                          ordered = paste0("u", 1:6),
##'                                          parameterization = "theta",
##'                                          ID.fac = "std.lv",
##'                                          ID.cat = "Wu.Estabrook.2016",
##'                                          group = "g",
##'                                          group.equal = group.equal,
##'                                          #longFacNames = longFacNames,
##'                                          #long.equal = long.equal,
##'                                          return.fit = TRUE)
##' }
##'
##' compareFit(meq.list)
##'
#TODO: add ternary example? or note to start with EQ thresholds?
##'
##' ## ---------------------
##' ## Multilevel Invariance
##' ## ---------------------
##'
##' ## To test invariance across levels in a MLSEM, specify syntax as though
##' ## you are fitting to 2 groups instead of 2 levels.
##'
##' mlsem <- ' f1 =~ y1 + y2 + y3
##'            f2 =~ y4 + y5 + y6 '
##' ## metric invariance
##' syntax.metric <- measEq.syntax(configural.model = mlsem, meanstructure = TRUE,
##'                                ID.fac = "std.lv", sample.nobs = c(1, 1),
##'                                group = "cluster", group.equal = "loadings")
##' ## by definition, Level-1 means must be zero, so fix them
##' syntax.metric <- update(syntax.metric,
##'                         change.syntax = paste0("y", 1:6, " ~ c(0, NA)*1"))
##' ## save as a character string
##' mod.metric <- as.character(syntax.metric, groups.as.blocks = TRUE)
##' ## convert from multigroup to multilevel
##' mod.metric <- gsub(pattern = "group:", replacement = "level:",
##'                    x = mod.metric, fixed = TRUE)
##' ## fit model to data
##' fit.metric <- lavaan(mod.metric, data = Demo.twolevel, cluster = "cluster")
##' summary(fit.metric)
##' }
##' @export
measEq.syntax <- function(configural.model, ..., ID.fac = "std.lv",
                          ID.cat = "Wu.Estabrook.2016", ID.thr = c(1L, 2L),
                          group = NULL, group.equal = "", group.partial = "",
                          longFacNames = list(), longIndNames = list(),
                          long.equal = "", long.partial = "", auto = "all",
                          warn = TRUE, debug = FALSE, return.fit = FALSE) {

  mc <- match.call(expand.dots = TRUE)
  ## evaluate promises that might change before being evaluated
  ## (e.g., for-loops or Monte Carlo studies)
  mc$ID.fac        <- eval(ID.fac)
  mc$ID.cat        <- eval(ID.cat)
  mc$ID.thr        <- eval(ID.thr)
  mc$group         <- eval(group)
  mc$group.equal   <- eval(group.equal)
  mc$group.partial <- eval(group.partial)
  mc$longFacNames  <- eval(longFacNames)
  mc$longIndNames  <- eval(longIndNames)
  mc$long.equal    <- eval(long.equal)
  mc$long.partial  <- eval(long.partial)
  mc$auto          <- eval(auto)

  ## -------------------------------
  ## Preliminary checks on arguments
  ## -------------------------------

  ## check identification arguments
  ID.fac <- tolower(as.character(ID.fac)[1])
  if (ID.fac %in% c("std.lv","unit.variance","uv",
                    "fixed.factor","fixed-factor")) {
    ID.fac <- "uv"
    mc$ID.fac <- "uv"
  } else if (ID.fac %in% c("auto.fix.first","unit.loading","ul","marker","ref",
                           "marker.variable","marker-variable","ref.indicator",
                           "reference.indicator","reference-indicator")) {
    ID.fac <- "ul"
    mc$ID.fac <- "ul"
  } else if (ID.fac %in% c("fx","ec","effects","effects.coding",
                           "effects-coding","effects.code","effects-code")) {
    ID.fac <- "fx"
    mc$ID.fac <- "fx"
  } else stop('Invalid choice for argument: ID.fac = "', ID.fac, '"')

  ID.cat <- tolower(as.character(ID.cat)[1])
  if (ID.cat %in% c("wu.estabrook.2016","wu.2016","wu.estabrook","wu","wu2016")) {
    ID.cat <- "wu"
    mc$ID.cat <- "wu"
  } else if (ID.cat %in% c("millsap","millsap.2004","millsap.tein.2004")) {
    ID.cat <- "millsap"
    mc$ID.cat <- "millsap"
  } else if (ID.cat %in% c("default","mplus","muthen")) {
    ID.cat <- "mplus"
    mc$ID.cat <- "mplus"
  } else if (ID.cat %in% c("joreskog","lisrel")) {
    ID.cat <- "lisrel"
    mc$ID.cat <- "lisrel"
  } else stop('Invalid choice for argument: ID.cat = "', ID.cat, '"')

  ## pass arguments to lavaan
  dots <- list(...)
  dots$debug <- debug
  dots$warn <- warn
  dots$group <- group
  ## check lavaan arguments
  if (!is.null(dots$model)) stop('A model should be specified only with the ',
                                 '"configural.model=" argument, not "model=".')
  if (is.null(dots$meanstructure)) {
    constrMeanStr <- c("intercepts","means") %in% c(group.equal, long.equal)
    if (is.null(dots$sample.mean) && is.null(dots$data) && !any(constrMeanStr)) {
      dots$meanstructure <- FALSE
      mc$meanstructure <- FALSE
    } else {
      dots$meanstructure <- TRUE
      mc$meanstructure <- TRUE
    }
  }

  ## lavaan template from configural model
  if (inherits(configural.model, c("lavaan","lavaanList"))) {
    lavTemplate <- configural.model
    ## check that first loading is not constrained unless ID.fac == "ul"
    if (ID.fac != "ul" && lavInspect(lavTemplate, "options")$auto.fix.first) {
      stop('The "configural.model" argument is a lavaan model fitted using ',
           'auto.fix.first=TRUE (or std.lv=FALSE), which conflicts with the ',
           'requested "ID.fac" method. To generate syntax using the fixed-',
           'factor or effects-coding method of identification, set std.lv=TRUE',
           ' to prevent initial loadings from being fixed to 1 in the syntax.')
    }
    ## check that if (!meanstructure), not set TRUE in call
    if (!is.null(mc$meanstructure)) {
      if (!lavInspect(lavTemplate, "options")$meanstructure && mc$meanstructure)
        stop('Request for meanstructure=TRUE requires configural.model to be ',
             'fitted with meanstructure=TRUE')
    }
  } else {
    lavArgs <- dots
    if (ID.fac != "ul") lavArgs$std.lv <- TRUE
    lavArgs$model <- configural.model # let lavaan() do its own checks
    lavArgs$do.fit <- FALSE
    lavTemplate <- do.call("cfa", lavArgs) #FIXME: violates NAMESPACE rules?  Import cfa()?
    mc$meanstructure <- lavInspect(lavTemplate, "options")$meanstructure # just in case
    mc$configural.model <- lavTemplate
  }


  ## prevent inconsistency
  if (lavInspect(lavTemplate, "options")$categorical &&
      ID.cat %in% c("wu","mplus") &&
      ID.fac != "uv") warning('For factors measured only by categorical ',
                              'indicators, constraints on intercepts are ',
                              'insufficient to identify latent means when the ',
                              'intercepts are already fixed to zero in order ',
                              'to identify latent item scales.  To prevent',
                              'underidentified models, it is recommended to ',
                              'instead set ID.fac = "std.lv".')


  ## convert *.partial strings to parTables
  if (is.character(group.partial)) {
    if (group.partial == "" && length(group.partial) == 1L) {
      group.partial <- data.frame(stringsAsFactors = FALSE, lhs = character(0),
                                  op = character(0), rhs = character(0))
    } else {
      group.partial <- lavaan::lavParseModelString(group.partial,
                                                   as.data.frame. = TRUE,
                                                   warn = warn, debug = debug)
    }
  } #TODO: else {extract information from a measEq.partial object}
  if (is.character(long.partial)) {
    if (long.partial == "" && length(long.partial) == 1L) {
      long.partial <- data.frame(stringsAsFactors = FALSE, lhs = character(0),
                                 op = character(0), rhs = character(0))
    } else {
      long.partial <- lavaan::lavParseModelString(long.partial,
                                                  as.data.frame. = TRUE,
                                                  warn = warn, debug = debug)
    }
  } #TODO: else {extract information from a measEq.partial object}


  ## only relevant when there are longitudinal factors
  if (length(longFacNames) > 0L) {
    if (!is.atomic(auto)) stop("'auto' must be a non-negative integer or the character 'all'.")
    if (is.logical(auto)) { if (auto) auto <- "all" else auto <- 0L}
    if (is.factor(auto)) auto <- as.character(auto)
    if (is.character(auto) && auto != "all")
      stop("'auto' must be a non-negative integer or the character 'all'.")
    if (is.numeric(auto)) {
      auto <- as.integer(auto[1]) # only the first integer
      if (auto < 1L) auto <- NULL
    }
    mc$auto <- auto
  }

  ## extract options and other information
  parameterization <- mc$parameterization
  if (is.null(parameterization)) {
    parameterization <- lavInspect(lavTemplate, "options")$parameterization
  }
  meanstructure <- mc$meanstructure
  if (is.null(meanstructure)) {
    meanstructure <- lavInspect(lavTemplate, "options")$meanstructure
  }
  nG <- lavInspect(lavTemplate, "ngroups")
  ## names of ordinal indicators, number of thresholds for each
  allOrdNames <- lavNames(lavTemplate, type = "ov.ord")
  if (length(allOrdNames)) {
    #TODO: add nThr= argument (named numeric vector?) so data= not required
    nThr <- table(sapply(strsplit(lavNames(lavTemplate, "th"),
                                  split = "|", fixed = TRUE),
                         "[", i = 1))
  } else nThr <- numeric(0)

  if (length(allOrdNames) && ID.cat == "millsap") {
    ## Check for ID.thr
    if (is.numeric(ID.thr)) {
      if (length(ID.thr) == 1L) ID.thr <- rep(ID.thr, 2)
      ID.thr <- sapply(allOrdNames, function(x) ID.thr[1:2], simplify = FALSE)
    } else if (is.list(ID.thr)) {
      if (length((setdiff(allOrdNames, names(ID.thr)))))
        stop('If the same thresholds will not be used for all ordered indicators,',
             ' then "ID.thr" must specify 2 integers per ordered indicator in ',
             'a named list (using names of ordered indicators).')
    }

    ## check identification methods
    ID.fac <- "ul"
    if (parameterization != "theta") stop('If ID.cat == "millsap", you must ',
                                          'use parameterization = "theta"')
  }

  if (length(allOrdNames) && ID.cat == "lisrel") {
    if (parameterization != "theta") stop('If ID.cat == "lisrel", you must ',
                                          'use parameterization = "theta"')
    ## thresholds must be constrained to equality
    if (!"thresholds" %in% group.equal) group.equal <- c("thresholds", group.equal)
    if (!"thresholds" %in% long.equal) long.equal <- c("thresholds", long.equal)
    ## so remove any thresholds from *.partial
    partial.thr <- group.partial$op == "|"
    if (any(partial.thr)) group.partial <- group.partial[!partial.thr, ]
    partial.thr <- long.partial$op == "|"
    if (any(partial.thr)) long.partial <- long.partial[!partial.thr, ]
  }

  if (!meanstructure) {
    ## make sure *.equal includes no mean-structure parameters
    eq.means <- which(group.equal %in% c("means","intercepts"))
    if (length(eq.means)) group.equal <- group.equal[-eq.means]
    eq.means <- which(long.equal %in% c("means","intercepts"))
    if (length(eq.means)) long.equal <- long.equal[-eq.means]

    ## make sure *.partial includes no mean-structure parameters
    partial.means <- group.partial$op == "~1"
    if (any(partial.means)) group.partial <- group.partial[!partial.means, ]
    partial.means <- long.partial$op == "~1"
    if (any(partial.means)) long.partial <- long.partial[!partial.means, ]
  }
  mc$group.partial <- group.partial[c("lhs","op","rhs")] #FIXME: any more? "block" for multilevel?
  mc$long.partial  <- long.partial[c("lhs","op","rhs")]

  ## check logic of constraints
  if (length(allOrdNames) && parameterization == "delta") {
    if ("residuals" %in% long.equal) {
      stop('Residual variances cannot be tested for invariance ',
           'across repeated measures when parameterization =  "delta". \n',
           'Please set parameterization = "theta". \n')
    }
    if ("residuals" %in% group.equal) {
      stop('Residual variances cannot be tested for invariance ',
           'across groups when parameterization = "delta". \n',
           'Please set parameterization = "theta". \n')
    }
  }
  if (warn) {
    if (any(c("lv.variances","lv.autocov") %in% long.equal) && !"loadings" %in% long.equal)
      warning('Latent (co)variances are not comparable over repeated measures ',
              'if their respective factor loadings are not equal ',
              'over repeated measures.')
    if (any(c("lv.variances","lv.covariances") %in% group.equal) && !"loadings" %in% group.equal)
      warning('Latent (co)variances are not comparable across groups ',
              'if their respective factor loadings are not equal across groups.')

    if ("intercepts" %in% long.equal && !"loadings" %in% long.equal)
      warning('Indicator intercepts are not comparable over repeated measures ',
              'if their respective factor loadings are not equal ',
              'over repeated measures.')
    if ("intercepts" %in% group.equal && !"loadings" %in% group.equal)
      warning('Indicator intercepts are not comparable over across groups ',
              'if their respective factor loadings are not equal across groups.')

    if ("means" %in% long.equal && !all(c("loadings","intercepts") %in% long.equal))
      warning('Latent means are not comparable over repeated measures if their ',
              'respective factor loadings and intercepts are not equal ',
              'over repeated measures.')
    if ("means" %in% group.equal && !all(c("loadings","intercepts") %in% group.equal))
      warning('Latent means are not comparable across groups if their ',
              'respective factor loadings and intercepts are not equal ',
              'across groups.')

    if ("resid.autocov" %in% long.equal && !"residuals" %in% long.equal)
      warning('Residual auto-covariances might not be comparable over repeated ',
              'measures if their respective residual variances are not equal ',
              'over repeated measures.')
    if ("residual.covariances" %in% group.equal && !"residuals" %in% group.equal)
      warning('Residual covariances might not be comparable across groups if ',
              'their respective residual variances are not equal across groups.')

    if ("lv.autocov" %in% long.equal && !"lv.variances" %in% long.equal)
      warning('Latent auto-covariances might not be comparable over repeated ',
              'measures if their respective latent variances are not equal ',
              'over repeated measures.')
    if ("lv.covariances" %in% group.equal && !"lv.variances" %in% group.equal)
      warning('Latent covariances might not be comparable across groups if ',
              'their respective latent variances are not equal across groups.')
  }



  ## ------------------
  ## Parameter Matrices
  ## ------------------

  ## Parameter matrices used for labels, fixed/free values, and whether to specify
  GLIST.free <- lavInspect(lavTemplate, "free")
  if (nG == 1L) GLIST.free   <- list(`1` = GLIST.free)
  ## only save relevant matrices to specify
  pmats <- intersect(c("tau","lambda","beta",
                       if (meanstructure) "nu" else NULL ,
                       "theta",
                       if (meanstructure) "alpha" else NULL ,
                       "psi",
                       if (length(allOrdNames) && parameterization == "delta") "delta" else NULL),
                     names(GLIST.free[[1]]))
  if ("beta" %in% pmats && ID.fac != "ul") {
    ID.fac <- "ul" #FIXME: could use effects-coding with relative ease?
    mc$ID.fac <- ID.fac
    message('Higher-order factors detected. ID.fac set to "ul".')
  }

  ## matrices with estimates depends on class of model
  if (inherits(lavTemplate, "lavaan")) {
    GLIST.est <- lavInspect(lavTemplate, "est")
    if (nG == 1L) GLIST.est <- list(`1` = GLIST.est)
  } else if (inherits(lavTemplate, "lavaanList")) {
    nn <- names(lavTemplate@Model@GLIST) #FIXME: will @Model continue to exist?

    GLIST.est <- list()
    for (g in 1:nG) {
      GLIST.est[[g]] <- list()
      for (p in pmats) {
        GLIST.est[[g]][[p]] <- lavTemplate@Model@GLIST[[ which(nn == p)[g] ]]
        ## add dimnames to matrices
        dimnames(GLIST.est[[g]][[p]]) <- dimnames(GLIST.free[[g]][[p]])
      }
    }

  }

  for (g in 1:nG) {
    GLIST.est[[g]] <- GLIST.est[[g]][pmats]
    GLIST.free[[g]] <- GLIST.free[[g]][pmats]
    if (g > 1L) {
      ## make sure all groups have the same observed & latent variables
      same.obs <- all(rownames(GLIST.free[[g]]$lambda) == rownames(GLIST.free[[1]]$lambda))
      same.lat <- all(colnames(GLIST.free[[g]]$lambda) == colnames(GLIST.free[[1]]$lambda))
      if (!same.obs) stop('Models contain different observed variables across ',
                          'groups/blocks.  Configural invariance impossible.')
      if (!same.lat) stop('Models contain different latent variables across ',
                          'groups/blocks.  Configural invariance impossible.')
    }
  }
  ## FIXME: check for others? (e.g., test invariance across multiple levels?)



  ## In general, specify if GLIST.free > 0 | (GLIST.free == 0 & GLIST.est != 0)
  ##  - tau   : specify all
  ##  - lambda: specify any nonzero in free + fixed-nonzero (e.g., auto.fix.first)
  ##  - beta  : treat as second-order lambda
  ##  - nu    : specify all
  ##  - theta : specify diagonal (unless delta?) + any nonzero off-diagonal
  ##  - delta : specify ONLY if parameterization == "delta"
  ##  - alpha : specify all
  ##  - psi   : specify all
  GLIST.specify <- sapply(names(GLIST.free), function(g) list())
  for (g in 1:nG) {
    for (p in pmats) {

      ## THRESHOLDS
      if (p == "tau") {
        GLIST.specify[[g]]$tau <- GLIST.free[[g]]$tau == 0
        GLIST.specify[[g]]$tau[ , 1] <- TRUE
      }
      ## LOADINGS
      if (p == "lambda") {
        free.loading <- GLIST.free[[g]]$lambda > 0L
        fixed.nonzero.loading <- GLIST.free[[g]]$lambda == 0L & GLIST.est[[g]]$lambda != 0
        GLIST.specify[[g]]$lambda <- free.loading | fixed.nonzero.loading
      }
      ## SECOND-ORDER LOADINGS
      if (p == "beta") {
        free.loading <- GLIST.free[[g]]$beta > 0L
        fixed.nonzero.loading <- GLIST.free[[g]]$beta == 0L & GLIST.est[[g]]$beta != 0
        GLIST.specify[[g]]$beta <- free.loading | fixed.nonzero.loading
      }
      ## INTERCEPTS
      if (p == "nu") {
        GLIST.specify[[g]]$nu <- GLIST.free[[g]]$nu == 0
        GLIST.specify[[g]]$nu[ , 1] <- TRUE
      }
      ## LATENT MEANS
      if (p == "alpha") {
        GLIST.specify[[g]]$alpha <- GLIST.free[[g]]$alpha == 0
        GLIST.specify[[g]]$alpha[ , 1] <- TRUE
      }
      ## LATENT (CO)VARIANCES
      if (p == "psi") {
        GLIST.specify[[g]]$psi <- matrix(TRUE, nrow = nrow(GLIST.free[[g]]$psi),
                                         ncol = ncol(GLIST.free[[g]]$psi),
                                         dimnames = dimnames(GLIST.free[[g]]$psi))
        ## only specify lower triangle
        GLIST.specify[[g]]$psi[upper.tri(GLIST.specify[[g]]$psi)] <- FALSE
      }
      ## RESIDUAL (CO)VARIANCES
      if (p == "theta") {
        free.var <- GLIST.free[[g]]$theta > 0L
        fixed.nonzero.var <- GLIST.free[[g]]$theta == 0L & GLIST.est[[g]]$theta != 0
        GLIST.specify[[g]]$theta <- free.var | fixed.nonzero.var
        ## can't specify for ordinal indicators using delta parameterization
        if (parameterization == "delta")
          diag(GLIST.specify[[g]]$theta)[allOrdNames] <- FALSE
        ## only specify lower triangle
        GLIST.specify[[g]]$theta[upper.tri(GLIST.specify[[g]]$theta)] <- FALSE
      }
      ## SCALING FACTORS (delta parameters for latent item-responses)
      if (p == "delta") {
        GLIST.specify[[g]]$delta <- GLIST.free[[g]]$delta == 1
        GLIST.specify[[g]]$delta[ , 1] <- parameterization == "delta"
      }
      ## end loops
    }
  }
  ## check for any cross-loadings
  #TODO: special check for bifactor models possible? Find factors whose indicators all cross-load...
  anyXload <- FALSE
  for (g in 1:nG) {
    if (any(apply(GLIST.specify[[g]]$lambda, 1, sum) > 1)) anyXload <- TRUE
  }
  ## can the effects-coding identification method be used?
  if (ID.fac == "fx" && anyXload) {
    stop('Effects-coding method of factor identification ',
         '("ID.fac") unavailable in models with cross-loadings.')
  }
  ## Warn about constraining intercepts but not means
  freeMeans <- ("intercepts" %in% group.equal && !("means" %in% group.equal)) ||
               ("intercepts" %in% long.equal  && !("means" %in% long.equal) )
  if (ID.fac == "uv" && anyXload && freeMeans) {
    warning('A factor\'s mean cannot be freed unless it has at least one ',
            'indicator without a cross-loading whose intercept is constrained ',
            'to equality. Use cat(as.character()) to check whether the syntax ',
            'returned by measEq.syntax() must be manually adapted to free the ',
            'necessary latent means.')
  }



  ## If it is estimated in the user's configural model, free it (NA).
  ## If it is specified as fixed but != 0, retain fixed value.
  GLIST.values <- sapply(names(GLIST.free), function(g) list())
  for (g in 1:nG) {
    GLIST.values[[g]] <- mapply(function(est, free) {
      est[free > 0L] <- NA
      est
    }, SIMPLIFY = FALSE, est = GLIST.est[[g]], free = GLIST.free[[g]])

    ## constrain first loadings to 1 and first indicators to 0?
    if (ID.fac == "ul") {

      ## matrix to store whether each indicator is a reference indicator
      lambda.first <- matrix(FALSE, nrow = nrow(GLIST.values[[g]]$lambda),
                             ncol = ncol(GLIST.values[[g]]$lambda),
                             dimnames = dimnames(GLIST.values[[g]]$lambda))
      if ("beta" %in% pmats) {
        beta.first <- matrix(FALSE, nrow = nrow(GLIST.values[[g]]$beta),
                             ncol = ncol(GLIST.values[[g]]$beta),
                             dimnames = dimnames(GLIST.values[[g]]$beta))
      }

      ## loop over factors to constrain loadings to 1
      for (f in colnames(GLIST.values[[g]]$lambda)) {
        ## if any loading(s) is(are) fixed to 1 already, no changes needed
        ones.lambda <- which(GLIST.values[[g]]$lambda[ , f] == 1L)
        if ("beta" %in% pmats) ones.beta <- which(GLIST.values[[g]]$beta[ , f] == 1L)
        any1.lambda <- length(ones.lambda) > 0L
        any1.beta <- if ("beta" %in% pmats) length(ones.beta) > 0L else FALSE
        if (!any1.lambda && !any1.beta) {
          ## If not already indicated, find the first indicator and fix it to 1.
          ## Prioritize latent indicators to be first (in case observed has cross-loading)
          if ("beta" %in% pmats) {
            indicators <- names(which(GLIST.specify[[g]]$beta[ , f]))
            if (length(indicators)) {
              first.indicator <- indicators[1]
            } else first.indicator <- NULL
          } else first.indicator <- NULL
          if (length(first.indicator)) {
            ## only true if ("beta" %in% pmats)
            GLIST.values[[g]]$beta[first.indicator, f] <- 1L
            beta.first[first.indicator, f] <- TRUE
          } else {
            ## no latent indicators, so look in lambda
            indicators <- names(which(GLIST.specify[[g]]$lambda[ , f]))
            first.indicator <- indicators[1] #FIXME: no chance of NA by now, right?
            GLIST.values[[g]]$lambda[first.indicator, f] <- 1L
            lambda.first[first.indicator, f] <- TRUE
          }
          ## otherwise, use first fixed == 1 indicator as the marker variable
        } else if (any1.beta) {
          beta.first[ones.beta[1], f] <- TRUE
        } else if (any1.lambda) {
          lambda.first[ones.lambda[1], f] <- TRUE
        }
      }

      ## loop over indicators to constrain intercepts to zero
      if (meanstructure) {

        ## manifest indicators
        for (i in rownames(GLIST.specify[[g]]$lambda)) {
          ## for the first indicator of a construct, constrain to zero
          first.indicator <- lambda.first[i, ]
          if (sum(first.indicator) > 1L)
            stop('The intercept of indicator "', i, '" can only be fixed to zero ',
                 'in order to identify one latent mean, but it is specified as ',
                 'the first indicator of the following factors:\n\t',
                 paste(names(which(first.indicator)), collapse = ", "), '\n',
                 'Please respecify the model so that each factor has a unique ',
                 'first indicator to use as a reference indicator.')
          if (any(first.indicator)) GLIST.values[[g]]$nu[i, 1] <- 0
        }

        ## latent indicators of higher-order constructs
        if ("beta" %in% pmats) for (i in rownames(GLIST.specify[[g]]$beta)) {
          ## for the first indicator of a construct, constrain to zero
          first.indicator <- beta.first[i, ]
          if (sum(first.indicator) > 1L)
            stop('The intercept of indicator "', i, '" can only be fixed to zero ',
                 'in order to identify one factor mean, but it is specified as ',
                 'the first indicator of the following factors:\n\t',
                 paste(names(which(first.indicator)), collapse = ", "), '\n',
                 'Please respecify the model so that each factor has a unique ',
                 'first indicator to use as a reference indicator.')
          if (any(first.indicator)) {
            GLIST.values[[g]]$alpha[i, 1] <- 0
          } else GLIST.values[[g]]$alpha[i, 1] <- NA
        }

      }

    }

  }


  ## Make labels
  GLIST.labels <- sapply(names(GLIST.free), function(g) list())
  for (g in 1:nG) {
    for (p in pmats) {

      if (p == "tau") {
        ## THRESHOLDS
        GLIST.labels[[g]]$tau <- cbind(gsub(x = rownames(GLIST.free[[g]]$tau),
                                            pattern = "|t", replacement = ".thr",
                                            fixed = TRUE))
        dimnames(GLIST.labels[[g]]$tau) <- dimnames(GLIST.free[[g]]$tau)
      } else {
        ## ANY OTHER PARAMETERS
        GLIST.labels[[g]][[p]] <- matrix("", nrow = nrow(GLIST.free[[g]][[p]]),
                                         ncol = ncol(GLIST.free[[g]][[p]]),
                                         dimnames = dimnames(GLIST.free[[g]][[p]]))
        for (RR in rownames(GLIST.free[[g]][[p]])) {
          for (CC in colnames(GLIST.free[[g]][[p]])) {
            GLIST.labels[[g]][[p]][RR, CC] <- getLabel(GLIST.labels[[g]],
                                                       parMat = p,
                                                       RR = RR, CC = CC)
          }
        }

      }
      ## end loops
    }
    ## no labels for scaling factors (cannot equate, not a measuremet parameter)
    GLIST.labels[[g]]$delta <- NULL
  }



  ## ------------------------------------
  ## Preliminary checks on model and data
  ## ------------------------------------

  ## check longitudinal factor names
  if (!is.list(longFacNames)) stop('"longFacNames" must be a list of character vectors.')
  ## check that no longitudinal factors are only at 1 occasion
  if (length(longFacNames)) longFacNames <- longFacNames[sapply(longFacNames, length) > 1L]
  ## also check longIndNames, and each non-NULL element
  if (!is.list(longIndNames)) stop('"longIndNames" must be a list of character vectors.')
  if (length(longIndNames)) {
    longIndList <- sapply(longIndNames, is.character)
    if (!all(longIndList)) stop('"longIndNames" must be a list of character vectors.')
    ## No problem if any(length == 1L). It just won't be constrained.
  }

  ## names of factors in syntax
  allFacNames <- lapply(GLIST.free, function(x) colnames(x$lambda))
  ## collapse names of longitudinal factors (plus non-longitudinal factors)
  # reducedFacNames <- c(names(longFacNames), setdiff(unlist(allFacNames),
  #                                                   unlist(longFacNames)))

  ## check for longitudinal indicator names, automatically generate if empty
  make.longIndNames <- length(longIndNames) == 0L
  for (f in names(longFacNames)) {
    ## time-specific factor names
    fs <- longFacNames[[f]]
    nT <- length(fs) # number of occasions
    ## get indicators of each
    indNames <- sapply(fs, function(ff) {
      names(which(GLIST.specify[[1]]$lambda[ , ff]))
    }, simplify = FALSE)
    if (make.longIndNames) {
      # check for same number of indicators, match across factors
      nInd <- length(indNames[[1]])
      if (!all(sapply(indNames, length) == nInd))
        stop('The number of indicators for longitudinal factor "', f,
             '" differs across measurement occasions. Please use the ',
             '"longIndNames" argument to specify which longitudinal indicators',
             ' are the same indicator on different occasions of measurement.')
      if (nInd > 0L) for (i in 1:nInd) {
        longIndNames[[paste0("._", f, "_.ind.", i)]] <- sapply(indNames, "[",
                                                               i = i,
                                                               USE.NAMES = FALSE)
      }
    } else {
      ## add unique indicators per factor (omitted from user-specified matches) ## NO LONGER NECESSARY
      # for (i in fs) {
      #   extraIndicators <- setdiff(indNames[[i]], unlist(longIndNames[[f]]))
      #   longIndNames[[f]][extraIndicators] <- extraIndicators
      # }
    }
  }
  ## check none have cross-loadings
  longIndTable <- table(unlist(longIndNames))
  if (any(longIndTable > 1L))
    stop('Some longitudinal indicators define more than one factor:\n  ',
         paste(names(longIndTable[longIndTable > 1L]), collapse = ", "), "\n  ",
         'The "longIndNames=" argument must be explicitly declared.')
  ## check equivalence of data type (ordinal vs. continuous) across time
  longOrdNames <- sapply(longIndNames, "%in%", table = allOrdNames, simplify = FALSE)
  someNotAll <- sapply(longOrdNames, function(i) any(i) & !all(i))
  if (any(someNotAll)) {
    stop('At least one longitudinal indicator is declared as "ordered" on',
         ' at least one, but not every, occasion: \n  ',
         paste(names(which(someNotAll)), collapse = ", "))
  }
  ## check number of thresholds/categories is equivalent across time
  allOrd <- sapply(longOrdNames, all)
  if (length(allOrd)) for (i in which(allOrd)) {
    checkThr <- nThr[ longIndNames[[ names(allOrd)[i] ]] ]
    if (!all(checkThr == checkThr[1]))
      stop('These "ordered" longitudinal indicators do not have the same ',
           'number of thresholds (endorsed categories) on every occasion: \n',
           paste(names(checkThr), collapse = ", "),
           "\nConsider collapsing rarely endorsed categories.")
  }
  ## create a backward-key for finding long(Fac/Ind)Names from variable names
  longFacKey <- rep(names(longFacNames), times = sapply(longFacNames, length))
  names(longFacKey) <- unlist(longFacNames)
  longIndKey <- rep(names(longIndNames), times = sapply(longIndNames, length))
  names(longIndKey) <- unlist(longIndNames)

  mc$longFacNames <- longFacNames
  mc$longIndNames <- longIndNames



  ## -----------------
  ## Apply constraints
  ## -----------------

  ## THRESHOLDS (+ intercept & variance ID constraints for allOrdNames)

  ## longitudinal constraints (one group at a time, but same across groups)
  for (g in 1:nG) {
    ## loop over ordinal indicators
    for (i in allOrdNames) {

      ## when other variables are this same indicator?
      longInds <- names(longIndKey)[ which(longIndKey == longIndKey[i]) ]
      if (length(longInds) == 0L) next

      ## keep track of how many thresholds for the i_th indicator have been
      ## constrained, in case identification constraints can be released
      nEqThr <- 0L

      ## loop over thresholds of the i_th ordinal indicator
      for (th in 1:(nThr[i])) {

        ## (ADD) constraints across repeated measures?
        equate.long <- "thresholds" %in% long.equal
        ## check whether not to equate because it is in long.partial
        partial.th <- long.partial$op == "|" & long.partial$rhs == paste0("t", th)
        if (equate.long && any(partial.th)) {
          partial.inds <- longIndNames[[ long.partial$lhs[which(partial.th)] ]]
          equate.long <- !i %in% partial.inds
        }

        ## check whether to equate for identification (overrides *.partial)
        if (ID.cat == "millsap") {
          ## always equate the first (or only, if binary)
          if (th == ID.thr[[i]][1]) {
            equate.long <- TRUE
            ## remove this from long.partial, if necessary
            rm.th <- which(long.partial$lhs == longIndKey[i] & partial.th)
            if (length(rm.th)) long.partial <- long.partial[-rm.th, ]
          }
          ## for the first indicator of a construct, equate the second
          fs <- which(GLIST.specify[[g]]$lambda[i, ])
          first.indicator <- sapply(fs, function(f) {
            lams <- GLIST.specify[[g]]$lambda[ , f]
            lam.eq.1 <- which(GLIST.values[[g]]$lambda[ , f] == 1)
            if (length(lam.eq.1)) return(names(lams[ lam.eq.1[1] ]) == i)
            names(which(lams))[1] == i
          })
          if (th == ID.thr[[i]][2] && any(first.indicator)) {
            equate.long <- TRUE
            if (length(fs) > 1L && warn)
              warning('Millsap & Tein`s (2004) identification constraints might ',
                      'not be optimal when the reference indicator ("', i,
                      '") has a cross-loading (on factors "',
                      paste0(names(fs), collapse = '", "'), '")')
            ## remove this from long.partial, if necessary
            rm.th <- which(long.partial$lhs == longIndKey[i] & partial.th)
            if (length(rm.th)) long.partial <- long.partial[-rm.th, ]
          }
        }

        ## apply longitudinal constraint?
        if (equate.long) {

          ## iterate count of constrained thresholds
          nEqThr <- nEqThr + 1L

          ## apply longitudinal constraint
          this.th <- paste0(i, "|t", th)
          first.th <- paste0(longInds[1], "|t", th)
          GLIST.labels[[g]]$tau[this.th, 1] <- GLIST.labels[[g]]$tau[first.th, 1]

        }
      } ## end loop over thresholds

      ## check whether enough thresholds were equated to free
      ## IDENTIFICATION CONSTRAINTS on intercepts & residuals
      equate.int <- "intercepts" %in% long.equal &&
        !any(long.partial$lhs == longIndKey[i] & long.partial$op == "~1")
      equate.resid <- "residuals" %in% long.equal &&
        !any(long.partial$lhs == longIndKey[i] &
             long.partial$rhs == longIndKey[i] &
             long.partial$op == "~~") #FIXME: leave resid==0 for reference indicators

      if (i == longInds[1]) {

        if (ID.cat == "lisrel") {
          ## always estimate intercepts, and variances unless binary
          GLIST.values[[g]]$nu[i, 1] <- NA
          diag(GLIST.values[[g]]$theta)[i] <- if (nThr[i] == 1L) 1 else NA
        } else {
          ## always set reference occasion's intercepts to 0 and variances to 1
          GLIST.values[[g]]$nu[i, 1] <- 0
          if (parameterization == "theta") {
            diag(GLIST.values[[g]]$theta)[i] <- 1
          } else {
            GLIST.values[[g]]$delta[i, 1] <- 1
          }
        }

      } else if (ID.cat == "wu") {

        ## priority to freeing intercepts
        if (nEqThr == 0L || equate.int) {
          GLIST.values[[g]]$nu[i, 1] <- 0
        } else GLIST.values[[g]]$nu[i, 1] <- NA

        if (nEqThr == 0L || (nEqThr < 2L && !equate.int) || equate.resid) {
          ## keep (residual) variances fixed
          if (parameterization == "theta") {
            diag(GLIST.values[[g]]$theta)[i] <- 1
          } else {
            GLIST.values[[g]]$delta[i, 1] <- 1
          }
        } else {
          ## free (residual) variances
          if (parameterization == "theta") {
            diag(GLIST.values[[g]]$theta)[i] <- NA
          } else {
            GLIST.values[[g]]$delta[i, 1] <- NA
          }
        }

      } else if (ID.cat %in% c("mplus","millsap")) {
        ## never free intercepts, only variances

        if (nEqThr == 0L || equate.resid) {
          ## keep (residual) variances fixed
          if (parameterization == "theta") {
            diag(GLIST.values[[g]]$theta)[i] <- 1
          } else {
            GLIST.values[[g]]$delta[i, 1] <- 1
          }
        } else {
          ## free (residual) variances
          if (parameterization == "theta") {
            diag(GLIST.values[[g]]$theta)[i] <- NA
          } else {
            GLIST.values[[g]]$delta[i, 1] <- NA
          }
        }

      } else if (ID.cat == "lisrel") {
        ## always estimate intercepts, and variances unless binary
        GLIST.values[[g]]$nu[i, 1] <- NA
        diag(GLIST.values[[g]]$theta)[i] <- if (nThr[i] == 1L) 1 else NA
      }

    }
  }
  ## group constraints
  if (nG == 1L && ID.cat == "lisrel") {
    ## Single-group model for repeated measures:
    ## Longitudinal loop above only places LISREL equality constraints on
    ## thresholds. Here, still neeed to fix the first 2 == {0, 1}.

    ## loop over ordinal indicators
    for (i in allOrdNames) {
      ## loop over thresholds of the i_th ordinal indicator
      for (th in 1:(nThr[i])) {

        ## always fix the first (or only, if binary) to zero
        if (th == 1L) GLIST.values[[g]]$tau[paste0(i, "|t", th), 1] <- 0
        ## always fix the second to one
        if (th == 2L) GLIST.values[[g]]$tau[paste0(i, "|t", th), 1] <- 1
        ## estimate any others
        if (th > 2L) GLIST.values[[g]]$tau[paste0(i, "|t", th), 1] <- NA

      } ## end loop over thresholds
    } ## end loop over ordinal indicators

  } else if (nG > 1L) for (g in 1:nG) {
    ## loop over ordinal indicators
    for (i in allOrdNames) {

      ## keep track of how many thresholds for the i_th indicator have
      ## constrained, in case identification constraints can be released
      nEqThr <- 0L

      ## loop over thresholds of the i_th ordinal indicator
      for (th in 1:(nThr[i])) {

        ## (REMOVE) constraints across groups?
        equate.group <- "thresholds" %in% group.equal
        ## check whether not to equate because it is in group.partial
        partial.th <- group.partial$lhs == i & group.partial$op == "|" &
          group.partial$rhs == paste0("t", th)
        if (equate.group) equate.group <- !any(partial.th)

        ## check whether to equate for identification (overrides *.partial)
        if (ID.cat == "millsap") {
          ## always equate the first (or only, if binary)
          if (th == ID.thr[[i]][1]) {
            equate.group <- TRUE
            ## remove this from long.partial, if necessary
            rm.th <- which(partial.th)
            if (length(rm.th)) group.partial <- group.partial[-rm.th, ]
          }
          ## for the first indicator of a construct, equate the second
          fs <- which(GLIST.specify[[g]]$lambda[i, ])
          first.indicator <- sapply(fs, function(f) {
            lams <- GLIST.specify[[g]]$lambda[ , f]
            lam.eq.1 <- which(GLIST.values[[g]]$lambda[ , f] == 1)
            if (length(lam.eq.1)) return(names(lams[ lam.eq.1[1] ]) == i)
            names(which(lams))[1] == i
          })
          if (th == ID.thr[[i]][2] && any(first.indicator)) {
            equate.group <- TRUE
            if (length(fs) > 1L && warn)
              warning('Millsap & Tein`s (2004) identification constraints might ',
                      'not be optimal when the reference indicator ("', i,
                      '") has a cross-loading (on factors "',
                      paste0(names(fs), collapse = '", "'), '")')
            ## remove this from long.partial, if necessary
            rm.th <- which(partial.th)
            if (length(rm.th)) group.partial <- group.partial[-rm.th, ]
          }
        } else if (ID.cat == "lisrel") {
          ## always fix the first (or only, if binary) to zero
          if (th == 1L) GLIST.values[[g]]$tau[paste0(i, "|t", th), 1] <- 0
          ## always fix the second to one
          if (th == 2L) GLIST.values[[g]]$tau[paste0(i, "|t", th), 1] <- 1
          ## estimate any others
          if (th > 2L) GLIST.values[[g]]$tau[paste0(i, "|t", th), 1] <- NA
        }

        ## apply group-specific label, unless constrained
        if (!equate.group) {
          ## row in GLIST
          RR <- paste0(i, "|t", th)
          GLIST.labels[[g]]$tau[RR, 1] <- paste0(GLIST.labels[[g]]$tau[RR, 1], ".g", g)
        } else nEqThr <- nEqThr + 1L # iterate count of constrained thresholds

      } ## end loop over thresholds

      ## check whether enough thresholds were equated to free
      ## IDENTIFICATION CONSTRAINTS on intercepts & residuals.
      ## Note: Group 1 constraints already set in longitudinal loop, ONLY if
      ##       there are repeated measures identified by longInds.
      ##       Section below only RELEASES constraints.
      ##       DON'T OVERWRITE FREED CONSTRAINTS AFTER TIME 1.
      equate.int <- "intercepts" %in% group.equal &&
        !any(group.partial$lhs == i & group.partial$op == "~1")
      equate.resid <- "residuals" %in% group.equal &&
        !any(group.partial$lhs == i & group.partial$rhs == i & group.partial$op == "~~")

      if (g > 1L && ID.cat == "wu") {

        ## priority to freeing intercepts
        #FIXME: binary indicators, latent mean arbitrarily freed, nesting problems
        if (nEqThr >= 1L && !equate.int) GLIST.values[[g]]$nu[i, 1] <- NA

        if ((nEqThr >= 2L || (nEqThr >= 1L && equate.int)) && !equate.resid) {
          ## free (residual) variances
          if (parameterization == "theta") {
            diag(GLIST.values[[g]]$theta)[i] <- NA
          } else {
            GLIST.values[[g]]$delta[i, 1] <- NA
          }
        }

      } else if (g > 1L && ID.cat %in% c("mplus","millsap")) {
        ## never free intercepts, only variances

        if (nEqThr >= 1L && !equate.resid) {
          ## free (residual) variances
          if (parameterization == "theta") {
            diag(GLIST.values[[g]]$theta)[i] <- NA
          } else {
            GLIST.values[[g]]$delta[i, 1] <- NA
          }
        }

      } else if (ID.cat == "lisrel") {
        ## always estimate intercepts, and variances unless binary
        GLIST.values[[g]]$nu[i, 1] <- NA
        diag(GLIST.values[[g]]$theta)[i] <- if (nThr[i] == 1L) 1 else NA
      }

    }
  }


  ## LATENT MEANS

  ## longitudinal constraints (one group at a time, but same across groups)
  if (meanstructure) for (g in 1:nG) {

    ## fix or free factor means?
    if (ID.fac == "uv") {
      GLIST.values[[g]]$alpha[ , 1] <- 0 # free below, if any loading is constrained
      ## freed when any loading is constrained to equality
    } else if ("beta" %in% pmats) {
      ## latent indicators of any higher-order factors already set to 0 or NA
      ## in GLIST.values loop above
    } else GLIST.values[[g]]$alpha[ , 1] <- NA

    ## loop over factors
    for (f in rownames(GLIST.labels[[g]]$alpha)) {

      ## which other variables are this same factor?
      longFacs <- names(longFacKey)[ which(longFacKey == longFacKey[f]) ]
      if (length(longFacs) == 0L) {
        ## not a longitudinal factor, set first group's mean to 0 for Millsap
        if (ID.cat == "millsap" && g == 1L) GLIST.values[[g]]$alpha[f, 1] <- 0
        next
      }
      ## first time a factor is measured, set first group's mean to 0 for Millsap
      if (ID.cat == "millsap" && g == 1L && longFacs[1] == f) {
        GLIST.values[[g]]$alpha[f, 1] <- 0
      }

      ## assign labels
      equate.means <- "means" %in% long.equal &&
        !any(long.partial$lhs == longFacKey[f] & long.partial$op  == "~1")
      if (equate.means) {
        GLIST.labels[[g]]$alpha[f, 1] <- GLIST.labels[[g]]$alpha[longFacs[1], 1]
      }

    }
  }
  ## group constraints
  if (meanstructure && nG > 1L) for (g in 1:nG) {

    ## loop over factors
    for (f in rownames(GLIST.labels[[g]]$alpha)) {

      ## assign labels
      equate.means <- "means" %in% group.equal &&
        !any(group.partial$lhs == f & group.partial$op  == "~1")
      if (!equate.means) {
        GLIST.labels[[g]]$alpha[f, 1] <- paste0(GLIST.labels[[g]]$alpha[f, 1], ".g", g)
      }

    }
  }


  ## LATENT VARIANCES

  ## longitudinal constraints (one group at a time, but same across groups)
  for (g in 1:nG) {

    ## fix or free factor variances?
    if (ID.fac == "uv") {
      diag(GLIST.values[[g]]$psi) <- 1 # free below, if any loading is constrained
      ## freed when any loading is constrained to equality
    } else diag(GLIST.values[[g]]$psi) <- NA

    ## loop over factors
    for (f in colnames(GLIST.labels[[g]]$lambda)) {

      ## which other variables are this same factor?
      longFacs <- names(longFacKey)[ which(longFacKey == longFacKey[f]) ]
      if (length(longFacs) == 0L) next

      ## assign labels
      equate.var <- "lv.variances" %in% long.equal &&
        !any(long.partial$lhs == longFacKey[f] &
             long.partial$op  == "~~" &
             long.partial$rhs == longFacKey[f])
      if (equate.var) {
        GLIST.labels[[g]]$psi[f, f] <- GLIST.labels[[g]]$psi[longFacs[1], longFacs[1]]
      }

    }
  }
  ## group constraints
  if (nG > 1L) for (g in 1:nG) {

    ## loop over factors
    for (f in colnames(GLIST.labels[[g]]$lambda)) {

      ## assign labels
      equate.var <- "lv.variances" %in% group.equal &&
        !any(group.partial$lhs == f &
             group.partial$op  == "~~" &
             group.partial$rhs == f)
      if (!equate.var) {
        GLIST.labels[[g]]$psi[f, f] <- paste0(GLIST.labels[[g]]$psi[f, f], ".g", g)
      }

    }
  }



  ## LOADINGS

  ## longitudinal constraints (one group at a time, but same across groups)
  for (g in 1:nG) {

    ## loop over factors
    for (f in colnames(GLIST.labels[[g]]$lambda)) {

      ## which other factors are this same factor?
      longFacs <- names(longFacKey)[ which(longFacKey == longFacKey[f]) ]
      if (length(longFacs) == 0L) next

      ## loop over any manifest indicators within each factor
      for (i in names(which(GLIST.specify[[g]]$lambda[ , f])) ) {

        ## which other variables are this same indicator?
        longInds <- names(longIndKey)[ which(longIndKey == longIndKey[i]) ]
        if (length(longInds) == 0L) next

        ## assign labels
        equate.load <- "loadings" %in% long.equal &&
          !any(long.partial$lhs == longFacKey[f] &
               long.partial$op  == "=~" &
               long.partial$rhs == longIndKey[i])
        if (equate.load) {
          GLIST.labels[[g]]$lambda[i, f] <- GLIST.labels[[g]]$lambda[longInds[1], longFacs[1]]

          ## free factor variance(s) after Time 1
          if (ID.fac == "uv" && f %in% longFacs[-1]) diag(GLIST.values[[g]]$psi)[f] <- NA
        }

      }

      ## loop over any latent indicators within each factor
      if ("beta" %in% pmats) for (i in names(which(GLIST.specify[[g]]$beta[ , f])) ) {

        ## which other factors are this same factor?
        longInds <- names(longFacKey)[ which(longFacKey == longFacKey[i]) ]
        if (length(longInds) == 0L) next

        ## assign labels
        equate.load <- "regressions" %in% long.equal &&
          !any(long.partial$lhs == longFacKey[f] &
               long.partial$op  == "=~" &
               long.partial$rhs == longFacKey[i])
        if (equate.load) {
          GLIST.labels[[g]]$beta[i, f] <- GLIST.labels[[g]]$beta[longInds[1], longFacs[1]]
        }
      }

    }
  }
  ## group constraints
  if (nG > 1L) for (g in 1:nG) {

    ## loop over factors
    for (f in colnames(GLIST.labels[[g]]$lambda)) {

      ## loop over any manifest indicators within each factor
      for (i in names(which(GLIST.specify[[g]]$lambda[ , f])) ) {

        ## assign labels
        equate.load <- "loadings" %in% group.equal &&
          !any(group.partial$lhs == f &
               group.partial$op  == "=~" &
               group.partial$rhs == i)
        if (!equate.load) {
          GLIST.labels[[g]]$lambda[i, f] <- paste0(GLIST.labels[[g]]$lambda[i, f],
                                                   ".g", g)
        } else if (ID.fac == "uv" && g > 1L) {
          ## free factor variance(s) in group(s) other than the first
          diag(GLIST.values[[g]]$psi)[f] <- NA
        }

      }

      ## loop over any latent indicators within each factor
      if ("beta" %in% pmats) for (i in names(which(GLIST.specify[[g]]$beta[ , f])) ) {
        ## assign labels
        equate.load <- "regressions" %in% group.equal &&
          !any(group.partial$lhs == f &
               group.partial$op  == "=~" &
               group.partial$rhs == i)
        if (!equate.load) {
          GLIST.labels[[g]]$beta[i, f] <- paste0(GLIST.labels[[g]]$beta[i, f],
                                                 ".g", g)
        }
      }

    }
  }


  ## INTERCEPTS

  ## longitudinal constraints (one group at a time, but same across groups)
  if (meanstructure) for (g in 1:nG) {
    ## loop over indicators
    for (i in lavNames(lavTemplate, "ov.ind", group = g)) {

      ## when other variables are this same indicator?
      longInds <- names(longIndKey)[ which(longIndKey == longIndKey[i]) ]
      if (length(longInds) == 0L) next

      ## assign labels
      equate.int <- "intercepts" %in% long.equal &&
        !any(long.partial$lhs == longIndKey[i] & long.partial$op == "~1")
      if (equate.int) {
        GLIST.labels[[g]]$nu[i, 1] <- GLIST.labels[[g]]$nu[longInds[1], 1]

        ## free factor mean(s) after Time 1 only if an indicator without a
        ## cross-loading has an equality-constrained intercept
        if (ID.fac == "uv") {
          ## factors this indicator measures
          fs <- colnames(GLIST.specify[[g]]$lambda)[ GLIST.specify[[g]]$lambda[i,] ]
          only.measures.1 <- length(fs) == 1L
          ## name(s) of longitudinal factor(s)
          LFN <- longFacKey[fs]
          not.time.1 <- fs[1] %in% names(which(longFacKey == LFN))[-1]

          if (only.measures.1 && not.time.1) GLIST.values[[g]]$alpha[fs, 1] <- NA
        }

      }

    }
  }
  ## group constraints
  if (meanstructure && nG > 1L) for (g in 1:nG) {
    ## loop over indicators
    for (i in lavNames(lavTemplate, "ov.ind", group = g)) {

      ## assign labels
      equate.int <- "intercepts" %in% group.equal &&
        !any(group.partial$lhs == i & group.partial$op == "~1")
      if (!equate.int) {
        GLIST.labels[[g]]$nu[i, 1] <- paste0(GLIST.labels[[g]]$nu[i, 1], ".g", g)
      } else if (ID.fac == "uv") {
        ## factors this indicator measures
        fs <- colnames(GLIST.specify[[g]]$lambda)[ GLIST.specify[[g]]$lambda[i,] ]
        only.measures.1 <- length(fs) == 1L

        ## free factor mean(s) other than group 1 only if an indicator without a
        ## cross-loading has an equality-constrained intercept
        if (only.measures.1 && g > 1L) GLIST.values[[g]]$alpha[fs, 1] <- NA
      }

    }
  }


  ## RESIDUAL VARIANCES

  ## longitudinal constraints (one group at a time, but same across groups)
  for (g in 1:nG) {
    ## loop over indicators
    for (i in lavNames(lavTemplate, "ov.ind", group = g)) {

      ## when other variables are this same indicator?
      longInds <- names(longIndKey)[ which(longIndKey == longIndKey[i]) ]
      if (length(longInds) == 0L) next

      ## assign labels
      equate.resid <- "residuals" %in% long.equal &&
        !any(long.partial$lhs == longIndKey[i] &
               long.partial$rhs == longIndKey[i] &
               long.partial$op == "~~")
      if (equate.resid) {
        diag(GLIST.labels[[g]]$theta)[i] <- diag(GLIST.labels[[g]]$theta)[ longInds[1] ]
      }

    }
  }
  ## group constraints
  if (nG > 1L) for (g in 1:nG) {
    ## loop over indicators
    for (i in lavNames(lavTemplate, "ov.ind", group = g)) {

      ## assign labels
      equate.resid <- "residuals" %in% group.equal &&
        !any(group.partial$lhs == i & group.partial$rhs == i & group.partial$op == "~~")
      if (!equate.resid) {
        diag(GLIST.labels[[g]]$theta)[i] <- paste0(diag(GLIST.labels[[g]]$theta)[i],
                                                   ".g", g)
      }

    }
  }


  ## RESIDUAL AUTO-COVARIANCES: longitudinal constraints only

  if (length(longIndNames) && !is.null(auto)) for (g in 1:nG) {
    ## loop over longitudinal indicators
    for (i in names(longIndNames)) {
      nn <- longIndNames[[i]]
      nT <- length(nn) # number repeated measures of indicator i
      auto.i <- suppressWarnings(as.integer(auto))[1] # nT can vary over i
      if (auto == "all" | is.na(auto.i)) auto.i <- nT - 1L # max lag
      if (auto.i >= nT | auto.i < 0L ) auto.i <- nT - 1L # max lag

      ## for each lag...
      for (lag in 1:auto.i) {
        for (tt in 1:(nT - lag)) {

          ## sort indices to ensure the lower.tri is always specified, in case
          ## order of longIndNames does not match order in syntax/theta
          nn.idx <- c(which(rownames(GLIST.specify[[g]]$theta) == nn[tt]),
                      which(rownames(GLIST.specify[[g]]$theta) == nn[tt + lag]))
          idx1 <- nn.idx[ which.max(nn.idx) ] # row index
          idx2 <- nn.idx[ which.min(nn.idx) ] # column index

          ## specify and set free
          GLIST.specify[[g]]$theta[idx1, idx2] <- TRUE
          GLIST.values[[g]]$theta[ idx1, idx2] <- NA

          ## constrain to equality across repeated measures?
          if ("resid.autocov" %in% long.equal && tt > 1L) {
            o.idx <- c(which(rownames(GLIST.specify[[g]]$theta) == nn[1]),
                       which(rownames(GLIST.specify[[g]]$theta) == nn[1 + lag]))
            o1 <- o.idx[ which.max(o.idx) ] # row index
            o2 <- o.idx[ which.min(o.idx) ] # column index
            first.label <- GLIST.labels[[g]]$theta[o1, o2]
            GLIST.labels[[g]]$theta[idx1, idx2] <- first.label
          }

        }
      }

    }
  }

  ## group constraints on any RESIDUAL COVARIANCES

  if (nG > 1) for (g in 1:nG) {
    ## add group-specific labels to any off-diagonal GLIST.specify?
    freeTheta <- which(GLIST.specify[[g]]$theta, arr.ind = TRUE)
    offDiag <- freeTheta[ , "row"] > freeTheta[ , "col"]
    if (sum(offDiag) == 0) break # nothing to do

    ## loop over elements that require action
    free.offDiag <- freeTheta[offDiag, , drop = FALSE]
    for (RR in 1:nrow(free.offDiag)) {
      i <- free.offDiag[RR, "row"]
      j <- free.offDiag[RR, "col"]

      ## check group.partial in both directions
      partial.ij <- any(group.partial$lhs == i & group.partial$rhs == j & group.partial$op == "~~")
      partial.ji <- any(group.partial$lhs == j & group.partial$rhs == i & group.partial$op == "~~")
      equate.rescov <- "residual.covariances" %in% group.equal && !any(partial.ij | partial.ji)

      ## assign group-specific labels?
      if (!equate.rescov) {
        GLIST.labels[[g]]$theta[i, j] <- paste0(GLIST.labels[[g]]$theta[i, j],
                                                ".g", g)
      }

    }

  }


  ## LATENT AUTO-COVARIANCES: longitudinal constraints only

  if (length(longFacNames)) for (g in 1:nG) {
    ## loop over longitudinal indicators
    for (i in names(longFacNames)) {
      nn <- longFacNames[[i]]
      nT <- length(nn) # number repeated measures of indicator i

      ## for each lag...
      for (lag in 1:(nT - 1) ) {
        for (tt in 1:(nT - lag) ) {

          ## specify and set free (overwrite possible "orthogonal=TRUE")
          GLIST.specify[[g]]$psi[ nn[tt + lag], nn[tt] ] <- TRUE
          GLIST.values[[g]]$psi[ nn[tt + lag], nn[tt] ] <- NA

          ## constrain to equality across repeated measures?
          if ("lv.autocov" %in% long.equal && tt > 1L) {
            first.label <- GLIST.labels[[g]]$psi[ nn[1 + lag], nn[1] ]
            GLIST.labels[[g]]$psi[ nn[tt + lag], nn[tt] ] <- first.label
          }

        }
      }

    }
  }


  ## group constraints on any LATENT COVARIANCES

  if (nG > 1) for (g in 1:nG) {
    ## add group-specific labels to any off-diagonal GLIST.specify?
    freePsi <- which(GLIST.specify[[g]]$psi, arr.ind = TRUE)
    offDiag <- freePsi[ , "row"] > freePsi[ , "col"]
    if (sum(offDiag) == 0) break # nothing to do

    ## loop over elements that require action
    free.offDiag <- freePsi[offDiag, , drop = FALSE]
    for (RR in 1:nrow(free.offDiag)) {
      i <- free.offDiag[RR, "row"]
      j <- free.offDiag[RR, "col"]

      ## check group.partial in both directions
      partial.ij <- any(group.partial$lhs == i & group.partial$rhs == j & group.partial$op == "~~")
      partial.ji <- any(group.partial$lhs == j & group.partial$rhs == i & group.partial$op == "~~")
      equate.latcov <- "lv.covariances" %in% group.equal && !any(partial.ij | partial.ji)

      ## assign group-specific labels?
      if (!equate.latcov) {
        GLIST.labels[[g]]$psi[i, j] <- paste0(GLIST.labels[[g]]$psi[i, j], ".g", g)
      }

    }

  }


  ## assemble parameter labels for effects-code identification constraints
  fxList <- character(0)
  if (ID.fac == "fx") {

    listLabels.L <- list()
    if (meanstructure) listLabels.I <- list()
    for (g in 1:nG) {

      ## loadings labels
      listLabels.L[[g]] <- sapply(colnames(GLIST.labels[[g]]$lambda), function(f) {
        GLIST.labels[[g]]$lambda[GLIST.specify[[g]]$lambda[ , f], f]
      }, simplify = FALSE)

      ## intercept labels
      if (meanstructure) {
        listLabels.I[[g]] <- sapply(colnames(GLIST.labels[[g]]$lambda), function(f) {
          GLIST.labels[[g]]$nu[GLIST.specify[[g]]$lambda[ , f], 1]
        }, simplify = FALSE)
      }

    }

    ## names of factors measured in each group
    gFacNames <- lapply(listLabels.L, names)

    ## loop over common-factor names
    for (f in unique(unlist(allFacNames))) {
      ## in which groups is this factor measured?
      groups.with.f <- which(sapply(gFacNames, function(gn) f %in% gn))
      ## get the labels used for indicators in each group
      allLabels.L <- lapply(listLabels.L[groups.with.f], "[[", i = f)
      if (meanstructure) allLabels.I <- lapply(listLabels.I[groups.with.f],
                                               "[[", i = f)


      ## one group, one time --> no checks necessary
      if (length(groups.with.f) == 1L && !f %in% names(longFacKey)) {
        fxList <- c(fxList, make.FX.constraint(allLabels.L[[1]], "loadings"))
        if (meanstructure) {
          fxList <- c(fxList, make.FX.constraint(allLabels.I[[1]], "intercepts"))
        }
      }


      ## one group, multiple times
      if (length(groups.with.f) == 1L && f %in% names(longFacKey)) {
        ## this factor's name on all occasions
        LFN <- names(which(longFacKey == longFacKey[f]))

        ## count constraints on loadings across time
        allConstrained <- which(table(unlist(listLabels.L[[1]][LFN])) == length(LFN))
        if (length(allConstrained)) {
          if (f == LFN[1]) {
            fxList <- c(fxList, make.FX.constraint(names(allConstrained), "loadings"))
          }
        } else {
          ## no constraints, each factor gets its own
          fxList <- c(fxList, make.FX.constraint(allLabels.L[[1]], "loadings"))
        }

        ## count constraints on intercepts across time
        if (meanstructure) {
          allConstrained <- which(table(unlist(listLabels.I[[1]][LFN])) == length(LFN))
          if (length(allConstrained)) {
            if (f == LFN[1]) {
              fxList <- c(fxList, make.FX.constraint(names(allConstrained), "intercepts"))
            }
          } else {
            ## no constraints, each factor gets its own
            fxList <- c(fxList, make.FX.constraint(allLabels.I[[1]], "intercepts"))
          }
        }

      }


      ## multiple groups, one time
      if (length(groups.with.f) > 1L && !f %in% names(longFacKey)) {

        ## count constraints on loadings across groups
        allConstrained <- which(table(unlist(allLabels.L)) == length(groups.with.f))
        if (length(allConstrained)) {
          fxList <- c(fxList, make.FX.constraint(names(allConstrained), "loadings"))
        } else {
          ## no constraints, each group gets its own
          for (g in groups.with.f) {
            fxList <- c(fxList, make.FX.constraint(allLabels.L[[g]], "loadings"))
          }
        }

        ## count constraints on intercepts across groups
        if (meanstructure) {
          allConstrained <- which(table(unlist(allLabels.I)) == length(groups.with.f))
          if (length(allConstrained)) {
            fxList <- c(fxList, make.FX.constraint(names(allConstrained), "intercepts"))
          } else {
            ## no constraints, each group gets its own
            for (g in groups.with.f) {
              fxList <- c(fxList, make.FX.constraint(allLabels.I[[g]], "loadings"))
            }
          }
        }

      }


      ## multiple groups, multiple times: Constrain across any/all dimensions?
      if (length(groups.with.f) > 1L && f %in% names(longFacKey)) {

        ## This factor's name on all occasions
        LFN <- names(which(longFacKey == longFacKey[f]))
        ## Number of dimensions (number of groups times number of occasions).
        ## Assumes each occasion was measured in each group.
        nGT <- length(LFN)*length(groups.with.f)

        ## count constraints on loadings across both dimensions
        all.GL.Labels.L <- lapply(LFN, function(ff) {
          lapply(listLabels.L[groups.with.f], "[[", i = ff)
        })
        all.GL.Constrained.L <- which(table(unlist(all.GL.Labels.L)) == nGT)
        if (length(all.GL.Constrained.L)) {
          if (f == LFN[1]) {
            fxList <- c(fxList, make.FX.constraint(names(all.GL.Constrained.L), "loadings"))
          }
        } else {
          if (f == LFN[1])
            warning('No indicators of longitudinal factor "', longFacKey[f],
                    '" have loadings constrained across all groups and all ',
                    'occasions, so the automatically generated syntax applies ',
                    'effects-code identification constraints separately for each',
                    ' occasion and group. If at least 1 loading is constrained ',
                    'across either groups or occasions, the user should save the',
                    ' syntax to manually reduce the number of identification ',
                    'constraints by applying them only to loadings constrained ',
                    'to equality across groups or occasions.') #TODO: update() method
          for (g in groups.with.f) {
            fxList <- c(fxList, make.FX.constraint(allLabels.L[[g]], "loadings"))
          }
        }

        ## count constraints on intercepts across both dimensions
        if (meanstructure) {

          all.GL.Labels.I <- lapply(LFN, function(ff) {
            lapply(listLabels.I[groups.with.f], "[[", i = ff)
          })
          all.GL.Constrained.I <- which(table(unlist(all.GL.Labels.I)) == nGT)
          if (length(all.GL.Constrained.I)) {
            if (f == LFN[1]) {
              fxList <- c(fxList, make.FX.constraint(names(all.GL.Constrained.I), "intercepts"))
            }
          } else {
            if (f == LFN[1])
              warning('No indicators of longitudinal factor "', longFacKey[f],
                      '" have intercepts constrained across all groups and all ',
                      'occasions, so the automatically generated syntax applies ',
                      'effects-code identification constraints separately for each',
                      ' occasion and group. If at least 1 loading is constrained ',
                      'across either groups or occasions, the user should save the',
                      ' syntax to manually reduce the number of identification ',
                      'constraints by applying them only to intercepts constrained ',
                      'to equality across groups or occasions.') #TODO: update() method
            for (g in groups.with.f) {
              fxList <- c(fxList, make.FX.constraint(allLabels.I[[g]], "intercepts"))
            }
          }

        }

      }


    } # end loop over common factors

    #TODO: Implement effects-coding constraints for thresholds?
    #      For each latent item-response, mean(thresholds) == 0, which
    #      identifies intercepts, resolving the problem of effects-coding with
    #      categorical indicators!
    #      (i.e., constraining intercepts that == 0 to average 0 is redundant)
  }



  ## -------------
  ## Return object
  ## -------------

  out <- new("measEq.syntax", package = "lavaan", model.type = "cfa", call = mc,
             meanstructure = meanstructure,
             numeric = lavNames(lavTemplate, "ov.num"),
             ordered = lavNames(lavTemplate, "ov.ord"),
             parameterization = parameterization,
             specify = GLIST.specify,
             values  = GLIST.values,
             labels  = GLIST.labels,
             constraints = fxList,
             updates = list(values = data.frame(NULL),
                            labels = data.frame(NULL)),
             ngroups = nG)

  if (return.fit) {
    if (inherits(configural.model, "lavaan")) {
      fit <- try(lavaan::update(configural.model,
                                model = as.character(out), ...),
                 silent = TRUE)
    } else if (inherits(configural.model, "lavaanList")) {
      configural.model@call$model <- as.character(out)
      configural.model@call$do.fit <- TRUE
      fit <- try(eval(configural.model@call, parent.frame()), silent = TRUE)
    } else {
      lavArgs$model <- as.character(out)
      lavArgs$do.fit <- TRUE
      fit <- try(do.call("cfa", lavArgs), silent = TRUE)
    }

    ## check whether the model could be fit
    if (inherits(fit, "try-error")) {
      warning('The generated model syntax was not successfully fit to the ',
              'data, and generated the following error message(s): \n\n',
              fit[1:length(fit)], "\n",
              "The measEq.syntax object was returned instead.")
    } else {
      fit@external$measEq.syntax <- out # save the syntax to the lavaan(.mi) object
      out <- fit # return the fitted lavaan(.mi) object
    }

  }

  out
}



## ----------------
## Hidden Functions
## ----------------


## function to label a parameter by its location in a parameter matrix
getLabel <- function(GLIST, parMat, RR, CC = 1L) {
  dn <- dimnames(GLIST[[parMat]])
  out <- paste(parMat, which(dn[[1]] == RR), sep = ".")
  if (!parMat %in% c("alpha","nu")) out <- paste(out, which(dn[[2]] == CC),
                                                 sep = "_")
  out
}

## function to assemble a model constraint for effects-code identification
make.FX.constraint <- function(parLabels, param) {
  nCon <- length(parLabels)
  conVal <- if (param == "loadings") nCon else 0 #TODO: algorithm for thresholds
  out <- paste0(parLabels[1], " == ", conVal)
  if (nCon > 1) out <- paste(c(out, parLabels[-1]), collapse = " - ")
  out
}

## function to generate a character vector of lines of syntax (for as.character)
write.lavaan.syntax <- function(pmat, specify, value, label) {
  nG <- length(specify)

  ## LOADINGS
  if (pmat == "lambda") {
    params <- "## LOADINGS:\n"

    for (fac in colnames(specify[[1]])) {
      for (ind in rownames(specify[[1]])) {

        if (!specify[[1]][ind, fac]) next

        if (nG > 1L) {
          params <- c(params,
                      paste0(fac, " =~ c(",
                             paste(sapply(value, "[", i = ind, j = fac),
                                   collapse = ", "),
                             ")*", ind, " + c(",
                             paste(sapply(label, "[", i = ind, j = fac),
                                   collapse = ", "),
                             ")*", ind))
        } else {
          params <- c(params,
                      paste0(fac, " =~ ", value[[1]][ind, fac], "*", ind,
                             " + ", label[[1]][ind, fac], "*", ind))
        }

      }
    }

    return(c(params, ""))
  }

  ## THRESHOLDS
  if (pmat == "tau") {
    params <- sapply(rownames(specify[[1]]), function(th) {
      th.names <- strsplit(th, split = "|", fixed = TRUE)[[1]]
      if (nG > 1L) {
        param <- paste0(th.names[1], " | c(",
                        paste(sapply(value, "[", i = th, j = 1),
                              collapse = ", "),
                        ")*", th.names[2], " + c(",
                        paste(sapply(label, "[", i = th, j = 1),
                              collapse = ", "),
                        ")*", th.names[2])
      } else {
        param <- paste0(th.names[1], " | ", value[[1]][th, 1], "*", th.names[2],
                        " + ", label[[1]][th, 1], "*", th.names[2])
      }
      param
    })

    return(c("## THRESHOLDS:\n", params, ""))
  }

  ## INTERCEPTS or LATENT MEANS
  if (pmat %in% c("nu","alpha")) {
    ## specify all, so no need to check
    params <- sapply(rownames(specify[[1]]), function(x) {
      if (nG > 1L) {
        param <- paste0(x, " ~ c(",
                        paste(sapply(value, "[", i = x, j = 1),
                              collapse = ", "),
                        ")*1 + c(",
                        paste(sapply(label, "[", i = x, j = 1),
                              collapse = ", "),
                        ")*1")
      } else {
        param <- paste0(x, " ~ ", value[[1]][x, 1], "*1 + ",
                        label[[1]][x, 1], "*1")
      }
      param
    })

    if (pmat == "nu") params <- c("## INTERCEPTS:\n", params)
    if (pmat == "alpha") params <- c("## LATENT MEANS/INTERCEPTS:\n", params)

    return(c(params, ""))
  }

  ## SCALING FACTORS (delta)
  if (pmat == "delta") {
    ## specify any?
    spec.delta <- which(specify[[1]][ , 1])
    if (length(spec.delta) == 0L) return(NULL)
    ## if so...
    params <- sapply(names(spec.delta), function(x) {
      if (nG > 1L) {
        param <- paste0(x, " ~*~ c(",
                        paste(sapply(value, "[", i = x, j = 1),
                              collapse = ", "),
                        ")*", x)
      } else {
        param <- paste0(x, " ~*~ ", value[[1]][x, 1], "*", x)
      }
      param
    })

    return(c("## SCALING FACTORS:\n", params, ""))
  }

  ## LATENT or RESIDUAL (CO)VARIANCES
  if (pmat %in% c("theta","psi")) {

    ## do diagonal first, then check for off-diagonal
    spec.vars <- which(diag(specify[[1]]))
    if (pmat == "psi") {
      params <- "## COMMON-FACTOR VARIANCES:\n"
    } else if (pmat == "theta" && length(spec.vars)) {
      params <- "## UNIQUE-FACTOR VARIANCES:\n"
    } else params <- character(0)

    ## variances
    if (length(spec.vars)) {
      params <- c(params,
                  sapply(names(spec.vars), function(x) {
                    if (nG > 1L) {
                      param <- paste0(x, " ~~ c(",
                                      paste(sapply(value, function(j) diag(j)[x]),
                                            collapse = ", "),
                                      ")*", x, " + c(",
                                      paste(sapply(label, function(j) diag(j)[x]),
                                            collapse = ", "),
                                      ")*", x)
                    } else {
                      param <- paste0(x, " ~~ ", diag(value[[1]])[x], "*", x,
                                      " + ", diag(label[[1]])[x], "*", x)
                    }
                    param
                  }))
    }

    ## covariances
    if (any(specify[[1]][lower.tri(specify[[1]], diag = FALSE)])) {
      if (pmat == "psi")   params <- c(params, "\n## COMMON-FACTOR COVARIANCES:\n")
      if (pmat == "theta") params <- c(params, "\n## UNIQUE-FACTOR COVARIANCES:\n")
    }
    nn <- rownames(specify[[1]])
    if (length(nn) > 1L) for (CC in 1:(length(nn) - 1)) {
      for (RR in (CC + 1):length(nn)) {

        if (!specify[[1]][RR, CC]) next

        if (nG > 1L) {
          params <- c(params,
                      paste0(nn[CC], " ~~ c(",
                             paste(sapply(value, "[", i = RR, j = CC),
                                   collapse = ", "),
                             ")*", nn[RR], " + c(",
                             paste(sapply(label, "[", i = RR, j = CC),
                                   collapse = ", "),
                             ")*", nn[RR]))
        } else {
          params <- c(params,
                      paste0(nn[CC], " ~~ ", value[[1]][RR, CC], "*", nn[RR],
                             " + ", label[[1]][RR, CC], "*", nn[RR]))
        }

      }
    }

    return(c(params, ""))
  }

  ## out of options, should never get this far
  invisible(NULL)
}
#TODO: adapt routine to write Mplus MODEL statements and OpenMx RAM commands
write.mplus.syntax <- function(object, group = 1, params = NULL) {
  out <- character()

  pmatList <- intersect(c("lambda","tau","nu", object@parameterization,
                          "alpha","psi"), names(object@specify[[group]]))
  names(pmatList) <- c("loadings","thresholds","intercepts",
                       ifelse(object@parameterization == "delta",
                              "scales", "residuals"),"means","lv.variances")
  ## selected parameter types?
  if (!is.null(params)) {
    requested <- intersect(names(pmatList), params)
    if (!length(requested)) stop('invalid choice: params = c("',
                                 paste(params, collapse = '", "'), '")\n',
                                 'Valid choices include: ',
                                 paste(names(pmatList), collapse = ", "))
    pmatList <- pmatList[requested]
  }


  ## concatenate all latent-variable definitions
  if ("beta" %in% names(object@specify[[group]])) {
    specify.lambda <- rbind(object@specify[[group]]$lambda,
                            object@specify[[group]]$beta)
    values.lambda  <- rbind(object@values[[group]]$lambda,
                            object@values[[group]]$beta)
    labels.lambda  <- rbind(object@labels[[group]]$lambda,
                            object@labels[[group]]$beta)
  } else {
    specify.lambda <- object@specify[[group]]$lambda
    values.lambda  <- object@values[[group]]$lambda
    labels.lambda  <- object@labels[[group]]$lambda
  }


  ## check for @ordered, define latent item-response factors,
  if (length(object@ordered)) {
    out <- c(out, "! Define LATENT ITEM-RESPONSES as factors",
             paste0("LIR", 1:length(object@ordered), " BY ", object@ordered,
                    "@1;  LIR", 1:length(object@ordered), "@0;"))

    for (i in seq_along(object@ordered)) {
      ## update rownames in Lambda
      #FIXME: update names in concatenated Lambda instead?
      idx <- which(rownames(object@specify[[group]]$lambda) == object@ordered[i])
      rownames(specify.lambda)[idx] <- paste0("LIR", i)
      rownames(values.lambda)[ idx] <- paste0("LIR", i)
      rownames(labels.lambda)[ idx] <- paste0("LIR", i)
      ## update rownames in Nu
      idx <- which(rownames(object@specify[[group]]$nu) == object@ordered[i])
      rownames(object@specify[[group]]$nu)[idx] <- paste0("LIR", i)
      rownames(object@values[[group]]$nu)[idx] <- paste0("LIR", i)
      rownames(object@labels[[group]]$nu)[idx] <- paste0("LIR", i)
    }

  }


  out <- c(out, "! FACTOR LOADINGS")
  ## shorten labels
  labels.lambda <- gsub(pattern = "lambda.", replacement = "L",
                        x = labels.lambda, fixed = TRUE)
  labels.lambda <- gsub(pattern = ".g", replacement = "_",
                        x = labels.lambda, fixed = TRUE)
  ## loop over factors
  for (fac in colnames(specify.lambda)) {
    out <- c(out, paste(fac, "BY"))

    ind <- names(which(specify.lambda[ , fac]))
    lastInd <- rev(ind)[1]

    for (i in ind) {
      val <- values.lambda[i, fac]

      out <- c(out, paste0("    ", i,
                           if (is.na(val)) "*" else paste0("@", val),
                           " (", labels.lambda[i, fac],
                           ")", if (i == lastInd) ";" else ""))
    }
  }


  if ("tau" %in% pmatList) {
    out <- c(out, "! THRESHOLDS")
    ## find unique names to shorten labels
    allThrNames <- unique(do.call(c, lapply(object@labels, "[[", i = "tau")))
    ## loop over ordinal indicators, specify set on a single line
    for (i in object@ordered) {
      iThr <- grep(i, rownames(object@labels[[group]]$tau))
      specify <- object@specify[[group]]$tau[iThr, 1] #NOTE: These are now vectors
      values  <- object@values[[ group]]$tau[iThr, 1]
      labels  <- object@labels[[ group]]$tau[iThr, 1]
      ## identify unique parameter number among thresholds (for short labels)
      idx <- integer()
      for (lab in labels) idx <- c(idx, which(allThrNames == lab))

      out <- c(out,
               paste0("[", i, "$", 1:length(iThr),
                      ifelse(is.na(values), "", paste("@", values)),
                      "] (T", idx, ");", collapse = " "))
    }
  }


  ## INDICATOR-LEVEL PARAMETERS

  hasInts <- object@meanstructure
  hasResid <- length(object@numeric) || object@parameterization == "theta"
  hasScales <- length(object@ordered) && object@parameterization == "delta"
  ## assemble comment for this section
  if (sum(hasInts, hasResid, hasScales) == 3L) {
    out <- c(out, "! INDICATOR INTERCEPTS, RESIDUAL VARIANCES, & SCALING FACTORS")
  } else {
    element.names <- c("INTERCEPTS","RESIDUAL VARIANCES","SCALING FACTORS")
    element.tests <- c(hasInts, hasResid, hasScales)
    out <- c(out, paste0("! INDICATOR ", paste(element.names[element.tests],
                                              collapse = " and ")))
  }

  i.nu <- character()
  i.var <- character()
  ## Loop over indicators
  for (i in 1:nrow(object@specify[[group]]$lambda)) {

    LIR <- rownames(specify.lambda)[i] # LIR names
    RR <- rownames(object@specify[[group]]$lambda)[i]

    if (object@meanstructure) {
      ## INTERCEPTS
      N.val  <- object@values[[group]]$nu[LIR, 1]
      ## shorten labels
      N.lab <- gsub(pattern = "nu.", replacement = "N",
                    x = object@labels[[group]]$nu[LIR, 1], fixed = TRUE)
      N.lab <- gsub(pattern = ".g", replacement = "_", x = N.lab, fixed = TRUE)

      i.nu <- c(i.nu, paste0("[", LIR, ifelse(is.na(N.val), yes = "*",
                                            no = paste0("@", N.val)),
                             "] (", N.lab, ");  "))
    }

    if (RR %in% object@ordered && object@parameterization == "delta") {
      ## SCALING FACTORS
      E.val  <- object@values[[group]]$delta[RR, 1]
      E.lab  <- ""

      i.var <- c(i.var, paste0("{", RR, ifelse(is.na(E.val), yes = "*",
                                              no = paste0("@", E.val)), "};"))
    } else {
      ## RESIDUAL VARIANCES
      E.val  <- object@values[[group]]$theta[RR, RR]
      ## shorten labels
      E.lab <- gsub(pattern = "theta.", replacement = "E",
                    x = object@labels[[group]]$theta[RR, RR], fixed = TRUE)
      E.lab <- gsub(pattern = ".g", replacement = "_", x = E.lab, fixed = TRUE)

      i.var <- c(i.var, paste0(RR, ifelse(is.na(E.val), yes = "*",
                                         no = paste0("@", E.val)),
                               " (", E.lab, ");"))
    }

  }
  out <- c(out, paste(i.nu, i.var))


  E.specify <- object@specify[[group]]$theta
  LT <- E.specify & lower.tri(E.specify, diag = FALSE)
  if (any(LT)) {
    out <- c(out, "! RESIDUAL COVARIANCES")

    E.values <- object@values[[group]]$theta
    ## shorten labels
    E.labels <- gsub(pattern = "theta.", replacement = "E",
                     x = object@labels[[group]]$theta, fixed = TRUE)
    E.labels <- gsub(pattern = ".g", replacement = "_",
                     x = E.labels, fixed = TRUE)

    for (CC in 1:(ncol(LT) - 1)) {
      if (!any(LT[ , CC])) next
      if (sum(LT[ , CC]) == 1L) {
        RR <- which(LT[ , CC])
        out <- c(out,
                 paste0(colnames(LT)[CC], " WITH ", rownames(LT)[RR],
                        ifelse(is.na(E.values[RR, CC]), yes = "",
                               no = paste("@", E.values[RR, CC])),
                        " (", E.labels[RR, CC], ");"))
        next
      }
      ## else, there are multiple covariates with LT[CC]
      out <- c(out, paste(colnames(LT)[CC], "WITH"))

      ind <- names(which(LT[ , CC]))
      lastInd <- rev(ind)[1]

      for (RR in ind) {
        val <- E.values[RR, CC]

        out <- c(out, paste0("    ", RR,
                             if (is.na(val)) "" else paste0("@", val),
                             " (", E.labels[RR, CC],
                             ")", if (RR == lastInd) ";" else ""))
      }

    }
  }


  ## FACTOR-LEVEL PARAMETERS
  out <- c(out, paste("! FACTOR",
                      if (object@meanstructure) "INTERCEPTS &" else NULL,
                      "(RESIDUAL) VARIANCES"))
  i.alpha <- character()
  i.psi <- character()
  ## Loop over factors
  for (i in rownames(object@specify[[group]]$psi)) {

    if (object@meanstructure) {
      ## INTERCEPTS
      A.val  <- object@values[[group]]$alpha[i, 1]
      ## shorten labels
      A.lab <- gsub(pattern = "alpha.", replacement = "A",
                    x = object@labels[[group]]$alpha[i, 1], fixed = TRUE)
      A.lab <- gsub(pattern = ".g", replacement = "_", x = A.lab, fixed = TRUE)

      i.alpha <- c(i.alpha, paste0("[", i, ifelse(is.na(A.val), yes = "*",
                                                  no = paste0("@", A.val)),
                                   "] (", A.lab, ");  "))
    }

    ## RESIDUAL VARIANCES
    P.val  <- object@values[[group]]$psi[i, i]
    ## shorten labels
    P.lab <- gsub(pattern = "psi.", replacement = "P",
                  x = object@labels[[group]]$psi[i, i], fixed = TRUE)
    P.lab <- gsub(pattern = ".g", replacement = "_", x = P.lab, fixed = TRUE)

    i.psi <- c(i.psi, paste0(i, ifelse(is.na(P.val), yes = "",
                                        no = paste0("@", P.val)),
                             " (", P.lab, ");"))
  }
  out <- c(out, paste(i.alpha, i.psi))



  P.specify <- object@specify[[group]]$psi
  LT <- P.specify & lower.tri(P.specify, diag = FALSE)
  if (any(LT)) {
    out <- c(out, "! FACTOR COVARIANCES")

    P.values <- object@values[[group]]$psi
    ## shorten labels
    P.labels <- gsub(pattern = "psi.", replacement = "P",
                     x = object@labels[[group]]$psi, fixed = TRUE)
    P.labels <- gsub(pattern = ".g", replacement = "_",
                     x = P.labels, fixed = TRUE)

    for (CC in 1:(ncol(LT) - 1)) {
      if (!any(LT[ , CC])) next
      if (sum(LT[ , CC]) == 1L) {
        RR <- which(LT[ , CC])
        out <- c(out,
                 paste0(colnames(LT)[CC], " WITH ", rownames(LT)[RR],
                        ifelse(is.na(P.values[RR, CC]), yes = "",
                               no = paste("@", P.values[RR, CC])),
                        " (", P.labels[RR, CC], ");"))
        next
      }
      ## else, there are multiple covariates with LT[CC]
      out <- c(out, paste(colnames(LT)[CC], "WITH"))

      ind <- names(which(LT[ , CC]))
      lastInd <- rev(ind)[1]

      for (RR in ind) {
        val <- P.values[RR, CC]

        out <- c(out, paste0("    ", RR,
                             if (is.na(val)) "" else paste0("@", val),
                             " (", P.labels[RR, CC],
                             ")", if (RR == lastInd) ";" else ""))
      }

    }
  }


  ## MODEL CONSTRAINTs
  if (length(object@constraints) && group == object@ngroups) {
    con <- object@constraints
    con <- gsub("lambda.", "L", con)
    con <- gsub("theta.", "E", con)
    con <- gsub("psi.", "P", con)
    if (length(object@ordered)) for (th in object@labels[[group]]$tau[ , 1]) {
      con <- gsub(th, paste0("T", which(allThrNames == th)), con)
    }
    if (object@meanstructure) {
      con <- gsub("nu.", "N", con)
      con <- gsub("alpha.", "A", con)
    }
    con <- gsub(".g", "_", con)
    con <- gsub("==", "=", con)

    out <- c(out, "\nMODEL CONSTRAINT:", paste0(con, ";"))
  }
  #TODO: gsub = for ==, add ";", anything else? object@constraints, "")

  paste(out, collapse = "\n")
}
# write.OpenMx.syntax <- function(pmat, specify, value, label) {}


## function to allow users to customize syntax with update(),
## so they don't necessarily have to copy/paste a script to adapt it.
override <- function(object, slotName = "values", group = 1L,
                     matName, row, col, replacement) {
  stopifnot(inherits(object, "measEq.syntax"))
  MM <- methods::slot(object, slotName)[[group]] # only "values" or "labels"

  ## check indices
  if (is.character(row)) {
    if (! row %in% rownames(MM[[matName]]))
      stop("'", row, "' not found in rownames(",
           deparse(substitute(object)), "@", slotName, "[[", group, "]]$",
           matName, ")")
  } else if (is.numeric(row)) {
    if (! as.integer(row) %in% 1:nrow(MM[[matName]]))
      stop(as.integer(row), "' is outside the number of nrow(",
           deparse(substitute(object)), "@", slotName, "[[", group, "]]$",
           matName, ")")
  } else stop('row argument must be numeric/character indices')
  ## repeat for col
  if (matName %in% c("nu","alpha","delta","tau")) col <- 1L else {
    if (is.character(col)) {
      if (! col %in% colnames(MM[[matName]]))
        stop("'", col, "' not found in colnames(",
             deparse(substitute(object)), "@", slotName, "[[", group, "]]$",
             matName, ")")
    } else if (is.numeric(col)) {
      if (! as.integer(col) %in% 1:ncol(MM[[matName]]))
        stop(as.integer(col), "' is outside the number of ncol(",
             deparse(substitute(object)), "@", slotName, "[[", group, "]]$",
             matName, ")")
    } else stop('col argument must be numeric/character indices')
  }

  newM <- MM[[matName]]
  newM[row, col] <- replacement
  if (matName %in% c("theta","psi")) newM[col, row] <- replacement
  newM
}

## function to assemble values/labels to update
char2update <- function(object, model, return.object = TRUE) {
  stopifnot(inherits(object, "measEq.syntax"))
  stopifnot(inherits(model, "character"))

  PT <- lavaan::lavParseModelString(model, as.data.frame. = TRUE)
  indNames <- lapply(object@values, function(x) rownames(x$lambda)) # per block
  facNames <- lapply(object@values, function(x) colnames(x$lambda))

  values <- PT$fixed
  labels <- PT$label
  ## check for multigroup specification of values/labels
  if (any(grepl(pattern = ";", x = values))) {
    values <- strsplit(values, split = ";")
    nValues <- sapply(values, length)
  } else nValues <- rep(1L, length(values))
  if (any(grepl(pattern = ";", x = labels))) {
    labels <- strsplit(labels, split = ";")
    nLabels <- sapply(labels, length)
  } else nLabels <- rep(1L, length(labels))
  nBlocks <- length(facNames)

  values.DF <- data.frame(NULL)
  labels.DF <- data.frame(NULL)
  for (RR in 1:nrow(PT)) {
    ## check whether numbers match
    if (nValues > 1L && nValues != nBlocks) {
      stop('Number of fixed/free values (', nValues[RR],
           ') specified for parameter "', PT$lhs[RR], PT$op[RR], PT$rhs[RR],
           '" does not match the number of groups (', nBlocks, ')')
    }
    if (nLabels > 1L && nLabels != nBlocks) {
      stop('Number of labels (', nLabels[RR],
           ') specified for parameter "', PT$lhs[RR], PT$op[RR], PT$rhs[RR],
           '" does not match the number of groups (', nBlocks, ')')
    }

    ## loop over blocks (currently only groups)
    for (BB in 1:nBlocks) {
      ## make template for values and labels, depending on parameter matrix

      ## INTERCEPTS
      if (PT$op[RR] == "~1" && PT$lhs[RR] %in% indNames[[BB]]) {
        DF <- data.frame(stringsAsFactors = FALSE, group = BB, matName = "nu",
                         row = PT$lhs[RR], col = "intercept")

        ## LATENT MEANS
      } else if (PT$op[RR] == "~1" && PT$lhs[RR] %in% facNames[[BB]]) {
        DF <- data.frame(stringsAsFactors = FALSE, group = BB, matName = "alpha",
                         row = PT$lhs[RR], col = "intercept")

        ## LOADINGS
      } else if (PT$op[RR] == "=~" && PT$rhs[RR] %in% indNames[[BB]]) {
        DF <- data.frame(stringsAsFactors = FALSE, group = BB, matName = "lambda",
                         row = PT$rhs[RR], col = PT$lhs[RR])

        ## SECOND-ORDER LOADINGS
      } else  if (PT$op[RR] == "=~" && PT$rhs[RR] %in% facNames[[BB]]) {
        DF <- data.frame(stringsAsFactors = FALSE, group = BB, matName = "beta",
                         row = PT$rhs[RR], col = PT$lhs[RR])

        ## LATENT (CO)VARIANCES
      } else if (PT$op[RR] == "~~" && PT$rhs[RR] %in% facNames[[BB]]) {
        DF <- data.frame(stringsAsFactors = FALSE, group = BB, matName = "psi",
                         # symmetry handled in override
                         row = PT$rhs[RR], col = PT$lhs[RR])

        ## RESIDUAL (CO)VARIANCES
      } else if (PT$op[RR] == "~~" && PT$rhs[RR] %in% indNames[[BB]]) {
        DF <- data.frame(stringsAsFactors = FALSE, group = BB, matName = "theta",
                         # symmetry handled in override
                         row = PT$rhs[RR], col = PT$lhs[RR])

        ## THRESHOLDS
      } else if (PT$op[RR] == "|") {
        if (!length(object@ordered)) {
          warning('Thresholds ignored when no indicators are declared as ordered')
        }
        DF <- data.frame(stringsAsFactors = FALSE, group = BB, matName = "tau",
                         row = paste0(PT$lhs[RR], "|", PT$rhs[RR]),
                         col = "threshold")

        ## SCALING FACTORS (delta parameters for latent item-responses)
      } else if (PT$op[RR] == "~*~") {
        if (!length(object@ordered)) {
          warning('Thresholds ignored when no indicators are declared as ordered')
        }
        if (object@parameterization == "theta") {
          warning('Latent-response scales (specified with the "~*~" operator) ',
                  'ignored when parameterization = "theta"')
        }
        if (PT$lhs[RR] != PT$rhs[RR]) {
          warning('Latent-response scales (specified with the "~*~" operator) ',
                  'ignored when left- and right-hand side do not match (',
                  PT$lhs[RR], '~*~', PT$rhs[RR], ')')
          next
        }
        DF <- data.frame(stringsAsFactors = FALSE, group = BB, matName = "delta",
                         row = PT$lhs[RR], col = "scales")
      }

      #FIXME? anything that does not match is simply ignored (no error messages)

      ## change labels?
      if (BB > 1L && nLabels[RR] == 1L) {

        if (labels[[RR]] != "") {
          labels.DF <- rbind(labels.DF, cbind(DF, stringsAsFactors = FALSE,
                                              replacement = labels[[RR]]))
        }

      } else if (labels[[RR]][BB] != "") {
          labels.DF <- rbind(labels.DF, cbind(DF, stringsAsFactors = FALSE,
                                              replacement = labels[[RR]][BB]))
      }

      ## change fixed/free values?
      if (BB > 1L && nValues[RR] == 1L) {

        if (values[[RR]] != "") {
          values.DF <- rbind(values.DF, cbind(DF, stringsAsFactors = FALSE,
                                              replacement = values[[RR]]))
        }

      } else if (values[[RR]][BB] != "") {
        values.DF <- rbind(values.DF, cbind(DF, stringsAsFactors = FALSE,
                                            replacement = values[[RR]][BB]))
      }

    } # end loop over blocks
  }   # end loop over parameters

  ## make sure values are stored as numeric, not character
  if (nrow(values.DF)) {
    suppressWarnings(values.DF$replacement <- as.numeric(values.DF$replacement))
  }

  if (return.object) {
    object@updates$values <- rbind(object@updates$values, values.DF)
    object@updates$labels <- rbind(object@updates$labels, labels.DF)
    return(object)
  }
  ## else return the list of data.frames with updates to make
  list(values = values.DF, labels = labels.DF)
}



