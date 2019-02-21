### Terrence D. Jorgensen
### Last updated: 21 February 2019


##' Random Allocation of Items to Parcels in a Structural Equation Model
##'
##' This function generates a given number of randomly generated item-to-parcel
##' allocations, fits a model to each allocation, and provides averaged results
##' over all allocations.
##'
##' This function implements the random item-to-parcel allocation procedure
##' described in Sterba (2011) and Sterba and MacCallum (2010). The function
##' takes a single data set with item-level data, randomly assigns items to
##' parcels, fits a structural equation model to the parceled data (using
##' \link[lavaan]{lavaanList}), and repeats this process for a user-specified
##' number of random allocations. Results from all fitted models are summarized
##' in the output. For further details on the benefits of randomly allocating
##' items to parcels, see Sterba (2011) and Sterba and MccCallum (2010).
##'
##' @importFrom stats sd qnorm
##' @importFrom lavaan parTable lavInspect lavaanList lavaanify lavNames
##'
##' @param model \code{\link[lavaan]{lavaan}} model syntax specifying the model
##'   fit to (at least some) parceled data. Note that there can be a mixture of
##'   items and parcels (even within the same factor), in case certain items
##'   should never be parceled. Can be a character string or parameter table.
##'   Also see \code{\link[lavaan]{lavaanify}} for more details.
##' @param data A \code{data.frame} containing all observed variables appearing
##'   in the \code{model}, as well as those in the \code{item.syntax} used to
##'   create parcels. If the data have missing values, multiple imputation
##'   before parceling is recommended: submit a stacked data set (with a variable
##'   for the imputation number, so they can be separateed later) and set
##'   \code{do.fit = FALSE} to return the list of \code{data.frame}s (one per
##'   allocation), each of which is a stacked, imputed data set with parcels.
##' @param parcel.names \code{character} vector containing names of all parcels
##' appearing as indicators in \code{model}.
##' @param item.syntax \link[lavaan]{lavaan} model syntax specifying the model
##'   that would be fit to all of the unparceled items, including items that
##'   should be randomly allocated to parcels appearing in \code{model}.
##' @param nAlloc The number of random items-to-parcels allocations to generate.
##' @param fun \code{character} string indicating the name of the
##'   \code{\link[lavaan]{lavaan}} function used to fit \code{model} to
##'   \code{data}. Can only take the values \code{"lavaan"}, \code{"sem"},
##'   \code{"cfa"}, or \code{"growth"}.
##' @param alpha Alpha level used as criterion for significance.
##' @param fit.measures \code{character} vector containing names of fit measures
##'   to request from each fitted \code{\link[lavaan]{lavaan}} model.  See the
##'   output of \code{\link[lavaan]{fitMeasures}} for a list of available measures.
##' @param \dots Additional arguments to be passed to
##'   \code{\link[lavaan]{lavaanList}}. See also \code{\link[lavaan]{lavOptions}}
##' @param show.progress If \code{TRUE}, show a \code{\link[utils]{txtProgressBar}}
##'   indicating how fast the model-fitting iterates over allocations.
##' @param iseed (Optional) Random seed used for parceling items. When the same
##'   random seed is specified and the program is re-run, the same allocations
##'   will be generated. Using the same \code{iseed} argument will ensure the
##'   any model is fit to the same parcel allocations. \emph{Note}: When using
##'   \pkg{parallel} options, you must first type \code{RNGkind("L'Ecuyer-CMRG")}
##'   into the R Console, so that the seed will be controlled across cores.
##' @param do.fit If \code{TRUE} (default), the \code{model} is fitted to each
##'   parceled data set, and the summary of results is returned (see the Value
##'   section below). If \code{FALSE}, the items are randomly parceled, but the
##'   model is not fit; instead, the \code{list} of \code{data.frame}s is
##'   returned (so assign it to an object).
##' @param return.fit If \code{TRUE}, a \code{\link[lavaan]{lavaanList}} object
##'   is returned with the \code{list} of results across allocations
##' @param warn Whether to print warnings when fitting \code{model} to each allocation
##'
##' @return
##'   \item{Estimates}{A \code{data.frame} containing results related to
##'     parameter estimates with columns corresponding to their names; average
##'     and standard deviation across allocations; minimum, maximum, and range
##'     across allocations; and the proportion of allocations in which each
##'     parameter estimate was significant.}
##'   \item{SE}{A \code{data.frame} containing results similar to
##'     \code{Estimates}, but related to the standard errors of parameter
##'     estimates.}
##'   \item{Fit}{A \code{data.frame} containing results related to model fit,
##'     with columns corresponding to fit index names; their average and
##'     standard deviation across allocations; the minimum, maximum, and range
##'     across allocations; and (if the test statistic or RMSEA is included in
##'     \code{fit.measures}) the proportion of allocations in which each
##'     test of (exact or close) fit was significant.}
##'   \item{Model}{A \code{\link[lavaan]{lavaanList}} object containing results
##'     of the \code{model} fitted to each parcel allocation. Only returned if
##'     \code{return.fit = TRUE}.}
##'
##' @author
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @seealso \code{\link{PAVranking}} for comparing 2 models,
##'   \code{\link{poolMAlloc}} for choosing the number of allocations
##'
##' @references
##'
##' Sterba, S. K. (2011). Implications of parcel-allocation
##' variability for comparing fit of item-solutions and parcel-solutions.
##' \emph{Structural Equation Modeling, 18}(4), 554--577.
##' doi:10.1080/10705511.2011.607073
##'
##' Sterba, S. K. & MacCallum, R. C. (2010). Variability in parameter estimates
##' and model fit across random allocations of items to parcels.
##' \emph{Multivariate Behavioral Research, 45}(2), 322--358.
##' doi:10.1080/00273171003680302
##'
##' Sterba, S. K., & Rights, J. D. (2016). Accounting for parcel-allocation
##' variability in practice: Combining sources of uncertainty and choosing the
##' number of allocations. \emph{Multivariate Behavioral Research, 51}(2--3),
##' 296--313. doi:10.1080/00273171.2016.1144502
##'
##' Sterba, S. K., & Rights, J. D. (2017). Effects of parceling on model
##' selection: Parcel-allocation variability in model ranking.
##' \emph{Psychological Methods, 22}(1), 47--68. doi:10.1037/met0000067
##'
##' @examples
##'
##' ## Fit 2-factor CFA to simulated data. Each factor has 9 indicators.
##'
##' ## Specify the item-level model (if NO parcels were created)
##' item.syntax <- c(paste0("f1 =~ f1item", 1:9),
##'                  paste0("f2 =~ f2item", 1:9))
##' cat(item.syntax, sep = "\n")
##' ## Below, we reduce the size of this same model by
##' ## applying different parceling schemes
##'
##'
##' ## 3-indicator parcels
##' mod.parcels <- '
##' f1 =~ par1 + par2 + par3
##' f2 =~ par4 + par5 + par6
##' '
##' ## names of parcels
##' (parcel.names <- paste0("par", 1:6))
##'
##' \dontrun{
##' ## override default random-number generator to use parallel options
##' RNGkind("L'Ecuyer-CMRG")
##'
##' parcelAllocation(mod.parcels, data = simParcel, nAlloc = 100,
##'                  parcel.names = parcel.names, item.syntax = item.syntax,
##'                  std.lv = TRUE,       # any addition lavaan arguments
##'                  parallel = "snow")   # parallel options
##'
##'
##'
##' ## POOL RESULTS by treating parcel allocations as multiple imputations
##' ## Details provided in Sterba & Rights (2016); see ?poolMAlloc.
##'
##' ## save list of data sets instead of fitting model yet
##' dataList <- parcelAllocation(mod.parcels, data = simParcel, nAlloc = 100,
##'                              parcel.names = parcel.names,
##'                              item.syntax = item.syntax,
##'                              do.fit = FALSE)
##' ## now fit the model to each data set
##' fit.parcels <- cfa.mi(mod.parcels, data = dataList, std.lv = TRUE)
##' summary(fit.parcels) # uses Rubin's rules
##' anova(fit.parcels)   # pooled test statistic
##' class?lavaan.mi      # find more methods for pooling results
##' }
##'
##'
##' ## multigroup example
##' simParcel$group <- 0:1 # arbitrary groups for example
##' mod.mg <- '
##' f1 =~ par1 + c(L2, L2)*par2 + par3
##' f2 =~ par4 + par5 + par6
##' '
##' ## names of parcels
##' (parcel.names <- paste0("par", 1:6))
##'
##' parcelAllocation(mod.mg, data = simParcel, parcel.names, item.syntax,
##'                  std.lv = TRUE, group = "group", group.equal = "loadings",
##'                  nAlloc = 20, show.progress = TRUE)
##'
##'
##'
##' ## parcels for first factor, items for second factor
##' mod.items <- '
##' f1 =~ par1 + par2 + par3
##' f2 =~ f2item2 + f2item7 + f2item8
##' '
##' ## names of parcels
##' (parcel.names <- paste0("par", 1:3))
##'
##' parcelAllocation(mod.items, data = simParcel, parcel.names, item.syntax,
##'                  nAlloc = 20, std.lv = TRUE)
##'
##'
##'
##' ## mixture of 1- and 3-indicator parcels for second factor
##' mod.mix <- '
##' f1 =~ par1 + par2 + par3
##' f2 =~ f2item2 + f2item7 + f2item8 + par4 + par5 + par6
##' '
##' ## names of parcels
##' (parcel.names <- paste0("par", 1:6))
##'
##' parcelAllocation(mod.mix, data = simParcel, parcel.names, item.syntax,
##'                  nAlloc = 20, std.lv = TRUE)
##'
##' @export
parcelAllocation <- function(model, data, parcel.names, item.syntax,
                             nAlloc = 100, fun = "sem", alpha = .05,
                             fit.measures = c("chisq","df","cfi",
                                              "tli","rmsea","srmr"), ...,
                             show.progress = FALSE, iseed = 12345,
                             do.fit = TRUE, return.fit = FALSE, warn = FALSE) {
  if (nAlloc < 2) stop("Minimum of two allocations required.")
  if (!fun %in% c("sem","cfa","growth","lavaan"))
    stop("'fun' argument must be either 'lavaan', 'cfa', 'sem', or 'growth'")

  lavArgs <- list(...)
  lavArgs$model <- item.syntax
  lavArgs$data <- data
  lavArgs$do.fit <- FALSE

  ## fit item-level model to data
  item.fit <- do.call(fun, lavArgs)
  item.PT <- parTable(item.fit)

  ## construct parameter table for parcel-level model
  if (is.character(model)) {
    ## default lavaanify arguments
    ptArgs <- formals(lavaanify)
    ## arguments passed to lavaan by user
    fitArgs <- lavInspect(item.fit, "call")[-1]
    ## overwrite defaults with user's values
    sameArgs <- intersect(names(ptArgs), names(fitArgs))
    ptArgs[sameArgs] <- fitArgs[sameArgs]
    ptArgs$model <- model
    if (is.null(ptArgs$model.type)) ptArgs$model.type <- "sem"
    if (ptArgs$model.type != "growth") ptArgs$model.type <- "sem"
    ptArgs$ngroups <- lavInspect(item.fit, "ngroups")
    PT <- do.call("lavaanify", ptArgs)
  } else if (is.data.frame(model)) {
    PT <- model
  } else stop("'model' argument must be a character string of lavaan model",
              " syntax or a lavaan parameter table.  See ?lavaanify help page.")

  ## check that both models specify the same factors
  factorNames <- lavNames(PT, type = "lv")
  if (!all(sort(lavNames(item.PT, type = "lv")) == sort(factorNames))) {
    stop("'model' and 'item.syntax' arguments specify different factors.\n",
         "'model' specifies: ", paste(sort(factorNames), collapse = ", "), "\n",
         "'item.syntax' specifies: ", paste(sort(lavNames(item.PT,
                                                                  type = "lv")),
                                            collapse = ", "))
  }

  ## for each factor, assign item sets to parcel sets
  assignments <- list()
  for (i in factorNames) {
    ## all indicators from parcel-level model
    parcels <- PT$rhs[PT$lhs == i & PT$op == "=~"]
    ## all indicators from item-level model
    items <- item.PT$rhs[item.PT$lhs == i & item.PT$op == "=~"]
    ## exclude observed indicators from parceling scheme if specified
    ## in parcel-level model
    assignments[[i]]$parcels <- setdiff(parcels, names(data))
    assignments[[i]]$items <- setdiff(items, parcels)

    ## Does this factor have parcels?  If not, omit this factor from next loop
    if (length(assignments[[i]]$parcels) == 0L) {
      factorNames <- factorNames[-which(factorNames == i)]
      next
    }

    ## how many items per parcel?
    nItems <- length(assignments[[i]]$items)
    nParcels <- length(assignments[[i]]$parcels)
    assignments[[i]]$nPerParcel <- rep(nItems %/% nParcels, nParcels)
    if (nItems %% nParcels > 0) for (j in 1:(nItems %% nParcels)) {
      assignments[[i]]$nPerParcel[j] <- assignments[[i]]$nPerParcel[j] + 1
    }
    names(assignments[[i]]$nPerParcel) <- assignments[[i]]$parcels
  }

  ## for each allocation, create parcels from items
  dataList <- list()
  for (i in 1:nAlloc) {
    dataList[[i]] <- data
    for (j in factorNames) {
      ## create a random assignment pattern
      ranAss <- sample(rep(names(assignments[[j]]$nPerParcel),
                           times = assignments[[j]]$nPerParcel))
      ## add each parcel to a copy of the original data set
      for (k in assignments[[j]]$parcels) {
        ## which items were selected for this parcel?
        ranVars <- assignments[[j]]$items[ranAss == k]
        ## calculate row means of those items, save as parcel
        dataList[[i]][ , k] <- rowMeans(data[ , ranVars])
      }
    }
  }
  if (!do.fit) return(dataList)

  ## fit parcel-level model to list of data sets
  set.seed(iseed) # in case not using parallel
  fitList <- lavaanList(model, dataList, cmd = fun, ..., warn = warn, iseed = iseed,
                        FUN = lavaan::fitMeasures, show.progress = show.progress)
  ## for which data sets did the model converge?
  conv <- fitList@meta$ok
  if (!any(conv)) stop("The model did not converge for any allocations.")
  if (!all(conv)) message("The model did not converge for the following ",
                          "allocations: ", paste(which(!conv), collapse = ", "))

  ## tools to extract output
  getOutput <- function(x, sig = FALSE) {
    c(Avg = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE),
      Min = min(x, na.rm = TRUE), Max = max(x, na.rm = TRUE),
      Range = max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  out <- list()
  myCols <- c("lhs","op","rhs","group", "block","label")
  template <- data.frame(fitList@ParTableList[[which(conv)[1]]][myCols])

  ## parameter estimates
  Est <- sapply(fitList@ParTableList[conv], function(x) x$est)
  out$Estimates <- cbind(template, t(apply(Est, 1, getOutput)))

  ## standard errors
  SE <- sapply(fitList@ParTableList[conv], function(x) x$se)
  ## Any for which SE could not be calculated?
  missingSE <- apply(SE, 2, function(x) any(is.na(x)))
  if (!all(missingSE)) {
    if (any(missingSE)) message("Standard errors could not be computed for ",
                                "the following allocations: ",
                                paste(which(missingSE), collapse = ", "))
    out$SE <- cbind(template, t(apply(SE[ , !missingSE], 1, getOutput)))

    ## add significance test results to $Estimates
    Sig <- abs(Est[, !missingSE] / SE[, !missingSE]) > qnorm(alpha / 2,
                                                             lower.tail = FALSE)
    out$Estimates$Percent_Sig <- rowMeans(Sig)
    out$Estimates$Percent_Sig[fitList@ParTableList[[which(conv)[1]]]$free == 0L] <- NA
  } else {
    message("Standard errors could not be calculated for any converged",
            " data sets, so no significance tests could be conducted.")
    out$SE <- NULL
  }

  ## fit measures
  Fit <- do.call(cbind, fitList@funList[conv])[fit.measures, ]
  out$Fit <- data.frame(t(apply(Fit, 1, getOutput)))
  if (any(grepl(pattern = "chisq", fit.measures))) {
    out$Fit$Percent_Sig <- NA
    if ("chisq" %in% fit.measures) {
      pvalues <- sapply(fitList@funList[conv], "[", i = "pvalue")
      out$Fit["chisq", "Percent_Sig"] <- mean(pvalues < alpha, na.rm = TRUE)
    }
    if ("chisq.scaled" %in% fit.measures) {
      pvalues <- sapply(fitList@funList[conv], "[", i = "pvalue.scaled")
      out$Fit["chisq.scaled", "Percent_Sig"] <- mean(pvalues < alpha, na.rm = TRUE)
    }
  }
  if (any(grepl(pattern = "rmsea", fit.measures))) {
    if (is.null(out$Fit$Percent_Sig)) out$Fit$Percent_Sig <- NA
    if ("rmsea" %in% fit.measures) {
      pvalues <- sapply(fitList@funList[conv], "[", i = "rmsea.pvalue")
      out$Fit["rmsea", "Percent_Sig"] <- mean(pvalues < alpha, na.rm = TRUE)
    }
    if ("rmsea.scaled" %in% fit.measures) {
      pvalues <- sapply(fitList@funList[conv], "[", i = "rmsea.pvalue.scaled")
      out$Fit["rmsea.scaled", "Percent_Sig"] <- mean(pvalues < alpha, na.rm = TRUE)
    }
  }
  ## check for robust test
  if (any(grepl(pattern = "scaled", names(fitList@funList[conv][[1]]))) &
      !any(grepl(pattern = "scaled", fit.measures))) {
    warning('Robust test requested, but "fit.measures" argument does not',
            ' include any scaled measures (e.g., "chisq.scaled", ',
            '"rmsea.scaled", or "rmsea.robust").')
  }

  ## remove rows that do not correspond to estimates
  out$Estimates <- out$Estimates[fitList@ParTableList[[which(conv)[1]]]$group > 0L, ]
  if (!is.null(out$SE)) out$SE <- out$SE[fitList@ParTableList[[which(conv)[1]]]$group > 0L, ]

  ## assign class for lavaan's print method
  class(out$Estimates) <- c("lavaan.data.frame","data.frame")
  if (!is.null(out$SE)) class(out$SE) <- c("lavaan.data.frame","data.frame")
  class(out$Fit) <- c("lavaan.data.frame","data.frame")

  ## return output
  if (return.fit) {
    out$Model <- fitList
    out$Model@external$dataList <- dataList
  }
  out
}



