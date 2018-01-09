
parcelAllocation <- function(model, data, parcel.names, item.syntax,
                             nAlloc = 100, fun = "sem", alpha = .05,
                             fit.measures = c("chisq","df","cfi",
                                              "tli","rmsea","srmr"), ...,
                             show.progress = FALSE) {
  if (nAlloc < 2) stop("Minimum of two allocations required.")
  if (!fun %in% c("sem","cfa","growth","lavaan"))
    stop("'fun' argument must be either 'lavaan', 'cfa', 'sem', or 'growth'")

  lavArgs <- list(...)
  lavArgs$model <- item.syntax
  lavArgs$data <- data
  lavArgs$do.fit <- FALSE

  ## fit item-level model to data
  item.fit <- do.call(fun, lavArgs)
  item.PT <- lavaan::parTable(item.fit)

  ## construct parameter table for parcel-level model
  if (is.character(model)) {
    ## default lavaanify arguments
    ptArgs <- formals(lavaanify)
    ## arguments passed to lavaan by user
    fitArgs <- lavaan::lavInspect(item.fit, "call")[-1]
    ## overwrite defaults with user's values
    sameArgs <- intersect(names(ptArgs), names(fitArgs))
    ptArgs[sameArgs] <- fitArgs[sameArgs]
    ptArgs$model <- model
    if (is.null(ptArgs$model.type)) ptArgs$model.type <- "sem"
    if (ptArgs$model.type != "growth") ptArgs$model.type <- "sem"
    ptArgs$ngroups <- lavaan::lavInspect(item.fit, "ngroups")
    PT <- do.call("lavaanify", ptArgs)
  } else if (is.data.frame(model)) {
    PT <- model
  } else stop("'model' argument must be a character string of lavaan model",
              " syntax or a lavaan parameter table.  See ?lavaanify help page.")

  ## check that both models specify the same factors
  factorNames <- lavaan::lavNames(PT, type = "lv")
  if (!all(sort(lavaan::lavNames(item.PT, type = "lv")) == sort(factorNames))) {
    stop("'model' and 'item.syntax' arguments specify different factors.\n",
         "'model' specifies: ", paste(sort(factorNames), collapse = ", "), "\n",
         "'item.syntax' specifies: ", paste(sort(lavaan::lavNames(item.PT,
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

  ## fit parcel-level model to list of data sets
  fitList <- lavaanList(model, dataList, cmd = fun, ...,
                        FUN = lavaan::fitMeasures, show.progress = show.progress)
  ## for which data sets did the model converge?
  conv <- fitList@meta$ok
  if (!any(conv)) stop("The model did not converge for any allocations.")
  if (!all(conv)) message("The model did not converge for the following ",
                          "allocations: ", paste(which(!conv), collapse = ", "))

  ## tools to extract output
  getOutput <- function(x) {
    c(Avg = mean(x), SD = sd(x), Min = min(x), Max = max(x))
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
    out$Estimates[ , "Percent_Sig"] <- rowMeans(Sig)
    out$Estimates[fitList@ParTableList[[which(conv)[1]]]$free == 0L, "Percent_Sig"] <- NA
  } else {
    message("Standard errors could not be calculated for any converged",
            " data sets, so no significance tests could be conducted.")
    out$SE <- NULL
  }

  ## fit measures
  Fit <- do.call(cbind, fitList@funList[conv])[fit.measures, ]
  out$Fit <- data.frame(t(apply(Fit, 1, getOutput)))

  ## remove rows that do not correspond to estimates
  out$Estimates <- out$Estimates[fitList@ParTableList[[which(conv)[1]]]$group > 0L, ]
  if (!is.null(out$SE)) out$SE <- out$SE[fitList@ParTableList[[which(conv)[1]]]$group > 0L, ]

  ## assign class for lavaan's print method
  class(out$Estimates) <- c("lavaan.data.frame","data.frame")
  if (!is.null(out$SE)) class(out$SE) <- c("lavaan.data.frame","data.frame")
  class(out$Fit) <- c("lavaan.data.frame","data.frame")

  ## return output
  out
}



## New version requires lavaanList, so depends on lavaan version 0.6 or higher:
# install.packages("lavaan", repos="http://www.da.ugent.be", type="source")
library(lavaan)
simParcel <- semTools::simParcel

## item-level model (if NO parcels were created)
item.syntax <- c(paste0("f1 =~ f1item", 1:9),
                 paste0("f2 =~ f2item", 1:9))
cat(item.syntax, sep = "\n")



## 3-indicator parcels
mod.parcels <- '
f1 =~ par1 + par2 + par3
f2 =~ par4 + par5 + par6
'
## names of parcels
(parcel.names <- paste0("par", 1:6))

parcelAllocation(mod.parcels, data = simParcel, parcel.names, item.syntax,
                 nAlloc = 20, std.lv = TRUE, parallel = "snow", iseed = 12345)



## multigroup example
simParcel$group <- 0:1 # arbitrary groups for example
mod.mg <- '
f1 =~ par1 + c(L2, L2)*par2 + par3
f2 =~ par4 + par5 + par6
'
## names of parcels
(parcel.names <- paste0("par", 1:6))

set.seed(12345)
parcelAllocation(mod.mg, data = simParcel, parcel.names, item.syntax,
                 std.lv = TRUE, group = "group", group.equal = "loadings",
                 nAlloc = 20, show.progress = TRUE)



## parcels for first factor, items for second factor
mod.items <- '
f1 =~ par1 + par2 + par3
f2 =~ f2item2 + f2item7 + f2item8
'
## names of parcels
(parcel.names <- paste0("par", 1:3))

parcelAllocation(mod.items, data = simParcel, parcel.names, item.syntax,
                 nAlloc = 20, std.lv = TRUE)



## mixture of 1- and 3-indicator parcels for second factor
mod.mix <- '
f1 =~ par1 + par2 + par3
f2 =~ f2item2 + f2item7 + f2item8 + par4 + par5 + par6
'
## names of parcels
(parcel.names <- paste0("par", 1:6))

parcelAllocation(mod.mix, data = simParcel, parcel.names, item.syntax,
                 nAlloc = 20, std.lv = TRUE)



