### Terrence D. Jorgensen
### Last updated: 11 June 2019
### Savalei & Yuan's (2009) model-based bootstrap for missing data


## ----------------------------
## "BootMiss" Class and Methods
## ----------------------------

##' Class For the Results of Bollen-Stine Bootstrap with Incomplete Data
##'
##' This class contains the results of Bollen-Stine bootstrap with missing data.
##'
##'
##' @name BootMiss-class
##' @aliases BootMiss-class show,BootMiss-method summary,BootMiss-method
##' hist,BootMiss-method
##' @docType class
##' @section Objects from the Class: Objects can be created via the
##' \code{\link{bsBootMiss}} function.
##' @slot time A list containing 2 \code{difftime} objects (\code{transform}
##'  and \code{fit}), indicating the time elapsed for data transformation and
##'  for fitting the model to bootstrap data sets, respectively.
##' @slot transData Transformed data
##' @slot bootDist The vector of \eqn{chi^2} values from bootstrap data sets
##'  fitted by the target model
##' @slot origChi The \eqn{chi^2} value from the original data set
##' @slot df The degree of freedom of the model
##' @slot bootP The \emph{p} value comparing the original \eqn{chi^2} with the
##'  bootstrap distribution
##' @author Terrence D. Jorgensen (University of Amsterdam;
##' \email{TJorgensen314@@gmail.com})
##' @seealso \code{\link{bsBootMiss}}
##' @examples
##'
##' # See the example from the bsBootMiss function
##'
setClass("BootMiss", representation(time = "list",
                                    transData = "data.frame",
                                    bootDist = "vector",
                                    origChi = "numeric",
                                    df = "numeric",
                                    bootP = "numeric"))

##' @rdname BootMiss-class
##' @aliases show,BootMiss-method
##' @importFrom stats pchisq
##' @export
setMethod("show", "BootMiss",
function(object) {
  cat("Chi-Squared = ", object@origChi, "\nDegrees of Freedom = ",
      object@df, "\nTheoretical p value = ",
      pchisq(object@origChi, object@df, lower.tail = FALSE),
      "\n    i.e., pchisq(", object@origChi, ", df = ",
                    object@df, ", lower.tail = FALSE)\n",
      "\nBootstrapped p value = ", object@bootP, "\n\n", sep = "")
  invisible(object)
})

##' @rdname BootMiss-class
##' @aliases summary,BootMiss-method
##' @importFrom stats var
##' @export
setMethod("summary", "BootMiss",
function(object) {
  cat("Time elapsed to transform the data:\n")
  print(object@time$transform)
  cat("\nTime elapsed to fit the model to", length(object@bootDist),
      "bootstrapped samples:\n")
  print(object@time$fit)
  cat("\nMean of Theoretical Distribution = DF =", object@df,
      "\nVariance of Theoretical Distribution = 2*DF =", 2*object@df,
      "\n\nMean of Bootstrap Distribution =", mean(object@bootDist),
      "\nVariance of Bootstrap Distribution =",
      var(object@bootDist), "\n\n")
  show(object)
  invisible(object)
})

##' @rdname BootMiss-class
##' @aliases hist,BootMiss-method
##' @importFrom stats qchisq dchisq quantile
##' @param object,x object of class \code{BootMiss}
##' @param ... Additional arguments to pass to \code{\link[graphics]{hist}}
##' @param alpha alpha level used to draw confidence limits
##' @param nd number of digits to display
##' @param printLegend \code{logical}. If \code{TRUE} (default), a legend will
##'  be printed with the histogram
##' @param legendArgs \code{list} of arguments passed to the
##'  \code{\link[graphics]{legend}} function.  The default argument is a list
##'  placing the legend at the top-left of the figure.
##' @return The \code{hist} method returns a list of \code{length == 2},
##'  containing the arguments for the call to \code{hist} and the arguments
##'  to the call for \code{legend}, respectively.
##' @export
setMethod("hist", "BootMiss",
function(x, ..., alpha = .05, nd = 2, printLegend = TRUE,
         legendArgs = list(x = "topleft")) {
  ChiSq <- x@origChi
  DF <- x@df
  bootDist <- x@bootDist
  bCrit <- round(quantile(bootDist, probs = 1 - alpha), nd)
  theoDist <- dchisq(seq(.1, max(c(ChiSq, bootDist)), by = .1), df = DF)
  Crit <- round(qchisq(p = alpha, df = DF, lower.tail = FALSE), nd)
  Lim <- c(0, max(c(ChiSq, bootDist, Crit)))
  if (ChiSq > Lim[2]) Lim[2] <- ChiSq

  histArgs <- list(...)
  histArgs$x <- bootDist
  histArgs$freq <- FALSE
  if (is.null(histArgs$col)) histArgs$col <- "grey69"
  if (is.null(histArgs$xlim)) histArgs$xlim <- Lim
  if (is.null(histArgs$main)) histArgs$main <- expression("Model-Based Bootstrap Distribution of" ~ chi^2)
  if (is.null(histArgs$ylab)) histArgs$ylab <- "Probability Density"
  if (is.null(histArgs$xlab)) histArgs$xlab <- expression(chi^2)

  if (printLegend) {
    if (nd < length(strsplit(as.character(1 / alpha), "")[[1]]) - 1) {
      warning(paste0("The number of digits argument (nd = ", nd ,
                     ") is too low to display your p value at the",
                     " same precision as your requested alpha level (alpha = ",
                     alpha, ")"))
    }

    if (x@bootP < (1 / 10^nd)) {
      pVal <- paste(c("< .", rep(0, nd - 1),"1"), collapse = "")
    } else {
      pVal <- paste("=", round(x@bootP, nd))
    }

    if (is.null(legendArgs$box.lty)) legendArgs$box.lty <- 0
    if (is.null(legendArgs$lty)) legendArgs$lty <- c(1, 2, 2, 1, 0, 0)
    if (is.null(legendArgs$lwd)) legendArgs$lwd <- c(2, 2, 2, 3, 0, 0)
    #if (is.null(legendArgs$cex)) legendArgs$cex <- c(1.1, 1, 1, 1, 1, 1)
    if (is.null(legendArgs$col)) legendArgs$col <- c("black","black","grey69","red","", "")
    legendArgs$legend <- c(bquote(chi[.(paste("df =", DF))]^2),
                           bquote(Critical ~ chi[alpha ~ .(paste(" =", alpha))]^2 == .(Crit)),
                           bquote(Bootstrap~Critical~chi[alpha ~ .(paste(" =", alpha))]^2 == .(bCrit)),
                           expression(Observed ~ chi^2),
                           bquote(.("")),
                           bquote(Bootstrap ~ italic(p) ~~ .(pVal)))
  }

  H <- do.call(hist, c(histArgs["x"], plot = FALSE))
  histArgs$ylim <- c(0, max(H$density, theoDist))
  suppressWarnings({
    do.call(hist, histArgs)
    lines(x = seq(.1, max(c(ChiSq, bootDist)), by = .1), y = theoDist, lwd = 2)
    abline(v = Crit, col = "black", lwd = 2, lty = 2)
    abline(v = bCrit, col = "grey69", lwd = 2, lty = 2)
    abline(v = ChiSq, col = "red", lwd = 3)
    if (printLegend) do.call(legend, legendArgs)
  })
  ## return arguments to create histogram (and optionally, legend)
  invisible(list(hist = histArgs, legend = legendArgs))
})



## --------------------
## Constructor Function
## --------------------

##' Bollen-Stine Bootstrap with the Existence of Missing Data
##'
##' Implement the Bollen and Stine's (1992) Bootstrap when missing observations
##' exist. The implemented method is proposed by Savalei and Yuan (2009). This
##' can be used in two ways. The first and easiest option is to fit the model to
##' incomplete data in \code{lavaan} using the FIML estimator, then pass that
##' \code{lavaan} object to \code{bsBootMiss}.
##'
##' The second is designed for users of other software packages (e.g., LISREL,
##' EQS, Amos, or Mplus). Users can import their data, \eqn{\chi^2} value, and
##' model-implied moments from another package, and they have the option of
##' saving (or writing to a file) either the transformed data or bootstrapped
##' samples of that data, which can be analyzed in other programs. In order to
##' analyze the bootstrapped samples and return a \emph{p} value, users of other
##' programs must still specify their model using lavaan syntax.
##'
##'
##' @importFrom lavaan lavInspect parTable
##' @param x A target \code{lavaan} object used in the Bollen-Stine bootstrap
##' @param transformation The transformation methods in Savalei and Yuan (2009).
##' There are three methods in the article, but only the first two are currently
##' implemented here.  Use \code{transformation = 1} when there are few missing
##' data patterns, each of which has a large size, such as in a
##' planned-missing-data design.  Use \code{transformation = 2} when there are
##' more missing data patterns. The currently unavailable
##' \code{transformation = 3} would be used when several missing data patterns
##' have n = 1.
##' @param nBoot The number of bootstrap samples.
##' @param model Optional. The target model if \code{x} is not provided.
##' @param rawData Optional. The target raw data set if \code{x} is not
##'  provided.
##' @param Sigma Optional. The model-implied covariance matrix if \code{x} is
##'  not provided.
##' @param Mu Optional. The model-implied mean vector if \code{x} is not
##'  provided.
##' @param group Optional character string specifying the name of the grouping
##'  variable in \code{rawData} if \code{x} is not provided.
##' @param ChiSquared Optional. The model's \eqn{\chi^2} test statistic if
##'  \code{x} is not provided.
##' @param EMcov Optional, if \code{x} is not provided. The EM (or Two-Stage ML)
##' estimated covariance matrix used to speed up Transformation 2 algorithm.
##' @param transDataOnly Logical. If \code{TRUE}, the result will provide the
##' transformed data only.
##' @param writeTransData Logical. If \code{TRUE}, the transformed data set is
##' written to a text file, \code{transDataOnly} is set to \code{TRUE}, and the
##' transformed data is returned invisibly.
##' @param bootSamplesOnly Logical. If \code{TRUE}, the result will provide
##' bootstrap data sets only.
##' @param writeBootData Logical. If \code{TRUE}, the stacked bootstrap data
##' sets are written to a text file, \code{bootSamplesOnly} is set to
##' \code{TRUE}, and the list of bootstrap data sets are returned invisibly.
##' @param writeArgs Optional \code{list}. If \code{writeBootData = TRUE} or
##' \code{writeBootData = TRUE}, user can pass arguments to the
##' \code{\link[utils]{write.table}} function as a list.  Some default values
##' are provided: \code{file} = "bootstrappedSamples.dat", \code{row.names} =
##' \code{FALSE}, and \code{na} = "-999", but the user can override all of these
##' by providing other values for those arguments in the \code{writeArgs} list.
##' @param seed The seed number used in randomly drawing bootstrap samples.
##' @param suppressWarn Logical. If \code{TRUE}, warnings from \code{lavaan}
##' function will be suppressed when fitting the model to each bootstrap sample.
##' @param showProgress Logical. Indicating whether to display a progress bar
##' while fitting models to bootstrap samples.
##' @param \dots The additional arguments in the \code{\link[lavaan]{lavaan}}
##' function. See also \code{\link[lavaan]{lavOptions}}
##' @return As a default, this function returns a \code{\linkS4class{BootMiss}}
##' object containing the results of the bootstrap samples. Use \code{show},
##' \code{summary}, or \code{hist} to examine the results. Optionally, the
##' transformed data set is returned if \code{transDataOnly = TRUE}. Optionally,
##' the bootstrap data sets are returned if \code{bootSamplesOnly = TRUE}.
##' @author Terrence D. Jorgensen (University of Amsterdam;
##' \email{TJorgensen314@@gmail.com})
##'
##' Syntax for transformations borrowed from http://www2.psych.ubc.ca/~vsavalei/
##' @seealso \code{\linkS4class{BootMiss}}
##' @references
##'
##' Bollen, K. A., & Stine, R. A. (1992). Bootstrapping goodness-of-fit measures
##' in structural equation models. \emph{Sociological Methods &
##' Research, 21}(2), 205--229. doi:10.1177/0049124192021002004
##'
##' Savalei, V., & Yuan, K.-H. (2009). On the model-based bootstrap with missing
##' data: Obtaining a p-value for a test of exact fit. \emph{Multivariate
##' Behavioral Research, 44}(6), 741--763. doi:10.1080/00273170903333590
##' @examples
##'
##' \dontrun{
##' dat1 <- HolzingerSwineford1939
##' dat1$x5 <- ifelse(dat1$x1 <= quantile(dat1$x1, .3), NA, dat1$x5)
##' dat1$x9 <- ifelse(is.na(dat1$x5), NA, dat1$x9)
##'
##' targetModel <- "
##' visual  =~ x1 + x2 + x3
##' textual =~ x4 + x5 + x6
##' speed   =~ x7 + x8 + x9
##' "
##' targetFit <- sem(targetModel, dat1, meanstructure = TRUE, std.lv = TRUE,
##'                  missing = "fiml", group = "school")
##' summary(targetFit, fit = TRUE, standardized = TRUE)
##'
##' # The number of bootstrap samples should be much higher.
##' temp <- bsBootMiss(targetFit, transformation = 1, nBoot = 10, seed = 31415)
##'
##' temp
##' summary(temp)
##' hist(temp)
##' hist(temp, printLegend = FALSE) # suppress the legend
##' ## user can specify alpha level (default: alpha = 0.05), and the number of
##' ## digits to display (default: nd = 2).  Pass other arguments to hist(...),
##' ## or a list of arguments to legend() via "legendArgs"
##' hist(temp, alpha = .01, nd = 3, xlab = "something else", breaks = 25,
##'      legendArgs = list("bottomleft", box.lty = 2))
##' }
##'
##' @export
bsBootMiss <- function(x, transformation = 2, nBoot = 500, model, rawData,
                       Sigma, Mu, group, ChiSquared, EMcov,
                       writeTransData = FALSE, transDataOnly = FALSE,
                       writeBootData = FALSE, bootSamplesOnly = FALSE,
                       writeArgs, seed = NULL, suppressWarn = TRUE,
                       showProgress = TRUE, ...) {
  if(writeTransData) transDataOnly <- TRUE
  if(writeBootData) bootSamplesOnly <- TRUE

  check.nBoot <- (!is.numeric(nBoot) | nBoot < 1L) & !transDataOnly
  if (check.nBoot) stop("The \"nBoot\" argument must be a positive integer.")

  ## Which transformation?
  if (!(transformation %in% 1:2)) stop("User must specify transformation 1 or 2.
       Consult Savalei & Yuan (2009) for advice.
       Transformation 3 is not currently available.")
  if (transformation == 2) SavaleiYuan <- trans2
  #if (transformation == 3) SavaleiYuan <- trans3

  ######################
  ## Data Preparation ##
  ######################

  ## If a lavaan object is supplied, the extracted values for rawData, Sigma, Mu,
  ## EMcov, and EMmeans will override any user-supplied arguments.
  if (hasArg(x)) {
    if (!lavInspect(x, "options")$meanstructure)
      stop('You must fit the lavaan model with the argument "meanstructure=TRUE".')
    nG <- lavInspect(x, "ngroups")
    if (nG == 1L) {
      rawData <- list(as.data.frame(lavInspect(x, "data")))
    } else rawData <- lapply(lavInspect(x, "data"), as.data.frame)
    for (g in seq_along(rawData)) {
      colnames(rawData[[g]]) <- lavaan::lavNames(x)
      checkAllMissing <- apply(rawData[[g]], 1, function(x) all(is.na(x)))
      if (any(checkAllMissing)) rawData[[g]] <- rawData[[g]][!checkAllMissing, ]
    }
    ChiSquared <- lavInspect(x, "fit")[c("chisq", "chisq.scaled")]
    ChiSquared <- ifelse(is.na(ChiSquared[2]), ChiSquared[1], ChiSquared[2])
    group <- lavInspect(x, "group")
    if (length(group) == 0) group <- "group"
    group.label <- lavInspect(x, "group.label")
    if (length(group.label) == 0) group.label <- 1
    Sigma <- lavInspect(x, "cov.ov")
    Mu <- lavInspect(x, "mean.ov")
    EMcov <- lavInspect(x, "sampstat")$cov
    if (nG == 1L) {
      Sigma <- list(Sigma)
      Mu <- list(Mu)
      EMcov <- list(EMcov)
    }
  } else {
  ## If no lavaan object is supplied, check that required arguments are.
    suppliedData <- c(hasArg(rawData), hasArg(Sigma), hasArg(Mu))
    if (!all(suppliedData)) {
      stop("Without a lavaan fitted object, user must supply raw data and",
           " model-implied moments.")
    }
    if (!hasArg(model) & !(transDataOnly | bootSamplesOnly)) {
      stop("Without model syntax or fitted lavaan object, user can only call",
           " this function to save transformed data or bootstrapped samples.")
    }
    if (!hasArg(ChiSquared) & !(transDataOnly | bootSamplesOnly)) {
      stop("Without a fitted lavaan object or ChiSquared argument, user can",
           " only call this function to save transformed data, bootstrapped",
           " samples, or bootstrapped chi-squared values.")
    }
    if (!any(c(transDataOnly, bootSamplesOnly))) {
      if (!is.numeric(ChiSquared)) stop("The \"ChiSquared\" argument must be numeric.")
    }

    ## If user supplies one-group data & moments, convert to lists.
    if (class(rawData) == "data.frame") rawData <- list(rawData)
    if (class(rawData) != "list") {
      stop("The \"rawData\" argument must be a data.frame or list of data frames.")
    } else {
      if (!all(sapply(rawData, is.data.frame))) stop("Every element of \"rawData\" must be a data.frame")
    }
    if (class(Sigma) == "matrix") Sigma <- list(Sigma)
    if (is.numeric(Mu)) Mu <- list(Mu)

    ## check whether EMcov was supplied for starting values in Trans2/Trans3
    if (!hasArg(EMcov)) {
      EMcov <- vector("list", length(Sigma))
    } else {
      if (class(EMcov) == "matrix") EMcov <- list(EMcov)
      ## check EMcov is symmetric and dimensions match Sigma
      for (g in seq_along(EMcov)) {
        if (!isSymmetric(EMcov[[g]])) stop("EMcov in group ", g, " not symmetric.")
        unequalDim <- !all(dim(EMcov[[g]]) == dim(Sigma[[g]]))
        if (unequalDim) stop("Unequal dimensions in Sigma and EMcov.")
      }
    }

    ## Check the number of groups by the size of the lists.
    unequalGroups <- !all(length(rawData) == c(length(Sigma), length(Mu)))
    if (unequalGroups) stop("Unequal number of groups in rawData, Sigma, Mu.
       For multiple-group models, rawData must be a list of data frames,
       NOT a single data frame with a \"group\" column.")
    nG <- length(Sigma)

    ## In each group, check Sigma is symmetric and dimensions match rawData and Mu.
    for (g in seq_along(rawData)) {
      if (!isSymmetric(Sigma[[g]])) stop("Sigma in group ", g, " not symmetric.")
      unequalDim <- !all(ncol(rawData[[g]]) == c(nrow(Sigma[[g]]), length(Mu[[g]])))
      if (unequalDim) stop("Unequal dimensions in rawData, Sigma, Mu.")
    }

    ## Check for names of group levels. If NULL, assign arbitrary ones.
    if (!hasArg(group)) group <- "group"
    if (!is.character(group)) stop("The \"group\" argument must be a character string.")
    if (is.null(names(rawData))) {
      group.label <- paste0("g", seq_along(rawData))
    } else {
      group.label <- names(rawData)
    }
  }

  ## save a copy as myTransDat, whose elements will be replaced iteratively by
  ## group and by missing data pattern within group.
  myTransDat <- rawData
  names(myTransDat) <- group.label
  output <- list()

  #########################
  ## Data Transformation ##
  #########################

  for (g in seq_along(group.label)) {
    if (transformation == 1) {
      ## get missing data patterns
      R <- ifelse(is.na(rawData[[g]]), 1, 0)
      rowMissPatt <- apply(R, 1, function(x) paste(x, collapse = ""))
      patt <- unique(rowMissPatt)
      myRows <- lapply(patt, function(x) which(rowMissPatt == x))

      ## for each pattern, apply transformation
      tStart <- Sys.time()
      transDatList <- lapply(patt, trans1, rowMissPatt = rowMissPatt,
                             dat = rawData[[g]], Sigma = Sigma[[g]], Mu = Mu[[g]])
      output$timeTrans <- Sys.time() - tStart
      for (i in seq_along(patt)) myTransDat[[g]][myRows[[i]], ] <- transDatList[[i]]
    } else {
      tStart <- Sys.time()
      myTransDat[[g]] <- SavaleiYuan(dat = rawData[[g]], Sigma = Sigma[[g]],
                                     Mu = Mu[[g]], EMcov = EMcov[[g]])
      output$timeTrans <- Sys.time() - tStart
    }
  }

  ## option to end function here
  if (transDataOnly) {
    for (g in seq_along(myTransDat)) myTransDat[[g]][ , group] <- group.label[g]
    ## option to write transformed data to a file
    if (writeTransData) {
      ## Set a few options, if the user didn't.
      if (!hasArg(writeArgs)) writeArgs <- list(file = "transformedData.dat",
                                                row.names = FALSE, na = "-999")
      if (!exists("file", where = writeArgs)) writeTransArgs$file <- "transformedData.dat"
      if (!exists("row.names", where = writeArgs)) writeArgs$row.names <- FALSE
      if (!exists("na", where = writeArgs)) writeArgs$na <- "-999"

      ## add grouping variable and bind together into one data frame
      for (g in seq_along(myTransDat)) myTransDat[[g]][ , group] <- group.label[g]
      writeArgs$x <- do.call("rbind", myTransDat)

      ## write to file, print details to screen
      do.call("write.table", writeArgs)
      cat("Transformed data was written to file \"", writeArgs$file, "\" in:\n\n",
          getwd(), "\n\nunless path specified by user in 'file' argument.\n", sep = "")
      return(invisible(writeArgs$x))
    }
    return(do.call("rbind", myTransDat))
  }

  #############################################
  ## Bootstrap distribution of fit statistic ##
  #############################################

  ## draw bootstrap samples
  if (!is.null(seed)) set.seed(seed)
  bootSamples <- lapply(1:nBoot, function(x) getBootSample(myTransDat, group, group.label))

  ## option to write bootstrapped samples to file(s)
  if (writeBootData) {
    ## Set a few options, if the user didn't.
    if (!hasArg(writeArgs)) writeArgs <- list(file = "bootstrappedSamples.dat",
                                              row.names = FALSE, na = "-999")
    if (!exists("file", where = writeArgs)) writeTransArgs$file <- "bootstrappedSamples.dat"
    if (!exists("row.names", where = writeArgs)) writeArgs$row.names <- FALSE
    if (!exists("na", where = writeArgs)) writeArgs$na <- "-999"

    ## add indicator for bootstrapped sample, bind together into one data frame
    for (b in seq_along(bootSamples)) bootSamples[[b]]$bootSample <- b
    writeArgs$x <- do.call("rbind", bootSamples)

    ## write to file, print details to screen
    do.call("write.table", writeArgs)
    cat("Bootstrapped samples written to file \"", writeArgs$file, "\" in:\n\n",
        getwd(), "\n\nunless path specified by user in 'file' argument.\n", sep = "")
    return(invisible(bootSamples))
  }

  ## option to end function here
  if (bootSamplesOnly) return(bootSamples)

  ## check for lavaan arguments in (...)
  lavaanArgs <- list(...)
  lavaanArgs$group <- group

  ## fit model to bootstrap samples, save distribution of chi-squared test stat
  if (hasArg(x)) {
    ## grab defaults from lavaan object "x"
    lavaanArgs$slotParTable <- as.list(parTable(x))
    lavaanArgs$slotModel <- x@Model
    lavaanArgs$slotOptions <- lavInspect(x, "options")
  } else {
    lavaanArgs$model <- model
    lavaanArgs$missing <- "fiml"
    ## set defaults that will be necessary for many models to run, that will
    ## probably not be specified explictly or included in lavaan syntax
    lavaanArgs$meanstructure <- TRUE
    if (!exists("auto.var", where = lavaanArgs)) lavaanArgs$auto.var <- TRUE
    if (!exists("auto.cov.y", where = lavaanArgs)) lavaanArgs$auto.cov.y <- TRUE
    if (!exists("auto.cov.lv.x", where = lavaanArgs)) lavaanArgs$auto.cov.lv.x <- TRUE
  }
  ## run bootstrap fits
  if (showProgress) {
    mypb <- utils::txtProgressBar(min = 1, max = nBoot, initial = 1, char = "=",
                                  width = 50, style = 3, file = "")
    bootFits <- numeric()
    tStart <- Sys.time()
    for (j in 1:nBoot) {
      bootFits[j] <- fitBootSample(bootSamples[[j]], args = lavaanArgs,
                                   suppress = suppressWarn)
      utils::setTxtProgressBar(mypb, j)
    }
    close(mypb)
    output$timeFit <- Sys.time() - tStart
  } else {
    tStart <- Sys.time()
    bootFits <- sapply(bootSamples, fitBootSample, args = lavaanArgs,
                       suppress = suppressWarn)
    output$timeFit <- Sys.time() - tStart
  }

  ## stack groups, save transformed data and distribution in output object
  for (g in seq_along(myTransDat)) myTransDat[[g]][ , group] <- group.label[g]
  output$Transformed.Data <- do.call("rbind", myTransDat)
  output$Bootstrapped.Distribution <- bootFits
  output$Original.ChiSquared <- ChiSquared
  if (hasArg(x)) {
    output$Degrees.Freedom <- lavInspect(x, "fit")["df"]
  } else {
    convSamp <- which(!is.na(bootFits))[1]
    lavaanArgs$data <- bootSamples[[convSamp]]
	lavaanlavaan <- function(...) { lavaan::lavaan(...) }
    output$Degrees.Freedom <- lavInspect(do.call(lavaanlavaan, lavaanArgs), "fit")["df"]
  }

  ## calculate bootstrapped p-value
  output$Bootstrapped.p.Value <- mean(bootFits >= ChiSquared, na.rm = TRUE)

  ## print warning if any models didn't converge
  if (any(is.na(bootFits))) {
    nonConvMessage <- paste("Model did not converge for the following bootstrapped samples",
                            paste(which(is.na(bootFits)), collapse = "\t"), sep = ":\n")
    warning(nonConvMessage)
  }

  finalResult <- new("BootMiss", time = list(transform = output$timeTrans, fit = output$timeFit), transData = output$Transformed.Data, bootDist = output$Bootstrapped.Distribution, origChi = output$Original.ChiSquared, df = output$Degrees.Freedom, bootP = output$Bootstrapped.p.Value)

  finalResult
}


## ----------------
## Hidden Functions
## ----------------

## Function to execute Transformation 1 on a single missing-data pattern
trans1 <- function(MDpattern, rowMissPatt, dat, Sigma, Mu) {
  myRows <- which(rowMissPatt == MDpattern)
  X <- apply(dat[myRows, ], 2, scale, scale = FALSE)
  observed <- !is.na(X[1, ])
  Xreduced <- X[ , observed]
  Mreduced <- as.numeric(Mu[observed])
  SigmaChol <- chol(Sigma[observed, observed])
  S <- t(Xreduced) %*% Xreduced / nrow(X)
  Areduced <- t(SigmaChol) %*% t(solve(chol(S)))
  Yreduced <- t(Areduced %*% t(Xreduced) + Mreduced)
  Y <- replace(X, !is.na(X), Yreduced)
  Y
}

## Function to execute Transformation 2 on a single group
trans2 <- function(dat, Sigma, Mu, EMcov) {
  ## Computing Function of A (eq. 12), whose root is desired
  eq12 <- function(A) {
    ga <- rep(0, pStar)
    for (j in 1:J) {
      Tj <- Mjs[[j]] %*% A %*% Hjs[[j]] %*% A %*% Mjs[[j]] - Mjs[[j]]
      ga <- ga + Njs[j] * Dupinv %*% c(Tj) # same as vech(Tj)
    }
    ga
  }

  ## Computing Derivative of Function of A (eq. 13)
  eq13 <- function(A) {
    deriv12 <- matrix(0, nrow = pStar, ncol = pStar)
    for (j in 1:J) {
      Tj1 <- Mjs[[j]] %*% A %*% Hjs[[j]]
      deriv12 <- deriv12 + 2*Njs[j]*Dupinv %*% kronecker(Tj1, Mjs[[j]]) %*% Dup
    }
    deriv12
  }

  ## get missing data patterns
  R <- ifelse(is.na(dat), 1, 0)
  rowMissPatt <- apply(R, 1, function(x) paste(x, collapse = ""))
  MDpattern <- unique(rowMissPatt)
  ## sample size within each MD pattern
  Njs <- sapply(MDpattern, function(patt) sum(rowMissPatt == patt))
  J <- length(MDpattern) # number of MD patterns
  p <- ncol(dat) # number of variables in model
  pStar <- p*(p + 1) / 2  # number of nonredundant covariance elements

  ## create empty lists for each MD pattern
  Xjs <- vector("list", J)
  Hjs <- vector("list", J)
  Mjs <- vector("list", J)

  ## Create Duplication Matrix and its inverse (Magnus & Neudecker, 1999)
  Dup <- lavaan::lav_matrix_duplication(p)
  Dupinv <- solve(t(Dup) %*% Dup) %*% t(Dup)

  ## step through each MD pattern, populate Hjs and Mjs
  for (j in 1:J) {
    Xjs[[j]] <- apply(dat[rowMissPatt == MDpattern[j], ], 2, scale, scale = FALSE)
    if (!is.matrix(Xjs[[j]])) Xjs[[j]] <- t(Xjs[[j]])
    observed <- !is.na(Xjs[[j]][1, ])
    Sj <- t(Xjs[[j]]) %*% Xjs[[j]] / Njs[j]
    Hjs[[j]] <- replace(Sj, is.na(Sj), 0)
    Mjs[[j]] <- replace(Sj, !is.na(Sj), solve(Sigma[observed, observed]))
    Mjs[[j]] <- replace(Mjs[[j]], is.na(Mjs[[j]]), 0)
  }

  ## Compute starting Values for A
  if (is.null(EMcov)) {
    A <- diag(p)
  } else {
    EMeig <- eigen(EMcov)
    EMrti <- EMeig$vectors %*% diag(1 / sqrt(EMeig$values)) %*% t(EMeig$vectors)
    Sigeig <- eigen(Sigma)
    Sigrt <- Sigeig$vectors %*% diag(sqrt(Sigeig$values)) %*% t(Sigeig$vectors)
    B <- Sigrt %*% EMrti
    A <- .5*(B + t(B))
  }

  ## Newton Algorithm for finding root (eq. 14)
  crit <- .1
  a <- c(A)
  fA <- eq12(A)
  while (crit > 1e-11) {
    dvecF <- eq13(A)
    a <- a - Dup %*% solve(dvecF) %*% fA
    A <- matrix(a, ncol = p)
    fA <- eq12(A)
    crit <- max(abs(fA))
  }

  ## Transform dataset X to dataset Y
  Yjs <- Xjs
  for (j in 1:J) {
    observed <- !is.na(Xjs[[j]][1, ])
    XjReduced <- Xjs[[j]][ , observed, drop = FALSE]
    Aj <- A[observed, observed, drop = FALSE]
    Mj <- as.numeric(Mu[observed])
    Yj <- t(Aj %*% t(XjReduced) + Mj)
    Yjs[[j]] <- replace(Yjs[[j]], !is.na(Yjs[[j]]), Yj)
  }
  Y <- as.data.frame(do.call("rbind", Yjs))
  colnames(Y) <- colnames(dat)
  Y
}


## Function to execute Transformation 3 on a single group -- TRANSFORMATION DOES NOT RETURN CH-SQ = 0
trans3 <- function(dat, Sigma, Mu, EMcov) {
  # Computing Saturated Means as a Function of A (eq. B1 in Appendix B)
  mut <- function(A) {
    M <- matrix(0, ncol = 1, nrow = p)
    for (j in 1:J) {
      M <- M + Njs[[j]] * Mjs[[j]] %*% A %*% Ybarjs[[j]]
    }
    Mjtoti %*% M
  }

  # Computing Function of A (eq. 18) whose root is desired
  eq18 <- function(A) {
    ga <- rep(0, pStar)
    mutilda <- mut(A)
    for (j in 1:J) {
      Tj <- Mjs[[j]] %*% A %*% Hjs[[j]] %*% A %*% Mjs[[j]] - Mjs[[j]]
      dif <- A %*% Ybarjs[[j]] - mutilda
      middle <- dif %*% t(dif)
      Tjnew <- Tj + Mjs[[j]] %*% middle %*% Mjs[[j]]
      ga <- ga + Njs[j] * Dupinv %*% c(Tjnew)
    }
    ga
  }

  # Computing Derivative of Function eq. 18
  deriv18 <- function(A) {
    d18 <- matrix(0, nrow = pStar, ncol = pStar)
    for (j in 1:J) {
      Tj1 <- Mjs[[j]] %*% A %*% Hjs[[j]]
      mutilda <- mut(A)
      dif <- A %*% Ybarjs[[j]] - mutilda
      Tj2 <- Mjs[[j]] %*% dif %*% t(Ybarjs[[j]])
      Tj3 <- kronecker(Mjs[[j]] %*% dif, Mjs[[j]]) %*% Mjtoti %*% Tj3add
      d18 <- d18 + 2*Njs[j]*Dupinv %*% ((kronecker((Tj1 + Tj2), Mjs[[j]])) - Tj3) %*% Dup
    }
    d18
  }

  ## get missing data patterns
  R <- ifelse(is.na(dat), 1, 0)
  rowMissPatt <- apply(R, 1, function(x) paste(x, collapse = ""))
  MDpattern <- unique(rowMissPatt)
  ## sample size within each MD pattern
  Njs <- sapply(MDpattern, function(patt) sum(rowMissPatt == patt))
  J <- length(MDpattern) # number of MD patterns
  p <- ncol(dat) # number of variables in model
  pStar <- p*(p + 1) / 2  # number of nonredundant covariance elements

  ## create empty lists for each MD pattern
  Xjs <- vector("list", J)
  Ybarjs <- vector("list", J)
  Hjs <- vector("list", J)
  Mjs <- vector("list", J)
  Mjtot <- matrix(0, ncol = p, nrow = p)
  Tj3add <- matrix(0, nrow = p, ncol = p * p)

  ## Create Duplication Matrix and its inverse (Magnus & Neudecker, 1999)
  Dup <- lavaan::duplicationMatrix(p)
  Dupinv <- solve(t(Dup) %*% Dup) %*% t(Dup)

  ## step through each MD pattern, populate Hjs and Mjs
  for (j in 1:J) {
    Xjs[[j]] <- apply(dat[rowMissPatt == MDpattern[j], ], 2, scale, scale = FALSE)
    if (!is.matrix(Xjs[[j]])) Xjs[[j]] <- t(Xjs[[j]])
    observed <- !is.na(Xjs[[j]][1, ])
    pj <- p - sum(observed)
    means <- colMeans(dat[rowMissPatt == MDpattern[j], ])
    Ybarjs[[j]] <- replace(means, is.na(means), 0)
    Sj <- t(Xjs[[j]]) %*% Xjs[[j]] / Njs[j]
    Hjs[[j]] <- replace(Sj, is.na(Sj), 0)
    Mjs[[j]] <- replace(Sj, !is.na(Sj), solve(Sigma[observed, observed]))
    Mjs[[j]] <- replace(Mjs[[j]], is.na(Mjs[[j]]), 0)
    Mjtot <- Mjtot + Njs[[j]] * Mjs[[j]]
    Tj3add <- Tj3add + Njs[[j]] * kronecker(t(Ybarjs[[j]]), Mjs[[j]])
  }
  Mjtoti <- solve(Mjtot)

  ## Compute starting Values for A
  if (is.null(EMcov)) {
    A <- diag(p)
  } else {
    EMeig <- eigen(EMcov)
    EMrti <- EMeig$vectors %*% diag(1 / sqrt(EMeig$values)) %*% t(EMeig$vectors)
    Sigeig <- eigen(Sigma)
    Sigrt <- Sigeig$vectors %*% diag(sqrt(Sigeig$values)) %*% t(Sigeig$vectors)
    B <- Sigrt %*% EMrti
    A <- .5*(B + t(B))
  }

  ## Newton Algorithm for finding root (eq. 14)
  crit <- .1
  a <- c(A)
  fA <- eq18(A)
  while (crit > 1e-11) {
    dvecF <- deriv18(A)
    a <- a - Dup %*% solve(dvecF) %*% fA
    A <- matrix(a, ncol = p)
    fA <- eq18(A)
    crit <- max(abs(fA))
  }

  ## Transform dataset X to dataset Y (Z in the paper, eqs. 15-16)
  Yjs <- Xjs
  for (j in 1:J) {
    observed <- !is.na(Xjs[[j]][1, ])
    XjReduced <- Xjs[[j]][ , observed, drop = FALSE]
    Aj <- A[observed, observed, drop = FALSE]
    Mj <- as.numeric((Mu - mut(A))[observed])
    Yj <- t(Aj %*% t(XjReduced) + Mj)
    Yjs[[j]] <- replace(Yjs[[j]], !is.na(Yjs[[j]]), Yj)
  }
  Y <- as.data.frame(do.call("rbind", Yjs))
  colnames(Y) <- colnames(dat)
  Y
}


## Get a single bootstrapped sample from the transformed data.  If there are
## multiple groups, bootstrapping occurs independently within each group, and
## a single data frame is returned.  A new column is added to indicate group
## membership, which will be ignored in a single-group analysis.
getBootSample <- function(groupDat, group, group.label) {
  bootSamp <- list()
  for (g in seq_along(groupDat)) {
    dat <- groupDat[[g]]
    dat[ , group] <- group.label[g]
    bootSamp[[g]] <- dat[sample(1:nrow(dat), nrow(dat), replace = TRUE), ]
  }
  do.call("rbind", bootSamp)
}

## fit the model to a single bootstrapped sample and return chi-squared
##' @importFrom lavaan lavInspect
fitBootSample <- function(dat, args, suppress) {
  args$data <- dat
  lavaanlavaan <- function(...) { lavaan::lavaan(...) }
  if (suppress) {
    fit <- suppressWarnings(do.call(lavaanlavaan, args))
  } else {
    fit <- do.call(lavaanlavaan, args)
  }
  if (!exists("fit")) return(c(chisq = NA))
  if (lavInspect(fit, "converged")) {
    chisq <- lavInspect(fit, "fit")[c("chisq", "chisq.scaled")]
  } else {
    chisq <- NA
  }
  if (is.na(chisq[2])) return(chisq[1]) else return(chisq[2])
}

