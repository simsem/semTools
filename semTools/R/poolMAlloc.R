### Authors:
### Jason D. Rights (Vanderbilt University; jason.d.rights@vanderbilt.edu)
### - based on research from/with Sonya Sterba
### - adapted from OLD parcelAllocation() by Corbin Quick and Alexander Schoemann
### - additional "indices" argument added by Terrence D. Jorgensen
### Last updated: 10 January 2021


##' Pooled estimates and standard errors across M parcel-allocations: Combining
##' sampling variability and parcel-allocation variability.
##'
##' This function employs an iterative algorithm to pick the number of random
##' item-to-parcel allocations needed to meet user-defined stability criteria
##' for a fitted structural equation model (SEM) (see \bold{Details} below for
##' more information). Pooled point and standard-error estimates from this SEM
##' can be outputted at this final selected number of allocations (however, it
##' is more efficient to save the allocations and treat them as multiple
##' imputations using \code{\link{runMI}}; see \bold{See Also} for links with
##' examples). Additionally, new indices (see Sterba & Rights, 2016) are
##' outputted for assessing the relative contributions of parcel-allocation
##' variability vs. sampling variability in each estimate. At each iteration,
##' this function generates a given number of random item-to-parcel allocations,
##' fits a SEM to each allocation, pools estimates across allocations from that
##' iteration, and then assesses whether stopping criteria are met. If stopping
##' criteria are not met, the algorithm increments the number of allocations
##' used (generating all new allocations).
##'
##' For further details on the benefits of the random allocation of items to
##' parcels, see Sterba (2011) and Sterba & MacCallum (2010).
##'
##' This function implements an algorithm for choosing the number of allocations
##' (\emph{M}; described in Sterba & Rights, 2016), pools point and
##' standard-error estimates across these \emph{M} allocations, and produces
##' indices for assessing the relative contributions of parcel-allocation
##' variability vs. sampling variability in each estimate.
##'
##' To obtain pooled test statistics for model fit or model comparison, the
##' \code{list} or parcel allocations can be passed to \code{\link{runMI}}
##' (find \bold{Examples} on the help pages for \code{\link{parcelAllocation}}
##' and \code{\link{PAVranking}}).
##'
##' This function randomly generates a given number (\code{nAllocStart}) of
##' item-to-parcel allocations, fits a SEM to each allocation, and then
##' increments the number of allocations used (by \code{nAllocAdd}) until the
##' pooled point and standard-error estimates fulfill stopping criteria
##' (\code{stopProp} and \code{stopValue}, defined above). A summary of results
##' from the model that was fit to the \emph{M} allocations are returned.
##'
##' Additionally, this function outputs the proportion of allocations with
##' solutions that converged (using a maximum likelihood estimator) as well as
##' the proportion of allocations with solutions that were converged and proper.
##' The converged and proper solutions among the final \emph{M} allocations are
##' used in computing pooled results.
##'
##' Additionally, after each iteration of the algorithm, information useful in
##' monitoring the algorithm is outputted. The number of allocations used at
##' that iteration, the proportion of pooled parameter estimates meeting
##' stopping criteria at the previous iteration, the proportion of pooled
##' standard errors meeting stopping criteria at the previous iteration, and the
##' runtime of that iteration are outputted. When stopping criteria are
##' satisfied, the full set of results are outputted.
##'
##' @importFrom stats sd pnorm pt qt runif pchisq
##' @importFrom lavaan lavInspect
##'
##' @param nPerPar A list in which each element is a vector, corresponding to
##' each factor, indicating sizes of parcels. If variables are left out of
##' parceling, they should not be accounted for here (i.e., there should not be
##' parcels of size "1").
##' @param facPlc A list of vectors, each corresponding to a factor, specifying
##' the item indicators of that factor (whether included in parceling or not).
##' Either variable names or column numbers. Variables not listed will not be
##' modeled or included in output datasets.
##' @param nAllocStart The number of random allocations of items to parcels to
##' generate in the first iteration of the algorithm.
##' @param nAllocAdd The number of allocations to add with each iteration of the
##' algorithm. Note that if only one iteration is desired, \code{nAllocAdd} can
##' be set to \eqn{0} and results will be output for \code{nAllocStart}
##'  allocationsonly.
##' @param syntax lavaan syntax that defines the model.
##' @param dataset Item-level dataset
##' @param parceloutput Optional \code{character}. Path (folder/directory) where
##' \emph{M} (the final selected number of allocations) parceled data sets will
##' be outputted from the iteration where the algorithm met stopping criteria.
##' Note for Windows users: file path must be specified using forward slashes
##' (\code{/}), not backslashes (\code{\\}). See \code{\link[base]{path.expand}}
##' for details.  If \code{NULL} (default), nothing is saved to disk.
##' @param stopProp Value used in defining stopping criteria of the algorithm
##' (\eqn{\delta_a} in Sterba & Rights, 2016). This is the minimum proportion of
##' change (in any pooled parameter or pooled standard error estimate listed in
##' \code{selectParam}) that is allowable from one iteration of the algorithm to
##' the next. That is, change in pooled estimates and pooled standard errors
##' from one iteration to the next must all be less than (\code{stopProp}) x
##' (value from former iteration). Note that \code{stopValue} can override this
##' criterion (see below). Also note that values less than .01 are unlikely to
##' lead to more substantively meaningful precision. Also note that if only
##' \code{stopValue} is a desired criterion, \code{stopProp} can be set to 0.
##' @param stopValue Value used in defining stopping criteria of the algorithm
##' (\eqn{\delta_b} in Sterba & Rights, 2016). \code{stopValue} is a minimum
##' allowable amount of absolute change (in any pooled parameter or pooled
##' standard error estimate listed in \code{selectParam}) from one iteration of
##' the algorithm to the next. For a given pooled estimate or pooled standard
##' error, \code{stopValue} is only invoked as a stopping criteria when the
##' minimum change required by \code{stopProp} is less than \code{stopValue}.
##' Note that values less than .01 are unlikely to lead to more substantively
##' meaningful precision. Also note that if only \code{stopProp} is a desired
##' criterion, \code{stopValue} can be set to 0.
##' @param selectParam (Optional) A list of the pooled parameters to be used in
##' defining stopping criteria (i.e., \code{stopProp} and \code{stopValue}).
##' These parameters should appear in the order they are listed in the lavaan
##' syntax. By default, all pooled parameters are used. Note that
##' \code{selectParam} should only contain freely-estimated parameters. In one
##' example from Sterba & Rights (2016) \code{selectParam} included all free
##' parameters except item intercepts and in another example \code{selectParam}
##' included only structural parameters.
##' @param indices Optional \code{character} vector indicating the names of
##' available \code{\link[lavaan]{fitMeasures}} to be included in the output.
##' The first and second elements should be a chi-squared test statistic and its
##' associated degrees of freedom, both of which will be added if missing. If
##' \code{"default"}, the indices will be \code{c("chisq", "df", "cfi", "tli",
##' "rmsea","srmr")}. If a robust test statistic is requested (see
##' \code{\link[lavaan]{lavOptions}}), \code{c("chisq","df")} will be replaced
##' by \code{c("chisq.scaled","df.scaled")}. For the output to include both the
##' naive and robust test statistics, \code{indices} should include both, but
##' put the scaled test statistics first, as in \code{indices =
##' c("chisq.scaled", "df.scaled", "chisq", "df")}
##' @param double (Optional) If set to \code{TRUE}, requires stopping criteria
##' (\code{stopProp} and \code{stopValue}) to be met for all parameters (in
##' \code{selectParam}) for two consecutive iterations of the algorithm. By
##' default, this is set to \code{FALSE}, meaning stopping criteria need only be
##' met at one iteration of the algorithm.
##' @param names (Optional) A character vector containing the names of parceled
##' variables.
##' @param leaveout (Optional) A vector of variables to be left out of
##' randomized parceling. Either variable names or column numbers are allowed.
##' @param useTotalAlloc (Optional) If set to \code{TRUE}, function will output
##' a separate set of results that uses all allocations created by the
##' algorithm, rather than \emph{M} allocations (see "Allocations needed for
##' stability" below). This distinction is further discussed in Sterba and
##' Rights (2016).
##' @param checkConv (Optional) If set to TRUE, function will output pooled
##' estimates and standard errors from 10 iterations post-convergence.
##' @param \dots Additional arguments to be passed to
##' \code{\link[lavaan]{lavaan}}. See also \code{\link[lavaan]{lavOptions}}
##'
##' @return
##' \item{Estimates}{A table containing pooled results across \emph{M}
##' allocations at the iteration where stopping criteria were met. Columns
##' correspond to individual parameter name, pooled estimate, pooled standard
##' error, \emph{p}-value for a \emph{z}-test of the parameter, \emph{z}-based
##' 95\% confidence interval, \emph{p}-value for a \emph{t}-test of the
##' parameter (using degrees of freedom described in Sterba & Rights, 2016), and
##' \emph{t}-based 95\% confidence interval for the parameter.}
##' \item{Fit}{A table containing results related to model fit from the \emph{M}
##' allocations at the iteration where stopping criteria were met. Columns
##' correspond to fit index names, the average of each index across allocations,
##' the standard deviation of each fit index across allocations, the maximum of
##' each fit index across allocations, the minimum of each fit index across
##' allocations, the range of each fit index across allocations, and the percent
##' of the \emph{M} allocations where the chi-square test of absolute fit was
##' significant.}
##' \item{Proportion of converged and proper allocations}{A table
##' containing the proportion of the final \emph{M} allocations that converged
##' (using a maximum likelihood estimator) and the proportion of allocations
##' that converged to proper solutions. Note that pooled estimates, pooled
##' standard errors, and other results are computed using only the converged,
##' proper allocations.}
##' \item{Allocations needed for stability (M)}{The number of allocations
##' (\emph{M}) at which the algorithm's stopping criteria (defined above) were
##' met.}
##' \item{Indices used to quantify uncertainty in estimates due to sample vs.
##' allocation variability}{A table containing individual parameter names, an
##' estimate of the proportion of total variance of a pooled parameter estimate
##' that is attributable to parcel-allocation variability (PPAV), and an estimate
##' of the ratio of the between-allocation variance of a pooled parameter
##' estimate to the within-allocation variance (RPAV). See Sterba & Rights (2016)
##' for more detail.}
##' \item{Total runtime (minutes)}{The total runtime of the function, in minutes.
##' Note that the total runtime will be greater when the specified model
##' encounters convergence problems for some allocations, as is the case with the
##' \code{\link{simParcel}} dataset used below.}
##'
##' @author
##' Jason D. Rights (Vanderbilt University; \email{jason.d.rights@@vanderbilt.edu})
##'
##' The author would also like to credit Corbin Quick and Alexander Schoemann
##' for providing the original parcelAllocation function on which this function
##' is based.
##'
##' @seealso
##'   \code{\link{runMI}} for treating allocations as multiple imputations to
##'   pool results across allocations. See \bold{Examples} on help pages for:
##'   \itemize{
##'     \item{\code{\link{parcelAllocation}} for fitting a single model}
##'     \item{\code{\link{PAVranking}} for comparing 2 models}
##'   }
##'
##' @references
##'
##' Sterba, S. K. (2011). Implications of parcel-allocation
##' variability for comparing fit of item-solutions and parcel-solutions.
##' \emph{Structural Equation Modeling, 18}(4), 554--577.
##' \doi{10.1080/10705511.2011.607073}
##'
##' Sterba, S. K., & MacCallum, R. C. (2010). Variability in parameter estimates
##' and model fit across random allocations of items to parcels.
##' \emph{Multivariate Behavioral Research, 45}(2), 322--358.
##' \doi{10.1080/00273171003680302}
##'
##' Sterba, S. K., & Rights, J. D. (2016). Accounting for parcel-allocation
##' variability in practice: Combining sources of uncertainty and choosing the
##' number of allocations. \emph{Multivariate Behavioral Research, 51}(2--3),
##' 296--313. \doi{10.1080/00273171.2016.1144502}
##'
##' Sterba, S. K., & Rights, J. D. (2017). Effects of parceling on model
##' selection: Parcel-allocation variability in model ranking.
##' \emph{Psychological Methods, 22}(1), 47--68. \doi{10.1037/met0000067}
##'
##' @examples
##'
##' \dontrun{
##' ## lavaan syntax: A 2 Correlated
##' ## factor CFA model to be fit to parceled data
##'
##' parmodel <- '
##'    f1 =~ NA*p1f1 + p2f1 + p3f1
##'    f2 =~ NA*p1f2 + p2f2 + p3f2
##'    p1f1 ~ 1
##'    p2f1 ~ 1
##'    p3f1 ~ 1
##'    p1f2 ~ 1
##'    p2f2 ~ 1
##'    p3f2 ~ 1
##'    p1f1 ~~ p1f1
##'    p2f1 ~~ p2f1
##'    p3f1 ~~ p3f1
##'    p1f2 ~~ p1f2
##'    p2f2 ~~ p2f2
##'    p3f2 ~~ p3f2
##'    f1 ~~ 1*f1
##'    f2 ~~ 1*f2
##'    f1 ~~ f2
##' '
##'
##' ## specify items for each factor
##' f1name <- colnames(simParcel)[1:9]
##' f2name <- colnames(simParcel)[10:18]
##'
##' ## run function
##' poolMAlloc(nPerPar = list(c(3,3,3), c(3,3,3)),
##'            facPlc = list(f1name, f2name), nAllocStart = 10, nAllocAdd = 10,
##'            syntax = parmodel, dataset = simParcel, stopProp = .03,
##'            stopValue = .03, selectParam = c(1:6, 13:18, 21),
##'            names = list("p1f1","p2f1","p3f1","p1f2","p2f2","p3f2"),
##'            double = FALSE, useTotalAlloc = FALSE)
##' }
##'
##' ## See examples on ?parcelAllocation and ?PAVranking for how to obtain
##' ## pooled test statistics and other pooled lavaan output.
##' ## Details provided in Sterba & Rights (2016).
##'
##' @export
poolMAlloc <- function(nPerPar, facPlc, nAllocStart, nAllocAdd = 0,
                       parceloutput = NULL, syntax, dataset, stopProp, stopValue,
                       selectParam = NULL, indices = "default", double = FALSE,
                       checkConv = FALSE, names = "default", leaveout = 0,
                       useTotalAlloc = FALSE, ...) {
  message('Note that more options for pooled results are available using the ',
          'runMI() function (see Examples on ?parcelAllocation and ?PAVranking)')
  if (!is.null(parceloutput)) {
    if (!dir.exists(parceloutput)) stop('invalid directory:\n',
                                        paste(parceloutput), "\n\n")
  }

  StartTimeFull <- proc.time()
  if (is.character(dataset)) dataset <- utils::read.csv(dataset)
  if (indices[1] == "default") indices <- c("chisq", "df", "cfi", "tli", "rmsea","srmr")
  ## make sure chi-squared and df are the first and second elements
  requestedChiSq <- grep(pattern = "chisq", indices, value = TRUE)
  if (length(requestedChiSq) == 0L) {
    indices <- unique(c("chisq", indices))
  } else {
    indices <- unique(c(requestedChiSq[1], indices))
  }
  requestedDF <- grep(pattern = "df", indices, value = TRUE)
  if (length(requestedDF) == 0L) {
    indices <- unique(c(indices[1], "df", indices[-1]))
  } else {
    indices <- unique(c(indices[1], requestedDF[1], indices[-1]))
  }

  isProperSolution <- function(object) {
    lavpartable <- object@ParTable
    lavfit <- object@Fit
    lavdata <- object@Data
    lavmodel <- object@Model
    var.idx <- which(lavpartable$op == "~~" & lavpartable$lhs == lavpartable$rhs)
    if (length(var.idx) > 0L && any(lavfit@est[var.idx] < 0)) return(FALSE)
    if (length(lavaan::lavaanNames(lavpartable, type = "lv.regular")) > 0L) {
      ETA <- list(lavInspect(object, "cov.lv"))
      for (g in 1:lavdata@ngroups) {
        eigvals <- eigen(ETA[[g]], symmetric = TRUE, only.values = TRUE)$values
        if (any(eigvals < -1 * .Machine$double.eps^(3/4))) return(FALSE)
      }
    }
    THETA <- list(lavInspect(object, "theta"))
    for (g in 1:lavdata@ngroups) {
      num.idx <- lavmodel@num.idx[[g]]
      if (length(num.idx) > 0L) {
        eigvals <- eigen(THETA[[g]][unlist(num.idx),
                                    unlist(num.idx), drop = FALSE], symmetric = TRUE,
                         only.values = TRUE)$values
        if (any(eigvals < -1 * .Machine$double.eps^(3/4))) return(FALSE)
      }
    }
    TRUE
  }

  nloop <- 0
  nAllocStarttemp <- nAllocStart
  options(max.print = 1e+06)
  BreakCounter <- NA
  repeat {
    StartTime <- proc.time()
    nloop <- nloop + 1
    if (double == TRUE & is.na(BreakCounter) == FALSE)
      BreakCounter <- BreakCounter + 1
    if (checkConv == TRUE & is.na(BreakCounter) == FALSE)
      BreakCounter <- BreakCounter + 1
    if (nloop > 1) {
      if (is.na(BreakCounter) == TRUE) {
        Parmn_revFinal <- Parmn_rev[[nloop - 1]]
        nConvergedOutput <- nConverged
        nConvergedProperOutput <- nConvergedProper
        PooledSEwithinvarFinal <- PooledSEwithinvar
        PooledSEbetweenvarFinal <- PooledSEbetweenvar
        PooledSEFinal <- PooledSE
        FitsumOutput <- Fitsum
        nAllocOutput <- nAllocStart - nAllocAdd
        AllocationsOutput <- Allocations
        #ParamFinal <- Param # defined, but never used
      }
      ParamPooledSE_temp <- ParamPooledSE
      ParamTest_temp <- ParamTest
      PooledSE_temp <- PooledSE
      ParamPoolSEdiffmin <- abs(ParamPooledSE_temp * stopProp)
      ParamPoolSEdiffmin[ParamPoolSEdiffmin < stopValue] <- stopValue
      ParamDiffMin <- abs(ParamTest * stopProp)
      ParamDiffMin[ParamDiffMin < stopValue] <- stopValue
      PooledSEmin <- abs(PooledSE * stopProp)
      PooledSEmin[PooledSEmin < stopValue] <- stopValue
    }
    dataset <- as.matrix(dataset)
    if (nAllocStart < 2) stop("Minimum of two allocations required.")
    if (is.list(facPlc)) {
      if (is.numeric(facPlc[[1]][1]) == FALSE) {
        facPlcb <- facPlc
        Namesv <- colnames(dataset)
        for (i in 1:length(facPlc)) {
          for (j in 1:length(facPlc[[i]])) {
            facPlcb[[i]][j] <- match(facPlc[[i]][j],
                                     Namesv)
          }
          facPlcb[[i]] <- as.numeric(facPlcb[[i]])
        }
        facPlc <- facPlcb
      }
      facPlc2 <- rep(0, ncol(dataset))
      for (i in 1:length(facPlc)) {
        for (j in 1:length(facPlc[[i]])) {
          facPlc2[facPlc[[i]][j]] <- i
        }
      }
      facPlc <- facPlc2
    }
    if (leaveout != 0) {
      if (is.numeric(leaveout) == FALSE) {
        leaveoutb <- rep(0, length(leaveout))
        Namesv <- colnames(dataset)
        for (i in 1:length(leaveout)) {
          leaveoutb[i] <- match(leaveout[i], Namesv)
        }
        leaveout <- as.numeric(leaveoutb)
      }
      k1 <- 0.001
      for (i in 1:length(leaveout)) {
        facPlc[leaveout[i]] <- facPlc[leaveout[i]] + k1
        k1 <- k1 + 0.001
      }
    }
    if (0 %in% facPlc == TRUE) {
      Zfreq <- sum(facPlc == 0)
      for (i in 1:Zfreq) {
        Zplc <- match(0, facPlc)
        dataset <- dataset[, -Zplc]
        facPlc <- facPlc[-Zplc]
      }
    }
    if (is.list(nPerPar)) {
      nPerPar2 <- c()
      for (i in 1:length(nPerPar)) {
        Onesp <- sum(facPlc > i & facPlc < i + 1)
        nPerPar2 <- c(nPerPar2, nPerPar[i], rep(1, Onesp), recursive = TRUE)
      }
      nPerPar <- nPerPar2
    }
    Npp <- c()
    for (i in 1:length(nPerPar)) {
      Npp <- c(Npp, rep(i, nPerPar[i]))
    }
    Locate <- sort(round(facPlc))
    Maxv <- max(Locate) - 1
    if (length(Locate) != length(Npp)) {
      stop("** ERROR! ** Parcels incorrectly specified. Check input!")
    }
    if (Maxv > 0) {
      for (i in 1:Maxv) {
        Mat <- match(i + 1, Locate)
        if (Npp[Mat] == Npp[Mat - 1]) {
          stop("** ERROR! ** Parcels incorrectly specified. Check input!")
        }
      }
    }
    Onevec <- facPlc - round(facPlc)
    NleaveA <- length(Onevec) - sum(Onevec == 0)
    NleaveP <- sum(nPerPar == 1)
    if (NleaveA < NleaveP) {
      warning("** WARNING! ** Single-variable parcels have been requested.",
              " Check input!")
    }
    if (NleaveA > NleaveP)
      warning("** WARNING! ** More non-parceled variables have been requested",
              " than provided for in parcel vector. Check input!")
    if (length(names) > 1) {
      if (length(names) != length(nPerPar)) {
        warning("** WARNING! ** Number of parcel names provided not equal to",
                " number of parcels requested. Check input!")
      }
    }
    Data <- c(1:ncol(dataset))
    # Nfactors <- max(facPlc) # defined but never used
    Nindicators <- length(Data)
    Npar <- length(nPerPar)
    Rmize <- runif(Nindicators, 1, Nindicators)
    Data <- rbind(facPlc, Rmize, Data)
    Results <- matrix(numeric(0), nAllocStart, Nindicators)
    Pin <- nPerPar[1]
    for (i in 2:length(nPerPar)) {
      Pin <- c(Pin, nPerPar[i] + Pin[i - 1])
    }
    for (i in 1:nAllocStart) {
      Data[2, ] <- runif(Nindicators, 1, Nindicators)
      Data <- Data[, order(Data[2, ])]
      Data <- Data[, order(Data[1, ])]
      Results[i, ] <- Data[3, ]
    }
    Alpha <- rbind(Results[1, ], dataset)
    Allocations <- list()
    for (i in 1:nAllocStart) {
      Ineff <- rep(NA, ncol(Results))
      Ineff2 <- c(1:ncol(Results))
      for (inefficient in 1:ncol(Results)) {
        Ineff[Results[i, inefficient]] <- Ineff2[inefficient]
      }
      Alpha[1, ] <- Ineff
      Beta <- Alpha[, order(Alpha[1, ])]
      Temp <- matrix(NA, nrow(dataset), Npar)
      TempAA <- if (length(1:Pin[1]) > 1) {
        Beta[2:nrow(Beta), 1:Pin[1]]
      } else cbind(Beta[2:nrow(Beta), 1:Pin[1]], Beta[2:nrow(Beta), 1:Pin[1]])
      Temp[, 1] <- rowMeans(TempAA, na.rm = TRUE)
      for (al in 2:Npar) {
        Plc <- Pin[al - 1] + 1
        TempBB <- if (length(Plc:Pin[al]) > 1) {
          Beta[2:nrow(Beta), Plc:Pin[al]]
        } else cbind(Beta[2:nrow(Beta), Plc:Pin[al]],
                     Beta[2:nrow(Beta), Plc:Pin[al]])
        Temp[, al] <- rowMeans(TempBB, na.rm = TRUE)
      }
      if (length(names) > 1) {
        colnames(Temp) <- names
      }
      Allocations[[i]] <- Temp
    }
    Param <- list()
    Fitind <- list()
    Converged <- list()
    ProperSolution <- list()
    ConvergedProper <- list()
    for (i in 1:(nAllocStart)) {
      data_parcel <- as.data.frame(Allocations[[i]], row.names = NULL, optional = FALSE)
      fit <- lavaan::sem(syntax, data = data_parcel, ...)
      ## if a robust estimator was requested, update fit indices accordingly
      requestedTest <- lavInspect(fit, "options")$test
      if (any(requestedTest %in% c("satorra.bentler","yuan.bentler",
                                   "yuan.bentler.mplus","scaled.shifted",
                                   "mean.var.adjusted","satterthwaite"))) {
        indices[1:2] <- c("chisq.scaled","df.scaled")
      } else indices[1:2] <- c("chisq","df")
      ## check convergence and solution
      if (lavInspect(fit, "converged") == TRUE) {
        Converged[[i]] <- 1
      } else Converged[[i]] <- 0
      Param[[i]] <- lavaan::parameterEstimates(fit)[,
                                                    c("lhs", "op", "rhs", "est", "se", "z", "pvalue",
                                                      "ci.lower", "ci.upper")]
      if (isProperSolution(fit) == TRUE & Converged[[i]] == 1) {
        ProperSolution[[i]] <- 1
      } else ProperSolution[[i]] <- 0
      if (any(is.na(Param[[i]][, 5] == TRUE)))
        ProperSolution[[i]] <- 0
      if (Converged[[i]] == 1 & ProperSolution[[i]] == 1) {
        ConvergedProper[[i]] <- 1
      } else ConvergedProper[[i]] <- 0
      if (ConvergedProper[[i]] == 0)
        Param[[i]][, 4:9] <- matrix(data = NA, nrow(Param[[i]]), 6)
      if (ConvergedProper[[i]] == 1) {
        Fitind[[i]] <- lavaan::fitMeasures(fit, indices)
        if (!all(indices %in% names(Fitind[[i]]))) {
          invalidIndices <- setdiff(indices, names(Fitind[[i]]))
          Fitind[[i]][invalidIndices] <- NA
        }
      } else Fitind[[i]] <- rep(NA, length(indices))
    }
    nConverged <- Reduce("+", Converged)
    nProperSolution <- Reduce("+", ProperSolution)
    nConvergedProper <- Reduce("+", ConvergedProper)
    if (nConvergedProper == 0) stop("All allocations failed to converge and/or",
                                    " yielded improper solutions for a given loop.")
    Parmn <- Param[[1]]
    if (is.null(selectParam))
      selectParam <- 1:nrow(Parmn)
    ParSE <- matrix(NA, nrow(Parmn), nAllocStart)
    ParSEmn <- Parmn[, 5]
    Parsd <- matrix(NA, nrow(Parmn), nAllocStart)
    Fitmn <- Fitind[[1]]
    Fitsd <- matrix(NA, length(Fitmn), nAllocStart)
    Sigp <- matrix(NA, nrow(Parmn), nAllocStart)
    Fitind <- data.frame(Fitind)
    ParamSEsquared <- list()
    for (i in 1:nAllocStart) {
      ParamSEsquared[[i]] <- cbind(Param[[i]][, 5], Param[[i]][, 5])
      if (any(is.na(ParamSEsquared[[i]]) == TRUE)) ParamSEsquared[[i]] <- 0
      ParamSEsquared[[i]] <- apply(as.data.frame(ParamSEsquared[[i]]), 1, prod)
      Parsd[, i] <- Param[[i]][, 4]
      ParSE[, i] <- Param[[i]][, 5]
      Sigp[, ncol(Sigp) - i + 1] <- Param[[i]][, 7]
      Fitsd[, i] <- Fitind[[i]]
    }
    Sigp <- Sigp + 0.45
    Sigp <- apply(Sigp, c(1, 2), round)
    Sigp <- 1 - as.vector(rowMeans(Sigp, na.rm = TRUE))
    Parsum <- cbind(apply(Parsd, 1, mean, na.rm = TRUE),
                    apply(Parsd, 1, sd, na.rm = TRUE),
                    apply(Parsd, 1, max, na.rm = TRUE),
                    apply(Parsd, 1, min, na.rm = TRUE),
                    apply(Parsd, 1, max, na.rm = TRUE) - apply(Parsd, 1, min, na.rm = TRUE),
                    Sigp)
    colnames(Parsum) <- c("Avg Est.", "S.D.", "MAX",
                          "MIN", "Range", "% Sig")
    ParSEmn <- Parmn[, 1:3]
    ParSEfn <- cbind(ParSEmn, apply(ParSE, 1, mean, na.rm = TRUE),
                     apply(ParSE, 1, sd, na.rm = TRUE),
                     apply(ParSE, 1, max, na.rm = TRUE),
                     apply(ParSE, 1, min, na.rm = TRUE),
                     apply(ParSE, 1, max, na.rm = TRUE) - apply(ParSE, 1, min, na.rm = TRUE))
    colnames(ParSEfn) <- c("lhs", "op", "rhs", "Avg SE",
                           "S.D.", "MAX", "MIN", "Range")
    Fitsum <- cbind(apply(Fitsd, 1, mean, na.rm = TRUE),
                    apply(Fitsd, 1, sd, na.rm = TRUE),
                    apply(Fitsd, 1, max, na.rm = TRUE),
                    apply(Fitsd, 1, min, na.rm = TRUE),
                    apply(Fitsd, 1, max, na.rm = TRUE) - apply(Fitsd, 1, min, na.rm = TRUE))
    rownames(Fitsum) <- indices
    Parmn[, 4:ncol(Parmn)] <- Parmn[, 4:ncol(Parmn)]/nConvergedProper
    Parmn <- Parmn[, 1:3]
    Parmn <- cbind(Parmn, Parsum)
    Fitmn <- Fitmn/nConvergedProper
    pChisq <- list()
    sigChisq <- list()
    for (i in 1:nAllocStart) {
      pChisq[[i]] <- (1 - pchisq(Fitsd[1, i], Fitsd[2, i]))
      if (is.na(pChisq[[i]]) == FALSE & pChisq[[i]] < 0.05) {
        sigChisq[[i]] <- 1
      }
      else sigChisq[[i]] <- 0
    }
    PerSigChisq <- (Reduce("+", sigChisq))/nConvergedProper * 100
    PerSigChisq <- round(PerSigChisq, 4)
    PerSigChisqCol <- c(PerSigChisq, # however many indices != chisq(.scaled)
                        rep("n/a", sum(!grepl(pattern = "chisq", x = indices))))
    options(stringsAsFactors = FALSE)
    Fitsum <- data.frame(Fitsum, PerSigChisqCol)
    colnames(Fitsum) <- c("Avg Ind", "S.D.", "MAX", "MIN",
                          "Range", "% Sig")
    options(stringsAsFactors = TRUE)
    PooledSEwithinvar <- Reduce("+", ParamSEsquared)/nConvergedProper
    PooledSEbetweenvar <- Parmn[, 5]^2
    PooledSE <- sqrt(PooledSEwithinvar + PooledSEbetweenvar + PooledSEbetweenvar/nConvergedProper)
    ParamPooledSE <- c(Parmn[, 4], PooledSE)
    ParamTest <- Parmn[, 4]
    if (nloop > 1) {
      ParamPoolSEdiff <- abs(ParamPooledSE_temp - ParamPooledSE)
      Paramdiff <- abs(ParamTest_temp - ParamTest)
      PooledSEdiff <- abs(PooledSE - PooledSE_temp)
      ParamPoolSEdifftest <- ParamPoolSEdiff - ParamPoolSEdiffmin
      ParamPoolSEdifftest[ParamPoolSEdifftest <= 0] <- 0
      ParamPoolSEdifftest[ParamPoolSEdifftest > 0] <- 1
      Paramdifftest <- Paramdiff - ParamDiffMin
      Paramdifftest[Paramdifftest <= 0] <- 0
      Paramdifftest[Paramdifftest > 0] <- 1
      PooledSEdifftest <- PooledSEdiff - PooledSEmin
      PooledSEdifftest[PooledSEdifftest <= 0] <- 0
      PooledSEdifftest[PooledSEdifftest > 0] <- 1
      if (nloop == 2) {
        ParamPoolSEdifftesttable <- cbind(ParamPoolSEdifftest)
        Paramdifftesttable <- cbind(Paramdifftest)
        PooledSEdifftesttable <- cbind(PooledSEdifftest)
      }
      if (nloop > 2) {
        ParamPoolSEdifftesttable <- cbind(ParamPoolSEdifftesttable,
                                          ParamPoolSEdifftest)
        Paramdifftesttable <- cbind(Paramdifftesttable,
                                    Paramdifftest)
        PooledSEdifftesttable <- cbind(PooledSEdifftesttable,
                                       PooledSEdifftest)
      }
      PropStopParam <- 1 - (Reduce("+", Paramdifftesttable[selectParam,
                                                           nloop - 1])/length(selectParam))
      PropStopPooled <- 1 - (Reduce("+", PooledSEdifftesttable[selectParam,
                                                               nloop - 1])/length(selectParam))
      PropStopParamPooled <- 1 - (Reduce("+", ParamPoolSEdifftesttable[c(selectParam, selectParam + nrow(Parmn)), nloop - 1]) /
                                    (2 * length(selectParam)))
      if (checkConv == TRUE & is.na(BreakCounter) == TRUE) {
        print(nAllocStart)
        print("Proportion of pooled estimates meeting stop criteria:")
        print(PropStopParam)
        print("Proportion of pooled SE meeting stop criteria:")
        print(PropStopPooled)
      }
      if (checkConv == FALSE) {
        print(nAllocStart)
        print("Proportion of pooled estimates meeting stop criteria:")
        print(PropStopParam)
        print("Proportion of pooled SE meeting stop criteria:")
        print(PropStopPooled)
      }
    }
    nAllocStart <- nAllocStart + nAllocAdd
    StopTime <- proc.time() - StartTime
    print("Runtime:")
    print(StopTime)
    Parmn_rev <- list()
    Parmn_rev[[nloop]] <- cbind(Parmn[, 1:4], PooledSE)
    Parmn_rev[[nloop]][, 4:5] <- sapply(Parmn_rev[[nloop]][,4:5], as.numeric)
    colnames(Parmn_rev[[nloop]]) <- c("lhs", "op", "rhs","Estimate", "Pooled SE")
    if (nloop == 1) {
      Param_revTemp <- cbind(Parmn[, 1:3], Parmn_rev[[nloop]][,4])
      Param_revTemp[, 4] <- as.numeric(Param_revTemp[,4])
      Param_revTotal <- cbind(Param_revTemp)
      PooledSE_revTemp <- cbind(Parmn[, 1:3], Parmn_rev[[nloop]][,5])
      PooledSE_revTemp[, 4] <- as.numeric(PooledSE_revTemp[,4])
      PooledSE_revTotal <- cbind(PooledSE_revTemp)
    }
    if (nloop > 1) {
      Param_revTemp <- cbind(Parmn_rev[[nloop]][, 4])
      Param_revTemp <- as.numeric(Param_revTemp)
      Param_revTotal <- cbind(Param_revTotal, Param_revTemp)
      PooledSE_revTemp <- cbind(Parmn_rev[[nloop]][,
                                                   5])
      PooledSE_revTemp <- as.numeric(PooledSE_revTemp)
      PooledSE_revTotal <- cbind(PooledSE_revTotal,
                                 PooledSE_revTemp)
    }
    if (nloop == 1) {
      ParamTotal <- Param
      FitindTotal <- Fitind
      AllocationsTotal <- Allocations
      nAllocTotal <- nAllocStart - nAllocAdd
      nConvergedTotal <- nConverged
      nProperSolutionTotal <- nProperSolution
      nConvergedProperTotal <- nConvergedProper
    }
    if (nloop > 1) {
      ParamTotal <- c(ParamTotal, Param)
      FitindTotal <- c(FitindTotal, Fitind)
      AllocationsTotal <- c(AllocationsTotal, Allocations)
      nAllocTotal <- nAllocTotal + nAllocStart - nAllocAdd
      nConvergedTotal <- nConverged + nConvergedTotal
      nProperSolution <- nProperSolution + nProperSolutionTotal
      nConvergedProperTotal <- nConvergedProper + nConvergedProperTotal
    }
    if (nloop > 1 & double == TRUE & is.na(BreakCounter) == FALSE & BreakCounter == 2) {
      if (Reduce("+", ParamPoolSEdifftesttable[c(selectParam,
                                                 selectParam + nrow(Parmn_rev[[nloop]])), nloop - 1]) == 0)
        break
    }
    if (nloop > 1 & double == TRUE) {
      if (Reduce("+", ParamPoolSEdifftesttable[c(selectParam,
                                                 selectParam + nrow(Parmn_rev[[nloop]])), nloop - 1]) == 0) {
        BreakCounter <- 1
      }
      else BreakCounter <- NA
    }
    if (nloop > 1 & checkConv == TRUE & is.na(BreakCounter) == TRUE) {
      if (Reduce("+", ParamPoolSEdifftesttable[c(selectParam,
                                                 selectParam + nrow(Parmn_rev[[nloop]])), nloop - 1]) == 0)
        BreakCounter <- 0
    }
    if (nloop > 1 & double == FALSE & checkConv == FALSE) {
      if (Reduce("+", ParamPoolSEdifftesttable[c(selectParam,
                                                 selectParam + nrow(Parmn_rev[[nloop]])), nloop - 1]) == 0)
        break
    }
    if (nAllocAdd == 0)
      break
    if (checkConv == TRUE & is.na(BreakCounter) == FALSE & BreakCounter == 9)
      break
  }
  if (nAllocAdd == 0) {
    Parmn_revFinal <- Parmn_rev[[nloop]]
    nConvergedOutput <- nConverged
    nConvergedProperOutput <- nConvergedProper
    PooledSEwithinvarFinal <- PooledSEwithinvar
    PooledSEbetweenvarFinal <- PooledSEbetweenvar
    PooledSEFinal <- PooledSE
    FitsumOutput <- Fitsum
    nAllocOutput <- nAllocStart - nAllocAdd
    AllocationsOutput <- Allocations
  }
  if (!is.null(parceloutput)) {
    replist <- matrix(NA, nAllocOutput, 1)
    for (i in 1:(nAllocOutput)) {
      colnames(AllocationsOutput[[i]]) <- names
      utils::write.table(AllocationsOutput[[i]],
                         file = paste(parceloutput, "/parcelruns", i, ".dat", sep = ""),
                         row.names = FALSE, col.names = TRUE)
      replist[i, 1] <- paste("parcelruns", i, ".dat", sep = "")
    }
    utils:: write.table(replist, paste(parceloutput, "/parcelrunsreplist.dat",
                                       sep = ""), quote = FALSE, row.names = FALSE,
                        col.names = FALSE)
  }
  if (useTotalAlloc == TRUE) {
    ParmnTotal <- ParamTotal[[1]]
    ParSETotal <- matrix(NA, nrow(ParmnTotal), nAllocTotal)
    ParSEmnTotal <- ParmnTotal[, 5]
    ParsdTotal <- matrix(NA, nrow(ParmnTotal), nAllocTotal)
    FitmnTotal <- FitindTotal[[1]]
    FitsdTotal <- matrix(NA, length(FitmnTotal), nAllocTotal)
    SigpTotal <- matrix(NA, nrow(ParmnTotal), nAllocTotal)
    FitindTotal <- data.frame(FitindTotal)
    ParamSEsquaredTotal <- list()
    for (i in 1:nAllocTotal) {
      ParamSEsquaredTotal[[i]] <- cbind(ParamTotal[[i]][,5], ParamTotal[[i]][, 5])
      if (any(is.na(ParamSEsquaredTotal[[i]]) == TRUE))
        ParamSEsquaredTotal[[i]] <- 0
      ParamSEsquaredTotal[[i]] <- apply(as.data.frame(ParamSEsquaredTotal[[i]]),1, prod)
      ParsdTotal[, i] <- ParamTotal[[i]][, 4]
      ParSETotal[, i] <- ParamTotal[[i]][, 5]
      SigpTotal[, ncol(Sigp) - i + 1] <- ParamTotal[[i]][,7]
      FitsdTotal[, i] <- FitindTotal[[i]]
    }
    SigpTotal <- SigpTotal + 0.45
    SigpTotal <- apply(SigpTotal, c(1, 2), round)
    SigpTotal <- 1 - as.vector(rowMeans(SigpTotal, na.rm = TRUE))
    ParsumTotal <- cbind(apply(ParsdTotal, 1, mean, na.rm = TRUE),
                         apply(ParsdTotal, 1, sd, na.rm = TRUE),
                         apply(ParsdTotal, 1, max, na.rm = TRUE),
                         apply(ParsdTotal, 1, min, na.rm = TRUE),
                         apply(ParsdTotal, 1, max, na.rm = TRUE) - apply(ParsdTotal, 1, min, na.rm = TRUE),
                         SigpTotal)
    colnames(ParsumTotal) <- c("Avg Est.", "S.D.", "MAX", "MIN", "Range", "% Sig")
    ParSEmnTotal <- ParmnTotal[, 1:3]
    ParSEfnTotal <- cbind(ParSEmnTotal,
                          apply(ParSETotal, 1, mean, na.rm = TRUE),
                          apply(ParSETotal, 1, sd, na.rm = TRUE),
                          apply(ParSETotal, 1, max, na.rm = TRUE),
                          apply(ParSETotal, 1, min, na.rm = TRUE),
                          apply(ParSETotal, 1, max, na.rm = TRUE) - apply(ParSETotal, 1, min, na.rm = TRUE))
    colnames(ParSEfnTotal) <- c("lhs", "op", "rhs", "Avg SE",
                                "S.D.", "MAX", "MIN", "Range")
    FitsumTotal <- cbind(apply(FitsdTotal, 1, mean, na.rm = TRUE),
                         apply(FitsdTotal, 1, sd, na.rm = TRUE),
                         apply(FitsdTotal, 1, max, na.rm = TRUE),
                         apply(FitsdTotal, 1, min, na.rm = TRUE),
                         apply(FitsdTotal, 1, max, na.rm = TRUE) - apply(FitsdTotal, 1, min, na.rm = TRUE))
    rownames(FitsumTotal) <- indices
    ParmnTotal[, 4:ncol(ParmnTotal)] <- ParmnTotal[,4:ncol(Parmn)]/nConvergedProperTotal
    ParmnTotal <- ParmnTotal[, 1:3]
    ParmnTotal <- cbind(ParmnTotal, ParsumTotal)
    FitmnTotal <- FitmnTotal/nConvergedProperTotal
    pChisqTotal <- list()
    sigChisqTotal <- list()
    for (i in 1:nAllocTotal) {
      pChisqTotal[[i]] <- (1 - pchisq(FitsdTotal[1,i], FitsdTotal[2, i]))
      if (is.na(pChisqTotal[[i]]) == FALSE & pChisqTotal[[i]] < 0.05) {
        sigChisqTotal[[i]] <- 1
      } else sigChisqTotal[[i]] <- 0
    }
    PerSigChisqTotal <- (Reduce("+", sigChisqTotal))/nConvergedProperTotal * 100
    PerSigChisqTotal <- round(PerSigChisqTotal, 4)
    PerSigChisqColTotal <- c(PerSigChisqTotal, "n/a", "n/a", "n/a", "n/a")
    options(stringsAsFactors = FALSE)
    FitsumTotal <- data.frame(FitsumTotal, PerSigChisqColTotal)
    colnames(FitsumTotal) <- c("Avg Ind", "S.D.", "MAX", "MIN", "Range", "% Sig")
    options(stringsAsFactors = TRUE)
    PooledSEwithinvarTotal <- Reduce("+", ParamSEsquaredTotal)/nConvergedProperTotal
    PooledSEbetweenvarTotal <- ParmnTotal[, 5]^2
    PooledSETotal <- sqrt(PooledSEwithinvarTotal + PooledSEbetweenvarTotal +
                            PooledSEbetweenvarTotal/nConvergedProperTotal)
    ParamPooledSETotal <- c(ParmnTotal[, 4], PooledSETotal)
    ParamTestTotal <- ParmnTotal[, 4]
    Parmn_revTotal <- cbind(ParmnTotal[, 1:4], PooledSETotal)
    Parmn_revTotal[, 4:5] <- sapply(Parmn_revTotal[,4:5], as.numeric)
    colnames(Parmn_revTotal) <- c("lhs", "op", "rhs",
                                  "Estimate", "Pooled SE")
    df_tTotal <- (nConvergedProperTotal - 1) *
      (1 + (nConvergedProperTotal * PooledSEwithinvarTotal)/(nConvergedProperTotal *
                                                               PooledSEbetweenvarTotal + PooledSEbetweenvarTotal))^2
    crit_tTotal <- abs(qt(0.05/2, df_tTotal))
    pval_zTotal <- 2 * (1 - pnorm(abs(Parmn_revTotal[, 4]/PooledSETotal)))
    pval_tTotal <- 2 * (1 - pt(abs(Parmn_revTotal[, 4]/PooledSETotal),
                               df = df_tTotal))
    CI95_Lower_zTotal <- Parmn_revTotal[, 4] - 1.959963985 * PooledSETotal
    CI95_Upper_zTotal <- Parmn_revTotal[, 4] + 1.959963985 * PooledSETotal
    CI95_Lower_tTotal <- Parmn_revTotal[, 4] - crit_tTotal * PooledSETotal
    CI95_Upper_tTotal <- Parmn_revTotal[, 4] + crit_tTotal * PooledSETotal
    Parmn_revTotal <- cbind(Parmn_revTotal, pval_zTotal,
                            CI95_Lower_zTotal, CI95_Upper_zTotal, pval_tTotal,
                            CI95_Lower_tTotal, CI95_Upper_tTotal)
    colnames(Parmn_revTotal) <- c("lhs", "op", "rhs",
                                  "Pooled Est", "Pooled SE", "pval_z", "CI95_LB_z",
                                  "CI95_UB_z", "pval_t", "CI95_LB_t", "CI95_UB_t")
    for (i in 1:nrow(Parmn_revTotal)) {
      if (Parmn_revTotal[i, 5] == 0)
        Parmn_revTotal[i, 6:11] <- NA
    }
    RPAVTotal <- (PooledSEbetweenvarTotal + (PooledSEbetweenvarTotal/(nConvergedProperTotal)))/PooledSEwithinvarTotal
    PPAVTotal <- (((nConvergedProperTotal + 1)/(nConvergedProperTotal)) *
                    PooledSEbetweenvarTotal)/(PooledSEwithinvarTotal +
                                                (((nConvergedProperTotal + 1)/(nConvergedProperTotal)) * PooledSEbetweenvarTotal))
    PAVtableTotal <- cbind(ParmnTotal[1:3], RPAVTotal, PPAVTotal)
    Parmn_revTotal[, 4:11] <- apply(Parmn_revTotal[, 4:11], 2, round, digits = 4)
    FitsumTotal[, 1:5] <- apply(FitsumTotal[, 1:5], 2, round, digits = 4)
    PAVtableTotal[, 4:5] <- apply(PAVtableTotal[, 4:5], 2, round, digits = 4)
    FitsumTotal[2, 2:5] <- c("n/a", "n/a", "n/a", "n/a")
    ConvergedProperSumTotal <- rbind((nConvergedTotal)/(nAllocTotal),
                                     (nConvergedProperTotal)/(nAllocTotal))
    rownames(ConvergedProperSumTotal) <- c("Converged", "Converged and Proper")
    colnames(ConvergedProperSumTotal) <- "Proportion of Allocations"
  }
  if (nAllocAdd != 0) {
    if (nloop == 2) {
      PropParamMet <- matrix(data = 1, nrow(Parmn), 1)
      PropPooledSEMet <- matrix(data = 1, nrow(Parmn), 1)
    }
    if (nloop != 2) {
      PropParamMet <- (1 - apply(Paramdifftesttable[, 1:nloop - 1], 1, mean)) * 100
      PropPooledSEMet <- (1 - apply(PooledSEdifftesttable[,1:nloop - 1], 1, mean)) * 100
    }
    FirstParamMet <- apply(Paramdifftesttable == 0, 1, which.max)
    FirstPooledSEMet <- apply(PooledSEdifftesttable == 0, 1, which.max)
  }
  if (nAllocAdd == 0) {
    PropParamMet <- matrix(data = NA, nrow(Parmn), 1)
    PropPooledSEMet <- matrix(data = NA, nrow(Parmn), 1)
    FirstParamMet <- matrix(data = NA, nrow(Parmn), 1)
    FirstPooledSEMet <- matrix(data = NA, nrow(Parmn), 1)
  }
  PerLoops <- cbind(Parmn[, 1:3], PropParamMet, PropPooledSEMet)
  colnames(PerLoops) <- c("lhs", "op", "rhs", "Param Criteria Met",
                          "PooledSE Criteria Met")
  FirstLoops <- cbind(Parmn[, 1:3], FirstParamMet, FirstPooledSEMet)
  colnames(FirstLoops) <- c("lhs", "op", "rhs", "Param Criteria Met",
                            "PooledSE Criteria Met")
  NumbAllocations <- cbind(Parmn[, 1:3],
                           (FirstParamMet - 1) * nAllocAdd + nAllocStarttemp,
                           (FirstPooledSEMet - 1) * nAllocAdd + nAllocStarttemp)
  colnames(NumbAllocations) <- c("lhs", "op", "rhs", "Param Criteria Met",
                                 "PooledSE Criteria Met")
  if (nAllocAdd != 0) {
    for (i in 1:nrow(Parmn)) {
      if ((i %in% selectParam) == FALSE)
        PerLoops[i, 4:5] <- NA
      if ((i %in% selectParam) == FALSE)
        FirstLoops[i, 4:5] <- NA
      if ((i %in% selectParam) == FALSE)
        NumbAllocations[i, 4:5] <- NA
    }
  }
  df_t <- (nConvergedProperOutput - 1) *
    (1 + (nConvergedProperOutput * PooledSEwithinvarFinal) /
       (nConvergedProperOutput * PooledSEbetweenvarFinal + PooledSEbetweenvarFinal))^2
  crit_t <- abs(qt(0.05/2, df_t))
  pval_z <- 2 * (1 - pnorm(abs(Parmn_revFinal[, 4]/PooledSEFinal)))
  pval_t <- 2 * (1 - pt(abs(Parmn_revFinal[, 4]/PooledSEFinal),
                        df = df_t))
  CI95_Lower_z <- Parmn_revFinal[, 4] - 1.959963985 * PooledSEFinal
  CI95_Upper_z <- Parmn_revFinal[, 4] + 1.959963985 * PooledSEFinal
  CI95_Lower_t <- Parmn_revFinal[, 4] - crit_t * PooledSEFinal
  CI95_Upper_t <- Parmn_revFinal[, 4] + crit_t * PooledSEFinal
  Parmn_revFinal <- cbind(Parmn_revFinal, pval_z, CI95_Lower_z,
                          CI95_Upper_z, pval_t, CI95_Lower_t, CI95_Upper_t)
  colnames(Parmn_revFinal) <- c("lhs", "op", "rhs", "Pooled Est",
                                "Pooled SE", "pval_z", "CI95_LB_z", "CI95_UB_z",
                                "pval_t", "CI95_LB_t", "CI95_UB_t")
  for (i in 1:nrow(Parmn_revFinal)) {
    if (Parmn_revFinal[i, 5] == 0 | is.na(Parmn_revFinal[i, 5]) == TRUE)
      Parmn_revFinal[i, 6:11] <- NA
  }
  RPAV <- (PooledSEbetweenvarFinal + (PooledSEbetweenvarFinal/(nConvergedProperOutput)))/PooledSEwithinvarFinal
  PPAV <- (((nConvergedProperOutput + 1)/(nConvergedProperOutput)) *
             PooledSEbetweenvarFinal)/(PooledSEwithinvarFinal +
                                         (((nConvergedProperOutput + 1)/(nConvergedProperOutput)) *
                                            PooledSEbetweenvarFinal))
  PAVtable <- cbind(Parmn[1:3], RPAV, PPAV)
  colnames(Param_revTotal) <- c("lhs", "op", "rhs", c(1:nloop))
  colnames(PooledSE_revTotal) <- c("lhs", "op", "rhs",
                                   c(1:nloop))
  Param_revTotal[, 4:(nloop + 3)] <- sapply(Param_revTotal[,
                                                           4:(nloop + 3)], as.numeric)
  PooledSE_revTotal[, 4:(nloop + 3)] <- sapply(PooledSE_revTotal[,
                                                                 4:(nloop + 3)], as.numeric)
  Parmn_revFinal[, 4:11] <- apply(Parmn_revFinal[, 4:11],
                                  2, round, digits = 4)
  FitsumOutput[, 1:5] <- apply(FitsumOutput[, 1:5], 2,
                               round, digits = 4)
  if (nAllocAdd != 0)
    Param_revTotal[, 4:(nloop + 3)] <- apply(Param_revTotal[,
                                                            4:(nloop + 3)], 2, round, digits = 8)
  if (nAllocAdd == 0)
    Param_revTotal[, 4] <- round(Param_revTotal[, 4],
                                 8)
  if (nAllocAdd != 0)
    PooledSE_revTotal[, 4:(nloop + 3)] <- apply(PooledSE_revTotal[,
                                                                  4:(nloop + 3)], 2, round, digits = 8)
  if (nAllocAdd == 0)
    PooledSE_revTotal[, 4] <- round(PooledSE_revTotal[,
                                                      4], 8)
  PAVtable[, 4:5] <- apply(PAVtable[, 4:5], 2, round, digits = 4)
  FitsumOutput[2, 2:5] <- c("n/a", "n/a", "n/a", "n/a")
  ConvergedProperSum <- rbind((nConvergedOutput)/(nAllocOutput),
                              (nConvergedProperOutput)/(nAllocOutput))
  rownames(ConvergedProperSum) <- c("Converged", "Converged and Proper")
  colnames(ConvergedProperSum) <- "Proportion of Allocations"
  StopTimeFull <- proc.time() - StartTimeFull
  if (useTotalAlloc == FALSE) {
    Output_mod <- list(Parmn_revFinal, FitsumOutput,
                       ConvergedProperSum, nAllocOutput, PAVtable, StopTimeFull[[3]]/60)
    names(Output_mod) <- c("Estimates", "Fit",
                           "Proportion of Converged and Proper Allocations",
                           "Allocations needed for stability (M)",
                           "Indices to quantify uncertainty in estimates due to sampling vs. allocation variability",
                           "Total runtime (minutes)")
  }
  if (useTotalAlloc == TRUE) {
    Output_mod <- list(Parmn_revFinal, FitsumOutput,
                       ConvergedProperSum, nAllocOutput, PAVtable, Parmn_revTotal,
                       FitsumTotal, ConvergedProperSumTotal, nAllocTotal,
                       PAVtableTotal, StopTimeFull[[3]]/60)
    names(Output_mod) <- c("Estimates (using M allocations)", "Fit (using M allocations)",
                           "Proportion of Converged and Proper Allocations (using M allocations)",
                           "Allocations needed for stability (M)",
                           "Indices to quantify uncertainty in estimates due to sampling vs. allocation variability (using M allocations)",
                           "Estimates (using all allocations)", "Fit (using all allocations)",
                           "Proportion of Converged and Proper Allocations (using all allocations)",
                           "Total Allocations used by algorithm",
                           "Indices to quantify uncertainty in estimates due to sampling vs. allocation variability (using all allocations)",
                           "Total runtime (minutes)")
  }

  if (exists("invalidIndices")) {
    if (length(invalidIndices)) message('\n\nInvalid fit indices requested: ',
                                        paste(invalidIndices, collapse = ", "),
                                        "\n\n")
  }

  return(Output_mod)
}

