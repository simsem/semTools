### Jason D. Rights
### Last updated: 9 March 2018


#' Parcel-Allocation Variability in Model Ranking
#'
#' This function quantifies and assesses the consequences of parcel-allocation
#' variability for model ranking of structural equation models (SEMs) that
#' differ in their structural specification but share the same parcel-level
#' measurement specification (see Sterba & Rights, 2016). This function is a
#' modified version of \code{\link{parcelAllocation}} which can be used with
#' only one SEM in isolation. The \code{PAVranking} function repeatedly
#' generates a specified number of random item-to-parcel allocations, and then
#' fits two models to each allocation. Output includes summary information
#' about the distribution of model selection results (including plots) and the
#' distribution of results for each model individually, across allocations
#' within-sample. Note that this function can be used when selecting among more
#' than two competing structural models as well (see instructions below
#' involving \code{seed}).
#'
#' This is a modified version of \code{\link{parcelAllocation}} which was, in
#' turn, based on the SAS macro \code{ParcelAlloc} (Sterba & MacCallum, 2010).
#' The \code{PAVranking} function produces results discussed in Sterba and
#' Rights (2016) relevant to the assessment of parcel-allocation variability in
#' model selection and model ranking. Specifically, the \code{PAVranking}
#' function first uses a modified version of parcelAllocation to generate a
#' given number (\code{nAlloc}) of item-to-parcel allocations. Then,
#' \code{PAVranking} provides the following new developments: specifying more
#' than one SEM and producing results for Model A and Model B separately that
#' summarize parcel allocation variability in estimates, standard errors, and
#' fit indices. \code{PAVranking} also newly produces results summarizing
#' parcel allocation variability in model selection index values and model
#' ranking between Models A and B. Additionally, \code{PAVranking} newly allows
#' for nonconverged solutions and outputs the proportion of allocations that
#' converged as well as the proportion of proper solutions (results are
#' summarized for converged and proper allocations only).
#'
#' For further details on the benefits of the random allocation of items to
#' parcels, see Sterba (2011) and Sterba and MacCallum (2010).
#'
#' \emph{Note}: This function requires the \code{lavaan} package. Missing data
#'  codeneeds to be \code{NA}. If function returns \code{"Error in plot.new() :
#' figure margins too large,"} user may need to increase size of the plot
#' window and rerun.
#'
#'
#' @importFrom stats sd runif pchisq
#' @importFrom lavaan lavInspect
#'
#' @param nPerPar A list in which each element is a vector, corresponding to
#' each factor, indicating sizes of parcels. If variables are left out of
#' parceling, they should not be accounted for here (i.e., there should not be
#' parcels of size "1").
#' @param facPlc A list of vectors, each corresponding to a factor, specifying
#' the item indicators of that factor (whether included in parceling or not).
#' Either variable names or column numbers. Variables not listed will not be
#' modeled or included in output datasets.
#' @param nAlloc The number of random allocations of items to parcels to
#' generate.
#' @param syntaxA lavaan syntax for Model A. Note that, for likelihood ratio
#' test (LRT) results to be interpreted, Model A should be nested within Model
#' B (though the function will still provide results when Models A and B are
#' nonnested).
#' @param syntaxB lavaan syntax for Model B. Note that, for likelihood ratio
#' test (LRT) results to be appropriate, Model A should be nested within Model
#' B (though the function will still provide results when Models A and B are
#' nonnested).
#' @param dataset Item-level dataset
#' @param parceloutput folder where parceled data sets will be outputted (note
#' for Windows users: file path must specified using forward slashes).
#' @param seed (Optional) Random seed used for parceling items. When the same
#' random seed is specified and the program is re-run, the same allocations
#' will be generated. The seed argument can be used to assess parcel-allocation
#' variability in model ranking when considering more than two models. For each
#' pair of models under comparison, the program should be rerun using the same
#' random seed. Doing so ensures that multiple model comparisons will employ
#' the same set of parcel datasets.
#' @param names (Optional) A character vector containing the names of parceled
#' variables.
#' @param leaveout (Optional) A vector of variables to be left out of
#' randomized parceling. Either variable names or column numbers are allowed.
#' @param \dots Additional arguments to be passed to
#' \code{\link[lavaan]{lavaan}}. See also \code{\link[lavaan]{lavOptions}}
#' @return
#' \item{Estimates_A, Estimates_B}{A table containing results related
#' to parameter estimates (in table Estimates_A for Model A and in table
#' Estimates_B for Model B) with columns corresponding to parameter name,
#' average parameter estimate across allocations, standard deviation of
#' parameter estimate across allocations, the maximum parameter estimate across
#' allocations, the minimum parameter estimate across allocations, the range of
#' parameter estimates across allocations, and the percent of allocations in
#' which the parameter estimate is significant.}
#' \item{SE_A, SE_B}{A table containing results related to standard errors (in
#' table SE_A for Model A and in table SE_B for Model B) with columns
#' corresponding to parameter name, average standard error across allocations,
#' the standard deviation of standard errors across allocations, the maximum
#' standard error across allocations, the minimum standard error across
#' allocations, and the range of standard errors across allocations.}
#' \item{Fit_A, Fit_B}{A table containing results related to model fit (in
#' table Fit_A for Model A and in table Fit_B for Model B) with columns
#' corresponding to fit index name, the average of the fit index across
#' allocations, the standard deviation of the fit index across allocations,
#' the maximum of the fit index across allocations, the minimum of the fit
#' index across allocations, the range of the fit index across allocations, and
#' the percent of allocations where the chi-square test of absolute fit was
#' significant.}
#' \item{LRT Summary, Model A vs. Model B}{A table with columns corresponding
#' to: average likelihood ratio test (LRT) statistic for comparing Model A vs.
#' Model B (null hypothesis is no difference in fit between Models A and B in
#' the population), degrees of freedom (i.e. difference in the number of free
#' parameters between Models A and B), as well as the standard deviation,
#' maximum, and minimum of LRT statistics across allocations, and the percent of
#' allocations where the LRT was significant (indicating preference for the more
#' complex Model B). }
#' \item{LRT Summary, Model A vs. Model B}{A table with columns corresponding
#' to: average likelihood ratio test (LRT) statistic for comparing Model A vs.
#' Model B (null hypothesis is no difference in fit between Models A and B in
#' the population), degrees of freedom (i.e. difference in the number of free
#' parameters between Models A and B), as well as the standard deviation,
#' maximum, and minimum of LRT statistics across allocations, and the percent
#' of allocations where the LRT was significant (indicating preference for the
#' more complex Model B). }
#' \item{Fit index differences}{A table containing percentage of allocations
#' where Model A is preferred over Model B according to BIC, AIC, RMSEA, CFI,
#' TLI and SRMR and where Model B is preferred over Model A according to the
#' same indices. Also includes the average amount by which the given model is
#' preferred (calculated only using allocations where it was preferred).}
#' \item{Fit index difference histograms}{Histograms are automatically outputted
#' showing the distribution of the differences (Model A - Model B) for each fit
#' index and for the p-value of the likelihood ratio difference test.}
#' \item{Percent of Allocations with | BIC Diff | > 10}{A table containing the
#' percentage of allocations with (BIC for Model A) - (BIC for Model B) < -10,
#' indicating "very strong evidence" to prefer Model A over Model B and the
#' percentage of allocations with (BIC for Model A) - (BIC for Model B) > 10,
#' indicating "very strong evidence" to prefer Model B over Model A (Raftery,
#' 1995).}
#' \item{Converged and proper}{A table containing the proportion of allocations
#' that converged for Model A, Model B, and both models, and the proportion of
#' allocations with converged and proper solutions for Model A, Model B, and
#' both models.}
#' @author
#' Jason D. Rights (Vanderbilt University; \email{jason.d.rights@@vanderbilt.edu})
#'
#' The author would also like to credit Corbin Quick and Alexander Schoemann
#' for providing the original parcelAllocation function on which this function
#' is based.
#' @seealso \code{\link{parcelAllocation}}, \code{\link{poolMAlloc}}
#' @references
#' Raftery, A. E. (1995). Bayesian model selection in social
#' research. \emph{Sociological Methodology, 25}, 111--163. doi:10.2307/271063
#'
#' Sterba, S. K. (2011). Implications of parcel-allocation variability for
#' comparing fit of item-solutions and parcel-solutions. \emph{Structural
#' Equation Modeling: A Multidisciplinary Journal, 18}(4), 554--577.
#' doi:10.1080/10705511.2011.607073
#'
#' Sterba, S. K., & MacCallum, R. C. (2010). Variability in parameter estimates
#' and model fit across repeated allocations of items to parcels.
#' \emph{Multivariate Behavioral Research, 45}(2), 322--358.
#' doi:10.1080/00273171003680302
#'
#' Sterba, S. K., & Rights, J. D. (2017). Effects of parceling on model
#' selection: Parcel-allocation variability in model ranking.
#' \emph{Psychological Methods, 22}(1), 47--68. doi:10.1037/met0000067
#' @examples
#'
#' \dontrun{
#' ## lavaan syntax for Model A: a 2 Uncorrelated
#' ## factor CFA model to be fit to parceled data
#'
#' parmodelA <- '
#'    f1 =~ NA*p1f1 + p2f1 + p3f1
#'    f2 =~ NA*p1f2 + p2f2 + p3f2
#'    p1f1 ~ 1
#'    p2f1 ~ 1
#'    p3f1 ~ 1
#'    p1f2 ~ 1
#'    p2f2 ~ 1
#'    p3f2 ~ 1
#'    p1f1 ~~ p1f1
#'    p2f1 ~~ p2f1
#'    p3f1 ~~ p3f1
#'    p1f2 ~~ p1f2
#'    p2f2 ~~ p2f2
#'    p3f2 ~~ p3f2
#'    f1 ~~ 1*f1
#'    f2 ~~ 1*f2
#'    f1 ~~ 0*f2
#' '
#'
#' ## lavaan syntax for Model B: a 2 Correlated
#' ## factor CFA model to be fit to parceled data
#'
#' parmodelB <- '
#'    f1 =~ NA*p1f1 + p2f1 + p3f1
#'    f2 =~ NA*p1f2 + p2f2 + p3f2
#'    p1f1 ~ 1
#'    p2f1 ~ 1
#'    p3f1 ~ 1
#'    p1f2 ~ 1
#'    p2f2 ~ 1
#'    p3f2 ~ 1
#'    p1f1 ~~ p1f1
#'    p2f1 ~~ p2f1
#'    p3f1 ~~ p3f1
#'    p1f2 ~~ p1f2
#'    p2f2 ~~ p2f2
#'    p3f2 ~~ p3f2
#'    f1 ~~ 1*f1
#'    f2 ~~ 1*f2
#'    f1 ~~ f2
#' '
#'
#' ## specify items for each factor
#' f1name <- colnames(simParcel)[1:9]
#' f2name <- colnames(simParcel)[10:18]
#'
#' ## run function
#' PAVranking(nPerPar = list(c(3,3,3), c(3,3,3)), facPlc = list(f1name,f2name),
#'            nAlloc = 100, parceloutput = 0, leaveout = 0,
#'            syntaxA = parmodelA, syntaxB = parmodelB, dataset = simParcel,
#'            names = list("p1f1","p2f1","p3f1","p1f2","p2f2","p3f2"))
#' }
#'
#' @export
PAVranking <- function(nPerPar, facPlc, nAlloc = 100, parceloutput = 0, syntaxA, syntaxB,
                       dataset, names = NULL, leaveout = 0, seed = NA, ...) {
  if (is.null(names))
    names <- matrix(NA, length(nPerPar), 1)
  ## set random seed if specified
  if (is.na(seed) == FALSE)
    set.seed(seed)
  ## allow many tables to be outputted
  options(max.print = 1e+06)

  ## Create parceled datasets
  if (is.character(dataset)) dataset <- utils::read.csv(dataset)
  dataset <- as.matrix(dataset)

  if (nAlloc < 2)
    stop("Minimum of two allocations required.")

  if (is.list(facPlc)) {
    if (is.numeric(facPlc[[1]][1]) == FALSE) {
      facPlcb <- facPlc
      Namesv <- colnames(dataset)

      for (i in 1:length(facPlc)) {
        for (j in 1:length(facPlc[[i]])) {
          facPlcb[[i]][j] <- match(facPlc[[i]][j], Namesv)
        }
        facPlcb[[i]] <- as.numeric(facPlcb[[i]])
      }
      facPlc <- facPlcb
    }

    # facPlc2 <- rep(0, sum(sapply(facPlc, length)))
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
    ## this allows for unused variables in dataset, which are specified by zeros, and
    ## deleted
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

  if (length(Locate) != length(Npp))
    stop("Parcels incorrectly specified.\",
         \" Check input!")

  if (Maxv > 0) {
    ## Bug was here. With 1 factor Maxv=0. Skip this with a single factor
    for (i in 1:Maxv) {
      Mat <- match(i + 1, Locate)
      if (Npp[Mat] == Npp[Mat - 1])
        stop("Parcels incorrectly specified.\",
             \" Check input!")
    }
  }
  ## warning message if parcel crosses into multiple factors vector, parcel to which
  ## each variable belongs vector, factor to which each variables belongs if
  ## variables are in the same parcel, but different factors error message given in
  ## output

  Onevec <- facPlc - round(facPlc)
  NleaveA <- length(Onevec) - sum(Onevec == 0)
  NleaveP <- sum(nPerPar == 1)

  if (NleaveA < NleaveP)
    warning("Single-variable parcels have been requested.\",
            \" Check input!")

  if (NleaveA > NleaveP)
    warning("More non-parceled variables have been", " requested than provided for in parcel",
            " vector. Check input!")

  if (length(names) > 1) {
    if (length(names) != length(nPerPar))
      warning("Number of parcel names provided not equal to number", " of parcels requested")
  }

  Data <- c(1:ncol(dataset))
  ## creates a vector of the number of indicators e.g. for three indicators, c(1, 2,
  ## 3)
  Nfactors <- max(facPlc)
  ## scalar, number of factors
  Nindicators <- length(Data)
  ## scalar, number of indicators
  Npar <- length(nPerPar)
  ## scalar, number of parcels
  Rmize <- runif(Nindicators, 1, Nindicators)
  ## create vector of randomly ordered numbers, length of number of indicators

  Data <- rbind(facPlc, Rmize, Data)
  ## 'Data' becomes object of three rows, consisting of 1) factor to which each
  ## indicator belongs (in order to preserve indicator/factor assignment during
  ## randomization) 2) randomly order numbers 3) indicator number

  Results <- matrix(numeric(0), nAlloc, Nindicators)
  ## create empty matrix for parcel allocation matrix

  Pin <- nPerPar[1]
  for (i in 2:length(nPerPar)) {
    Pin <- c(Pin, nPerPar[i] + Pin[i - 1])
    ## creates vector which indicates the range of columns (endpoints) in each parcel
  }

  for (i in 1:nAlloc) {
    Data[2, ] <- runif(Nindicators, 1, Nindicators)
    ## Replace second row with newly randomly ordered numbers

    Data <- Data[, order(Data[2, ])]
    ## Order the columns according to the values of the second row

    Data <- Data[, order(Data[1, ])]
    ## Order the columns according to the values of the first row in order to preserve
    ## factor assignment

    Results[i, ] <- Data[3, ]
    ## assign result to allocation matrix
  }

  Alpha <- rbind(Results[1, ], dataset)
  ## bind first random allocation to dataset 'Alpha'

  Allocations <- list()
  ## create empty list for allocation data matrices

  for (i in 1:nAlloc) {
    Ineff <- rep(NA, ncol(Results))
    Ineff2 <- c(1:ncol(Results))
    for (inefficient in 1:ncol(Results)) {
      Ineff[Results[i, inefficient]] <- Ineff2[inefficient]
    }

    Alpha[1, ] <- Ineff
    ## replace first row of dataset matrix with row 'i' from allocation matrix

    Beta <- Alpha[, order(Alpha[1, ])]
    ## arrangle dataset columns by values of first row assign to temporary matrix
    ## 'Beta'

    Temp <- matrix(NA, nrow(dataset), Npar)
    ## create empty matrix for averaged parcel variables

    TempAA <- if (length(1:Pin[1]) > 1)
      Beta[2:nrow(Beta), 1:Pin[1]] else cbind(Beta[2:nrow(Beta), 1:Pin[1]], Beta[2:nrow(Beta), 1:Pin[1]])
    Temp[, 1] <- rowMeans(TempAA, na.rm = TRUE)
    ## fill first column with averages from assigned indicators
    for (al in 2:Npar) {
      Plc <- Pin[al - 1] + 1
      ## placeholder variable for determining parcel width
      TempBB <- if (length(Plc:Pin[al]) > 1)
        Beta[2:nrow(Beta), Plc:Pin[al]] else cbind(Beta[2:nrow(Beta), Plc:Pin[al]], Beta[2:nrow(Beta), Plc:Pin[al]])
      Temp[, al] <- rowMeans(TempBB, na.rm = TRUE)
      ## fill remaining columns with averages from assigned indicators
    }
    if (length(names) > 1)
      colnames(Temp) <- names
    Allocations[[i]] <- Temp
    ## assign result to list of parcel datasets
  }

  ## Write parceled datasets
  if (as.vector(regexpr("/", parceloutput)) != -1) {
    replist <- matrix(NA, nAlloc, 1)
    for (i in 1:nAlloc) {
      ## if (is.na(names)==TRUE) names <- matrix(NA,nrow(
      colnames(Allocations[[i]]) <- names
      utils::write.table(Allocations[[i]], paste(parceloutput, "/parcelruns", i,
                                                 ".dat", sep = ""),
                         row.names = FALSE, col.names = TRUE)
      replist[i, 1] <- paste("parcelruns", i, ".dat", sep = "")
    }
    utils::write.table(replist, paste(parceloutput, "/parcelrunsreplist.dat",
                                      sep = ""),
                       quote = FALSE, row.names = FALSE, col.names = FALSE)
  }


  ## Model A estimation

  {
    Param_A <- list()
    ## list for parameter estimated for each imputation
    Fitind_A <- list()
    ## list for fit indices estimated for each imputation
    Converged_A <- list()
    ## list for whether or not each allocation converged
    ProperSolution_A <- list()
    ## list for whether or not each allocation has proper solutions
    ConvergedProper_A <- list()
    ## list for whether or not each allocation converged and has proper solutions

    for (i in 1:nAlloc) {
      data_A <- as.data.frame(Allocations[[i]], row.names = NULL, optional = FALSE)
      ## convert allocation matrix to dataframe for model estimation
      fit_A <- lavaan::sem(syntaxA, data = data_A, ...)
      ## estimate model in lavaan
      if (lavInspect(fit_A, "converged") == TRUE) {
        Converged_A[[i]] <- 1
      } else Converged_A[[i]] <- 0
      ## determine whether or not each allocation converged
      Param_A[[i]] <- lavaan::parameterEstimates(fit_A)[, c("lhs", "op", "rhs",
                                                            "est", "se", "z", "pvalue", "ci.lower", "ci.upper")]
      ## assign allocation parameter estimates to list
      if (lavInspect(fit_A, "post.check") == TRUE & Converged_A[[i]] == 1) {
        ProperSolution_A[[i]] <- 1
      } else ProperSolution_A[[i]] <- 0
      ## determine whether or not each allocation has proper solutions
      if (any(is.na(Param_A[[i]][, 5] == TRUE)))
        ProperSolution_A[[i]] <- 0
      ## make sure each allocation has existing SE
      if (Converged_A[[i]] == 1 & ProperSolution_A[[i]] == 1) {
        ConvergedProper_A[[i]] <- 1
      } else ConvergedProper_A[[i]] <- 0
      ## determine whether or not each allocation converged and has proper solutions

      if (ConvergedProper_A[[i]] == 0)
        Param_A[[i]][, 4:9] <- matrix(data = NA, nrow(Param_A[[i]]), 6)
      ## make parameter estimates null for nonconverged, improper solutions

      if (ConvergedProper_A[[i]] == 1) {
        Fitind_A[[i]] <- lavaan::fitMeasures(fit_A, c("chisq", "df", "cfi",
                                                      "tli", "rmsea", "srmr", "logl", "bic", "aic"))
      } else Fitind_A[[i]] <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA)
      ### assign allocation parameter estimates to list

    }


    nConverged_A <- Reduce("+", Converged_A)
    ## count number of converged allocations

    nProperSolution_A <- Reduce("+", ProperSolution_A)
    ## count number of allocations with proper solutions

    nConvergedProper_A <- Reduce("+", ConvergedProper_A)
    ## count number of allocations with proper solutions

    if (nConvergedProper_A == 0)
      stop("All allocations failed to converge and/or yielded improper solutions for Model A and/or B.")
    ## stop program if no allocations converge

    Parmn_A <- Param_A[[1]]
    ## assign first parameter estimates to mean dataframe

    ParSE_A <- matrix(NA, nrow(Parmn_A), nAlloc)
    ParSEmn_A <- Parmn_A[, 5]

    Parsd_A <- matrix(NA, nrow(Parmn_A), nAlloc)
    ## assign parameter estimates for S.D. calculation

    Fitmn_A <- Fitind_A[[1]]
    ## assign first fit indices to mean dataframe

    Fitsd_A <- matrix(NA, length(Fitmn_A), nAlloc)
    ## assign fit indices for S.D. calculation

    Sigp_A <- matrix(NA, nrow(Parmn_A), nAlloc)
    ## assign p-values to calculate percentage significant

    Fitind_A <- data.frame(Fitind_A)
    ### convert fit index table to data frame

    for (i in 1:nAlloc) {

      Parsd_A[, i] <- Param_A[[i]][, 4]
      ## assign parameter estimates for S.D. estimation

      ParSE_A[, i] <- Param_A[[i]][, 5]

      if (i > 1) {
        ParSEmn_A <- rowSums(cbind(ParSEmn_A, Param_A[[i]][, 5]), na.rm = TRUE)
      }

      Sigp_A[, ncol(Sigp_A) - i + 1] <- Param_A[[i]][, 7]
      ## assign p-values to calculate percentage significant

      Fitsd_A[, i] <- Fitind_A[[i]]
      ## assign fit indices for S.D. estimation

      if (i > 1) {
        Parmn_A[, 4:ncol(Parmn_A)] <- rowSums(cbind(Parmn_A[, 4:ncol(Parmn_A)],
                                                    Param_A[[i]][, 4:ncol(Parmn_A)]), na.rm = TRUE)
      }
      ## add together all parameter estimates

      if (i > 1)
        Fitmn_A <- rowSums(cbind(Fitmn_A, Fitind_A[[i]]), na.rm = TRUE)
      ## add together all fit indices

    }


    Sigp_A <- Sigp_A + 0.45
    Sigp_A <- apply(Sigp_A, c(1, 2), round)
    Sigp_A <- 1 - as.vector(rowMeans(Sigp_A, na.rm = TRUE))
    ## calculate percentage significant parameters

    Parsum_A <- cbind(apply(Parsd_A, 1, mean, na.rm = TRUE), apply(Parsd_A, 1,
                                                                   sd, na.rm = TRUE), apply(Parsd_A, 1, max, na.rm = TRUE), apply(Parsd_A,
                                                                                                                                  1, min, na.rm = TRUE), apply(Parsd_A, 1, max, na.rm = TRUE) - apply(Parsd_A,
                                                                                                                                                                                                      1, min, na.rm = TRUE), Sigp_A * 100)
    colnames(Parsum_A) <- c("Avg Est.", "S.D.", "MAX", "MIN", "Range", "% Sig")
    ## calculate parameter S.D., minimum, maximum, range, bind to percentage
    ## significant

    ParSEmn_A <- Parmn_A[, 1:3]
    ParSEfn_A <- cbind(ParSEmn_A, apply(ParSE_A, 1, mean, na.rm = TRUE), apply(ParSE_A,
                                                                               1, sd, na.rm = TRUE), apply(ParSE_A, 1, max, na.rm = TRUE), apply(ParSE_A,
                                                                                                                                                 1, min, na.rm = TRUE), apply(ParSE_A, 1, max, na.rm = TRUE) - apply(ParSE_A,
                                                                                                                                                                                                                     1, min, na.rm = TRUE))
    colnames(ParSEfn_A) <- c("lhs", "op", "rhs", "Avg SE", "S.D.", "MAX", "MIN",
                             "Range")

    Fitsum_A <- cbind(apply(Fitsd_A, 1, mean, na.rm = TRUE), apply(Fitsd_A, 1,
                                                                   sd, na.rm = TRUE), apply(Fitsd_A, 1, max, na.rm = TRUE), apply(Fitsd_A,
                                                                                                                                  1, min, na.rm = TRUE), apply(Fitsd_A, 1, max, na.rm = TRUE) - apply(Fitsd_A,
                                                                                                                                                                                                      1, min, na.rm = TRUE))
    rownames(Fitsum_A) <- c("chisq", "df", "cfi", "tli", "rmsea", "srmr", "logl",
                            "bic", "aic")
    ## calculate fit S.D., minimum, maximum, range

    Parmn_A[, 4:ncol(Parmn_A)] <- Parmn_A[, 4:ncol(Parmn_A)]/nConvergedProper_A
    ## divide totalled parameter estimates by number converged allocations
    Parmn_A <- Parmn_A[, 1:3]
    ## remove confidence intervals from output
    Parmn_A <- cbind(Parmn_A, Parsum_A)
    ## bind parameter average estimates to cross-allocation information
    Fitmn_A <- Fitmn_A/nConvergedProper_A
    ## divide totalled fit indices by number converged allocations

    pChisq_A <- list()
    ## create empty list for Chi-square p-values
    sigChisq_A <- list()
    ## create empty list for Chi-square significance

    for (i in 1:nAlloc) {
      pChisq_A[[i]] <- (1 - pchisq(Fitsd_A[1, i], Fitsd_A[2, i]))
      ## calculate p-value for each Chi-square
      if (is.na(pChisq_A[[i]]) == FALSE & pChisq_A[[i]] < 0.05) {
        sigChisq_A[[i]] <- 1
      } else sigChisq_A[[i]] <- 0
    }
    ## count number of allocations with significant chi-square

    PerSigChisq_A <- (Reduce("+", sigChisq_A))/nConvergedProper_A * 100
    PerSigChisq_A <- round(PerSigChisq_A, 3)
    ## calculate percent of allocations with significant chi-square

    PerSigChisqCol_A <- c(PerSigChisq_A, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a",
                          "n/a", "n/a")
    ## create list of Chi-square Percent Significant and 'n/a' (used for fit summary
    ## table)

    options(stringsAsFactors = FALSE)
    ## set default option to allow strings into dataframe without converting to factors

    Fitsum_A <- data.frame(Fitsum_A, PerSigChisqCol_A)
    colnames(Fitsum_A) <- c("Avg Ind", "S.D.", "MAX", "MIN", "Range", "% Sig")
    ### bind to fit averages

    options(stringsAsFactors = TRUE)
    ## unset option to allow strings into dataframe without converting to factors

    ParSEfn_A[, 4:8] <- apply(ParSEfn_A[, 4:8], 2, round, digits = 3)
    Parmn_A[, 4:9] <- apply(Parmn_A[, 4:9], 2, round, digits = 3)
    Fitsum_A[, 1:5] <- apply(Fitsum_A[, 1:5], 2, round, digits = 3)
    ## round output to three digits

    Fitsum_A[2, 2:5] <- c("n/a", "n/a", "n/a", "n/a")
    ## Change df row to 'n/a' for sd, max, min, and range

    Output_A <- list(Parmn_A, ParSEfn_A, Fitsum_A)
    names(Output_A) <- c("Estimates_A", "SE_A", "Fit_A")
    ## output summary for model A

  }

  ## Model B estimation

  {
    Param <- list()
    ## list for parameter estimated for each imputation
    Fitind <- list()
    ## list for fit indices estimated for each imputation
    Converged <- list()
    ## list for whether or not each allocation converged
    ProperSolution <- list()
    ## list for whether or not each allocation has proper solutions
    ConvergedProper <- list()
    ## list for whether or not each allocation is converged and proper

    for (i in 1:nAlloc) {
      data <- as.data.frame(Allocations[[i]], row.names = NULL, optional = FALSE)
      ## convert allocation matrix to dataframe for model estimation
      fit <- lavaan::sem(syntaxB, data = data, ...)
      ## estimate model in lavaan
      if (lavInspect(fit, "converged") == TRUE) {
        Converged[[i]] <- 1
      } else Converged[[i]] <- 0
      ## determine whether or not each allocation converged
      Param[[i]] <- lavaan::parameterEstimates(fit)[, c("lhs", "op", "rhs",
                                                        "est", "se", "z", "pvalue", "ci.lower", "ci.upper")]
      ## assign allocation parameter estimates to list
      if (lavInspect(fit, "post.check") == TRUE & Converged[[i]] == 1) {
        ProperSolution[[i]] <- 1
      } else ProperSolution[[i]] <- 0
      ## determine whether or not each allocation has proper solutions
      if (any(is.na(Param[[i]][, 5] == TRUE)))
        ProperSolution[[i]] <- 0
      ## make sure each allocation has existing SE
      if (Converged[[i]] == 1 & ProperSolution[[i]] == 1) {
        ConvergedProper[[i]] <- 1
      } else ConvergedProper[[i]] <- 0
      ## determine whether or not each allocation converged and has proper solutions

      if (ConvergedProper[[i]] == 0)
        Param[[i]] <- matrix(data = NA, nrow(Param[[i]]), ncol(Param[[i]]))
      ## make parameter estimates null for nonconverged, improper solutions

      if (ConvergedProper[[i]] == 1) {
        Fitind[[i]] <- lavaan::fitMeasures(fit, c("chisq", "df", "cfi", "tli",
                                                  "rmsea", "srmr", "logl", "bic", "aic"))
      } else Fitind[[i]] <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA)
      ### assign allocation parameter estimates to list


    }




    nConverged <- Reduce("+", Converged)
    ## count number of converged allocations

    nProperSolution <- Reduce("+", ProperSolution)
    ## count number of allocations with proper solutions

    nConvergedProper <- Reduce("+", ConvergedProper)
    ## count number of allocations with proper solutions

    if (nConvergedProper == 0)
      stop("All allocations failed to converge", " and/or yielded improper solutions for",
           " Model A and/or B.")
    ## stop program if no allocations converge

    Parmn <- Param[[1]]
    ## assign first parameter estimates to mean dataframe

    ParSE <- matrix(NA, nrow(Parmn), nAlloc)
    ParSEmn <- Parmn[, 5]

    Parsd <- matrix(NA, nrow(Parmn), nAlloc)
    ## assign parameter estimates for S.D. calculation

    Fitmn <- Fitind[[1]]
    ## assign first fit indices to mean dataframe

    Fitsd <- matrix(NA, length(Fitmn), nAlloc)
    ## assign fit indices for S.D. calculation

    Sigp <- matrix(NA, nrow(Parmn), nAlloc)
    ## assign p-values to calculate percentage significant

    Fitind <- data.frame(Fitind)
    ### convert fit index table to dataframe


    for (i in 1:nAlloc) {

      Parsd[, i] <- Param[[i]][, 4]
      ## assign parameter estimates for S.D. estimation

      ParSE[, i] <- Param[[i]][, 5]

      if (i > 1)
        ParSEmn <- rowSums(cbind(ParSEmn, Param[[i]][, 5]), na.rm = TRUE)

      Sigp[, ncol(Sigp) - i + 1] <- Param[[i]][, 7]
      ## assign p-values to calculate percentage significant


      Fitsd[, i] <- Fitind[[i]]
      ## assign fit indices for S.D. estimation

      if (i > 1) {
        Parmn[, 4:ncol(Parmn)] <- rowSums(cbind(Parmn[, 4:ncol(Parmn)], Param[[i]][,
                                                                                   4:ncol(Parmn)]), na.rm = TRUE)
      }
      ## add together all parameter estimates

      if (i > 1)
        Fitmn <- rowSums(cbind(Fitmn, Fitind[[i]]), na.rm = TRUE)
      ## add together all fit indices

    }


    Sigp <- Sigp + 0.45
    Sigp <- apply(Sigp, c(1, 2), round)
    Sigp <- 1 - as.vector(rowMeans(Sigp, na.rm = TRUE))
    ## calculate percentage significant parameters

    Parsum <- cbind(apply(Parsd, 1, mean, na.rm = TRUE), apply(Parsd, 1, sd, na.rm = TRUE),
                    apply(Parsd, 1, max, na.rm = TRUE), apply(Parsd, 1, min, na.rm = TRUE),
                    apply(Parsd, 1, max, na.rm = TRUE) - apply(Parsd, 1, min, na.rm = TRUE),
                    Sigp * 100)
    colnames(Parsum) <- c("Avg Est", "S.D.", "MAX", "MIN", "Range", "% Sig")
    ## calculate parameter S.D., minimum, maximum, range, bind to percentage
    ## significant

    ParSEmn <- Parmn[, 1:3]
    ParSEfn <- cbind(ParSEmn, apply(ParSE, 1, mean, na.rm = TRUE), apply(ParSE,
                                                                         1, sd, na.rm = TRUE), apply(ParSE, 1, max, na.rm = TRUE), apply(ParSE,
                                                                                                                                         1, min, na.rm = TRUE), apply(ParSE, 1, max, na.rm = TRUE) - apply(ParSE,
                                                                                                                                                                                                           1, min, na.rm = TRUE))
    colnames(ParSEfn) <- c("lhs", "op", "rhs", "Avg SE", "S.D.", "MAX", "MIN",
                           "Range")

    Fitsum <- cbind(apply(Fitsd, 1, mean, na.rm = TRUE), apply(Fitsd, 1, sd, na.rm = TRUE),
                    apply(Fitsd, 1, max, na.rm = TRUE), apply(Fitsd, 1, min, na.rm = TRUE),
                    apply(Fitsd, 1, max, na.rm = TRUE) - apply(Fitsd, 1, min, na.rm = TRUE))
    rownames(Fitsum) <- c("chisq", "df", "cfi", "tli", "rmsea", "srmr", "logl",
                          "bic", "aic")
    ## calculate fit S.D., minimum, maximum, range

    Parmn[, 4:ncol(Parmn)] <- Parmn[, 4:ncol(Parmn)]/nConvergedProper
    ## divide totalled parameter estimates by number converged allocations
    Parmn <- Parmn[, 1:3]
    ## remove confidence intervals from output
    Parmn <- cbind(Parmn, Parsum)
    ## bind parameter average estimates to cross-allocation information
    Fitmn <- as.numeric(Fitmn)
    ## make fit index values numeric
    Fitmn <- Fitmn/nConvergedProper
    ## divide totalled fit indices by number converged allocations

    pChisq <- list()
    ## create empty list for Chi-square p-values
    sigChisq <- list()
    ## create empty list for Chi-square significance

    for (i in 1:nAlloc) {

      pChisq[[i]] <- (1 - pchisq(Fitsd[1, i], Fitsd[2, i]))
      ## calculate p-value for each Chi-square

      if (is.na(pChisq[[i]]) == FALSE & pChisq[[i]] < 0.05) {
        sigChisq[[i]] <- 1
      } else sigChisq[[i]] <- 0
    }
    ## count number of allocations with significant chi-square

    PerSigChisq <- (Reduce("+", sigChisq))/nConvergedProper * 100
    PerSigChisq <- round(PerSigChisq, 3)
    ## calculate percent of allocations with significant chi-square

    PerSigChisqCol <- c(PerSigChisq, "n/a", "n/a", "n/a", "n/a", "n/a", "n/a",
                        "n/a", "n/a")
    ## create list of Chi-square Percent Significant and 'n/a' (used for fit summary
    ## table)

    options(stringsAsFactors = FALSE)
    ## set default option to allow strings into dataframe without converting to factors

    Fitsum <- data.frame(Fitsum, PerSigChisqCol)
    colnames(Fitsum) <- c("Avg Ind", "S.D.", "MAX", "MIN", "Range", "% Sig")
    ### bind to fit averages

    options(stringsAsFactors = TRUE)
    ## unset option to allow strings into dataframe without converting to factors

    ParSEfn[, 4:8] <- apply(ParSEfn[, 4:8], 2, round, digits = 3)
    Parmn[, 4:9] <- apply(Parmn[, 4:9], 2, round, digits = 3)
    Fitsum[, 1:5] <- apply(Fitsum[, 1:5], 2, round, digits = 3)
    ## round output to three digits

    Fitsum[2, 2:5] <- c("n/a", "n/a", "n/a", "n/a")
    ## Change df row to 'n/a' for sd, max, min, and range

    Output_B <- list(Parmn, ParSEfn, Fitsum)
    names(Output_B) <- c("Estimates_B", "SE_B", "Fit_B")
    ## output summary for model A

  }

  ## Model Comparison (everything in this section is new)

  {
    Converged_AB <- list()
    ## create list of convergence comparison for each allocation
    ProperSolution_AB <- list()
    ## create list of proper solution comparison for each allocation
    ConvergedProper_AB <- list()
    ## create list of convergence and proper solution comparison for each allocation
    lrtest_AB <- list()
    ## create list for likelihood ratio test for each allocation
    lrchisq_AB <- list()
    ## create list for likelihood ratio chi square value
    lrchisqp_AB <- list()
    ## create list for likelihood ratio test p-value
    lrsig_AB <- list()
    ## create list for likelihood ratio test significance

    for (i in 1:nAlloc) {
      if (Converged_A[[i]] == 1 & Converged[[i]] == 1) {
        Converged_AB[[i]] <- 1
      } else Converged_AB[[i]] <- 0
      ## compare convergence

      if (ProperSolution_A[[i]] == 1 & ProperSolution[[i]] == 1) {
        ProperSolution_AB[[i]] <- 1
      } else ProperSolution_AB[[i]] <- 0
      ## compare existence of proper solutions

      if (ConvergedProper_A[[i]] == 1 & ConvergedProper[[i]] == 1) {
        ConvergedProper_AB[[i]] <- 1
      } else ConvergedProper_AB[[i]] <- 0
      ## compare existence of proper solutions and convergence



      if (ConvergedProper_AB[[i]] == 1) {

        data <- as.data.frame(Allocations[[i]], row.names = NULL, optional = FALSE)
        ## convert allocation matrix to dataframe for model estimation
        fit_A <- lavaan::sem(syntaxA, data = data, ...)
        ## estimate model A in lavaan
        fit <- lavaan::sem(syntaxB, data = data, ...)
        ## estimate model B in lavaan
        lrtest_AB[[i]] <- lavaan::lavTestLRT(fit_A, fit)
        ## likelihood ratio test comparing A and B
        lrtestd_AB <- as.data.frame(lrtest_AB[[i]], row.names = NULL, optional = FALSE)
        ## convert lrtest results to dataframe
        lrchisq_AB[[i]] <- lrtestd_AB[2, 5]
        ## write lrtest chisq as single numeric variable
        lrchisqp_AB[[i]] <- lrtestd_AB[2, 7]
        ## write lrtest p-value as single numeric variable
        if (lrchisqp_AB[[i]] < 0.05) {
          lrsig_AB[[i]] <- 1
        } else {
          lrsig_AB[[i]] <- 0
        }
        ## determine statistical significance of lrtest

      }
    }

    lrchisqp_AB <- unlist(lrchisqp_AB, recursive = TRUE, use.names = TRUE)
    ## convert lrchisqp_AB from list to vector
    lrchisqp_AB <- as.numeric(lrchisqp_AB)
    ## make lrchisqp_AB numeric
    lrsig_AB <- unlist(lrsig_AB, recursive = TRUE, use.names = TRUE)
    ## convert lrsig_AB from list to vector
    lrsig_AB <- as.numeric(lrsig_AB)
    ### make lrsig_AB numeric


    nConverged_AB <- Reduce("+", Converged_AB)
    ## count number of allocations that converged for both A and B
    nProperSolution_AB <- Reduce("+", ProperSolution_AB)
    ## count number of allocations with proper solutions for both A and B
    nConvergedProper_AB <- Reduce("+", ConvergedProper_AB)
    ## count number of allocations that converged and have proper solutions for both A
    ## and B
    ProConverged_AB <- (nConverged_AB/nAlloc) * 100
    ## calc proportion of allocations that converged for both A and B
    nlrsig_AB <- Reduce("+", lrsig_AB)
    ## count number of allocations with significant lrtest between A and B
    Prolrsig_AB <- (nlrsig_AB/nConvergedProper_AB) * 100
    ## calc proportion of allocations with significant lrtest between A and B
    lrchisq_AB <- unlist(lrchisq_AB, recursive = TRUE, use.names = TRUE)
    ### convert lrchisq_AB from list to vector
    lrchisq_AB <- as.numeric(lrchisq_AB)
    ### make lrchisq_AB numeric
    AvgLRT_AB <- (Reduce("+", lrchisq_AB))/nConvergedProper_AB
    ## calc average LRT

    LRTsum <- cbind(AvgLRT_AB, lrtestd_AB[2, 3], sd(lrchisq_AB, na.rm = TRUE),
                    max(lrchisq_AB), min(lrchisq_AB),
                    max(lrchisq_AB) - min(lrchisq_AB), Prolrsig_AB)
    colnames(LRTsum) <- c("Avg LRT", "df", "S.D.", "MAX", "MIN", "Range", "% Sig")
    ## calculate LRT distribution statistics

    FitDiff_AB <- Fitsd_A - Fitsd
    ## compute fit index difference matrix

    for (i in 1:nAlloc) {
      if (ConvergedProper_AB[[i]] != 1)
        FitDiff_AB[1:9, i] <- 0
    }
    ### make fit differences zero for each non-converged allocation

    BICDiff_AB <- list()
    AICDiff_AB <- list()
    RMSEADiff_AB <- list()
    CFIDiff_AB <- list()
    TLIDiff_AB <- list()
    SRMRDiff_AB <- list()
    BICDiffGT10_AB <- list()
    ## create list noting each allocation in which A is preferred over B

    BICDiff_BA <- list()
    AICDiff_BA <- list()
    RMSEADiff_BA <- list()
    CFIDiff_BA <- list()
    TLIDiff_BA <- list()
    SRMRDiff_BA <- list()
    BICDiffGT10_BA <- list()
    ## create list noting each allocation in which B is preferred over A

    for (i in 1:nAlloc) {
      if (FitDiff_AB[8, i] < 0) {
        BICDiff_AB[[i]] <- 1
      } else BICDiff_AB[[i]] <- 0
      if (FitDiff_AB[9, i] < 0) {
        AICDiff_AB[[i]] <- 1
      } else AICDiff_AB[[i]] <- 0
      if (FitDiff_AB[5, i] < 0) {
        RMSEADiff_AB[[i]] <- 1
      } else RMSEADiff_AB[[i]] <- 0
      if (FitDiff_AB[3, i] > 0) {
        CFIDiff_AB[[i]] <- 1
      } else CFIDiff_AB[[i]] <- 0
      if (FitDiff_AB[4, i] > 0) {
        TLIDiff_AB[[i]] <- 1
      } else TLIDiff_AB[[i]] <- 0
      if (FitDiff_AB[6, i] < 0) {
        SRMRDiff_AB[[i]] <- 1
      } else SRMRDiff_AB[[i]] <- 0
      if (FitDiff_AB[8, i] < (-10)) {
        BICDiffGT10_AB[[i]] <- 1
      } else BICDiffGT10_AB[[i]] <- 0
    }
    nBIC_AoverB <- Reduce("+", BICDiff_AB)
    nAIC_AoverB <- Reduce("+", AICDiff_AB)
    nRMSEA_AoverB <- Reduce("+", RMSEADiff_AB)
    nCFI_AoverB <- Reduce("+", CFIDiff_AB)
    nTLI_AoverB <- Reduce("+", TLIDiff_AB)
    nSRMR_AoverB <- Reduce("+", SRMRDiff_AB)
    nBICDiffGT10_AoverB <- Reduce("+", BICDiffGT10_AB)
    ## compute number of 'A preferred over B' for each fit index

    for (i in 1:nAlloc) {
      if (FitDiff_AB[8, i] > 0) {
        BICDiff_BA[[i]] <- 1
      } else BICDiff_BA[[i]] <- 0
      if (FitDiff_AB[9, i] > 0) {
        AICDiff_BA[[i]] <- 1
      } else AICDiff_BA[[i]] <- 0
      if (FitDiff_AB[5, i] > 0) {
        RMSEADiff_BA[[i]] <- 1
      } else RMSEADiff_BA[[i]] <- 0
      if (FitDiff_AB[3, i] < 0) {
        CFIDiff_BA[[i]] <- 1
      } else CFIDiff_BA[[i]] <- 0
      if (FitDiff_AB[4, i] < 0) {
        TLIDiff_BA[[i]] <- 1
      } else TLIDiff_BA[[i]] <- 0
      if (FitDiff_AB[6, i] > 0) {
        SRMRDiff_BA[[i]] <- 1
      } else SRMRDiff_BA[[i]] <- 0
      if (FitDiff_AB[8, i] > (10)) {
        BICDiffGT10_BA[[i]] <- 1
      } else BICDiffGT10_BA[[i]] <- 0
    }
    nBIC_BoverA <- Reduce("+", BICDiff_BA)
    nAIC_BoverA <- Reduce("+", AICDiff_BA)
    nRMSEA_BoverA <- Reduce("+", RMSEADiff_BA)
    nCFI_BoverA <- Reduce("+", CFIDiff_BA)
    nTLI_BoverA <- Reduce("+", TLIDiff_BA)
    nSRMR_BoverA <- Reduce("+", SRMRDiff_BA)
    nBICDiffGT10_BoverA <- Reduce("+", BICDiffGT10_BA)
    ## compute number of 'B preferred over A' for each fit index

    BICDiffAvgtemp <- list()
    AICDiffAvgtemp <- list()
    RMSEADiffAvgtemp <- list()
    CFIDiffAvgtemp <- list()
    TLIDiffAvgtemp <- list()
    SRMRDiffAvgtemp <- list()
    BICgt10DiffAvgtemp <- list()
    ## create empty list for average fit index differences

    for (i in 1:nAlloc) {
      if (BICDiff_AB[[i]] != 1) {
        BICDiffAvgtemp[[i]] <- 0
      } else BICDiffAvgtemp[[i]] <- FitDiff_AB[8, i]
      if (AICDiff_AB[[i]] != 1) {
        AICDiffAvgtemp[[i]] <- 0
      } else AICDiffAvgtemp[[i]] <- FitDiff_AB[9, i]
      if (RMSEADiff_AB[[i]] != 1) {
        RMSEADiffAvgtemp[[i]] <- 0
      } else RMSEADiffAvgtemp[[i]] <- FitDiff_AB[5, i]
      if (CFIDiff_AB[[i]] != 1) {
        CFIDiffAvgtemp[[i]] <- 0
      } else CFIDiffAvgtemp[[i]] <- FitDiff_AB[3, i]
      if (TLIDiff_AB[[i]] != 1) {
        TLIDiffAvgtemp[[i]] <- 0
      } else TLIDiffAvgtemp[[i]] <- FitDiff_AB[4, i]
      if (SRMRDiff_AB[[i]] != 1) {
        SRMRDiffAvgtemp[[i]] <- 0
      } else SRMRDiffAvgtemp[[i]] <- FitDiff_AB[6, i]
      if (BICDiffGT10_AB[[i]] != 1) {
        BICgt10DiffAvgtemp[[i]] <- 0
      } else BICgt10DiffAvgtemp[[i]] <- FitDiff_AB[8, i]
    }
    ## make average fit index difference list composed solely of values where A is
    ## preferred over B

    BICDiffAvg_AB <- Reduce("+", BICDiffAvgtemp)/nBIC_AoverB * (-1)
    AICDiffAvg_AB <- Reduce("+", AICDiffAvgtemp)/nAIC_AoverB * (-1)
    RMSEADiffAvg_AB <- Reduce("+", RMSEADiffAvgtemp)/nRMSEA_AoverB * (-1)
    CFIDiffAvg_AB <- Reduce("+", CFIDiffAvgtemp)/nCFI_AoverB
    TLIDiffAvg_AB <- Reduce("+", TLIDiffAvgtemp)/nTLI_AoverB
    SRMRDiffAvg_AB <- Reduce("+", SRMRDiffAvgtemp)/nSRMR_AoverB * (-1)
    BICgt10DiffAvg_AB <- Reduce("+", BICgt10DiffAvgtemp)/nBICDiffGT10_AoverB *
      (-1)
    ## calc average fit index difference when A is preferred over B

    FitDiffAvg_AoverB <- list(BICDiffAvg_AB, AICDiffAvg_AB, RMSEADiffAvg_AB, CFIDiffAvg_AB,
                              TLIDiffAvg_AB, SRMRDiffAvg_AB)
    ## create list of all fit index differences when A is preferred over B

    FitDiffAvg_AoverB <- unlist(FitDiffAvg_AoverB, recursive = TRUE, use.names = TRUE)
    ### convert from list to vector

    for (i in 1:nAlloc) {
      if (BICDiff_BA[[i]] != 1) {
        BICDiffAvgtemp[[i]] <- 0
      } else BICDiffAvgtemp[[i]] <- FitDiff_AB[8, i]
      if (AICDiff_BA[[i]] != 1) {
        AICDiffAvgtemp[[i]] <- 0
      } else AICDiffAvgtemp[[i]] <- FitDiff_AB[9, i]
      if (RMSEADiff_BA[[i]] != 1) {
        RMSEADiffAvgtemp[[i]] <- 0
      } else RMSEADiffAvgtemp[[i]] <- FitDiff_AB[5, i]
      if (CFIDiff_BA[[i]] != 1) {
        CFIDiffAvgtemp[[i]] <- 0
      } else CFIDiffAvgtemp[[i]] <- FitDiff_AB[3, i]
      if (TLIDiff_BA[[i]] != 1) {
        TLIDiffAvgtemp[[i]] <- 0
      } else TLIDiffAvgtemp[[i]] <- FitDiff_AB[4, i]
      if (SRMRDiff_BA[[i]] != 1) {
        SRMRDiffAvgtemp[[i]] <- 0
      } else SRMRDiffAvgtemp[[i]] <- FitDiff_AB[6, i]
      if (BICDiffGT10_BA[[i]] != 1) {
        BICgt10DiffAvgtemp[[i]] <- 0
      } else BICgt10DiffAvgtemp[[i]] <- FitDiff_AB[8, i]
    }
    ## make average fit index difference list composed solely of values where B is
    ## preferred over A

    BICDiffAvg_BA <- Reduce("+", BICDiffAvgtemp)/nBIC_BoverA
    AICDiffAvg_BA <- Reduce("+", AICDiffAvgtemp)/nAIC_BoverA
    RMSEADiffAvg_BA <- Reduce("+", RMSEADiffAvgtemp)/nRMSEA_BoverA
    CFIDiffAvg_BA <- Reduce("+", CFIDiffAvgtemp)/nCFI_BoverA * (-1)
    TLIDiffAvg_BA <- Reduce("+", TLIDiffAvgtemp)/nTLI_BoverA * (-1)
    SRMRDiffAvg_BA <- Reduce("+", SRMRDiffAvgtemp)/nSRMR_BoverA
    BICgt10DiffAvg_BA <- Reduce("+", BICgt10DiffAvgtemp)/nBICDiffGT10_BoverA
    ## calc average fit index difference when B is preferred over A

    FitDiffAvg_BoverA <- list(BICDiffAvg_BA, AICDiffAvg_BA, RMSEADiffAvg_BA, CFIDiffAvg_BA,
                              TLIDiffAvg_BA, SRMRDiffAvg_BA)
    ## create list of all fit index differences when B is preferred over A

    FitDiffAvg_BoverA <- unlist(FitDiffAvg_BoverA, recursive = TRUE, use.names = TRUE)
    ### convert from list to vector

    FitDiffBICgt10_AoverB <- nBICDiffGT10_AoverB/nConvergedProper_AB * 100
    ### calculate portion of allocations where A strongly preferred over B

    FitDiffBICgt10_BoverA <- nBICDiffGT10_BoverA/nConvergedProper_AB * 100
    ### calculate portion of allocations where B strongly preferred over A

    FitDiffBICgt10 <- rbind(FitDiffBICgt10_AoverB, FitDiffBICgt10_BoverA)
    rownames(FitDiffBICgt10) <- c("Very Strong evidence for A>B", "Very Strong evidence for B>A")
    colnames(FitDiffBICgt10) <- "% Allocations"
    ### create table of proportions of 'A strongly preferred over B' and 'B strongly
    ### preferred over A'

    FitDiff_AoverB <- list(nBIC_AoverB/nConvergedProper_AB * 100, nAIC_AoverB/nConvergedProper_AB *
                             100, nRMSEA_AoverB/nConvergedProper_AB * 100, nCFI_AoverB/nConvergedProper_AB *
                             100, nTLI_AoverB/nConvergedProper_AB * 100, nSRMR_AoverB/nConvergedProper_AB *
                             100)
    ### create list of all proportions of 'A preferred over B'
    FitDiff_BoverA <- list(nBIC_BoverA/nConvergedProper_AB * 100, nAIC_BoverA/nConvergedProper_AB *
                             100, nRMSEA_BoverA/nConvergedProper_AB * 100, nCFI_BoverA/nConvergedProper_AB *
                             100, nTLI_BoverA/nConvergedProper_AB * 100, nSRMR_BoverA/nConvergedProper_AB *
                             100)
    ### create list of all proportions of 'B preferred over A'

    FitDiff_AoverB <- unlist(FitDiff_AoverB, recursive = TRUE, use.names = TRUE)
    ### convert from list to vector

    FitDiff_BoverA <- unlist(FitDiff_BoverA, recursive = TRUE, use.names = TRUE)
    ### convert from list to vector

    FitDiffSum_AB <- cbind(FitDiff_AoverB, FitDiffAvg_AoverB, FitDiff_BoverA,
                           FitDiffAvg_BoverA)
    colnames(FitDiffSum_AB) <- c("% A>B", "Avg Amount A>B", "% B>A", "Avg Amount B>A")
    rownames(FitDiffSum_AB) <- c("bic", "aic", "rmsea", "cfi", "tli", "srmr")
    ## create table showing number of allocations in which A>B and B>A as well as
    ## average difference values

    for (i in 1:nAlloc) {
      is.na(FitDiff_AB[1:9, i]) <- ConvergedProper_AB[[i]] != 1
    }
    ### make fit differences missing for each non-converged allocation

    LRThistMax <- max(hist(lrchisqp_AB, plot = FALSE)$counts)
    BIChistMax <- max(hist(FitDiff_AB[8, 1:nAlloc], plot = FALSE)$counts)
    AIChistMax <- max(hist(FitDiff_AB[9, 1:nAlloc], plot = FALSE)$counts)
    RMSEAhistMax <- max(hist(FitDiff_AB[5, 1:nAlloc], plot = FALSE)$counts)
    CFIhistMax <- max(hist(FitDiff_AB[3, 1:nAlloc], plot = FALSE)$counts)
    TLIhistMax <- max(hist(FitDiff_AB[4, 1:nAlloc], plot = FALSE)$counts)
    ### calculate y-axis height for each histogram

    LRThist <- hist(lrchisqp_AB, ylim = c(0, LRThistMax), xlab = "p-value", main = "LRT p-values")
    ## plot histogram of LRT p-values

    BIChist <- hist(FitDiff_AB[8, 1:nAlloc], ylim = c(0, BIChistMax), xlab = "BIC_modA - BIC_modB",
                    main = "BIC Diff")
    AIChist <- hist(FitDiff_AB[9, 1:nAlloc], ylim = c(0, AIChistMax), xlab = "AIC_modA - AIC_modB",
                    main = "AIC Diff")
    RMSEAhist <- hist(FitDiff_AB[5, 1:nAlloc], ylim = c(0, RMSEAhistMax), xlab = "RMSEA_modA - RMSEA_modB",
                      main = "RMSEA Diff")
    CFIhist <- hist(FitDiff_AB[3, 1:nAlloc], ylim = c(0, CFIhistMax), xlab = "CFI_modA - CFI_modB",
                    main = "CFI Diff")
    TLIhist <- hist(FitDiff_AB[4, 1:nAlloc], ylim = c(0, TLIhistMax), xlab = "TLI_modA - TLI_modB",
                    main = "TLI Diff")
    ### plot histograms for each index_modA - index_modB
    BIChist
    AIChist
    RMSEAhist
    CFIhist
    TLIhist

    ConvergedProperSum <- rbind(nConverged_A/nAlloc, nConverged/nAlloc, nConverged_AB/nAlloc,
                                nConvergedProper_A/nAlloc, nConvergedProper/nAlloc, nConvergedProper_AB/nAlloc)
    rownames(ConvergedProperSum) <- c("Converged_A", "Converged_B", "Converged_AB",
                                      "ConvergedProper_A", "ConvergedProper_B", "ConvergedProper_AB")
    colnames(ConvergedProperSum) <- "Proportion of Allocations"
    ### create table summarizing proportions of converged allocations and allocations
    ### with proper solutions

    Output_AB <- list(round(LRTsum, 3), "LRT results are interpretable specifically for nested models",
                      round(FitDiffSum_AB, 3), round(FitDiffBICgt10, 3), ConvergedProperSum)
    names(Output_AB) <- c("LRT Summary, Model A vs. Model B", "Note:", "Fit Index Differences",
                          "Percent of Allocations with |BIC Diff| > 10", "Converged and Proper Solutions Summary")
    ### output for model comparison

  }

  return(list(Output_A, Output_B, Output_AB))
  ## returns output for model A, model B, and the comparison of these
}

