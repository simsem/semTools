### Terrence D. Jorgensen & Yves Rosseel
### Last updated: 12 March 2025
### Pooled score test (= Lagrange Multiplier test) for multiple imputations
### Borrowed source code from lavaan/R/lav_test_score.R

### DEPRECATED: 16 June 2024
### supplanted by lavaan.mi package

## this function can run two modes:
## MODE 1: 'add'
##   add new parameters that are currently not included in de model
##   (aka fixed to zero), but should be released
## MODE 2: 'release' (the default)
##   release existing "==" constraints



##' Score Test for Multiple Imputations
##'
##' Score test (or "Lagrange multiplier" test) for lavaan models fitted to
##' multiple imputed data sets. Statistics for releasing one or more
##' fixed or constrained parameters in model can be calculated by pooling
##' the gradient and information matrices pooled across imputed data sets in a
##' method proposed by Mansolf, Jorgensen, & Enders (2020)---analogous to
##' the "D1" Wald test proposed by Li, Meng, Raghunathan, & Rubin's (1991)---or
##' by pooling the complete-data score-test statistics across imputed data sets
##' (i.e., "D2"; Li et al., 1991).
##'
## @aliases lavTestScore.mi
##' @importFrom lavaan lavListInspect parTable
##' @importFrom stats cov pchisq pf
##' @importFrom methods getMethod
##'
##' @param object An object of class [OLDlavaan.mi-class].
##' @param add Either a `character` string (typically between single
##'   quotes) or a parameter table containing additional (currently
##'   fixed-to-zero) parameters for which the score test must be computed.
##' @param release Vector of `integer`s. The indices of the *equality*
##'  constraints that should be released. The indices correspond to the order of
##'  the equality constraints as they appear in the parameter table.
##' @param test `character` indicating which pooling method to use.
##'   `"D1"` requests Mansolf, Jorgensen, & Enders' (2020) proposed
##'   Wald-like test for pooling the gradient and information, which are then
##'   used to calculate score-test statistics in the usual manner. `"D2"`
##'   (default because it is less computationall intensive) requests to pool the
##'   complete-data score-test statistics from each imputed data set, then pool
##'   them across imputations, described by Li et al. (1991) and Enders (2010).
##' @param scale.W `logical`. If `FALSE`, the pooled
##'   information matrix is calculated as the weighted sum of the
##'   within-imputation and between-imputation components. Otherwise, the pooled
##'   information is calculated by scaling the within-imputation component by
##'   the average relative increase in variance (ARIV; Enders, 2010, p. 235),
##'   which is *only* consistent when requesting the *F* test (i.e.,
##'   `asymptotic = FALSE`.  Ignored (irrelevant) if `test = "D2"`.
##' @param omit.imps `character` vector specifying criteria for omitting
##'   imputations from pooled results.  Can include any of
##'   `c("no.conv", "no.se", "no.npd")`, the first 2 of which are the
##'   default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (`"no.npd"`) would exclude any imputations which
##'   yielded a nonpositive definite covariance matrix for observed or
##'   latent variables, which would include any "improper solutions" such
##'   as Heywood cases. Specific imputation numbers can also be included in this
##'   argument, in case users want to  apply their own custom omission criteria
##'   (or simulations can use different numbers of imputations without
##'   redundantly refitting the model).
##' @param asymptotic `logical`. If `FALSE` (default when using
##'   `add` to test adding fixed parameters to the model), the pooled test
##'   will be returned as an *F*-distributed variable with numerator
##'   (`df1`) and denominator (`df2`) degrees of freedom.
##'   If `TRUE`, the pooled *F* statistic will be multiplied by its
##'   `df1` on the assumption that its `df2` is sufficiently large
##'   enough that the statistic will be asymptotically \eqn{\chi^2} distributed
##'   with `df1`. When using the `release` argument, `asymptotic`
##'   will be set to `TRUE` because (A)RIV can only be calculated for
##'   `add`ed parameters.
##' @param univariate `logical`. If `TRUE`, compute the univariate
##'   score statistics, one for each constraint.
##' @param cumulative `logical`. If `TRUE`, order the univariate score
##'   statistics from large to small, and compute a series of multivariate
##'   score statistics, each time including an additional constraint in the test.
##' @param epc `logical`. If `TRUE`, and we are releasing existing
##'   constraints, compute the expected parameter changes for the existing
##'   (free) parameters (and any specified with `add`), if all constraints
##'   were released. For EPCs associated with a particular (1-*df*)
##'   constraint, only specify one parameter in `add` or one constraint in
##'   `release`.
##' @param standardized If `TRUE`, two extra columns (`sepc.lv` and
##'   `sepc.all`) in the `$epc` table will contain standardized values
##'   for the EPCs. See [lavaan::lavTestScore()].
##' @param cov.std `logical`. See [lavaan::standardizedSolution()].
##' @param verbose `logical`. Not used for now.
##' @param warn `logical`. If `TRUE`, print warnings if they occur.
##' @param information `character` indicating the type of information
##'   matrix to use (check [lavaan::lavInspect()] for available options).
##'   `"expected"` information is the default, which provides better
##'   control of Type I errors.
##'
##' @return
##'  A list containing at least one `data.frame`:
##'  \itemize{
##'    \item{`$test`: The total score test, with columns for the score
##'      test statistic (`X2`), its degrees of freedom (`df`), its
##'      *p* value under the \eqn{\chi^2} distribution (`p.value`),
##'      and if `asymptotic=FALSE`, the average relative invrease in
##'      variance (ARIV) used to calculate the denominator *df* is also
##'      returned as a missing-data diagnostic, along with the fraction missing
##'      information (FMI = ARIV / (1 + ARIV)).}
##'    \item{`$uni`: Optional (if `univariate=TRUE`).
##'      Each 1-*df* score test, equivalent to modification indices. Also
##'      includes EPCs if `epc=TRUE`, and RIV and FMI if
##'      `asymptotic=FALSE`.}
##'    \item{`$cumulative`: Optional (if `cumulative=TRUE`).
##'      Cumulative score tests, with ARIV and FMI if `asymptotic=FALSE`.}
##'    \item{`$epc`: Optional (if `epc=TRUE`). Parameter estimates,
##'      expected parameter changes, and expected parameter values if ALL
##'      the tested constraints were freed.}
##'  }
##' See [lavaan::lavTestScore()] for details.
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' Adapted from \pkg{lavaan} source code, written by
##'   Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
##'
##' `test = "D1"` method proposed by
##'   Maxwell Mansolf (University of California, Los Angeles;
##'   \email{mamansolf@@gmail.com})
##'
##' @references
##'   Bentler, P. M., & Chou, C.-P. (1992). Some new covariance structure model
##'   improvement statistics. *Sociological Methods & Research, 21*(2),
##'   259--282. \doi{10.1177/0049124192021002006}
##'
##'   Enders, C. K. (2010). *Applied missing data analysis*.
##'   New York, NY: Guilford.
##'
##'   Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
##'   Significance levels from repeated *p*-values with multiply-imputed
##'   data. *Statistica Sinica, 1*(1), 65--92. Retrieved from
##'   <https://www.jstor.org/stable/24303994>
##'
##'   Mansolf, M., Jorgensen, T. D., & Enders, C. K. (2020). A multiple
##'   imputation score test for model modification in structural equation
##'   models. *Psychological Methods, 25*(4), 393--411.
##'   \doi{10.1037/met0000243}
##'
##' @seealso [lavaan::lavTestScore()]
##'
##' @examples
##'
##' ## See the new lavaan.mi package
##'
##' @name lavTestScore.mi-deprecated
##' @usage
##' lavTestScore.mi(object, add = NULL, release = NULL,
##'                 test = c("D2","D1"), scale.W = !asymptotic,
##'                 omit.imps = c("no.conv","no.se"),
##'                 asymptotic = is.null(add),
##'                 univariate = TRUE, cumulative = FALSE,
##'                 epc = FALSE, standardized = epc, cov.std = epc,
##'                 verbose = FALSE, warn = TRUE, information = "expected")
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL


##' @rdname semTools-deprecated
##' @export
lavTestScore.mi <- function(object, add = NULL, release = NULL,
                            test = c("D2","D1"), scale.W = !asymptotic,
                            omit.imps = c("no.conv","no.se"),
                            asymptotic = is.null(add), # as F or chi-squared
                            univariate = TRUE, cumulative = FALSE,
                            epc = FALSE, standardized = epc, cov.std = epc,
                            verbose = FALSE, warn = TRUE, information = "expected") {

  .Deprecated(msg = c("\nThe runMI() and related lavaan.mi functions have been",
                      " deprecated and will cease to be included in future ",
                      "versions of semTools.\n\nSupport is still provided for ",
                      "analyzing lavaan.mi-class objects (e.g., compRelSEM() ",
                      "can estimate reliability using pooled results), which ",
                      "can now be created using the lavaan.mi package.\n\nThe ",
                      "deprecated runMI() function now creates an object of ",
                      "class OLDlavaan.mi, which can be analyzed using the ",
                      "deprecated functions in semTools, like lavTestLRT.mi(),",
                      " that have been updated and improved in the lavaan.mi ",
                      "package.\n\nFind more details help('semTools-deprecated)"))

  stopifnot(inherits(object, "OLDlavaan.mi"))
  lavoptions <- object@Options

  useImps <- rep(TRUE, length(object@DataList))
  if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
  if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
  if ("no.npd" %in% omit.imps) {
    Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
    Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
    useImps <- useImps & !(Heywood.lv | Heywood.ov)
  }
  ## custom removal by imputation number
  rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
  if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
  ## whatever is left
  m <- sum(useImps)
  if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
  useImps <- which(useImps)

  test <- toupper(test[1])
  if (!test %in% c("D2","D1")) stop('Invalid choice of "test" argument.')

  ## check if model has converged
  if (m == 0L) stop("No models converged. Score tests unavailable.")

  # check for inequality constraints
  PT <- parTable(object)
  if (any(PT$op == ">" | PT$op == "<")) {
    stop("lavTestScore.mi() does not handle inequality constraints (yet)")
  }

  # check arguments
  if (cumulative) univariate <- TRUE
  if (sum(is.null(release), is.null(add)) == 0) {
    stop("`add' and `release' arguments cannot be used together.\n",
         "Fixed parameters can instead be labeled in the model syntax ",
         "and those labels can be constrained to fixed values, so that ",
         "the constraints can be tested using the `release' argument along ",
         "with other released constraints.")
  }

  oldCall <- object@lavListCall
  #oldCall$model <- parTable(object) # FIXME: necessary?

  if (test == "D2") {
    if (!is.null(oldCall$parallel)) {
      if (oldCall$parallel == "snow") {
        oldCall$parallel <- "no"
        oldCall$ncpus <- 1L
        if (warn) warning("Unable to pass lavaan::lavTestScore() arguments ",
                          "when parallel='snow'. Switching to parallel='no'.",
                          " Unless using Windows, parallel='multicore' works.")
      }
    }

    ## call lavaanList() again to run lavTestScore() on each imputation
    oldCall$FUN <- function(obj) {
      out <- try(lavaan::lavTestScore(obj, add = add, release = release,
                                      univariate = univariate, epc = epc,
                                      cumulative = cumulative, cov.std = cov.std,
                                      standardized = standardized,
                                      warn = FALSE, information = information),
                 silent = TRUE)
      if (inherits(out, "try-error")) return(NULL)
      out
    }
    FIT <- eval(as.call(oldCall))
    ## check if there are any results
    noScores <- sapply(FIT@funList, is.null)
    if (all(noScores)) stop("No success using lavTestScore() on any imputations.")

    ## template to fill in pooled values
    OUT <- FIT@funList[[ intersect(useImps, which(!noScores))[1] ]]

    ## at a minimum, pool the total score test
    chiList <- sapply(FIT@funList[intersect(useImps, which(!noScores))],
                      function(x) x$test$X2)
    chiPooled <- calculate.D2(chiList, DF = OUT$test$df, asymptotic)
    OUT$test$X2 <- chiPooled[1]
    if (!asymptotic) {
      names(OUT$test)[names(OUT$test) == "X2"] <- "F"
      names(OUT$test)[names(OUT$test) == "df"] <- "df1"
      OUT$test$df2 <- chiPooled[["df2"]]
      OUT$test$p.value <- NULL # so it appears after "df2" column
    }
    OUT$test$p.value <- chiPooled[["pvalue"]]
    OUT$test$ariv <- chiPooled[["ariv"]]
    OUT$test$fmi <- chiPooled[["fmi"]]

    ## univariate?
    if (univariate) {
      if (!asymptotic) {
        names(OUT$uni)[names(OUT$uni) == "X2"] <- "F"
        OUT$uni$p.value <- NULL # so it appears after "df2" column
        OUT$uni$df2 <- NA
        OUT$uni$p.value <- NA
        OUT$uni$riv <- NA
        OUT$uni$fmi <- NA
        if ("epc" %in% colnames(OUT$uni)) {
          OUT$uni$epc <- NULL # so it appears after "p.value"
          OUT$uni$epc <- NA # fill in below
        }
      }
      for (i in 1:nrow(OUT$uni)) {
        chiList <- sapply(FIT@funList[intersect(useImps, which(!noScores))],
                          function(x) x$uni$X2[i] )
        chiPooled <- calculate.D2(chiList, DF = OUT$uni$df[i], asymptotic)
        if (!asymptotic) {
          OUT$uni$F[i] <- chiPooled[[1]]
          OUT$uni$df2[i] <- chiPooled[["df2"]]
        } else OUT$uni$X2[i] <- chiPooled[[1]]
        OUT$uni$p.value[i] <- chiPooled[["pvalue"]]
        OUT$uni$riv[i] <- chiPooled[["ariv"]]
        OUT$uni$fmi[i]  <- chiPooled[["fmi"]]
        if ("epc" %in% colnames(OUT$uni)) {
          epcUniList <- sapply(FIT@funList[intersect(useImps, which(!noScores))],
                               function(x) x$uni$epc[i] )
          OUT$uni$epc[i] <- mean(epcUniList)
        }
      }
      if (!asymptotic) names(OUT$uni)[names(OUT$uni) == "df"] <- "df1"
    }

    ## cumulative?
    if (cumulative) {
      if (!asymptotic) {
        names(OUT$cumulative)[names(OUT$cumulative) == "X2"] <- "F"
        OUT$cumulative$p.value <- NULL # so it appears after "df2" column
        OUT$cumulative$df2 <- NA
        OUT$cumulative$p.value <- NA
        OUT$cumulative$ariv <- NA
        OUT$cumulative$fmi <- NA
      }
      for (i in 1:nrow(OUT$cumulative)) {
        chiList <- sapply(FIT@funList[intersect(useImps, which(!noScores))],
                          function(x) x$cumulative$X2[i] )
        chiPooled <- calculate.D2(chiList, DF = OUT$cumulative$df[i], asymptotic)
        if (!asymptotic) {
          OUT$cumulative$F[i] <- chiPooled[[1]]
          OUT$cumulative$df2[i] <- chiPooled[["df2"]]
        } else OUT$cumulative$X2[i] <- chiPooled[[1]]
        OUT$cumulative$p.value[i] <- chiPooled[["pvalue"]]
        OUT$cumulative$ariv[i] <- chiPooled[["ariv"]]
        OUT$cumulative$fmi[i]  <- chiPooled[["fmi"]]
      }
      if (!asymptotic) names(OUT$cumulative)[names(OUT$cumulative) == "df"] <- "df1"
    }

    ## EPCs?
    if (epc) {
      estList <- lapply(FIT@funList[intersect(useImps, which(!noScores))],
                        function(x) x$epc$est)
      OUT$epc$est <- rowMeans(do.call(cbind, estList))

      epcList <- lapply(FIT@funList[intersect(useImps, which(!noScores))],
                        function(x) x$epc$epc)
      OUT$epc$epc <- rowMeans(do.call(cbind, epcList))

      OUT$epc$epv <- OUT$epc$est + OUT$epc$epc
      if (standardized) {
        sepcList <- lapply(FIT@funList[intersect(useImps, which(!noScores))],
                           function(x) x$epc$sepc.lv)
        OUT$epc$sepc.lv <- rowMeans(do.call(cbind, sepcList))
        sepcList <- lapply(FIT@funList[intersect(useImps, which(!noScores))],
                           function(x) x$epc$sepc.all)
        OUT$epc$sepc.all <- rowMeans(do.call(cbind, sepcList))
        if ("sepc.nox" %in% colnames(FIT@funList[intersect(useImps, which(!noScores))][[1]])) {
          sepcList <- lapply(FIT@funList[intersect(useImps, which(!noScores))],
                             function(x) x$epc$sepc.nox)
          OUT$epc$sepc.nox <- rowMeans(do.call(cbind, sepcList))
        }
      }
    }

    return(OUT)
  } # else test == "D1", making 'scale.W=' relevant

  ## number of free parameters (regardless of whether they are constrained)
  npar <- object@Model@nx.free
  ## sample size
  N <- lavListInspect(object, "ntotal")
  if (lavoptions$mimic == "EQS") N <- N - 1

  # Mode 1: ADDING new parameters
  if (!is.null(add) && nchar(add) > 0L) {
    ## turn off SNOW cluster (can't past arguments)
    if (!is.null(oldCall$parallel)) {
      if (oldCall$parallel == "snow") {
        oldCall$parallel <- "no"
        oldCall$ncpus <- 1L
        if (warn) warning("Unable to pass lavaan::lavTestScore() arguments ",
                          "when parallel='snow'. Switching to parallel='no'.",
                          " Unless using Windows, parallel='multicore' works.")
      }
    }

    ## call lavaanList() to fit augmented model (do.fit = FALSE)
    oldCall$FUN <- function(obj) {
      ngroups <- lavaan::lavInspect(obj, "ngroups")
      nlevels <- lavaan::lavInspect(obj, "nlevels")

      ## --------------------------------------
      ## borrowed code from lav_object_extend()
      ## --------------------------------------

      ## select columns that should always be included below
      myCols <- c("lhs","op","rhs")
      if (ngroups > 1L) myCols <- c(myCols,"block","group")
      if (nlevels > 1L) myCols <- c(myCols,"block","level")
      myCols <- unique(myCols)

      # partable original model
      oldPT <- lavaan::parTable(obj)[c(myCols, "free","label","plabel")]
      # replace 'start' column, since lav_model will fill these in in GLIST
      oldPT$start <- lavaan::parameterEstimates(obj, remove.system.eq = FALSE,
                                                remove.def = FALSE,
                                                remove.eq = FALSE,
                                                remove.ineq = FALSE)$est

      # add new parameters, extend model
      # ADD <- lavaan::modindices(obj, standardized = FALSE)[myCols]
      if (is.list(add)) {
        stopifnot(!is.null(add$lhs),
                  !is.null(add$op),
                  !is.null(add$rhs))
        ADD <- as.data.frame(add)
      } else if (is.character(add)) {
        ADD <- lavaan::lavaanify(add, ngroups = ngroups)
        ADD <- ADD[ , c(myCols, "user","label")]
        remove.idx <- which(ADD$user == 0)
        if (length(remove.idx) > 0L) {
          ADD <- ADD[-remove.idx,]
        }
        ADD$start <- rep( 0, nrow(ADD))
        ADD$free  <- rep( 1, nrow(ADD))
        ADD$user  <- rep(10, nrow(ADD))
      } else stop("'add' must be lavaan model syntax or a parameter table.")
      # nR <- try(nrow(ADD), silent = TRUE)
      # if (inherits(nR, "try-error") || is.null(nR)) return(list(gradient = NULL,
      #                                                           information = NULL))
      # ADD$free <- rep(1L, nR)
      # ADD$user <- rep(10L, nR)

      # merge
      LIST <- lavaan::lav_partable_merge(oldPT, ADD, remove.duplicated = TRUE, warn = FALSE)
      # redo 'free'
      free.idx <- which(LIST$free > 0)
      LIST$free[free.idx] <- 1:length(free.idx)
      # adapt options
      lavoptions <- obj@Options
      if (any(LIST$op == "~1")) lavoptions$meanstructure <- TRUE
      lavoptions$do.fit <- FALSE

      obj2 <- lavaan::lavaan(LIST,
                             slotOptions     = lavoptions,
                             slotSampleStats = obj@SampleStats,
                             slotData        = obj@Data,
                             slotCache       = obj@Cache,
                             sloth1          = obj@h1)
      ## ---------------------------------
      list(gradient = lavaan::lavInspect(obj2, "gradient"),
           information = lavaan::lavInspect(obj2, paste("information",
                                                        information, sep = ".")),
           #TODO: Max wants to calculate EPCs as averages across imputations.
           #      Run lavTestScore(epc=TRUE) here?  or be consistent...
           nadd = nrow(ADD), parTable = lavaan::parTable(obj2))
    }
    FIT <- eval(as.call(oldCall))

    ## obtain list of inverted Jacobians: within-impuation covariance matrices
    R.model <- object@Model@con.jac[,,drop = FALSE]
    nadd <- FIT@funList[[ useImps[1] ]]$nadd

    ## pool gradients and information matrices
    gradList <- lapply(FIT@funList[useImps], "[[", i = "gradient")
    infoList <- lapply(FIT@funList[useImps], "[[", i = "information")
    score <- colMeans(do.call(rbind, gradList))  # pooled point estimates
    B <- cov(do.call(rbind, gradList) * sqrt(N)) # between-imputation UNIT information
    W <- Reduce("+", infoList) / m               # within-imputation UNIT information
    inv.W <- try(solve(W), silent = TRUE)
    if (inherits(inv.W, "try-error")) {
      if (warn && scale.W) warning("Could not invert W for total score test, ",
                                   "perhaps due to constraints on estimated ",
                                   "parameters. Generalized inverse used instead.\n",
                                   "If the model does not have equality constraints, ",
                                   "it may be safer to set `scale.W = FALSE'.")
      inv.W <- MASS::ginv(W)
    }
    ## relative increase in variance due to missing data
    ariv <- (1 + 1/m)/nadd * sum(diag(B %*% inv.W)) # ONLY relevant for scaling full INFO matrix

    if (scale.W) {
      information <- (1 + ariv) * W  # Enders (2010, p. 235) eqs. 8.20-21
    } else {
      ## less reliable, but constraints prevent inversion of W
      information <- W + B + (1/m)*B  # Enders (2010, p. 235) eq. 8.19
    }

    if (nrow(R.model) > 0L) {
      R.model <- cbind(R.model, matrix(0, nrow(R.model), ncol = nadd))
      R.add   <- cbind(matrix(0, nrow = nadd, ncol = npar), diag(nadd))
      R       <- rbind(R.model, R.add)

      Z <- cbind(rbind(information, R.model),
                 rbind(t(R.model),matrix(0,nrow(R.model),nrow(R.model))))
      Z.plus <- MASS::ginv(Z)
      J.inv  <- Z.plus[ 1:nrow(information), 1:nrow(information) ]

      r.idx <- seq_len(nadd) + nrow(R.model)
    } else {
      R <- cbind(matrix(0, nrow = nadd, ncol = npar), diag(nadd))
      J.inv <- MASS::ginv(information)

      r.idx <- seq_len(nadd)
    }

    PT <- FIT@funList[[ useImps[1] ]]$parTable
    if (is.null(PT$group)) PT$group <- PT$block
    # lhs/rhs
    lhs <- lavaan::lav_partable_labels(PT)[ PT$user == 10L ]
    op <- rep("==", nadd)
    rhs <- rep("0", nadd)
    Table <- data.frame(lhs = lhs, op = op, rhs = rhs)
    class(Table) <- c("lavaan.data.frame", "data.frame")
  } else {
    # MODE 2: releasing constraints
    if (!asymptotic) {
      message('The average relative increase in variance (ARIV) cannot be ',
              'calculated for releasing estimated constraints, preventing the ',
              'denominator degrees of freedom from being calculated for the F ',
              'test, so the "asymptotic" argument was switched to TRUE.' )
      asymptotic <- TRUE
      scale.W <- FALSE
    }
    if (is.character(release)) stop("not implemented yet") #FIXME: moved up to save time
    R <- object@Model@con.jac[,,drop = FALSE]
    if (nrow(R) == 0L) stop("No equality constraints found in the model.")


    ## use lavaanList() to get gradient/information from each imputation
    oldCall$FUN <- function(obj) {
      list(gradient = lavaan::lavInspect(obj, "gradient"),
           information = lavaan::lavInspect(obj, paste("information",
                                                       information, sep = ".")))
    }
    FIT <- eval(as.call(oldCall))
    ## pool gradients and information matrices
    gradList <- lapply(FIT@funList[useImps], "[[", i = "gradient")
    infoList <- lapply(FIT@funList[useImps], "[[", i = "information")
    score <- colMeans(do.call(rbind, gradList))  # pooled point estimates
    B <- cov(do.call(rbind, gradList) * sqrt(N)) # between-imputation UNIT information
    W <- Reduce("+", infoList) / m               # within-imputation UNIT information
    inv.W <- try(solve(W), silent = TRUE)
    if (inherits(inv.W, "try-error")) {
      if (warn && scale.W) warning("Could not invert W for total score test, ",
                                   "perhaps due to constraints on estimated ",
                                   "parameters. Generalized inverse used instead.\n",
                                   "If the model does not have equality constraints, ",
                                   "it may be safer to set `scale.W = FALSE'.")
      inv.W <- MASS::ginv(W)
    }
    ## relative increase in variance due to missing data
    k <- length(release)
    if (k == 0) k <- nrow(R)
    ariv <- (1 + 1/m)/k * sum(diag(B %*% inv.W)) #FIXME: Invalid extrapolation!
    if (scale.W) {
      #TODO: This option is disabled, kept only to update with future research
      information <- (1 + ariv) * W  # Enders (2010, p. 235) eqs. 8.20-21
    } else {
      ## less reliable, but constraints prevent inversion of W
      information <- W + B + (1/m)*B  # Enders (2010, p. 235) eq. 8.19
    }

    if (is.null(release)) {
      # ALL constraints
      r.idx <- seq_len( nrow(R) )
      J.inv <- MASS::ginv(information) #FIXME? Yves has this above if(is.null(release))
    } else if (is.numeric(release)) {
      r.idx <- release
      if(max(r.idx) > nrow(R)) {
        stop("lavaan ERROR: maximum constraint number (", max(r.idx),
             ") is larger than number of constraints (", nrow(R), ")")
      }

      # neutralize the non-needed constraints
      R1 <- R[-r.idx, , drop = FALSE]
      Z1 <- cbind( rbind(information, R1),
                   rbind(t(R1), matrix(0, nrow(R1), nrow(R1))) )
      Z1.plus <- MASS::ginv(Z1)
      J.inv <- Z1.plus[ 1:nrow(information), 1:nrow(information) ]
    } else if (is.character(release)) {
      stop("not implemented yet")
    }


    # lhs/rhs
    eq.idx <- which(object@ParTable$op == "==")
    if (length(eq.idx) > 0L) {
      lhs <- object@ParTable$lhs[eq.idx][r.idx]
      op <- rep("==", length(r.idx))
      rhs <- object@ParTable$rhs[eq.idx][r.idx]
    }
    Table <- data.frame(lhs = lhs, op = op, rhs = rhs)
    class(Table) <- c("lavaan.data.frame", "data.frame")
  }

  if (lavoptions$se == "standard") {
    stat <- as.numeric(N * score %*% J.inv %*% score)
  } else {
    # generalized score test
    if (warn) warning("se is not `standard'. Robust test not implemented yet. ",
                      "Falling back to ordinary score test.")
    # NOTE!!!
    # we can NOT use VCOV here, because it reflects the constraints,
    # and the whole point is to test for these constraints...

    stat <- as.numeric(N * score %*% J.inv %*% score)
  }

  # compute df, taking into account that some of the constraints may
  # be needed to identify the model (and hence information is singular)
  # information.plus <- information + crossprod(R)
  #df <- qr(R[r.idx,,drop = FALSE])$rank +
  #          ( qr(information)$rank - qr(information.plus)$rank )
  DF <- nrow( R[r.idx, , drop = FALSE] )
  if (asymptotic) {
    TEST <- data.frame(test = "score", X2 = stat, df = DF,
                       p.value = pchisq(stat, df = DF, lower.tail = FALSE))
  } else {
    ## calculate denominator DF for F statistic
    myDims <- 1:nadd + npar #TODO: not available in release mode, unless calculating Jacobian of constraints, like Wald test?
    ARIV <- (1 + 1/m)/nadd * sum(diag(B[myDims, myDims, drop = FALSE] %*% inv.W[myDims, myDims, drop = FALSE]))
    a <- DF*(m - 1)
    if (a > 4) {
      df2 <- 4 + (a - 4) * (1 + (1 - 2/a)*(1 / ARIV))^2 # Enders (eq. 8.24)
    } else {
      df2 <- a*(1 + 1/DF) * (1 + 1/ARIV)^2 / 2 # Enders (eq. 8.25)
    }
    TEST <- data.frame(test = "score", "F" = stat / DF, df1 = DF, df2 = df2,
                       p.value = pf(stat / DF, df1 = DF, df2 = df2, lower.tail = FALSE),
                       ariv = ARIV, fmi = ARIV / (1 + ARIV))
  }
  class(TEST) <- c("lavaan.data.frame", "data.frame")
  attr(TEST, "header") <- "total score test:"
  OUT <- list(test = TEST)

  if (univariate) {
    TS <- numeric( nrow(R) )
    EPC.uni <- numeric( nrow(R) ) #FIXME: to add univariate EPCs for added parameters
    for (r in r.idx) {
      R1 <- R[-r, , drop = FALSE]
      Z1 <- cbind( rbind(information, R1),
                   rbind(t(R1), matrix(0, nrow(R1), nrow(R1))) )
      Z1.plus <- MASS::ginv(Z1)
      Z1.plus1 <- Z1.plus[ 1:nrow(information), 1:nrow(information) ]
      TS[r] <- as.numeric(N * t(score) %*%  Z1.plus1 %*% score)

      ## FIXME: experimentally add univariate EPCs for added parameters, as would accompany modification indices
      if (epc && !is.null(add)) EPC.uni[r] <- 1 * utils::tail(as.numeric(score %*%  Z1.plus1), n = nrow(R))[r]
    }

    Table2 <- Table
    DF <- rep(1, length(r.idx))

    if (asymptotic) {
      Table2$X2 <- TS[r.idx]
      Table2$df <- DF
      Table2$p.value <- pchisq(Table2$X2, df = DF, lower.tail = FALSE)
    } else {
      Table2$F <- TS[r.idx] / DF
      Table2$df1 <- DF
      ## calculate denominator DF for F statistic using RIV per 1-df test (Enders eq. 8.10)
      myDims <- 1:nadd + npar
      RIVs <- diag((1 + 1/m) * B[myDims, myDims, drop = FALSE]) / diag(W[myDims, myDims, drop = FALSE])
      Table2$df2 <- sapply(RIVs, function(riv) {
        DF1 <- 1L # Univariate tests
        a <- DF1*(m - 1)
        DF2 <- if (a > 4) {
          4 + (a - 4) * (1 + (1 - 2/a)*(1 / riv))^2 # Enders (eq. 8.24)
        } else a*(1 + 1/DF1) * (1 + 1/riv)^2 / 2 # Enders (eq. 8.25)
        DF2
      })
      Table2$p.value <- pf(Table2$F, df1 = DF, df2 = Table2$df2, lower.tail = FALSE)
      Table2$riv <- RIVs
      Table2$fmi <- RIVs / (1 + RIVs)
    }

    ## add univariate EPCs, equivalent to modindices() output
    if (epc && !is.null(add)) Table2$epc <- EPC.uni[r.idx]

    attr(Table2, "header") <- "univariate score tests:"
    OUT$uni <- Table2
  }

  if (cumulative) {
    TS.order <- sort.int(TS, index.return = TRUE, decreasing = TRUE)$ix
    TS <- numeric( length(r.idx) )
    if (!asymptotic) ARIVs <- numeric( length(r.idx) )
    for (r in 1:length(r.idx)) {
      rcumul.idx <- TS.order[1:r]

      R1 <- R[-rcumul.idx, , drop = FALSE]
      Z1 <- cbind( rbind(information, R1),
                   rbind(t(R1), matrix(0, nrow(R1), nrow(R1))) )
      Z1.plus <- MASS::ginv(Z1)
      Z1.plus1 <- Z1.plus[ 1:nrow(information), 1:nrow(information) ]
      TS[r] <- as.numeric(N * t(score) %*%  Z1.plus1 %*% score)
      if (!asymptotic) {
        myDims <- rcumul.idx + npar
        ARIVs[r] <- (1 + 1/m)/length(myDims) * sum(diag(B[myDims, myDims, drop = FALSE] %*% inv.W[myDims, myDims, drop = FALSE]))
      }
    }

    Table3 <- Table[TS.order,]
    DF <- seq_len( length(TS) )
    if (asymptotic) {
      Table3$X2 <- TS
      Table3$df <- DF
      Table3$p.value <- pchisq(Table3$X2, df = DF, lower.tail = FALSE)
    } else {
      Table3$F <- TS / DF
      Table3$df1 <- DF
      ## calculate denominator DF for F statistic
      Table3$df2 <- mapply(FUN = function(DF1, ariv) {
        a <- DF1*(m - 1)
        DF2 <- if (a > 4) {
          4 + (a - 4) * (1 + (1 - 2/a)*(1 / ariv))^2 # Enders (eq. 8.24)
        } else a*(1 + 1/DF1) * (1 + 1/ariv)^2 / 2 # Enders (eq. 8.25)
        DF2
      }, DF1 = DF, ariv = ARIVs)
      Table3$p.value = pf(Table3$F, df1 = DF, df2 = Table3$df2, lower.tail = FALSE)
      Table3$riv <- ARIVs
      Table3$fmi <- ARIVs / (1 + ARIVs)
    }
    attr(Table3, "header") <- "cumulative score tests:"
    OUT$cumulative <- Table3
  }

  if (epc) {
    ngroups <- lavaan::lavInspect(object, "ngroups")
    nlevels <- lavaan::lavInspect(object, "nlevels")

    ################# source code Yves commented out.
    ################# Calculates 1 EPC-vector per constraint.
    ################# Better to call lavTestScore() multiple times?  Ugh...
    # EPC <- vector("list", length = length(r.idx))
    # for (i in 1:length(r.idx)) {
    #     r <- r.idx[i]
    #     R1 <- R[-r,,drop = FALSE]
    #     Z1 <- cbind( rbind(information, R1),
    #                  rbind(t(R1), matrix(0,nrow(R1),nrow(R1))) )
    #     Z1.plus <- MASS::ginv(Z1)
    #     Z1.plus1 <- Z1.plus[ 1:nrow(information), 1:nrow(information) ]
    #     EPC[[i]] <- -1 * as.numeric(score %*%  Z1.plus1)
    # }
    # OUT$EPC <- EPC

    # EPCs when freeing all constraints together (total test)
    R1 <- R[-r.idx, , drop = FALSE]
    Z1 <- cbind( rbind(information, R1),
                 rbind(t(R1), matrix(0, nrow(R1), nrow(R1))) )
    Z1.plus <- MASS::ginv(Z1)
    Z1.plus1 <- Z1.plus[ 1:nrow(information), 1:nrow(information) ]
    EPC.all <- 1 * as.numeric(score %*%  Z1.plus1)

    # create epc table for the 'free' parameters
    myCoefs <- getMethod("coef","OLDlavaan.mi")(object, omit.imps = omit.imps)
    myCols <- c("lhs","op","rhs","user")
    if (ngroups > 1L) myCols <- c(myCols, "block","group")
    if (nlevels > 1L) myCols <- c(myCols, "block","level")
    myCols <- c(unique(myCols), "free","exo","label","plabel")
    LIST <- if (!is.null(add) && nchar(add) > 0L) {
      PT[ , myCols]
    } else parTable(object)[ , myCols]

    nonpar.idx <- which(LIST$op %in% c("==", ":=", "<", ">"))
    if (length(nonpar.idx) > 0L) LIST <- LIST[ -nonpar.idx , ]

    LIST$est[ LIST$free > 0 & LIST$user != 10 ] <- myCoefs
    LIST$est[ LIST$user == 10L ] <- 0
    LIST$epc <- rep(as.numeric(NA), length(LIST$lhs))
    LIST$epc[ LIST$free > 0 ] <- EPC.all
    LIST$epv <- LIST$est + LIST$epc
    #TODO: add SEPCs
    if (standardized) {

      EPC <- LIST$epc

      if (cov.std) {
        # replace epc values for variances by est values
        var.idx <- which(LIST$op == "~~" & LIST$lhs == LIST$rhs & LIST$exo == 0L)
        EPC[ var.idx ] <- LIST$est[ var.idx ]
      }

      # two problems:
      #   - EPC of variances can be negative, and that is perfectly legal
      #   - EPC (of variances) can be tiny (near-zero), and we should
      #     not divide by tiny variables
      small.idx <- which(LIST$op == "~~" &
                         LIST$lhs == LIST$rhs &
                         abs(EPC) < sqrt( .Machine$double.eps ) )
      if (length(small.idx) > 0L) EPC[ small.idx ] <- as.numeric(NA)

      # get the sign
      EPC.sign <- sign(LIST$epc)

      ## pooled estimates for standardizedSolution()
      pooledest <- getMethod("coef", "OLDlavaan.mi")(object,
                                                     omit.imps = omit.imps)
      ## update @Model@GLIST for standardizedSolution(..., GLIST=)
      object@Model <- lavaan::lav_model_set_parameters(object@Model, x = pooledest)

      LIST$sepc.lv <- EPC.sign * lavaan::standardizedSolution(object, se = FALSE,
                                                              type = "std.lv",
                                                              cov.std = cov.std,
                                                              partable = LIST,
                                                              GLIST = object@Model@GLIST,
                                                              est = abs(EPC))$est.std
      LIST$sepc.all <- EPC.sign * lavaan::standardizedSolution(object, se = FALSE,
                                                               type = "std.all",
                                                               cov.std = cov.std,
                                                               partable = LIST,
                                                               GLIST = object@Model@GLIST,
                                                               est = abs(EPC))$est.std
      fixed.x <- lavListInspect(object, "options")$fixed.x && length(lavNames(object, "ov.x"))
      if (fixed.x) {
        LIST$sepc.nox <- EPC.sign * lavaan::standardizedSolution(object, se = FALSE,
                                                                 type = "std.nox",
                                                                 cov.std = cov.std,
                                                                 partable = LIST,
                                                                 GLIST = object@Model@GLIST,
                                                                 est = abs(EPC))$est.std
      }

      if (length(small.idx) > 0L) {
        LIST$sepc.lv[small.idx] <- 0
        LIST$sepc.all[small.idx] <- 0
        if (fixed.x) LIST$sepc.nox[small.idx] <- 0
      }

    }

    LIST$free[ LIST$user == 10L ] <- 0
    LIST$user <- NULL
    LIST$exo <- NULL

    DF <- if (asymptotic) OUT$test$df else OUT$test$df1
    attr(LIST, "header") <- paste0("expected parameter changes (epc) and ",
                                   "expected parameter values (epv)",
                                   if (DF < 2) ":" else {
                  " if ALL constraints in 'add' or 'release' were freed:" })

    OUT$epc <- LIST
  }

  OUT
}
