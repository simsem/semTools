### Terrence D. Jorgensen & Yves Rosseel
### Last updated: 12 March 2025
### adaptation of lavaan::modindices() for OLDlavaan.mi-class objects

### DEPRECATED: 16 June 2024
### supplanted by lavaan.mi package


##' Modification Indices for Multiple Imputations
##'
##' Modification indices (1-*df* Lagrange multiplier tests) from a
##' latent variable model fitted to multiple imputed data sets. Statistics
##' for releasing one or more fixed or constrained parameters in model can
##' be calculated by pooling the gradient and information matrices
##' across imputed data sets in a method proposed by Mansolf, Jorgensen, &
##' Enders (2020)---analogous to the "D1" Wald test proposed by Li, Meng,
##' Raghunathan, & Rubin (1991)---or by pooling the complete-data score-test
##' statistics across imputed data sets (i.e., "D2"; Li et al., 1991).
##'
## @aliases modificationIndices.mi modificationindices.mi modindices.mi
##' @importFrom lavaan lavInspect lavListInspect lavNames
##' @importFrom methods getMethod
##' @importFrom stats cov pchisq qchisq
##'
##' @param object An object of class [OLDlavaan.mi-class]
##' @param test `character` indicating which pooling method to use.
##'  `"D1"` requests Mansolf, Jorgensen, & Enders' (2020) proposed
##'  Wald-like test for pooling the gradient and information, which are then
##'  used to calculate score-test statistics in the usual manner. `"D2"`
##'  (default because it is less computationall intensive) requests to pool the
##'  complete-data score-test statistics from each imputed data set, then pool
##'  them across imputations, described by Li et al. (1991) and Enders (2010).
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
##' @param standardized `logical`. If `TRUE`, two extra columns
##'   (`$sepc.lv` and `$sepc.all`) will contain standardized values
##'   for the EPCs. In the first column (`$sepc.lv`), standardizization is
##'   based on the variances of the (continuous) latent variables. In the second
##'   column (`$sepc.all`), standardization is based on both the variances
##'   of both (continuous) observed and latent variables. (Residual) covariances
##'   are standardized using (residual) variances.
##' @param cov.std `logical`. `TRUE` if `test == "D2"`.
##'   If `TRUE` (default), the (residual)
##'   observed covariances are scaled by the square-root of the diagonal elements
##'   of the \eqn{\Theta} matrix, and the (residual) latent covariances are
##'   scaled by the square-root of the diagonal elements of the \eqn{\Psi}
##'   matrix. If `FALSE`, the (residual) observed covariances are scaled by
##'   the square-root of the diagonal elements of the model-implied covariance
##'   matrix of observed variables (\eqn{\Sigma}), and the (residual) latent
##'   covariances are scaled by the square-root of the diagonal elements of the
##'   model-implied covariance matrix of the latent variables.
##' @param information `character` indicating the type of information
##'   matrix to use (check [lavaan::lavInspect()] for available options).
##'   `"expected"` information is the default, which provides better
##'   control of Type I errors.
##' @param power `logical`. If `TRUE`, the (post-hoc) power is
##'   computed for each modification index, using the values of `delta`
##'   and `alpha`.
##' @param delta The value of the effect size, as used in the post-hoc power
##'   computation, currently using the unstandardized metric of the `$epc`
##'   column.
##' @param alpha The significance level used for deciding if the modification
##'   index is statistically significant or not.
##' @param high.power If the computed power is higher than this cutoff value,
##'   the power is considered 'high'. If not, the power is considered 'low'.
##'   This affects the values in the `$decision` column in the output.
##' @param sort. `logical`. If `TRUE`, sort the output using the
##'   values of the modification index values. Higher values appear first.
##' @param minimum.value `numeric`. Filter output and only show rows with a
##'   modification index value equal or higher than this minimum value.
##' @param maximum.number `integer`. Filter output and only show the first
##'   maximum number rows. Most useful when combined with the `sort.` option.
##' @param na.remove `logical`. If `TRUE` (default), filter output by
##'   removing all rows with `NA` values for the modification indices.
##' @param op `character` string. Filter the output by selecting only those
##'   rows with operator `op`.
##'
##' @note When `test = "D2"`, each (S)EPC will be pooled by taking its
##'   average across imputations. When `test = "D1"`, EPCs will be
##'   calculated in the standard way using the pooled gradient and information,
##'   and SEPCs will be calculated by standardizing the EPCs using model-implied
##'   (residual) variances.
##'
##' @return A `data.frame` containing modification indices and (S)EPCs.
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##'   Adapted from \pkg{lavaan} source code, written by
##'   Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
##'
##' `test = "D1"` method proposed by
##'   Maxwell Mansolf (University of California, Los Angeles;
##'   \email{mamansolf@@gmail.com})
##'
##' @references
##'   Enders, C. K. (2010). *Applied missing data analysis*.
##'   New York, NY: Guilford.
##'
##'   Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
##'   Significance levels from repeated *p*-values with multiply-imputed
##'   data.*Statistica Sinica, 1*(1), 65--92. Retrieved from
##'   <https://www.jstor.org/stable/24303994>
##'
##'   Mansolf, M., Jorgensen, T. D., & Enders, C. K. (2020). A multiple
##'   imputation score test for model modification in structural equation
##'   models. *Psychological Methods, 25*(4), 393--411.
##'   \doi{10.1037/met0000243}
##'
##' @examples
##'
##' ## See the new lavaan.mi package
##'
##' @name modindices.mi-deprecated
##' @usage
##' modindices.mi(object, test = c("D2","D1"), omit.imps = c("no.conv","no.se"),
##'               standardized = TRUE, cov.std = TRUE, information = "expected",
##'               power = FALSE, delta = 0.1, alpha = 0.05, high.power = 0.75,
##'               sort. = FALSE, minimum.value = 0, maximum.number = nrow(LIST),
##'               na.remove = TRUE, op = NULL)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL


##' @rdname semTools-deprecated
##' @export
modindices.mi <- function(object,
                          test = c("D2","D1"),
                          omit.imps = c("no.conv","no.se"),

                          standardized = TRUE,
                          cov.std = TRUE,
                          information = "expected",

                          # power statistics?
                          power = FALSE,
                          delta = 0.1,
                          alpha = 0.05,
                          high.power = 0.75,

                          # customize output
                          sort. = FALSE,
                          minimum.value = 0.0,
                          maximum.number = nrow(LIST),
                          na.remove = TRUE,
                          op = NULL) {

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

  test <- tolower(test[1])
  N <- lavListInspect(object, "ntotal")
  #FIXME: if (lavoptions$mimic == "EQS") N <- N - 1 # not in lavaan, why?

  # not ready for estimator = "PML"
  if (object@Options$estimator == "PML") {
      stop("Modification indices not yet implemented for estimator PML.")
  }

  # sanity check
  if (power) standardized <- TRUE

  ## use first available modification indices as template to store pooled results
  ngroups <- lavListInspect(object, "ngroups")
  nlevels <- lavListInspect(object, "nlevels")
  myCols <- c("lhs","op","rhs")
  if (ngroups > 1L) myCols <- c(myCols,"block","group")
  if (nlevels > 1L) myCols <- c(myCols,"block","level")
  myCols <- unique(myCols)

  for (i in useImps) {
    LIST <- object@miList[[i]][myCols]
    nR <- try(nrow(LIST), silent = TRUE)
    if (inherits(nR, "try-error") || is.null(nR)) {
      if (i == max(useImps)) {
        stop("No modification indices could be computed for any imputations.")
      } else next
    } else break
  }



  ## D2 pooling method
  if (test == "d2") {
    chiList <- lapply(object@miList[useImps], "[[", i = "mi")
    ## imputations in columns, parameters in rows
    pooledList <- apply(do.call(cbind, chiList), 1, function(x) {
      calculate.D2(x, DF = 1, asymptotic = TRUE)
    })
    LIST$mi <- pooledList[1, ] # could be "F" or "chisq"
    ## diagnostics
    LIST$riv <- pooledList["ariv", ]
    LIST$fmi <- pooledList["fmi", ]
    ## also take average of epc & sepc.all
    epcList <- lapply(object@miList[useImps], "[[", i = "epc")
    LIST$epc <- rowMeans(do.call(cbind, epcList))
    if (standardized) {
      sepcList <- lapply(object@miList[useImps], "[[", i = "sepc.lv")
      LIST$sepc.lv <- rowMeans(do.call(cbind, sepcList))
      sepcList <- lapply(object@miList[useImps], "[[", i = "sepc.all")
      LIST$sepc.all <- rowMeans(do.call(cbind, sepcList))
      fixed.x <- lavListInspect(object, "options")$fixed.x && length(lavNames(object, "ov.x"))
      if (fixed.x && "sepc.nox" %in% colnames(object@miList[useImps][[1]])) {
        sepcList <- lapply(object@miList[useImps], "[[", i = "sepc.nox")
        LIST$sepc.nox <- rowMeans(do.call(cbind, sepcList))
      }
    }

  } else {

    scoreOut <- lavTestScore.mi(object, add = cbind(LIST, user = 10L,
                                                    free = 1, start = 0),
                                test = "d1", omit.imps = omit.imps,
                                epc = TRUE, scale.W = FALSE, asymptotic = TRUE,
                                information = information)$uni
    LIST$mi <- scoreOut$X2
    LIST$riv <- scoreOut$riv
    LIST$fmi <- scoreOut$fmi
    LIST$epc <- scoreOut$epc #FIXME: use average across imputations?

    # standardize?
    if (standardized) {
      ## Need full parameter table for lavaan::standardizedSolution()
      ## Merge parameterEstimates() with modindices()
      oldPE <- getMethod("summary","OLDlavaan.mi")(object, se = FALSE,
                                                   output = "data.frame",
                                                   omit.imps = omit.imps)
      PE <- lavaan::lav_partable_merge(oldPE, cbind(LIST, est = 0),
                                       remove.duplicated = TRUE, warn = FALSE)
      ## merge EPCs, using parameter labels (unavailable for estimates)
      rownames(LIST) <- paste0(LIST$lhs, LIST$op, LIST$rhs, ".g", LIST$group) #FIXME: multilevel?
      rownames(PE) <- paste0(PE$lhs, PE$op, PE$rhs, ".g", PE$group)
      PE[rownames(LIST), "epc"] <- LIST$epc
      ## need "exo" column?
      PT <- parTable(object)
      if ("exo" %in% names(PT)) {
        rownames(PT) <- paste0(PT$lhs, PT$op, PT$rhs, ".g", PT$group)
        PE[rownames(PT), "exo"] <- PT$exo
      } else PE$exo <- 0L
      rownames(LIST) <- NULL
      rownames(PE) <- NULL

      EPC <- PE$epc

      if (cov.std) {
        # replace epc values for variances by est values
        var.idx <- which(PE$op == "~~" & PE$lhs == PE$rhs)
        EPC[ var.idx ] <- PE$est[ var.idx ]
      }

      # two problems:
      #   - EPC of variances can be negative, and that is perfectly legal
      #   - EPC (of variances) can be tiny (near-zero), and we should
      #     not divide by tiny variables
      small.idx <- which(PE$op == "~~" &
                           PE$lhs == PE$rhs &
                           abs(EPC) < sqrt( .Machine$double.eps ) )
      if (length(small.idx) > 0L) EPC[small.idx] <- as.numeric(NA)

      # get the sign
      EPC.sign <- sign(PE$epc)

      ## pooled estimates for standardizedSolution()
      pooledest <- getMethod("coef", "OLDlavaan.mi")(object,
                                                     omit.imps = omit.imps)
      ## update @Model@GLIST for standardizedSolution(..., GLIST=)
      object@Model <- lavaan::lav_model_set_parameters(object@Model, x = pooledest)

      PE$sepc.lv <- EPC.sign * lavaan::standardizedSolution(object, se = FALSE,
                                                            type = "std.lv",
                                                            cov.std = cov.std,
                                                            partable = PE,
                                                            GLIST = object@Model@GLIST,
                                                            est = abs(EPC))$est.std
      PE$sepc.all <- EPC.sign * lavaan::standardizedSolution(object, se = FALSE,
                                                             type = "std.all",
                                                             cov.std = cov.std,
                                                             partable = PE,
                                                             GLIST = object@Model@GLIST,
                                                             est = abs(EPC))$est.std
      fixed.x <- lavListInspect(object, "options")$fixed.x && length(lavNames(object, "ov.x"))
      if (fixed.x) {
        PE$sepc.nox <- EPC.sign * lavaan::standardizedSolution(object, se = FALSE,
                                                               type = "std.nox",
                                                               cov.std = cov.std,
                                                               partable = PE,
                                                               GLIST = object@Model@GLIST,
                                                               est = abs(EPC))$est.std
      }

      if (length(small.idx) > 0L) {
        PE$sepc.lv[small.idx] <- 0
        PE$sepc.all[small.idx] <- 0
        if (fixed.x) PE$sepc.nox[small.idx] <- 0
      }
      ## remove unnecessary columns, then merge
      if (is.null(LIST$block)) PE$block <- NULL
      PE$est <- NULL
      PE$mi <- NULL
      PE$epc <- NULL
      PE$exo <- NULL
      LIST <- merge(LIST, PE, sort = FALSE)
      class(LIST) <- c("lavaan.data.frame","data.frame")
    }
  }

  # power?
  if (power) {
    LIST$sepc.lv <- NULL
    LIST$delta <- delta
    # FIXME: this is using epc in unstandardized metric
    #        this would be much more useful in standardized metric
    #        we need a standardize.est.all.reverse function...
    LIST$ncp <- (LIST$mi / (LIST$epc*LIST$epc)) * (delta*delta)
    LIST$power <- 1 - pchisq(qchisq((1.0 - alpha), df=1),
                             df=1, ncp=LIST$ncp)
    LIST$decision <- character( length(LIST$power) )

    # five possibilities (Table 6 in Saris, Satorra, van der Veld, 2009)
    mi.significant <- ifelse( 1 - pchisq(LIST$mi, df=1) < alpha,
                              TRUE, FALSE )
    high.power <- LIST$power > high.power
    # FIXME: sepc.all or epc??
    #epc.high <- LIST$sepc.all > LIST$delta
    epc.high <- LIST$epc > LIST$delta

    LIST$decision[ which(!mi.significant & !high.power)] <- "(i)"
    LIST$decision[ which( mi.significant & !high.power)] <- "**(m)**"
    LIST$decision[ which(!mi.significant &  high.power)] <- "(nm)"
    LIST$decision[ which( mi.significant &  high.power &
                            !epc.high)] <-  "epc:nm"
    LIST$decision[ which( mi.significant &  high.power &
                            epc.high)] <-  "*epc:m*"

    #LIST$decision[ which(mi.significant &  high.power) ] <- "epc"
    #LIST$decision[ which(mi.significant & !high.power) ] <- "***"
    #LIST$decision[ which(!mi.significant & !high.power) ] <- "(i)"
  }

  # sort?
  if (sort.) {
    LIST <- LIST[order(LIST$mi, decreasing = TRUE),]
  }
  if (minimum.value > 0.0) {
    LIST <- LIST[!is.na(LIST$mi) & LIST$mi > minimum.value,]
  }
  if (maximum.number < nrow(LIST)) {
    LIST <- LIST[seq_len(maximum.number),]
  }
  if (na.remove) {
    idx <- which(is.na(LIST$mi))
    if (length(idx) > 0)  LIST <- LIST[-idx,]
  }
  if (!is.null(op)) {
    idx <- LIST$op %in% op
    if (length(idx) > 0) LIST <- LIST[idx,]
  }

  # add header
  # TODO: small explanation of the columns in the header?
#    attr(LIST, "header") <-
# c("modification indices for newly added parameters only; to\n",
#   "see the effects of releasing equality constraints, use the\n",
#   "lavTestScore() function")

  LIST
}

##' @name modindices.mi-deprecated
##' @usage
##' modificationIndices.mi(object, test = c("D2","D1"),
##'                        omit.imps = c("no.conv","no.se"),
##'                        standardized = TRUE, cov.std = TRUE,
##'                        information = "expected", power = FALSE, delta = 0.1,
##'                        alpha = 0.05, high.power = 0.75, sort. = FALSE,
##'                        minimum.value = 0, maximum.number = nrow(LIST),
##'                        na.remove = TRUE, op = NULL)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL

##' @rdname semTools-deprecated
##' @export
modificationIndices.mi <- modindices.mi


##' @name modindices.mi-deprecated
##' @usage
##' modificationindices.mi(object, test = c("D2","D1"),
##'                        omit.imps = c("no.conv","no.se"),
##'                        standardized = TRUE, cov.std = TRUE,
##'                        information = "expected", power = FALSE, delta = 0.1,
##'                        alpha = 0.05, high.power = 0.75, sort. = FALSE,
##'                        minimum.value = 0, maximum.number = nrow(LIST),
##'                        na.remove = TRUE, op = NULL)
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL

##' @rdname semTools-deprecated
##' @export
modificationindices.mi <- modindices.mi
