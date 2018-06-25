### Terrence D. Jorgensen & Yves rosseel
### Last updated: 25 June 2018
### adaptation of lavaan::modindices() for lavaan.mi-class objects


#' Modification Indices for Multiple Imputations
#'
#' Modification indices (1-\emph{df} Lagrange multiplier tests) from a
#' latent variable model fitted to multiple imputed data sets. Statistics
#' for releasing one or more fixed or constrained parameters in model can
#' be calculated by pooling the gradient and information matrices
#' across imputed data sets using Rubin's (1987) rules, or by pooling the
#' test statistics across imputed data sets (Li, Meng, Raghunathan, &
#' Rubin, 1991).
#'
#' @aliases modificationIndices.mi modificationindices.mi modindices.mi
#' @importFrom lavaan lavInspect lavListInspect
#' @importFrom methods getMethod
#' @importFrom stats cov pchisq qchisq
#'
#' @param object An object of class \code{\linkS4class{lavaan.mi}}
#' @param type \code{character} indicating which pooling method to use.
#'  \code{type = "D2"} (default), \code{"LMRR"}, or \code{"Li.et.al"} indicates
#'  that modification indices that were calculated within each imputed data set
#'  will be pooled across imputations, as described in Li, Meng, Raghunathan,
#'  & Rubin (1991) and Enders (2010).
#'  \code{"Rubin"} indicates Rubin's (1987) rules will be applied to the
#'  gradient and information, and those pooled values will be used to
#'  calculate modification indices in the usual manner.
#' @param standardized \code{logical}. If \code{TRUE}, two extra columns
#'  (\code{$sepc.lv} and \code{$sepc.all}) will contain standardized values for
#'  the EPCs. In the first column (\code{$sepc.lv}), standardizization is based
#'  on the variances of the (continuous) latent variables. In the second column
#'  (\code{$sepc.all}), standardization is based on both the variances of both
#'  (continuous) observed and latent variables. (Residual) covariances are
#'  standardized using (residual) variances.
#' @param cov.std \code{logical}. \code{TRUE} if \code{type == "D2"}.
#'  If \code{TRUE} (default), the (residual)
#'  observed covariances are scaled by the square-root of the diagonal elements
#'  of the \eqn{\Theta} matrix, and the (residual) latent covariances are
#'  scaled by the square-root of the diagonal elements of the \eqn{\Psi}
#'  matrix. If \code{FALSE}, the (residual) observed covariances are scaled by
#'  the square-root of the diagonal elements of the model-implied covariance
#'  matrix of observed variables (\eqn{\Sigma}), and the (residual) latent
#'  covariances are scaled by the square-root of the diagonal elements of the
#'  model-implied covariance matrix of the latent variables.
#' @param power \code{logical}. If \code{TRUE}, the (post-hoc) power is
#'  computed for each modification index, using the values of \code{delta}
#'  and \code{alpha}.
#' @param delta The value of the effect size, as used in the post-hoc power
#'  computation, currently using the unstandardized metric of the \code{$epc}
#'  column.
#' @param alpha The significance level used for deciding if the modification
#'  index is statistically significant or not.
#' @param high.power If the computed power is higher than this cutoff value,
#'  the power is considered 'high'. If not, the power is considered 'low'.
#'  This affects the values in the \code{$decision} column in the output.
#' @param sort. \code{logical}. If \code{TRUE}, sort the output using the
#'  values of the modification index values. Higher values appear first.
#' @param minimum.value \code{numeric}. Filter output and only show rows with a
#'  modification index value equal or higher than this minimum value.
#' @param maximum.number \code{integer}. Filter output and only show the first
#'  maximum number rows. Most useful when combined with the \code{sort.} option.
#' @param na.remove \code{logical}. If \code{TRUE} (default), filter output by
#'  removing all rows with \code{NA} values for the modification indices.
#' @param op \code{character} string. Filter the output by selecting only those
#'  rows with operator \code{op}.
#'
#' @note When \code{type = "D2"}, each (S)EPC will be pooled by taking its
#'  average across imputations. When \code{type = "Rubin"}, EPCs will be
#'  calculated in the standard way using the pooled gradient and information,
#'  and SEPCs will be calculated by standardizing the EPCs using model-implied
#'  (residual) variances.
#'
#' @return A \code{data.frame} containing modification indices and (S)EPCs.
#'
#' @author
#'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
#'
#' Adapted from \pkg{lavaan} source code, written by
#'   Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
#'
#' \code{type = "Rubin"} method proposed by
#'   Maxwell Mansolf (University of California, Los Angeles;
#'   \email{mamansolf@@gmail.com})
#'
#' @references
#' Enders, C. K. (2010). \emph{Applied missing data analysis}.
#' New York, NY: Guilford.
#'
#' Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
#' Significance levels from repeated \emph{p}-values with multiply-imputed data.
#' \emph{Statistica Sinica, 1}(1), 65--92. Retrieved from
#' \url{http://www.jstor.org/stable/24303994}
#'
#' Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}.
#' New York, NY: Wiley.
#'
#' @examples
#'  \dontrun{
#' ## impose missing data for example
#' HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
#'                                       "ageyr","agemo","school")]
#' set.seed(12345)
#' HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
#' age <- HSMiss$ageyr + HSMiss$agemo/12
#' HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
#'
#' ## impute missing data
#' library(Amelia)
#' set.seed(12345)
#' HS.amelia <- amelia(HSMiss, m = 20, noms = "school", p2s = FALSE)
#' imps <- HS.amelia$imputations
#'
#' ## specify CFA model from lavaan's ?cfa help page
#' HS.model <- '
#'   visual  =~ x1 + x2 + x3
#'   textual =~ x4 + x5 + x6
#'   speed   =~ x7 + x8 + x9
#' '
#'
#' out <- cfa.mi(HS.model, data = imps)
#'
#' modindices.mi(out) # default: Li et al.'s (1991) "D2" method
#' modindices.mi(out, type = "Rubin") # Rubin's rules
#'
#' }
#'
#' @export
modindices.mi <- function(object,
                          type = c("D2","Rubin"),

                          standardized = TRUE,
                          cov.std = TRUE,

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
  stopifnot(inherits(object, "lavaan.mi"))
  useSE <- sapply(object@convergence, "[[", i = "SE")
  useSE[is.na(useSE)] <- FALSE
  useImps <- useSE & sapply(object@convergence, "[[", i = "converged")
  m <- sum(useImps)
  type <- tolower(type[1])
  N <- lavListInspect(object, "ntotal")
  #FIXME: if (lavoptions$mimic == "EQS") N <- N - 1 # not in lavaan, why?

  ## check if model has converged
  if (m == 0L) stop("No models converged. Modification indices unavailable.")

  # not ready for estimator = "PML"
  if (object@Options$estimator == "PML") {
      stop("Modification indices not yet implemented for estimator PML.")
  }

  # sanity check
  if (power) standardized <- TRUE

  ## use first available modification indices as template to store pooled results
  myCols <- c("lhs","op","rhs")
  #FIXME: add "level" column?  how to check for multilevel data?
  if (lavListInspect(object, "ngroups") > 1L) myCols <- c(myCols,"block","group")
  for (i in which(useImps)) {
    LIST <- object@miList[[ which(useImps)[i] ]][myCols]
    nR <- try(nrow(LIST), silent = TRUE)
    if (class(nR) == "try-error" || is.null(nR)) {
      if (i == max(which(useImps))) {
        stop("No modification indices could be computed for any imputations.")
      } else next
    } else break
  }



  ## D2 pooling method
  if (type == "d2") {
    chiList <- lapply(object@miList[useImps], "[[", i = "mi")
    ## imputations in columns, parameters in rows
    LIST$mi <- apply(do.call(cbind, chiList), 1, function(x) {
      calculate.D2(x, DF = 1, asymptotic = TRUE)[1]
    })
    ## also take average of epc & sepc.all
    epcList <- lapply(object@miList[useImps], "[[", i = "epc")
    LIST$epc <- rowMeans(do.call(cbind, epcList))
    if (standardized) {
      sepcList <- lapply(object@miList[useImps], "[[", i = "sepc.lv")
      LIST$sepc.lv <- rowMeans(do.call(cbind, sepcList))
      sepcList <- lapply(object@miList[useImps], "[[", i = "sepc.all")
      LIST$sepc.all <- rowMeans(do.call(cbind, sepcList))
    }
  } else {

    scoreOut <- lavTestScore.mi(object, add = cbind(LIST, user = 10L,
                                                    free = 1, start = 0),
                                type = "Rubin", scale.W = FALSE,
                                epc = TRUE, asymptotic = TRUE)$uni
    LIST$mi <- scoreOut$X2
    LIST$epc <- scoreOut$epc

    # standardize?
    if (standardized) {
      ## Need full parameter table for lavaan::standardizedSolution()
      ## Merge parameterEstimates() with modindices()
      oldPE <- getMethod("summary","lavaan.mi")(object, se = FALSE,
                                                add.attributes = FALSE)
      PE <- lavaan::lav_partable_merge(oldPE, cbind(LIST, est = 0),
                                       remove.duplicated = TRUE, warn = FALSE)
      ## merge EPCs, using parameter labels (unavailable for estimates)
      rownames(LIST) <- paste0(LIST$lhs, LIST$op, LIST$rhs, ".g", LIST$group)
      rownames(PE) <- paste0(PE$lhs, PE$op, PE$rhs, ".g", PE$group)
      PE[rownames(LIST), "epc"] <- LIST$epc
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

      PE$sepc.lv <- EPC.sign * lavaan::standardizedSolution(object, se = FALSE,
                                                            type = "std.lv",
                                                            cov.std = cov.std,
                                                            partable = PE,
                                                            GLIST = object@GLIST,
                                                            est = abs(EPC))$est.std
      PE$sepc.all <- EPC.sign * lavaan::standardizedSolution(object, se = FALSE,
                                                             type = "std.all",
                                                             cov.std = cov.std,
                                                             partable = PE,
                                                             GLIST = object@GLIST,
                                                             est = abs(EPC))$est.std
      if (length(small.idx) > 0L) {
        PE$sepc.lv[small.idx] <- 0
        PE$sepc.all[small.idx] <- 0
      }
      ## remove unnecessary columns, then merge
      if (is.null(LIST$block)) PE$block <- NULL
      PE$est <- NULL
      PE$mi <- NULL
      PE$epc <- NULL
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

# aliases
modificationIndices.mi <- modificationindices.mi <- modindices.mi
