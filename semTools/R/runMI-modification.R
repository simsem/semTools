### Terrence D. Jorgensen & Yves rosseel
### Last updated: 2 April 2018
### adaptation of lavaan::modindices() for lavaan.mi-class objects


#' Modification Indices for Multiple Imputations
#'
#' Modification indices (1-\emph{df} Lagrange Multiplier tests) from a
#' latent variable model fitted to multiple imputed data sets. Statistics
#' for releasing one or more fixed or constrained parameters in model can
#' be calculated by pooling the gradient and information matrices pooled
#' across imputed data sets using Rubin's (1987) rules, or by pooling the
#' modification indices across imputed data sets (Li, Meng, Raghunathan, &
#' Rubin, 1991).
#'
#' @aliases modificationIndices.mi modificationindices.mi modindices.mi
#' @importFrom lavaan lavInspect lavListInspect parTable
#' @importFrom methods getMethod
#' @importFrom stats cov pchisq qchisq
#'
#' @param object An object of class \code{\linkS4class{lavaan.mi}}
#' @param type \code{character} indicating which pooling method to use.
#' \code{"Rubin"} indicates Rubin's (1987) rules will be applied to the
#' gradient and information, and those pooled values will be used to
#' calculate modification indices in the usual manner. \code{"D2"} (default),
#' \code{"LMRR"}, or \code{"Li.et.al"} indicate that modification indices
#' calculated from each imputed data set will be pooled across imputations,
#' as described in Li, Meng, Raghunathan, & Rubin (1991) and Enders (2010).
#' @param cov.std \code{logical}. Ignored if \code{type == "D2"}.
#'   If \code{TRUE} (default), the (residual)
#' observed covariances are scaled by the square-root of the diagonal elements
#' of the \eqn{\Theta} matrix, and the (residual) latent covariances are
#' scaled by the square-root of the diagonal elements of the \eqn{\Psi}
#' matrix. If \code{FALSE}, the (residual) observed covariances are scaled by
#' the square-root of the diagonal elements of the model-implied covariance
#' matrix of observed variables (\eqn{\Sigma}), and the (residual) latent
#' covariances are scaled by the square-root of the diagonal elements of the
#' model-implied covariance matrix of the latent variables.
#' @inheritParams lavaan::modificationIndices
#'
#' @note When \code{type = "D2"}, (S)EPCs will be pooled by taking the average
#' across imputations. When \code{type = "Rubin"}, EPCs will be calculated
#' in the standard way using the pooled gradient and information, and SEPCs
#' will be calculated by standardizing the EPCs using model-implied (residual)
#' variances.
#'
#' @return A \code{data.frame} containing modification indices and (S)EPCs.
#'
#' @author
#'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
#'
#'   Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
#'
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
  PT <- parTable(object)
  npar <- max(PT$free) - sum(PT$op == "==")
  useImps <- sapply(object@convergence, "[[", i = "converged")
  m <- sum(useImps)
  type <- tolower(type[1])
  #FIXME: cut useSE?
  # useSE <- sapply(object@convergence, "[[", i = "SE")
  # useSE[is.na(useSE)] <- FALSE

  ## check if model has converged
  if (m == 0L) stop("No models converged. Modification indices unavailable.")

  # not ready for estimator = "PML"
  if (object@Options$estimator == "PML") {
      stop("Modification indices not yet implemented for estimator PML.")
  }

  # sanity check
  if (power) standardized <- TRUE

  ## D2 pooling method
  if (type == "d2") {
    myCols <- c("lhs","op","rhs")
    if (lavListInspect(object, "ngroups") > 1L) myCols <- c(myCols,"block","group")
    LIST <- object@miList[[ which(useImps)[1] ]][myCols]
    nR <- try(nrow(LIST), silent = TRUE)
    if (class(nR) == "try-error" || is.null(nR)) stop("No modification indices were computed.")
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
    ## apply Rubin's rules to pool gradient and augmented information matrix
    oldCall <- object@lavListCall
    #oldCall$model <- parTable(object) #FIXME: necessary?
    lav_object_extend.mi <- function(obj) {
      ## --------------------------------------
      ## borrowed code from lav_object_extend()
      ## --------------------------------------

      # partable original model
      oldPT <- lavaan::parTable(obj)[c("lhs","op","rhs","block","group",
                                       "free","label","plabel")]
      oldPT$user <- rep(1L, length(oldPT$lhs))
      non.free.idx <- which(oldPT$free == 0L & !(oldPT$op %in% c("==",":=","<",">")))
      oldPT$free[ non.free.idx ] <- 1L
      oldPT$user[ non.free.idx ] <- 10L

      # replace 'start' column, since lav_model will fill these in in GLIST
      oldPT$start <- lavaan::parameterEstimates(obj, remove.system.eq = FALSE,
                                                remove.def = FALSE,
                                                remove.eq = FALSE,
                                                remove.ineq = FALSE)$est

      # add new parameters, extend model
      myCols <- c("lhs","op","rhs")
      if (lavaan::lavInspect(obj, "ngroups") > 1L) myCols <- c(myCols,"block","group")
      ADD <- lavaan::modindices(obj, standardized = FALSE)[myCols]
      nR <- try(nrow(ADD), silent = TRUE)
      if (class(nR) == "try-error" || is.null(nR)) return(list(gradient = NULL,
                                                               information = NULL))
      ADD$free <- rep(1L, nR)
      ADD$user <- rep(10L, nR)

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
      ## -------------------------------
      ## borrowed code from modindices()
      ## -------------------------------
      information <- lavaan::lavInspect(obj2, "information")
      LIST <- lavaan::parTable(obj2)
      model.idx <- LIST$free[ LIST$free > 0L & LIST$user != 10L ]
      extra.idx <- LIST$free[ LIST$free > 0L & LIST$user == 10L ]
      # partition
      I11 <- information[extra.idx, extra.idx, drop = FALSE]
      I12 <- information[extra.idx, model.idx, drop = FALSE]
      I21 <- information[model.idx, extra.idx, drop = FALSE]
      I22 <- information[model.idx, model.idx, drop = FALSE]
      I22.inv <- try(lavaan::lavInspect(obj2, "inverted.information"), silent = TRUE)
      # just in case...
      if (inherits(I22.inv, "try-error")) {
        I22.inv <- try(solve(I22), silent = TRUE)
        if (inherits(I22.inv, "try-error")) {
          I22.inv <- MASS::ginv(I22)
          warning("Could not invert I22; modification indices computed using generalized inverse.")
        }
      } else I22.inv <- I22.inv[model.idx, model.idx, drop = FALSE]
      list(gradient = lavaan::lavInspect(obj2, "gradient")[extra.idx],
           information = I11 - I12 %*% I22.inv %*% I21,
           parTable = LIST[LIST$free > 0L & LIST$user == 10L, ])
    }
    oldCall$FUN <- lav_object_extend.mi
    FIT <- eval(as.call(oldCall))
    ## pool gradients and information matrices
    gradList <- lapply(FIT@funList[useImps], "[[", i = "gradient")
    infoList <- lapply(FIT@funList[useImps], "[[", i = "information")
    score <- colMeans(do.call(rbind, gradList))
    B <- cov(do.call(rbind, gradList))
    W <- Reduce("+", infoList) / m
    ## check whether equality constraints prevent inversion of W
    inv.W <- try(solve(W), silent = TRUE)
    if (!inherits(inv.W, "try-error")) {
      ## relative increase in variance due to missing data
      r <- (1 + 1/m)/npar * sum(diag(B %*% inv.W)) # Enders (2010, p. 235) eqs. 8.20-21
      V <- (1 + r) * W
    } else {
      # warning("Could not invert W; information matrices pooled without ARIV")
      ## less reliable, but constraints prevent inversion of W
      V <- W + B + (1/m)*B ## Enders (2010, p. 235) eq. 8.19
    }

    ## template to save output
    myCols <- c("lhs","op","rhs")
    if (lavListInspect(object, "ngroups") > 1L) myCols <- c(myCols,"block","group")
    LIST <- FIT@funList[[ which(useImps)[1] ]]$parTable[myCols]
    ## remove rows with (in)equality constraints and user-defined parameters
    LIST <- LIST[!(LIST$op %in% c("==",":=","<",">")), ]
    rownames(LIST) <- NULL
    nR <- try(nrow(LIST), silent = TRUE)
    if (class(nR) == "try-error" || is.null(nR)) stop("No modification indices were computed.")

    V.diag <- diag(V)
    # dirty hack: catch very small or negative values in diag(V)
    # this is needed eg when parameters are not identified if freed-up;
    idx <- which(V.diag < sqrt(.Machine$double.eps))
    if (length(idx) > 0L) V.diag[idx] <- as.numeric(NA)

    # create and fill in mi
    N <- lavListInspect(object, "ntotal")
    LIST$mi <- N * (score*score) / V.diag

    # handle equality constraints (if any)
    #eq.idx <- which(LIST$op == "==")
    #if(length(eq.idx) > 0L) {
    #    OUT <- lavTestScore(object, warn = FALSE)
    #    LIST$mi[ eq.idx ] <- OUT$uni$X2
    #}

    # scaled?
    #if(length(object@test) > 1L) {
    #    LIST$mi.scaled <- LIST$mi / object@test[[2]]$scaling.factor
    #}

    # EPC
    d <- (-1 * N) * score
    # needed? probably not; just in case
    d[which(abs(d) < 1e-15)] <- 1.0
    LIST$epc <- LIST$mi / d

    # standardize?
    if (standardized) {
      ## Need full parameter table for lavaan::standardizedSolution()
      ## Merge parameterEstimates() with modindices()
      oldPE <- getMethod("summary","lavaan.mi")(object, se = FALSE,
                                                add.attributes = FALSE)
      PE <- lavaan::lav_partable_merge(oldPE, cbind(LIST, est = 0),
                                       remove.duplicated = TRUE, warn = FALSE)
      ## EPCs of user-specified fixed parameters were lost, so replace them.
      ## FIXME? irrelevant if type == "d2" ?  If so, harmless
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
                                                            partable = PE,
                                                            GLIST = object@GLIST,
                                                            est = abs(EPC))$est.std
      PE$sepc.all <- EPC.sign * lavaan::standardizedSolution(object, se = FALSE,
                                                             type = "std.all",
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