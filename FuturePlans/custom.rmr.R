### Terrence D. Jorgensen
### Last updated: 9 September 2019
### functions to propose new residuals-based indices of local misfit


## ------------------------
## Custom (S/C)RMR function
## ------------------------


##' Find Indices to Calculate Custom (S/C)RMRs
##'
##' The root mean-square residual (RMR), and its standardized counterparts
##' (e.g., the CRMR or SRMR; Maydeu-Olivares, 2017), have long been evaluated
##' using all residuals to summarize misfit of the entire model. However,
##' summaries can be calculated for mean and covariance structures to evaluate
##' each aspect of the model separately. More focused summaries can be used
##' to evaluate the local fit of particulat aspects of the model (e.g., a
##' separate (S/C)RMR for each factor in models with simple structure) or
##' to evaluate the practical equivalence of freeing a particular parameter
##' (e.g., to accompany modification indices and (S)EPCs). This function
##' facilitates automatically generating the list of indices necessary to a few
##' types of custom summaries.
##'
##' @importFrom lavaan lavInspect lavNames
##'
##' @param object An object of class \code{\linkS4class{lavaan}}
##' @param rmr.type \code{character} vector indicating the desired type(s) of
##'   summary. Can contain any of c("cross-loadings","latent.paths","per.factor")
##' @param resid.type \code{character} requesting a type of standardized
##'   residual (or \code{"raw"}). Passed to
##'   \code{\link[lavaan]{lavResiduals}(object, type = resid.type)}
##'
##' @return A named \code{list} containing summaries, one per (S/C)RMR. For
##'   each summary, there is a list (often only \code{$cov}) of indices locating
##'   which sample statistics' residuals (as returned by
##'   \code{\link[lavaan]{lavResiduals}}) should be included in the (S/C)RMR.
##'   This \code{list} is designed be passed to the \code{custom.rmr} argument
##'   of \code{\link[lavaan]{lavResiduals}}.
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Maydeu-Olivares, A. (2017). Assessing the size of model misfit in
##'   structural equation models. \emph{Psychometrika, 82}(3), 533--558.
##'   doi:10.1007/s11336-016-9552-7
##'
##' @seealso \code{\link[lavaan]{lavResiduals}}
##'
##' @examples
################ COMING
## example from pull request 146?  https://github.com/yrosseel/lavaan/pull/146
customRMRs <- function(object, rmr.type = "", resid.type = "cor.bentler") {
  ## both classes include @Model information (includes blavaan and lavaan.mi)
  stopifnot(inherits(object, c("lavaan","lavaanList")))
  lavoptions <- lavInspect(object, "options")


  rmr.type <- tolower(rmr.type)
  ## BORROWED FROM lav_residuals.R:  check type
  resid.type <- tolower(resid.type)[1]
  if (!resid.type %in% c("raw", "rmr",
                         "crmr", "cor", "cor.bollen",
                         "srmr", "cor.bentler", "cor.eqs",
                         "normalized", "standardized", "standardized.mplus")) {
    stop("lavaan ERROR: unknown argument for resid.type: ", dQuote(resid.type))
  }
  # if cor, choose 'default'
  if (resid.type == "cor") {
    if (lavoptions$mimic == "EQS") {
      resid.type <- "cor.bentler"
    } else {
      resid.type <- "cor.bollen"
    }
  }
  if (resid.type == "cor.eqs") resid.type <- "cor.bentler"
  if (resid.type == "rmr")     resid.type <- "raw"
  if (resid.type == "srmr")    resid.type <- "cor.bentler"
  if (resid.type == "crmr")    resid.type <- "cor.bollen"


  ## Retrieve basic model information
  meanstructure <- lavoptions$meanstructure
  categorical <- lavInspect(object, "categorical")
  ov.names <- lavNames(fit, "ov.model")
  lv.names <- lavNames(fit, "lv.regular")
  FREE <- lavInspect(fit, "free")
  EST <- lavInspect(fit, "est")
  if (lavInspect(fit, "ngroups") > 1) {
    FREE <- FREE[[1]]
    EST <- EST[[1]]
  }
  fix0lambda <- FREE$lambda == 0 & EST$lambda == 0
  fix0psi <- FREE$psi == 0 & EST$psi == 0
  if (!is.null(FREE$beta)) fix0psi <- fix0psi & FREE$beta == 0 & EST$beta == 0
  diag(fix0psi) <- FALSE # diagonal never relevant


  ## start with empty list
  customList <- list()
  ## returned if !any(rmr.type == c("cross-loadings","latent.cor","per.factor"))


  if (any(c("cross-loadings","cross.loadings","cross-loading","cross.loading",
            "cross.load","xload","cl","xl","lambda") %in% rmr.type)) {
    lambda.idx <- which(fix0lambda, arr.ind = TRUE)
    for (RR in seq_len(nrow(lambda.idx))) {
      ## indicator ov could potentially cross-load on lv
      ov <- rownames(fix0lambda)[ lambda.idx[RR, "row"] ]
      lv <- colnames(fix0lambda)[ lambda.idx[RR, "col"] ]
      if (!lv %in% lv.names) next

      ## indicators of lv
      ind <- names(which(!fix0lambda[ , lv]))
      if (length(ind) < 2L) next

      temp <- matrix(FALSE, nrow = length(ov.names), ncol = length(ov.names),
                     dimnames = list(ov.names, ov.names))
      temp[ov, ind] <- TRUE
      customList[[paste0(lv, "=~", ov)]]$cov <- which(temp, arr.ind = TRUE)
      rm(temp)
    }
  }


  if (any(c("latent.cor","factor.cor","latent.path","latent.paths","psi") %in% rmr.type)) {
    psi.idx <- which(fix0psi, arr.ind = TRUE)
    for (RR in seq_len(nrow(psi.idx))) {
      ## factors with no (un)directed paths
      lv1 <- rownames(fix0psi)[ psi.idx[RR, "row"] ]
      lv2 <- colnames(fix0psi)[ psi.idx[RR, "col"] ]
      if (!all(c(lv1, lv2) %in% lv.names)) next

      ## indicators of each factor
      ov1 <- names(which(!fix0lambda[ , lv1]))
      ov2 <- names(which(!fix0lambda[ , lv2]))
      ## remove shared indicators
      ov1 <- setdiff(ov1, ov2)
      ov2 <- setdiff(ov2, ov1)

      temp <- matrix(FALSE, nrow = length(ov.names), ncol = length(ov.names),
                     dimnames = list(ov.names, ov.names))
      temp[ov1, ov2] <- TRUE
      customList[[paste0(lv1, "~~", lv2)]]$cov <- which(temp, arr.ind = TRUE)
      rm(temp)
    }
  }


  if (any(c("factor","per.factor","per-factor","by.factor","by-factor") %in% rmr.type)) {
    for (lv in lv.names) {
      ## this factor's indicators
      ov <- names(which(!fix0lambda[ , lv]))
      if (length(ov) < 2L) next

      temp <- matrix(FALSE, nrow = length(ov.names), ncol = length(ov.names),
                     dimnames = list(ov.names, ov.names))
      temp[ov, ov] <- TRUE
      customList[[lv]]$cov <- which(temp, arr.ind = TRUE)
      rm(temp)

      if (meanstructure) {
        temp <- rep(FALSE, length(ov.names))
        names(temp) <- ov.names
        temp[ov] <- TRUE
        customList[[lv]]$mean <- which(temp)
        rm(temp)
      }

      if (categorical && type == "raw") {
        #TODO: add thresholds
        rm(temp)
      }

    }
  }



  #TODO: when lavResiduals() allows separate by block, write routines for
  #      common equality constraints (e.g., group invariance, add arguments
  #      from measEq.syntax for long.equal, eventually level.equal)
  #TODO: routine for LGCMs (fixed but nonzero slope loadings)

  customList
}



## --------------------------------------------------------
## Function to create a new "fitted" object implied by EPCs
## --------------------------------------------------------


##' Update a lavaan object's parameters using EPCs
##'
##' Expected parameter changes (EPCs) are a useful by-product of modification
##' indices, which are the expected chi-squared difference if a specified (set
##' of) constraint(s) were released. Based on a constrained model, EPCs estimate
##' the degree to which any estimated model parameter would change if a
##' specified (set of) constraint(s) were released. EPCs are added to current
##' estimates (and fixed value(s), if that is the constraint being tested) to
##' create expected parameter values. Those values are used to update the lavaan
##' object, so that expected local and global fit can be calculated based on
##' those EPCs.
##'
##' @importFrom lavaan parTable lavInspect lavaan
##'
##' @param object An object of class \code{\linkS4class{lavaan}}
##' @param ... Arguments passed to \code{\link[lavaan]{lavTestScore}}
##'
##' @return An object of  class \code{\linkS4class{lavaan}}. \bold{NOT} meant
##'   to be used to draw inferences, simply to estimate projected improvements
##'   in various aspects of global/local model fit.
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
## @examples
################ I think this should remain a hidden function
expectedFit <- function(object, ...) {
  ## score tests available for both classes
  stopifnot(inherits(object, "lavaan"))

  mc <- match.call(expand.dots = TRUE)
  add <- if (is.null(mc$add)) NULL else list(...)$add
  release <- if (is.null(mc$release)) NULL else list(...)$release
  mc[[1]] <- lavaan::lavTestScore
  mc$epc <- TRUE
  mc$univariate <- FALSE
  mc$cumulative <- FALSE
  mc$standardized <- FALSE
  if (is.null(mc$warn)) mc$warn <- FALSE
  scoreOut <- eval(mc)
  ## use EPC table as new model
  model <- scoreOut$epc
  ## set free $epv (or fixed $est) as $start
  model$start <- ifelse(is.na(model$epv), model$est, model$epv)
  model$est <- NULL
  model$epc <- NULL
  model$epv <- NULL


  ## rbind() any constraints from parTable() not in release=
  ## (same labels don't enforce constraints without those rows)
  PT <- parTable(object)
  PT <- PT[ PT$op %in% c("==", ":=", "<", ">") ,
            intersect(names(model), names(PT)) ]
  ## release mode?
  if (is.null(add) || nchar(add) == 0L) {
    if (is.null(release)) release <- 1:sum(PT$op == "==")
    PT <- PT[ which(PT$op == "==")[-release] , ]
  } else {
    #FIXME: make added parameter $free (to change degrees of freedom)?
  }
  ## any remaining constraints (or user-defined parameters)?
  if (nrow(PT)) model <- rbind(model, PT)


  ## pass new model as model without optimizing
  lavoptions1 <- lavInspect(object, "options")
  lavoptions1$verbose <- FALSE
  # lavoptions1$se <- "none"
  lavoptions1$optim.method <- "none"
  # lavmodel1 <- lavaan::lav_model_set_parameters(object@Model, x = EPV)
  # if ("control" %in% slotNames(lavmodel1)) {
  #   lavmodel1@control <- list(optim.method = "none")
  # }

  ## return newly "fitted" model
  lavaan(model, slotOptions = lavoptions1, #slotModel = lavmodel1, # updated slots
         slotSampleStats = object@SampleStats, slotData = object@Data,
         slotCache = object@Cache)
}



## -----------------------------------------------------------------
## Function to calculate expected change in (S/C)RMR implied by EPCs
## -----------------------------------------------------------------


##' Calculate expected change in (S/C)RMR implied by EPCs
##'
##' DESCRPITION
##'
##' @importFrom lavaan lavResiduals
##'
##' @param object An object of class \code{\linkS4class{lavaan}}
##    TODO: adapt for lavaan.mi
##' @param type \code{character} requesting a type of residual. Passed to
##'   \code{\link[lavaan]{lavResiduals}}.
##' @param unbiased \code{logical} indicating whether to use the standard
##'   (default) or unbiased estimate of a (S/C)RMR (see Maydeu-Olivares, 2017).
##' @param ... Arguments passed to \code{\link[lavaan]{modindices}}
##'
##' @return \code{data.frame} returned by \code{\link[lavaan]{modindices}},
##'   appended with expected changes in the requested \code{type} of RMR if
##'   each parameter were added to the model.
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Maydeu-Olivares, A. (2017). Assessing the size of model misfit in
##'   structural equation models. \emph{Psychometrika, 82}(3), 533--558.
##'   doi:10.1007/s11336-016-9552-7
##'
##' @seealso \code{\link[lavaan]{lavResiduals}}, \code{\link[lavaan]{modindices}}
##'
##' @examples
################ COMING
ec.rmr <- function(object, type = "cor.bentler", unbiased = FALSE, ...) {
  stopifnot(inherits(object, "lavaan"))

  idx <- if (unbiased) 5L else 1L

  ## grab (S/C)RMR using type=
  oldsum <- lavResiduals(object, #custom.rmr = custom.rmr,
                         type = type)$summary
  oldRMR <- oldsum[nrow(oldsum), idx]

  ## save modindices() output
  MI <- lavaan::modindices(object, ...)

  ## extract columns to pass to lavTestScore(object, add = )
  keepnames <- c("lhs","op","rhs")
  if ("block" %in% names(MI)) keepnames <- c(keepnames, "block")
  if ("group" %in% names(MI)) keepnames <- c(keepnames, "group")
  if ("level" %in% names(MI)) keepnames <- c(keepnames, "level")
  ADD <- MI[ , keepnames]
  ADD$user <- 10L
  ADD$free <- 1
  ADD$start <- 0

  ## Loop over rows to pass as add=
  for (pp in 1:nrow(ADD)) {
    efArgs <- list(object = object, add = ADD[pp, ])
    newfit <- do.call(expectedFit, efArgs)
    ## new (S/C)RMR
    newsum <- lavResiduals(newfit, #custom.rmr = custom.rmr,
                           type = type)$summary
    newRMR <- newsum[nrow(newsum), idx]

    ## append table with EC.RMR column
    MI[pp, paste0("ec.", colnames(newsum)[idx]) ] <- newRMR - oldRMR
  }

  #class(MI) <- c("lavaan.data.frame", "data.frame")
  MI
}



## --------------------
## Frequentist examples
## --------------------

library(lavaan)

## example from ?lavTestScore, with additional b4 label on each 3rd indicator
HS.model <- '
    visual  =~ x1 + b1*x2 + b4*x3
    textual =~ x4 + b2*x5 + b4*x6
    speed   =~ x7 + b3*x8 + b4*x9

    b1 == b2
    b2 == b3
  # visual ~ ageyr       # uncomment to check with fixed.x
'
fit <- cfa(HS.model, data = HolzingerSwineford1939)


## LOCAL RMRs
customList <- customRMRs(fit, rmr.type = c("xload","per.factor"))
lavResiduals(fit, custom.rmr = customList)
## compare to:
modindices(fit, op = "=~")


## Expected change in GLOBAL RMR
ec.rmr(fit) # release all (fit2) # srmr
ec.rmr(fit, type = "cor.bollen") # crmr
ec.rmr(fit, type = "raw")        # unstandardized



## -----------------
## Bayesian examples
## -----------------


library(blavaan)

## example from ?lavTestScore, with additional b4 label on each 3rd indicator
HS.model <- '
    visual  =~ x1 + b1*x2 + b4*x3
    textual =~ x4 + b1*x5 + b4*x6
    speed   =~ x7 + b1*x8 + b4*x9

  # b1 == b2
  # b2 == b3
  # visual ~ ageyr       # uncomment to check with fixed.x
'
bfit <- bcfa(HS.model, data = HolzingerSwineford1939, cp = "fa",
             #target = "jags", bcontrol = list(method = "rjparallel"),
             n.chains = 2, burnin = 1000, sample = 500)



## wrap 3*2*2 = 12 indices into discFUN for ppmc()
## - Raw, Correlation, and Standardized RMRs
## - standard and unbiased estimates
## - local RMRs and expected change in global RMRs
discFUN <- list(local.rmr = function(object) {
                  customList <- customRMRs(object, rmr.type = "xload")
                  summ <- lavResiduals(object, type = "raw",
                                       custom.rmr = customList)$summary
                  rmr <- summ$rmr
                  names(rmr) <- rownames(summ)
                  rmr
                },
                local.urmr = function(object) {
                  customList <- customRMRs(object, rmr.type = "xload")
                  summ <- lavResiduals(object, type = "raw",
                                       custom.rmr = customList)$summary
                  urmr <- summ$urmr
                  names(urmr) <- rownames(summ)
                  urmr
                },
                local.crmr = function(object) {
                  customList <- customRMRs(object, rmr.type = "xload")
                  summ <- lavResiduals(object, type = "crmr",
                                       custom.rmr = customList)$summary
                  crmr <- summ$crmr
                  names(crmr) <- rownames(summ)
                  crmr
                },
                local.ucrmr = function(object) {
                  customList <- customRMRs(object, rmr.type = "xload")
                  summ <- lavResiduals(object, type = "crmr",
                                       custom.rmr = customList)$summary
                  ucrmr <- summ$ucrmr
                  names(ucrmr) <- rownames(summ)
                  ucrmr
                },
                local.srmr = function(object) {
                  customList <- customRMRs(object, rmr.type = "xload")
                  summ <- lavResiduals(object, type = "srmr",
                                       custom.rmr = customList)$summary
                  srmr <- summ$srmr
                  names(srmr) <- rownames(summ)
                  srmr
                },
                local.usrmr = function(object) {
                  customList <- customRMRs(object, rmr.type = "xload")
                  summ <- lavResiduals(object, type = "srmr",
                                       custom.rmr = customList)$summary
                  usrmr <- summ$usrmr
                  names(usrmr) <- rownames(summ)
                  usrmr
                },
                global.rmr  = function(object) {
                  MI <- ec.rmr(object, op = "=~", type = "raw")
                  out <- MI$ec # partial matching (always starts with "ec.")
                  names(out) <- paste0(MI$lhs, MI$op, MI$rhs)
                  out
                },
                global.urmr  = function(object) {
                  MI <- ec.rmr(object, op = "=~", type = "raw", unbiased = TRUE)
                  out <- MI$ec # partial matching (always starts with "ec.")
                  names(out) <- paste0(MI$lhs, MI$op, MI$rhs)
                  out
                },
                global.crmr  = function(object) {
                  MI <- ec.rmr(object, op = "=~", type = "crmr")
                  out <- MI$ec # partial matching (always starts with "ec.")
                  names(out) <- paste0(MI$lhs, MI$op, MI$rhs)
                  out
                },
                global.ucrmr  = function(object) {
                  MI <- ec.rmr(object, op = "=~", type = "crmr", unbiased = TRUE)
                  out <- MI$ec # partial matching (always starts with "ec.")
                  names(out) <- paste0(MI$lhs, MI$op, MI$rhs)
                  out
                },
                global.srmr = function(object) {
                  MI <- ec.rmr(object, op = "=~", type = "srmr")
                  out <- MI$ec # partial matching (always starts with "ec.")
                  names(out) <- paste0(MI$lhs, MI$op, MI$rhs)
                  out
                },
                global.usrmr = function(object) {
                  MI <- ec.rmr(object, op = "=~", type = "srmr", unbiased = TRUE)
                  out <- MI$ec # partial matching (always starts with "ec.")
                  names(out) <- paste0(MI$lhs, MI$op, MI$rhs)
                  out
                })

## PPMC for each
out <- ppmc(bfit, thin = 10, discFUN = discFUN)
## THIS TAKES A LONG TIME. Try commenting-out some options, or thin=50 for examples
# saveRDS(out, file = "FuturePlans/ppmcOutput.rds")   # Save PPMC output
# readRDS("FuturePlans/ppmcOutput.rds")               # Load PPMC output

summary(out, discFUN = "local.crmr")
summary(out, discFUN = "global.crmr")


summary(out, discFUN = "local.usrmr")
summary(out, discFUN = "global.usrmr")

