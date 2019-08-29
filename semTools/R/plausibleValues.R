### Terrence D. Jorgensen
### Last updated: 29 August 2019
### function to draw plausible values of factor scores from lavPredict


# library(blavaan)
# bfit <- bcfa(HS.model, data=HolzingerSwineford1939, save.lvs = TRUE,
#              bcontrol=list(method="rjparallel"), group = "school",
#              #target = "stan", control=list(cores = 4, seed = 123),
#              burnin = 4000, sample = 30, n.chains = 2)
# bFS <- do.call(rbind, blavInspect(bfit, "lvs"))
# do.call()


## -------------
## Main function
## -------------

##' Plausible-Values Imputation of Factor Scores Estimated from a lavaan Model
##'
##' Draw plausible values of factor scores estimated from a fitted
##' \code{\link[lavaan]{lavaan}} model, then treat them as multiple imputations
##' of missing data using \code{\link{runMI}}.
##'
##'
##' Because latent variables are unobserved, they can be considered as missing
##' data, which can be imputed using Monte Carlo methods.  This may be of
##' interest to researchers with sample sizes too small to fit their complex
##' structural models.  Fitting a factor model as a first step,
##' \code{\link[lavaan]{lavPredict}} provides factor-score estimates, which can
##' be treated as observed values in a path analysis (Step 2).  However, the
##' resulting standard errors and test statistics could not be trusted because
##' the Step-2 analysis would not take into account the uncertainty about the
##' estimated factor scores.  Using the asymptotic sampling covariance matrix
##' of the factor scores provided by \code{\link[lavaan]{lavPredict}},
##' \code{plausibleValues} draws a set of \code{nDraws} imputations from the
##' sampling distribution of each factor score, returning a list of data sets
##' that can be treated like multiple imputations of incomplete data.  If the
##' data were already imputed to handle missing data, \code{plausibleValues}
##' also accepts an object of class \code{\linkS4class{lavaan.mi}}, and will
##' draw \code{nDraws} plausible values from each imputation.  Step 2 would
##' then take into account uncertainty about both missing values and factor
##' scores.  Bayesian methods can also be used to generate factor scores, as
##' available with the \pkg{blavaan} package, in which case plausible
##' values are simply saved parameters from the posterior distribution. See
##' Asparouhov and Muthen (2010) for further technical details and references.
##'
##' Each returned \code{data.frame} includes a \code{case.idx} column that
##' indicates the corresponding rows in the data set to which the model was
##' originally fitted (unless the user requests only Level-2 variables).  This
##' can be used to merge the plausible values with the original observed data,
##' but users should note that including any new variables in a Step-2 model
##' might not accurately account for their relationship(s) with factor scores
##' because they were not accounted for in the Step-1 model from which factor
##' scores were estimated.
##'
##' If \code{object} is a multilevel \code{lavaan} model, users can request
##' plausible values for latent variables at particular levels of analysis by
##' setting the \code{\link[lavaan]{lavPredict}} argument \code{level=1} or
##' \code{level=2}.  If the \code{level} argument is not passed via \dots,
##' then both levels are returned in a single merged data set per draw.  For
##' multilevel models, each returned \code{data.frame} also includes a column
##' indicating to which cluster each row belongs (unless the user requests only
##' Level-2 variables).
##'
##'
##' @importFrom lavaan lavInspect lavPredict
##'
##' @param object A fitted model of class \code{\linkS4class{lavaan}},
##'   \code{\link[blavaan]{blavaan}}, or \code{\linkS4class{lavaan.mi}}
##' @param nDraws \code{integer} specifying the number of draws, analogous to
##'   the number of imputed data sets. If \code{object} is of class
##'   \code{\linkS4class{lavaan.mi}}, this will be the number of draws taken
##'   \emph{per imputation}.  Ignored if \code{object} is of class
##'   \code{\link[blavaan]{blavaan}}, in which case the number of draws is the
##'   number of MCMC samples from the posterior.
##' @param seed \code{integer} passed to \code{\link{set.seed}()}.  Ignored if
##'   \code{object} is of class \code{\link[blavaan]{blavaan}},
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'   imputations when \code{object} is of class \code{\linkS4class{lavaan.mi}}.
##'   Can include any of \code{c("no.conv", "no.se", "no.npd")}.
##' @param ... Optional arguments to pass to \code{\link[lavaan]{lavPredict}}.
##'   \code{assemble} will be ignored because multiple groups are always
##'   assembled into a single \code{data.frame} per draw. \code{type} will be
##'   ignored because it is set internally to \code{type="lv"}.
##'
##' @return A \code{list} of length \code{nDraws}, each of which is a
##'   \code{data.frame} containing plausible values, which can be treated as
##'   a \code{list} of imputed data sets to be passed to \code{\link{runMI}}
##'   (see \bold{Examples}). If \code{object} is of class
##'   \code{\linkS4class{lavaan.mi}}, the \code{list} will be of length
##'   \code{nDraws*m}, where \code{m} is the number of imputations.
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Asparouhov, T. & Muthen, B. O. (2010). \emph{Plausible values for latent
##'   variables using M}plus. Technical Report. Retrieved from
##'   www.statmodel.com/download/Plausible.pdf
##'
##' @seealso \code{\link{runMI}}, \code{\linkS4class{lavaan.mi}}
##'
##' @examples
##'
##' ## example from ?cfa and ?lavPredict help pages
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ x7 + x8 + x9 '
##'
##' fit1 <- cfa(HS.model, data = HolzingerSwineford1939)
##' fs1 <- plausibleValues(fit1, nDraws = 3,
##'                        ## lavPredict() can add only the modeled data
##'                        append.data = TRUE)
##' lapply(fs1, head)
##'
##' ## To merge factor scores to original data.frame (not just modeled data)
##' fs1 <- plausibleValues(fit1, nDraws = 3)
##' idx <- lavInspect(fit1, "case.idx")      # row index for each case
##' if (is.list(idx)) idx <- do.call(c, idx) # for multigroup models
##' data(HolzingerSwineford1939)             # copy data to workspace
##' HolzingerSwineford1939$case.idx <- idx   # add row index as variable
##' ## loop over draws to merge original data with factor scores
##' for (i in seq_along(fs1)) {
##'   fs1[[i]] <- merge(fs1[[i]], HolzingerSwineford1939, by = "case.idx")
##' }
##' lapply(fs1, head)
##'
##'
##' ## multiple-group analysis, in 2 steps
##' step1 <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
##'             group.equal = c("loadings","intercepts"))
##' PV.list <- plausibleValues(step1)
##'
##' ## subsequent path analysis
##' path.model <- ' visual ~ c(t1, t2)*textual + c(s1, s2)*speed '
##' \dontrun{
##' step2 <- sem.mi(path.model, data = PV.list, group = "school")
##' ## test equivalence of both slopes across groups
##' lavTestWald.mi(step2, constraints = 't1 == t2 ; s1 == s2')
##' }
##'
##'
##' ## multilevel example from ?Demo.twolevel help page
##' model <- '
##'   level: 1
##'     fw =~ y1 + y2 + y3
##'     fw ~ x1 + x2 + x3
##'   level: 2
##'     fb =~ y1 + y2 + y3
##'     fb ~ w1 + w2
##' '
##' msem <- sem(model, data = Demo.twolevel, cluster = "cluster")
##' mlPVs <- plausibleValues(msem, nDraws = 3) # both levels by default
##' lapply(mlPVs, head, n = 10)
##' ## only Level 1
##' mlPV1 <- plausibleValues(msem, nDraws = 3, level = 1)
##' lapply(mlPV1, head)
##' ## only Level 2
##' mlPV2 <- plausibleValues(msem, nDraws = 3, level = 2)
##' lapply(mlPV2, head)
##'
##' @export
plausibleValues <- function(object, nDraws = 20L, seed = 12345,
                            omit.imps = c("no.conv","no.se"), ...) {

  if (class(object) == "lavaan") {
    ## generate vector of seeds
    set.seed(seed)
    seeds <- sample(100000:9999999, size = nDraws, replace = FALSE)

    PV <- lapply(seeds, plaus.lavaan, object = object, ...)

  } else if (class(object) == "lavaan.mi") {
    ## generate vector of seeds
    set.seed(seed)
    seeds <- sample(100000:9999999, size = nDraws, replace = FALSE)

    PV <- plaus.mi(object, seeds = seeds, omit.imps = omit.imps, ...)

  } else if (class(object) == "blavaan") {
    ## requireNamespace("blavaan")
    ## blavaan::blavInspect(object, "lvs")
    PV <- plaus.blavaan(object)

  } else stop("object's class not valid: ", class(object))

  PV
}



## ----------------
## Hidden functions
## ----------------

## draw 1 set of plausible values from a lavaan object
##' @importFrom lavaan lavInspect lavPredict lavNames
plaus.lavaan <- function(seed = 1, object, ...) {
  stopifnot(inherits(object, "lavaan"))
  if (lavInspect(object, "categorical")) {
    stop("Plausible values not available (yet) for categorical data")
  }
  if (lavInspect(object, "options")$missing %in% c("ml", "ml.x")) {
    stop("Plausible values not available (yet) for missing data + fiml.\n",
         "       Multiple imputations can be used via lavaan.mi()")
  }
  #FIXME?  https://github.com/yrosseel/lavaan/issues/156

  set.seed(seed)
  cluster <- lavInspect(object, "cluster")
  group <- lavInspect(object, "group")
  group.label <- lavInspect(object, "group.label")
  nG <- lavInspect(object, "ngroups")
  nL <- lavInspect(object, "nlevels")
  l.names <- o.names <- list()
  if (nG == 1L && nL == 1L) {
    ## single block
    l.names <- list(lavNames(object, "lv"))
    o.names <- list(lavNames(object, "ov"))
  } else if (nG == 1L && nL > 1L) {
    ## multilevel
    for (BB in 1:nL) {
      l.names <- c(l.names, list(lavNames(object, "lv", block = BB)))
      o.names <- c(o.names, list(lavNames(object, "ov", block = BB)))
    }
  } else if (nG > 1L && nL == 1L) {
    ## multigroup
    for (BB in 1:nG) {
      l.names <- c(l.names, list(lavNames(object, "lv", block = BB)))
      o.names <- c(o.names, list(lavNames(object, "ov", block = BB)))
    }
  } else {
    ## multilevel + multigroup
    for (BB in 1:(nG*nL)) { #FIXME: lavInspect(object, "nblocks")
      l.names <- c(l.names, list(lavNames(object, "lv", block = BB)))
      o.names <- c(o.names, list(lavNames(object, "ov", block = BB)))
    }
  }



  ## extract factor scores + covariance matrix
  fsArgs <- list(...)
  fsArgs$type <- "lv"
  fsArgs$assemble <- FALSE # assemble after drawing
  append.data <- fsArgs$append.data
  if (is.null(append.data)) append.data <- FALSE # default in lavPredict()
  only.L2 <- fsArgs$level == 2L
  if (length(only.L2) == 0L) only.L2 <- FALSE
  if (only.L2) fsArgs$append.data <- append.data <- FALSE #FIXME: how will Yves handle lavPredict(fit, append=T, level=2)?
  bothLevels <- nL > 1L && is.null(fsArgs$level)
  fsArgs$object <- object
  fsArgs$acov <- "standard" #FIXME: update if other options become available
  FS <- do.call(lavPredict, fsArgs) #FIXME: breaks when multigroup MLSEM: https://github.com/yrosseel/lavaan/issues/157
  ## also draw Level 2, if multilevel and no specific level requested
  if (bothLevels) {
    fsArgs$level <- 2L
    fsArgs$append.data <- FALSE #FIXME: how will Yves handle lavPredict(fit, append=T, level=2)?
    FS2 <- do.call(lavPredict, fsArgs)
  }

  ## draw plausible values, if factor scores exist
  if (nG == 1L) {

    if (ncol(FS) == 0L) {
      PV <- FS
    } else {
      ACOV <- attr(FS, "acov")[[1]]
      v.idx <- if (only.L2) 2L else 1L
      PV <- apply(FS[ , l.names[[v.idx]], drop = FALSE], 1, function(mu) {
        MASS::mvrnorm(n = 1, mu = mu, Sigma = ACOV)
      })
      if (is.null(dim(PV))) {
        PV <- as.matrix(PV)
        colnames(PV) <- l.names[[v.idx]]
      } else PV <- t(PV)
      if (append.data) {
        PV <- cbind(FS[ , o.names[[v.idx]], drop = FALSE], PV)
      }
    }

    ## add Level 2 if multilevel and no specific level requested
    if (bothLevels) {
      if (ncol(FS2) == 0L) {
        PV2 <- FS2
      } else {
        ACOV2 <- attr(FS2, "acov")[[1]]
        #FIXME: how will Yves handle lavPredict(fit, append=T, level=2)?
        PV2 <- apply(FS2, 1, function(mu) {
          out <- MASS::mvrnorm(n = 1, mu = mu, Sigma = ACOV2)
        })
        if (is.null(dim(PV2))) {
          PV2 <- as.matrix(PV2)
          colnames(PV2) <- l.names[[2]]
        } else PV2 <- t(PV2)
      }
    }

  } else {
    ACOV <- list()
    PV <- list()
    for (gg in 1:nG) {

      if (ncol(FS[[gg]]) == 0L) {
        PV[[gg]] <- FS[[gg]]
      } else {
        ACOV[[gg]] <- attr(FS, "acov")[[gg]]
        v.idx <- if (only.L2) (2L + (gg - 1L)*nL) else (1L + (gg - 1L)*nL)
        PV[[gg]] <- apply(FS[[gg]][ , l.names[[v.idx]], drop = FALSE], 1, function(mu) {
          MASS::mvrnorm(n = 1, mu = mu, Sigma = ACOV[[gg]])
        })
        if (is.null(dim(PV[[gg]]))) {
          PV[[gg]] <- as.matrix(PV[[gg]])
          colnames(PV[[gg]]) <- l.names[[v.idx]]
        } else PV[[gg]] <- t(PV[[gg]])
      }
      if (append.data) {
        PV[[gg]] <- cbind(FS[[gg]][ , o.names[[v.idx]], drop = FALSE], PV[[gg]])
      }
    }

    ## add Level 2 if multilevel and no specific level requested
    if (bothLevels) {
      ACOV2 <- list()
      PV2 <- list()
      for (gg in 1:nG) {

        if (ncol(FS2[[gg]]) == 0L) {
          PV2[[gg]] <- FS2[[gg]]
        } else {
          ACOV2[[gg]] <- attr(FS2, "acov")[[gg]]
          #FIXME: how will Yves handle lavPredict(fit, append=T, level=2)?
          PV2[[gg]] <- apply(FS2[[gg]], 1, function(mu) {
            MASS::mvrnorm(n = 1, mu = mu, Sigma = ACOV2[[gg]])
          })
          if (is.null(dim(PV2[[gg]]))) {
            PV2[[gg]] <- as.matrix(PV2[[gg]])
            colnames(PV2[[gg]]) <- colnames(FS2[[gg]])
          } else PV2[[gg]] <- t(PV2[[gg]])
        }

      }
    }

  }

  ## save as data.frame
  if (nG > 1L) {

    temp <- lapply(1:nG, function(gg) {
      dd <- data.frame(PV[[gg]])
      ## add groups if multiple
      dd[ , group] <- group.label[gg]
      dd <- dd[ , c(group, setdiff(names(dd), group)), drop = FALSE ]
      ## attach row indices from original data for optional merging
      if (only.L2) {
        dd[ , cluster] <- lavInspect(object, "cluster.id")[[gg]]
        dd <- dd[ , c(cluster, setdiff(names(dd), cluster)), drop = FALSE ]
      } else {
        dd <- cbind(case.idx = lavInspect(object, "case.idx")[[gg]], dd)
      }

      ## attach cluster IDs, if multilevel and no level requested
      if (bothLevels) {
        dd[ , cluster] <- lavInspect(object, "cluster.label")[[gg]]

        d2 <- data.frame(PV2[[gg]])
        d2[ , group] <- group.label[gg]
        d2[ , cluster] <- lavInspect(object, "cluster.id")[[gg]]

        dd <- merge(dd, d2, by = c(group, cluster), all = TRUE)
      }

      dd
    })
    PV <- do.call(rbind, temp)

  } else {

    PV <- data.frame(PV)
    ## attach row indices from original data for optional merging
    if (only.L2) {
      PV[ , cluster] <- lavInspect(object, "cluster.id")
      PV <- PV[ , c(cluster, setdiff(names(PV), cluster)), drop = FALSE ]
    } else {
      PV <- cbind(case.idx = lavInspect(object, "case.idx"), PV)
    }
    ## attach cluster IDs, if multilevel and no level requested
    if (bothLevels) {
      PV[ , cluster] <- lavInspect(object, "cluster.label")

      PV2 <- data.frame(PV2)
      PV2[ , cluster] <- lavInspect(object, "cluster.id")

      PV <- merge(PV, PV2, by = cluster, all = TRUE)
    }

  }

  PV
}


## draw plausible values from a lavaan.mi object
##' @importFrom lavaan lavInspect lavPredict
plaus.mi <- function(object, seeds = 1:5, omit.imps = c("no.conv","no.se"), ...) {
  stopifnot(inherits(object, "lavaan.mi"))

  useImps <- rep(TRUE, length(object@DataList))
  if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
  if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
  if ("no.npd" %in% omit.imps) {
    Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
    Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
    useImps <- useImps & !(Heywood.lv | Heywood.ov)
  }
  m <- sum(useImps)
  useImps <- which(useImps)
  ## check if model has converged
  if (m == 0L) stop("No models converged. Score tests unavailable.")

  oldCall <- object@lavListCall
  if (!is.null(oldCall$parallel)) {
    if (oldCall$parallel == "snow") {
      oldCall$parallel <- "no"
      oldCall$ncpus <- 1L
      message("Unable to pass lavaan::lavPredict() arguments ",
              "when parallel='snow'. Switching to parallel='no'.",
              " Unless using Windows, parallel='multicore' should work.")
    }
  }
  ## call lavaanList() again to run lavTestScore() on each imputation
  oldCall$dataList <- object@DataList[useImps]
  oldCall$FUN <- function(obj) lapply(seeds, plaus.lavaan, object = obj, ...)
  FIT <- eval(as.call(oldCall))

  ## check if there are any results
  noFS <- sapply(FIT@funList, is.null)
  if (all(noFS)) stop("No success drawing plausible values for any imputations.")

  do.call(c, FIT@funList) # concatenate lists
}



## draw 1 set of plausible values from a blavaan object
##' @importFrom lavaan lavNames lavInspect
plaus.blavaan <- function(object) {
  stopifnot(inherits(object, "blavaan"))
  requireNamespace("blavaan")
  if (!"package:blavaan" %in% search()) attachNamespace("blavaan")

  # cluster <- lavInspect(object, "cluster")
  group <- lavInspect(object, "group")
  group.label <- lavInspect(object, "group.label")
  nG <- lavInspect(object, "ngroups")
  # nL <- lavInspect(object, "nlevels")
  case.idx <- lavInspect(object, "case.idx")

  ## stack factor scores from each chain (one row per PV)
  FS <- do.call(rbind, blavaan::blavInspect(object, "lvs"))
  ## column names contain indices to store PVs in matrix
  eta.idx <- colnames(FS)
  ## N and latent variable names, to know dimensions of PV
  N <- lavInspect(object, "ntotal")
  etas <- lavNames(object, "lv") #FIXME: assumes same model in both groups
  PV <- list()
  ## loop over rows (draws), assign columns to eta matrix, save in PV list
  for (i in 1:nrow(FS)) {

    eta <- matrix(NA, nrow = N, ncol = length(etas), dimnames = list(NULL, etas))
    for (j in eta.idx) eval(parse(text = paste(j, "<-", FS[i, j]) ))
    PV[[i]] <- data.frame(eta)

    ## add case indices, and groups (if applicable)
    if (nG == 1L) PV[[i]]$case.idx <- case.idx else {
      PV[[i]]$case.idx <- do.call(c, case.idx)
      PV[[i]][ , group] <- rep(group.label, times = lavInspect(object, "nobs"))
    }

  }

  PV
}



## ------
## Checks
## ------


# HS.model <- ' visual  =~ x1 + x2 + x3
#               textual =~ x4 + x5 + x6
#               speed   =~ x7 + x8 + x9 '
#
# fit1 <- cfa(HS.model, data = HolzingerSwineford1939)
# fs1 <- plausibleValues(fit1, nDraws = 3, append.data = T)
# lapply(fs1, head)
#
#
#
# step1 <- cfa(HS.model, data = HolzingerSwineford1939, group = "school",
#              group.equal = c("loadings","intercepts"))
# PV.list <- plausibleValues(step1, append.data = T)
# lapply(PV.list[1:3], head)
#
#
# model <- '
#   level: 1
#     fw =~ y1 + y2 + y3
#     fw ~ x1 + x2 + x3
#   level: 2
#     fb =~ y1 + y2 + y3
#     fb ~ w1 + w2
# '
# msem <- sem(model, data = Demo.twolevel, cluster = "cluster")
# mlPVs <- plausibleValues(msem, nDraws = 3, append.data = T) # both levels by default
# lapply(mlPVs, head, n = 10)
# ## only Level 1
# mlPV1 <- plausibleValues(msem, nDraws = 3, level = 1, append.data = T)
# lapply(mlPV1, head)
# ## only Level 2
# mlPV2 <- plausibleValues(msem, nDraws = 3, level = 2, append.data = T)
# lapply(mlPV2, head)
#
#
#
# data(Demo.twolevel)
# Demo.twolevel$g <- ifelse(Demo.twolevel$cluster %% 2L, "foo", "bar") # arbitrary groups
# table(Demo.twolevel$g)
# model2 <- ' group: foo
# level: within
#   fw =~ y1 + L2*y2 + L3*y3
#   fw ~ x1 + x2 + x3
# level: between
#   fb =~ y1 + L2*y2 + L3*y3
#   fb ~ w1 + w2
#
# group: bar
#
# level: within
#   fw =~ y1 + L2*y2 + L3*y3
#   fw ~ x1 + x2 + x3
# level: between
#   fb =~ y1 + L2*y2 + L3*y3
#   fb ~ w1 + w2
# '
# msem2 <- sem(model2, data = Demo.twolevel, cluster = "cluster", group = "g")
# ml2PVs <- plausibleValues(msem2, nDraws = 3, append.data = T) # both levels by default
# lapply(ml2PVs, head, n = 10)
# ## only Level 1
# ml2PV1 <- plausibleValues(msem2, nDraws = 3, level = 1, append.data = T)
# lapply(ml2PV1, head)
# ## only Level 2
# ml2PV2 <- plausibleValues(msem2, nDraws = 3, level = 2, append.data = T)
# lapply(ml2PV2, head)



