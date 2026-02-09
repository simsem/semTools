### Terrence D. Jorgensen
### Last updated: 9 February 2026
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
##' [lavaan::lavaan()] model, then treat them as multiple imputations
##' of missing data using [lavaan.mi::lavaan.mi()].
##'
##'
##' Because latent variables are unobserved, they can be considered as missing
##' data, which can be imputed using Monte Carlo methods.  This may be of
##' interest to researchers with sample sizes too small to fit their complex
##' structural models.  Fitting a factor model as a first step,
##' [lavaan::lavPredict()] provides factor-score estimates, which can
##' be treated as observed values in a path analysis (Step 2).  However, the
##' resulting standard errors and test statistics could not be trusted because
##' the Step-2 analysis would not take into account the uncertainty about the
##' estimated factor scores.  Using the asymptotic sampling covariance matrix
##' of the factor scores provided by [lavaan::lavPredict()],
##' `plausibleValues` draws a set of `nDraws` imputations from the
##' sampling distribution of each factor score, returning a list of data sets
##' that can be treated like multiple imputations of incomplete data.  If the
##' data were already imputed to handle missing data, `plausibleValues`
##' also accepts an object of class [lavaan.mi::lavaan.mi-class], and will
##' draw `nDraws` plausible values from each imputation.  Step 2 would
##' then take into account uncertainty about both missing values and factor
##' scores.  Bayesian methods can also be used to generate factor scores, as
##' available with the \pkg{blavaan} package, in which case plausible
##' values are simply saved parameters from the posterior distribution. See
##' Asparouhov and Muthen (2010) for further technical details and references.
##'
##' Each returned `data.frame` includes a `case.idx` column that
##' indicates the corresponding rows in the data set to which the model was
##' originally fitted (unless the user requests only Level-2 variables).  This
##' can be used to merge the plausible values with the original observed data,
##' but users should note that including any new variables in a Step-2 model
##' might not accurately account for their relationship(s) with factor scores
##' because they were not accounted for in the Step-1 model from which factor
##' scores were estimated.
##'
##' If `object` is a multilevel `lavaan` model, users can request
##' plausible values for latent variables at particular levels of analysis by
##' setting the [lavaan::lavPredict()] argument `level=1` or
##' `level=2`.  If the `level` argument is not passed via \dots,
##' then both levels are returned in a single merged data set per draw.  For
##' multilevel models, each returned `data.frame` also includes a column
##' indicating to which cluster each row belongs (unless the user requests only
##' Level-2 variables).
##'
##'
##' @importFrom lavaan lavInspect lavPredict
##'
##' @param object A fitted model of class [lavaan::lavaan-class],
##'   [blavaan::blavaan-class], or [lavaan.mi::lavaan.mi-class]
##' @param nDraws `integer` specifying the number of draws, analogous to
##'   the number of imputed data sets. If `object` is of class
##'   [lavaan.mi::lavaan.mi-class], this will be the number of draws taken
##'   *per imputation*.  If `object` is of class
##'   [blavaan::blavaan-class], `nDraws` cannot exceed
##'   `blavInspect(object, "niter") * blavInspect(bfitc, "n.chains")`
##'   (number of MCMC samples from the posterior). The drawn samples will be
##'   evenly spaced (after permutation for `target="stan"`), using
##'   [ceiling()] to resolve decimals.
##' @param seed `integer` passed to [set.seed()].
##' @param omit.imps `character` vector specifying criteria for omitting
##'   imputations when `object` is of class [lavaan.mi::lavaan.mi-class].
##'   Can include any of `c("no.conv", "no.se", "no.npd")`.
##' @param ... Optional arguments to pass to [lavaan::lavPredict()].
##'   `assemble` will be ignored because multiple groups are always
##'   assembled into a single `data.frame` per draw. `type` will be
##'   ignored because it is set internally to `type="lv"`.
##'
##' @return A `list` of length `nDraws`, each of which is a
##'   `data.frame` containing plausible values, which can be treated as
##'   a `list` of imputed data sets to be passed to the `lavaan.mi` package
##'   (see **Examples**). If `object` is of class
##'   [lavaan.mi::lavaan.mi-class], the `list` will be of length
##'   `nDraws*m`, where `m` is the number of imputations.
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Asparouhov, T. & Muthen, B. O. (2010). *Plausible values for latent
##'   variables using M*plus. Technical Report. Retrieved from
##'   www.statmodel.com/download/Plausible.pdf
##'
##' @seealso [lavaan.mi::lavaan.mi()], [lavaan.mi::lavaan.mi-class]
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
##' \donttest{
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
##' if(requireNamespace("lavaan.mi")){
##'   library(lavaan.mi)
##'   step2 <- sem.mi(path.model, data = PV.list, group = "school")
##'   ## test equivalence of both slopes across groups
##'   lavTestWald.mi(step2, constraints = 't1 == t2 ; s1 == s2')
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
##'
##'
##' ## example with 20 multiple imputations of missing data:
##' nPVs <- 5
##' nImps <- 20
##'
##' if (requireNamespace("lavaan.mi")) {
##'   data(HS20imps, package = "lavaan.mi")
##'
##'   ## specify CFA model from lavaan's ?cfa help page
##'   HS.model <- '
##'     visual  =~ x1 + x2 + x3
##'     textual =~ x4 + x5 + x6
##'     speed   =~ x7 + x8 + x9
##'   '
##'   out2 <- cfa.mi(HS.model, data = HS20imps)
##'   PVs <- plausibleValues(out2, nDraws = nPVs)
##'
##'   idx <- out2@@Data@@case.idx # can't use lavInspect() on lavaan.mi
##'   ## empty list to hold expanded imputations
##'   impPVs <- list()
##'   for (m in 1:nImps) {
##'     HS20imps[[m]]["case.idx"] <- idx
##'     for (i in 1:nPVs) {
##'       impPVs[[ nPVs*(m - 1) + i ]] <- merge(HS20imps[[m]],
##'                                             PVs[[ nPVs*(m - 1) + i ]],
##'                                             by = "case.idx")
##'     }
##'   }
##'   lapply(impPVs, head)
##' }
##'
##' }
##'
##' @export
plausibleValues <- function(object, nDraws = 20L, seed = 12345,
                            omit.imps = c("no.conv","no.se"), ...) {

  if (inherits(object, "lavaan")) {
    ## generate vector of seeds
    set.seed(seed)
    seeds <- sample(100000:9999999, size = nDraws, replace = FALSE)

    PV <- lapply(seeds, plaus.lavaan, object = object, ...)

  } else if (inherits(object, "lavaan.mi")) {
    ## generate vector of seeds
    set.seed(seed)
    seeds <- sample(100000:9999999, size = nDraws, replace = FALSE)

    PV <- plaus.mi(object, seeds = seeds, omit.imps = omit.imps, ...)

  } else if (inherits(object, "blavaan")) {
    PV <- plaus.blavaan(object, nDraws = nDraws, seed = seed, ...)
    #TODO: pass nDraws to sample() iterations?

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
    #TODO: verify this:
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
        mnormt::rmnorm(n = 1, mean = mu, varcov = ACOV)
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
          out <- mnormt::rmnorm(n = 1, mean = mu, varcov = ACOV2)
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
          mnormt::rmnorm(n = 1, mean = mu, varcov = ACOV[[gg]])
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
            mnormt::rmnorm(n = 1, mean = mu, varcov = ACOV2[[gg]])
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
  if (!"package:lavaan.mi" %in% search()) attachNamespace("lavaan.mi")

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
plaus.blavaan <- function(object, nDraws = 20L, seed = 12345, ...) {
  stopifnot(inherits(object, "blavaan"))
  if (!"package:blavaan" %in% search()) attachNamespace("blavaan")

  # cluster <- lavInspect(object, "cluster")
  group <- lavInspect(object, "group")
  group.label <- lavInspect(object, "group.label")
  nG <- lavInspect(object, "ngroups")
  # nL <- lavInspect(object, "nlevels")
  case.idx <- lavInspect(object, "case.idx")

  ## plausible values of what? (could be latent item responses)
  dots <- list(...)
  if (is.null(dots$type)) dots$type <- "lv" # default to factor scores
  dots$object <- object

  ## stack factor scores from each chain (one row per PV)
  FS <- do.call(blavaan::blavPredict, dots)
  #NOTE: might be latent responses

  ## only save nDraws from posterior
  if (nDraws >= length(FS)) {
    ## why would anyone want this many?  Or sample so few during estimation?
    message('nDraws cannot exceed number of iterations in `object=`. \nSet to ',
            'nDraws = blavInspect(object, "niter") * blavInspect(bfitc, "n.chains")')
    nDraws <- length(FS)
  }
  set.seed(seed)
  idx.sample <- ceiling(1:nDraws * length(FS)/nDraws)

  #FIXME: if Ed accepts pull request, format will be the same as c("yhat","ypred")
  if (dots$type == "lv" && utils::compareVersion(utils::packageDescription('blavaan')$Version, '0.4-2.949') >= 0L) {
    ## column names contain indices to store PVs in matrix
    eta.idx <- colnames(FS)
    ## N and latent variable names, to know dimensions of PV
    N <- lavInspect(object, "ntotal")
    etas <- lavNames(object, "lv") #FIXME: assumes same model in both groups
    PV <- list()
    ## loop over nDraws rows, assign columns to eta matrix, save in PV list
    set.seed(seed)
    idx.sample <- ceiling(1:nDraws * length(FS)/nDraws)
    for (i in idx.sample) {

      eta <- matrix(NA, nrow = N, ncol = length(etas), dimnames = list(NULL, etas))
      for (j in eta.idx) eval(parse(text = paste(j, "<-", FS[i, j]) ))
      PV[[i]] <- data.frame(eta)

      ## add case indices, and groups (if applicable)
      if (nG == 1L) PV[[i]]$case.idx <- case.idx else {
        PV[[i]]$case.idx <- do.call(c, case.idx)
        PV[[i]][ , group] <- rep(group.label, times = lavInspect(object, "nobs"))
      }

    }
  } else {
    ## latent responses, already a list
    PV <- list()
    for (i in idx.sample) {
      ## convert matrix to data.frame
      PV[[i]] <- data.frame(FS[[i]])
      ## add case indices, and groups (if applicable)
      if (nG == 1L) PV[[i]]$case.idx <- case.idx else {
        PV[[i]]$case.idx <- do.call(c, case.idx)
        PV[[i]][ , group] <- rep(group.label, times = lavInspect(object, "nobs"))
      }
    }
  } #else stop('Not implemented for blavPredict(object, type="', dots$type,'")')


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


## ordered-categorical data
# data(datCat)
#
# modc <- ' ## Set thresholds equal across groups
#   ## thresholds at Time 1
#     u1 | c(tau1.1, tau1.1)*t1 + c(tau1.2, tau1.2)*t2 + c(tau1.3, tau1.3)*t3 + c(tau1.4, tau1.4)*t4
#     u2 | c(tau2.1, tau2.1)*t1 + c(tau2.2, tau2.2)*t2 + c(tau2.3, tau2.3)*t3 + c(tau2.4, tau2.4)*t4
#     u3 | c(tau3.1, tau3.1)*t1 + c(tau3.2, tau3.2)*t2 + c(tau3.3, tau3.3)*t3 + c(tau3.4, tau3.4)*t4
#     u4 | c(tau4.1, tau4.1)*t1 + c(tau4.2, tau4.2)*t2 + c(tau4.3, tau4.3)*t3 + c(tau4.4, tau4.4)*t4
#   ## thresholds at Time 2 equal to Time 1
#     u5 | c(tau1.1, tau1.1)*t1 + c(tau1.2, tau1.2)*t2 + c(tau1.3, tau1.3)*t3 + c(tau1.4, tau1.4)*t4
#     u6 | c(tau2.1, tau2.1)*t1 + c(tau2.2, tau2.2)*t2 + c(tau2.3, tau2.3)*t3 + c(tau2.4, tau2.4)*t4
#     u7 | c(tau3.1, tau3.1)*t1 + c(tau3.2, tau3.2)*t2 + c(tau3.3, tau3.3)*t3 + c(tau3.4, tau3.4)*t4
#     u8 | c(tau4.1, tau4.1)*t1 + c(tau4.2, tau4.2)*t2 + c(tau4.3, tau4.3)*t3 + c(tau4.4, tau4.4)*t4
#   ## define latent responses as single-indicator factors (resid. var = 0)
#     y1 =~ 1*u1   ;   u1 ~~ c(0, 0)*u1
#     y2 =~ 1*u2   ;   u2 ~~ c(0, 0)*u2
#     y3 =~ 1*u3   ;   u3 ~~ c(0, 0)*u3
#     y4 =~ 1*u4   ;   u4 ~~ c(0, 0)*u4
#     y5 =~ 1*u5   ;   u5 ~~ c(0, 0)*u5
#     y6 =~ 1*u6   ;   u6 ~~ c(0, 0)*u6
#     y7 =~ 1*u7   ;   u7 ~~ c(0, 0)*u7
#     y8 =~ 1*u8   ;   u8 ~~ c(0, 0)*u8
#   ## only fix mean=0 in first group/occasion
#     y1 + y2 + y3 + y4 ~ c( 0, NA)*1
#     y5 + y6 + y7 + y8 ~ c(NA, NA)*1
#   ## only fix variance=1 in first groop/occasion
#     y1 ~~ c( 1, NA)*y1
#     y2 ~~ c( 1, NA)*y2
#     y3 ~~ c( 1, NA)*y3
#     y4 ~~ c( 1, NA)*y4
#     y5 ~~ c(NA, NA)*y5
#     y6 ~~ c(NA, NA)*y6
#     y7 ~~ c(NA, NA)*y7
#     y8 ~~ c(NA, NA)*y8
#   ## estimate all covariances
#     y1 ~~ y2 + y3 + y4 + y5 + y6 + y7 + y8
#     y2 ~~ y3 + y4 + y5 + y6 + y7 + y8
#     y3 ~~ y4 + y5 + y6 + y7 + y8
#     y4 ~~ y5 + y6 + y7 + y8
#     y5 ~~ y6 + y7 + y8
#     y6 ~~ y7 + y8
#     y7 ~~ y8
# '
# ## fit in lavaan for sanity check
# fitc <- lavaan(modc, data = datCat, group = "g", parameterization = "theta")
# summary(fitc)
#
# ## impose 5% MCAR
# set.seed(123)
# for (i in 1:8) datCat[sample(1:nrow(datCat), size = .05*nrow(datCat)), i] <- NA
#
# ## try with pairwise deletion
# fitcm <- lavaan(modc, data = datCat, group = "g", parameterization = "theta",
#                 missing = "pairwise")
# summary(fitcm)
#
# ## try blavaan
# data(datCat) # doesn't yet work with missing ordinal data
# bfitc <- blavaan(modc, data = datCat, group = "g", ordered = TRUE, # why is this needed when they are already ordered?
#                  n.chains = 2, burnin = 500, sample = 101, seed = 123,
#                  bcontrol = list(cores = 2),
#                  save.lvs = TRUE)
# summary(bfitc)
# LIRs <- blavPredict(bfitc, type = "ypred")
# yhats <- blavPredict(bfitc, type = "yhat") # same when resid. var = 0?
# fscores <- blavPredict(bfitc, type = "lv") # same as LIRs?
#
#
# length(LIRs)  # list: 1 N*p matrix per chain
# length(yhats) # list: 1 N*p matrix per chain
# length(fscores) # matrix: 1 row per chain, 1 column per [N, fs]
# ## all are basically interchangeable
# ch <- 2
# cor(cbind(yhats = as.numeric(yhats[[ch]]),
#           LIRs = as.numeric(LIRs[[ch]]),
#           fscores = fscores[ch,]))
# ## compare means (by group)
# aggregate(yhats[[ch]], by = datCat["g"], FUN = mean)
# aggregate( LIRs[[ch]], by = datCat["g"], FUN = mean)
# aggregate(matrix(fscores[ch,], ncol = 8), by = datCat["g"], FUN = mean)
# ## compare to lavaan
# do.call(rbind, sapply(lavInspect(fitc, "est"),
#                       function(i) i$alpha[,1], simplify = FALSE))
#
#
#
# ## now a CFA (so there are latent item responses and common factors)
# mod <- ' FU1 =~ u1 + u2 + u3 + u4
#          FU2 =~ u5 + u6 + u7 + u8 '
# fit  <-  cfa(mod, data = datCat, std.lv = TRUE)
# bfit <- bcfa(mod, data = datCat, std.lv = TRUE, ordered = TRUE,
#              n.chains = 2, burnin = 100, sample = 101, seed = 123,
#              save.lvs = TRUE)
# summary(bfit)
# fscores <- blavPredict(bfit, type = "lv")
# LIRs <- blavPredict(bfit, type = "ypred")



# FS <- plausibleValues(bfitc)
# LIRs <- plausibleValues(bfitc, type = "ypred")
