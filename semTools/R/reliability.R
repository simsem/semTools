### Sunthud Pornprasertmanit, Terrence D. Jorgensen, Yves Rosseel
### Last updated: 9 May 2022



## -----
## AVE()
## -----

##' Calculate average variance extracted
##'
##' Calculate average variance extracted (AVE) per factor from `lavaan` object
##'
##' The average variance extracted (AVE) can be calculated by
##'
##' \deqn{ AVE = \frac{\bold{1}^\prime
##' \textrm{diag}\left(\Lambda\Psi\Lambda^\prime\right)\bold{1}}{\bold{1}^\prime
##' \textrm{diag}\left(\hat{\Sigma}\right) \bold{1}}, }
##'
##' Note that this formula is modified from Fornell & Larcker (1981) in the case
##' that factor variances are not 1. The proposed formula from Fornell & Larcker
##' (1981) assumes that the factor variances are 1. Note that AVE will not be
##' provided for factors consisting of items with dual loadings. AVE is the
##' property of items but not the property of factors. AVE is calculated with
##' polychoric correlations when ordinal indicators are used.
##'
##' @importFrom lavaan lavInspect
##' @importFrom methods getMethod
##'
##' @param object A \code{\linkS4class{lavaan}} or
##'   \code{\linkS4class{lavaan.mi}} object, expected to contain only
##'   exogenous common factors (i.e., a CFA model). Cross-loadings are not
##'   allowed and will result in \code{NA} for any factor with indicator(s)
##'   that cross-load.
##' @param obs.var \code{logical} indicating whether to compute AVE using
##'   observed variances in the denominator. Setting \code{FALSE} triggers
##'   using model-implied variances in the denominator.
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'   imputations from pooled results.  Can include any of
##'   \code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
##'   default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (\code{"no.npd"}) would exclude any imputations which
##'   yielded a nonpositive definite covariance matrix for observed or
##'   latent variables, which would include any "improper solutions" such
##'   as Heywood cases.  NPD solutions are not excluded by default because
##'   they are likely to occur due to sampling error, especially in small
##'   samples.  However, gross model misspecification could also cause
##'   NPD solutions, users can compare pooled results with and without
##'   this setting as a sensitivity analysis to see whether some
##'   imputations warrant further investigation.
##' @param omit.factors \code{character} vector naming any common factors
##'   modeled in \code{object} whose indicators' AVE is not of interest.
##' @param dropSingle \code{logical} indicating whether to exclude factors
##'   defined by a single indicator from the returned results. If \code{TRUE}
##'   (default), single indicators will still be included in the \code{total}
##'   column when \code{return.total = TRUE}.
##' @param return.df \code{logical} indicating whether to return reliability
##'   coefficients in a \code{data.frame} (one row per group/level), which is
##'   possible when every model block includes the same factors (after excluding
##'   those in \code{omit.factors} and applying \code{dropSingle}).
##'
##' @return \code{numeric} vector of average variance extracted from indicators
##'   per factor.  For models with multiple "blocks" (any combination of groups
##'   and levels), vectors may be returned as columns in a \code{data.frame}
##'   with additional columns indicating the group/level (see \code{return.df=}
##'   argument description for caveat).
##'
##' @author
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##' Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models
##' with unobservable variables and measurement errors. \emph{Journal of
##' Marketing Research, 18}(1), 39--50. \doi{10.2307/3151312}
##'
##' @seealso \code{\link{compRelSEM}} for composite reliability estimates
##'
##' @examples
##' data(HolzingerSwineford1939)
##' HS9 <- HolzingerSwineford1939[ , c("x7","x8","x9")]
##' HSbinary <- as.data.frame( lapply(HS9, cut, 2, labels=FALSE) )
##' names(HSbinary) <- c("y7","y8","y9")
##' HS <- cbind(HolzingerSwineford1939, HSbinary)
##'
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ y7 + y8 + y9 '
##'
##' fit <- cfa(HS.model, data = HS, ordered = c("y7","y8","y9"), std.lv = TRUE)
##'
##' ## works for factors with exclusively continuous OR categorical indicators
##' AVE(fit) # uses observed (or unconstrained polychoric/polyserial) by default
##' AVE(fit, obs.var = FALSE)
##'
##'
##' ## works for multigroup models and for multilevel models (and both)
##' data(Demo.twolevel)
##' ## assign clusters to arbitrary groups
##' Demo.twolevel$g <- ifelse(Demo.twolevel$cluster %% 2L, "type1", "type2")
##' model2 <- ' group: type1
##'   level: within
##'     fac =~ y1 + L2*y2 + L3*y3
##'   level: between
##'     fac =~ y1 + L2*y2 + L3*y3
##'
##' group: type2
##'   level: within
##'     fac =~ y1 + L2*y2 + L3*y3
##'   level: between
##'     fac =~ y1 + L2*y2 + L3*y3
##' '
##' fit2 <- sem(model2, data = Demo.twolevel, cluster = "cluster", group = "g")
##' AVE(fit2)
##'
##'@export
AVE <- function(object, obs.var = TRUE, omit.imps = c("no.conv","no.se"),
                omit.factors = character(0), dropSingle = TRUE,
                return.df = TRUE) {
  ## numbers of blocks
  ngroups <- lavInspect(object, "ngroups")
  nLevels <- lavInspect(object, "nlevels")
  nblocks <- ngroups*nLevels #FIXME: always true?

  ## labels for groups
  if (ngroups > 1L) {
    group.label <- lavInspect(object, "group.label")
    blk.g.lab <- if (!length(group.label)) paste0("g", 1:ngroups) else group.label
  } else {
    group.label <- blk.g.lab <- NULL
  }
  ## labels for clusters
  if (nLevels > 1L) {
    #FIXME? lavInspect(object, "level.label") is always ==
    #       c("within", lavInspect(object, "cluster"))
    PT <- parTable(object)
    clus.label <- unique(PT$level)
    clus.label <- clus.label[which(clus.label != "")]
    clus.label <- clus.label[which(clus.label != 0)]
    blk.clus.lab <- if (is.numeric(clus.label)) {
      c("within", lavInspect(object, "cluster"))
    } else clus.label
  } else clus.label <- blk.clus.lab <- NULL
  ## labels for blocks
  if (nblocks > 1L) {
    block.label <- paste(rep(blk.g.lab, each = nLevels), blk.clus.lab,
                         sep = if (ngroups > 1L && nLevels > 1L) "_" else "")
  } else block.label <- NULL

  ## check for categorical
  anyCategorical <- lavInspect(object, "categorical")

  if (inherits(object, "lavaan")) {
    PHI   <- lavInspect(object, "cov.lv") # common-factor variance
    EST   <- lavInspect(object, "est")    # to extract loadings
    SIGMA <- lavInspect(object,           # total variance
                        ifelse(obs.var, "sampstat", "fitted"))
    if (nblocks == 1L) {
      PHI    <- list(PHI)
      LAMBDA <- list(EST$lambda)
      SIGMA  <- list(SIGMA$cov)
    } else {
      LAMBDA <- sapply(EST, "[[", i = "lambda", simplify = FALSE)
      SIGMA <- sapply(SIGMA, "[[", i = "cov", simplify = FALSE)
    }

  } else if (inherits(object, "lavaan.mi")) {
    useImps <- rep(TRUE, length(object@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

    ## common-factor variance
    phiList <- object@phiList[useImps]
    if (nblocks == 1L) for (i in 1:m) phiList[[i]] <- list(phiList[[i]])
    PHI <- list()
    for (b in 1:nblocks) {
      PHI[[ block.label[b] ]] <- Reduce("+", lapply(phiList, "[[", i = b) ) / m
    }

    ## loadings
    if (nblocks == 1L) {
      lamList <- lapply(object@coefList[useImps], "[[", i = "lambda")
      LAMBDA <- Reduce("+", lamList) / length(lamList)
    } else {
      LAMBDA <- vector("list", nblocks)
      for (b in 1:nblocks) {
        lamList <- lapply(object@coefList[useImps], function(i) i[[b]]$lambda)
        LAMBDA[[b]] <- Reduce("+", lamList) / length(lamList)
      }
    }

    ## total variance
    if (obs.var) {
      SIGMA <- vector("list", nblocks)
      ## loop over blocks to pool saturated-model (observed) matrices
      for (b in 1:nblocks) {
        covList <- lapply(object@h1List[useImps], function(i) i$implied$cov[[b]])
        SIGMA[[ block.label[b] ]] <- Reduce("+", covList) / m
      }
    } else {
      ## pooled model-implied matrices
      if (nblocks == 1L) {
        SIGMA <- getMethod("fitted", class(object))(object)["cov"] # retain list format
      } else {
        SIGMA <- sapply(getMethod("fitted", class(object))(object),
                        "[[", "cov", simplify = FALSE)
      }
    }

  } # end lavaan vs. lavaan.mi conditional

  ## scale polychoric/polyserial to modeled LRV scale
  if (anyCategorical) {
    SDs <- sapply(getScales(object, omit.imps = omit.imps),
                  FUN = function(x) diag(1 / as.numeric(x)),
                  simplify = FALSE)

    for (b in 1:nblocks) {
      dimnames(SDs[[b]]) <- dimnames(SIGMA[[b]])
      SIGMA[[b]] <- SDs[[b]] %*% SIGMA[[b]] %*% SDs[[b]]
    }
  }

  avevar <- list()
  for (b in 1:nblocks) {
    ## extract factor and indicator names
    LY <- LAMBDA[[b]]
    allIndNames <- rownames(LY)
    allFacNames <- colnames(LY)
    myFacNames <- setdiff(allFacNames, omit.factors)
    if (dropSingle) {
      multInd <- sapply(myFacNames, function(fn) sum(LY[,fn] != 0) > 1L)
      myFacNames <- myFacNames[multInd]
    }
    subLY <- LY[ , myFacNames, drop = FALSE]
    myIndNames <- rownames(subLY)[apply(subLY, 1L, function(x) any(x != 0))]

    ## check for cross-loadings
    Xload <- apply(subLY, 1L, function(x) sum(round(x, 5) != 0) > 1L)

    avevar[[b]] <- setNames(rep(NA, length(myFacNames)), nm = myFacNames)

    ## loop over factors
    for (fn in myFacNames) {
      idx <- which(subLY[,fn] != 0)
      if (any(Xload[idx])) next # cross-loading violates AVE definition
      commonVar <- sum(subLY[idx, fn]^2) * PHI[[b]][fn, fn]
      avevar[[b]][fn] <- commonVar / sum(diag(SIGMA[[b]])[ myIndNames[idx] ])
    }

  }

  ## drop list structure?
  if (nblocks == 1L) {
    avevar <- avevar[[1]]
    class(avevar) <- c("lavaan.vector","numeric")
    return(avevar)

  } else {
    facList <- lapply(avevar, names)
    sameNames <- all(sapply(2:nblocks, function(i) {
      isTRUE(all.equal(facList[[1]], facList[[i]]))
    } ))
    if (!(sameNames && return.df)) {
      ## can't simplify, return as a list
      for (i in seq_along(avevar)) class(avevar[[i]]) <- c("lavaan.vector","numeric")
      names(avevar) <- block.label
      return(avevar)
    }
  }

  ## concatenate each factor's AVE across blocks
  facRel <- sapply(facList[[1]], simplify = FALSE, FUN = function(nn) {
    sapply(avevar, "[[", i = nn, USE.NAMES = FALSE) # extract AVE for factor i
  })
  if (ngroups > 1L && nLevels > 1L) {
    out <- data.frame(group = rep(blk.g.lab, each = nLevels),
                      level = rep(blk.clus.lab, times = ngroups),
                      facRel)
  } else if (ngroups > 1L) {
    out <- data.frame(group = blk.g.lab, facRel)
  } else if (nLevels > 1L) {
    out <- data.frame(level = blk.clus.lab, facRel)
  }
  class(out) <- c("lavaan.data.frame","data.frame")

  out
}



## ------------
## compRelSEM()
## ------------

##' Composite Reliability using SEM
##'
##' Calculate composite reliability from estimated factor-model parameters
##'
##' Several coefficients for factor-analysis reliability have been termed
##' "omega", which Cho (2021) argues is a misleading misnomer and argues for
##' using \eqn{\rho} to represent them all, differentiated by descriptive
##' subscripts.  In our package, we number \eqn{\omega} based on commonly
##' applied calculations.
##'
##' Bentler (1968) first introduced factor-analysis reliability for a
##' unidimensional factor model with congeneric indicators, labeling the
##' coeficients \eqn{\alpha}.  McDonald (1999) later referred to this
##' \emph{and other reliability coefficients}, first as \eqn{\theta} (in 1970),
##' then as \eqn{\omega}, which is a source of confusion when reporting
##' coefficients (Cho, 2021).  Coefficients based on factor models were later
##' generalized to account for multidimenisionality (possibly with
##' cross-loadings) and correlated errors. The general \eqn{\omega} formula
##' implemented in this function is:
##'
##' \deqn{ \omega = \frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
##' Var\left( \psi \right)}{\bold{1}^\prime \hat{\Sigma} \bold{1}}, }
##'
##' where \eqn{\hat{\Sigma}} can be the model-implied covariance matrix from
##' either the saturated model (i.e., the "observed" covariance matrix, used by
##' default) or from the hypothesized CFA model, controlled by the
##' \code{obs.var} argument. A \eqn{k}-dimensional vector \eqn{\bold{1}} is used
##' to sum elements in the matrix. Note that if the model includes any directed
##' effects (latent regression slopes), all coefficients are calculated
##' from \bold{total} factor variances:
##' \code{\link[lavaan]{lavInspect}(object, "cov.lv")}.
##'
##' Assuming (essential) tau-equivalence (\code{tau.eq=TRUE}) makes \eqn{\omega}
##' equivalent to coefficient \eqn{\alpha} from classical test theory
##' (Cronbach, 1951):
##'
##' \deqn{ \alpha = \frac{k}{k - 1}\left[ 1 - \frac{\sum^{k}_{i = 1}
##' \sigma_{ii}}{\sum^{k}_{i = 1} \sigma_{ii} + 2\sum_{i < j} \sigma_{ij}}
##' \right],}
##'
##' where \eqn{k} is the number of items in a factor's composite,
##' \eqn{\sigma_{ii}} signifies item \emph{i}'s variance, and \eqn{\sigma_{ij}}
##' signifies the covariance between items \emph{i} and \emph{j}. Again, the
##' \code{obs.var} argument controls whether \eqn{\alpha} is calculated using
##' the observed or model-implied covariance matrix.
##'
##' By setting \code{return.total=TRUE}, one can estimate reliability for a
##' single composite calculated using all indicators in a multidimensional
##' CFA (Bentler, 1972, 2009). Setting \code{return.total = -1} will return
##' \bold{only} the total-composite reliability (not per factor).
##'
##' \bold{Higher-Order Factors}:
##' The reliability of a composite that represents a higher-order construct
##' requires partitioning the model-implied factor covariance matrix \eqn{\Phi}
##' in order to isolate the common-factor variance associated only with the
##' higher-order factor. Using a second-order factor model, the model-implied
##' covariance matrix of observed indicators \eqn{\hat{\Sigma}} can be
##' partitioned into 3 sources:
##' \enumerate{
##'   \item the second-order common-factor (co)variance:
##'         \eqn{\Lambda \bold{B} \Phi_2 \bold{B}^{\prime} \Lambda^{\prime}}
##'   \item the residual variance of the first-order common factors (i.e., not
##'         accounted for by the second-order factor):
##'         \eqn{\Lambda \Psi_{u} \Lambda^{\prime}}
##'   \item the measurement error of observed indicators: \eqn{\Theta}
##' }
##'
##' where \eqn{\Lambda} contains first-order factor loadings, \eqn{\bold{B}}
##' contains second-order factor loadings, \eqn{\Phi_2} is the model-implied
##' covariance matrix of the second-order factor(s), and \eqn{\Psi_{u}} is the
##' covariance matrix of first-order factor disturbances. In practice, we can
##' use the full \eqn{\bold{B}} matrix and full model-implied \eqn{\Phi} matrix
##' (i.e., including all latent factors) because the zeros in \eqn{\bold{B}}
##' will cancel out unwanted components of \eqn{\Phi}. Thus, we can calculate
##' the proportion of variance of a composite score calculated from the observed
##' indicators (e.g., a total score or scale mean) that is attributable to the
##' second-order factor (i.e., coefficient \eqn{\omega}):
##'
##' \deqn{\omega_{L1}=\frac{\bold{1}^{\prime} \Lambda \bold{B} \Phi \bold{B}^{\prime}
##'   \Lambda^{\prime} \bold{1} }{ \bold{1}^{\prime} \hat{\Sigma} \bold{1}}, }
##'
##' where \eqn{\bold{1}} is the \emph{k}-dimensional vector of 1s and \emph{k}
##' is the number of observed indicators in the composite. Note that a
##' higher-order factor can also have observed indicators.
##'
##' \bold{Categorical Indicators}:
##' When all indicators (per composite) are ordinal, the \code{ord.scale}
##' argument controls whether the coefficient is calculated on the
##' latent-response scale (\code{FALSE}) or on the observed ordinal scale
##' (\code{TRUE}, the default).  For \eqn{\omega}-type coefficients
##' (\code{tau.eq=FALSE}), Green and Yang's (2009, formula 21) approach is used
##' to transform factor-model results back to the ordinal response scale.
##' When \code{ord.scale=TRUE}, coefficient \eqn{\alpha} is calculated using the
##' covariance matrix calculated from the integer-valued numeric weights for
##' ordinal categories, consistent with its definition (Chalmers, 2018) and the
##' \code{alpha} function in the \code{psych} package; this implies
##' \code{obs.var=TRUE}, so \code{obs.var=FALSE} will be ignored.  When
##' \code{ord.scale=FALSE}, the standard \eqn{\alpha} formula is applied to the
##' polychoric correlation matrix ("ordinal \eqn{\alpha}"; Zumbo et al., 2007),
##' estimated from the saturated or hypothesized model (see \code{obs.var}),
##' and \eqn{\omega} is calculated from CFA results without applying Green and
##' Yang's (2009) correction (see Zumbo & Kroc's, 2019, for a rationalization).
##' No method has been proposed for calculating reliability with a mixture of
##' categorical and continuous indicators, so an error is returned if
##' \code{object} includes factors with a mixture of indicator types (unless
##' omitted using \code{omit.factors}). If categorical indicators load on a
##' different factor(s) than continuous indicators, then reliability will still
##' be calculated separately for those factors, but \code{return.total} must be
##' \code{FALSE} (unless \code{omit.factors} is used to isolate factors with
##' indicators of the same type).
##'
##' \bold{Multilevel Measurement Models}:
##' Under the default settings, \code{compRelSEM()} will apply the same formula
##' in each "block" (group and/or level of analysis). In the case of multilevel
##' SEMs, this yields "reliability" for latent within- and between-level
##' components, as proposed by Geldhof et al. (2014).  This is not recommended
##' because the coefficients do not correspond to actual composites that would
##' be calculated from the observed data.  Lai (2021) proposed coefficients for
##' reliability of actual composites, depending on the type of construct, which
##' requires specifying the names of constructs for which reliability is desired
##' (or multiple constructs whose indicators would compose a multidimensional
##' composite). Configural (\code{config=}) and/or \code{shared=} constructs
##' can be specified; the same construct can be specified in both arguments, so
##' that overall scale-reliability can be estimated for a shared construct by
##' including it in \code{config}.  Instead of organizing the output by block
##' (the default), specifying \code{config=} and/or \code{shared=} will prompt
##' organizing the output by \code{$config} and/or \code{$shared}.
##'
##' \itemize{
##'   \item The overall (\code{_2L}) scale reliability for \code{config}ural
##'   constructs is returned, along with the reliability of a purely
##'   individual-level composite (\code{_W}, calculated by cluster-mean
##'   centering).
##'   \item The reliability for a \code{shared} construct quantifies
##'   generalizability across both indicators and raters (i.e., subjects rating
##'   their cluster's construct).  Lüdtke et al. (2011) refer to these as
##'   measurement error and sampling error, respectively.  An interrater
##'   reliability (IRR) coefficient is also returned, quantifying
##'   generalizability across rater/sampling-error only. To obtain a
##'   scale-reliability coefficient (quantifying a shared construct's
##'   generalizability across indicator/measurement-error only), include the
##'   same factor name in \code{config=}.  Jak et al. (2021) recommended
##'   modeling components of the same construct at both levels, but users may
##'   also saturate the within-level model (Lai, 2021).
##' }
##'
##' Be careful about including Level-2 variables in the model, especially
##' whether it makes sense to include them in a total composite for a Level-2
##' construct.  \code{dropSingle=TRUE} only prevents estimating reliability for
##' a single-indicator construct, not from including such an indicator in a
##' total composite.  It is permissible for \code{shared=} constructs to have
##' indicators at Level-2 only.  If it is necessary to model other Level-2
##' variables (e.g., to justify the missing-at-random assumption when using
##' \code{missing = "FIML" estimation}), they should be placed in the
##' \code{omit.indicators=} argument to exclude them from total composites.
##'
##'
##' @importFrom lavaan lavInspect lavNames parTable
##' @importFrom methods getMethod
##'
##' @param object A \code{\linkS4class{lavaan}} or
##'   \code{\linkS4class{lavaan.mi}} object, expected to contain only
##'   exogenous common factors (i.e., a CFA model).
##' @param obs.var \code{logical} indicating whether to compute AVE using
##'   observed variances in the denominator. Setting \code{FALSE} triggers
##'   using model-implied variances in the denominator.
##' @param tau.eq \code{logical} indicating whether to assume (essential)
##'   tau-equivalence, yielding a coefficient analogous to \eqn{\alpha}.
##'   Setting \code{FALSE} yields an \eqn{\omega}-type coefficient.
##' @param ord.scale \code{logical} indicating whether to apply Green and Yang's
##'   (2009, formula 21) correction, so that reliability is calculated for the
##'   actual ordinal response scale (ignored for factors with continuous
##'   indicators).  Setting \code{FALSE} yields coefficients that are
##'   only applicable to the continuous latent-response scale.
##' @param config \code{character} vector naming any configural constructs in
##'   a multilevel CFA. For these constructs (and optional total composite),
##'   Lai's (2021) coefficients \eqn{\omega^\textrm{W}} and \eqn{\omega^\textrm{2L}}
##'   are returned (or corresponding \eqn{\alpha} coefficients when
##'   \code{tau.eq=TRUE}), rather than Geldhof et al.'s (2014) coefficients for
##'   hypothetical composites of latent components (although the same formula
##'   is used for \eqn{\omega^\textrm{W}} in either case). Note that the same name
##'   must be used for the factor component represented at each level of the
##'   model.
##' @param shared \code{character} vector naming any shared constructs in
##'   a multilevel CFA. For these constructs (and optional total composite),
##'   Lai's (2021) coefficient \eqn{\omega^\textrm{B}} or \eqn{\alpha^\textrm{B}} is
##'   returned, rather than Geldhof et al.'s (2014) between-level coefficient
##'   for hypothetical composites of latent cluster means. Lai's (2021)
##'   coefficient quantifies reliability relative to error associated with both
##'   indicators (measurement error) and subjects (sampling error), like a
##'   generalizability coefficient.  Given that subjects can be considered as
##'   raters of their cluster's shared construct, an interrater reliability
##'   (IRR) coefficient is also returned, quantifying reliability relative to
##'   rater/sampling error alone.  To quantify reliability relative to
##'   indicator/measurement error alone (i.e., \eqn{\omega^\textrm{2L}}), the
##'   \code{shared=} construct name(s) can additionally be included in
##'   \code{config=} argument.
##' @param higher \code{character} vector naming any higher-order constructs in
##'   \code{object} for which composite reliability should be calculated.
##'   Ignored when \code{tau.eq=TRUE} because alpha is not based on a CFA model;
##'   instead, users must fit a CFA with tau-equivalence constraints.
##'   To obtain Lai's (2021) multilevel composite-reliability indices for a
##'   higher-order factor, do not use this argument; instead, specify the
##'   higher-order factor(s) using the \code{shared=} or \code{config=} argument
##'   (\code{compRelSEM} will automatically check whether it includes latent
##'   indicators and apply the appropriate formula).
##' @param return.total \code{logical} indicating whether to return a final
##'   column containing the reliability of a composite of all indicators (not
##'   listed in \code{omit.indicators}) of first-order factors not listed in
##'   \code{omit.factors}.  Ignored in 1-factor models, and should only be set
##'   \code{TRUE} if all factors represent scale dimensions that could be
##'   meaningfully collapsed to a single composite (scale sum or scale mean).
##'   Setting a negative value (e.g., \code{-1} returns \bold{only} the
##'   total-composite reliability (excluding coefficients per factor), except
##'   when requesting Lai's (2021) coefficients for multilevel \code{config}ural
##'   or \code{shared=} constructs.
##' @param dropSingle \code{logical} indicating whether to exclude factors
##'   defined by a single indicator from the returned results. If \code{TRUE}
##'   (default), single indicators will still be included in the \code{total}
##'   column when \code{return.total = TRUE}.
##' @param omit.factors \code{character} vector naming any common factors
##'   modeled in \code{object} whose composite reliability is not of
##'   interest. For example, higher-order or method factors. Note that
##'   \code{\link{reliabilityL2}()} should be used to calculate composite
##'   reliability of a higher-order factor.
##' @param omit.indicators \code{character} vector naming any observed variables
##'   that should be ignored when calculating composite reliability. This can
##'   be useful, for example, to estimate reliability when an indicator is
##'   removed.
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'   imputations from pooled results.  Can include any of
##'   \code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
##'   default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (\code{"no.npd"}) would exclude any imputations which
##'   yielded a nonpositive definite covariance matrix for observed or
##'   latent variables, which would include any "improper solutions" such
##'   as Heywood cases.  NPD solutions are not excluded by default because
##'   they are likely to occur due to sampling error, especially in small
##'   samples.  However, gross model misspecification could also cause
##'   NPD solutions, users can compare pooled results with and without
##'   this setting as a sensitivity analysis to see whether some
##'   imputations warrant further investigation.
##' @param return.df \code{logical} indicating whether to return reliability
##'   coefficients in a \code{data.frame} (one row per group/level), which is
##'   possible when every model block includes the same factors (after excluding
##'   those in \code{omit.factors} and applying \code{dropSingle}).
##'
##' @return A \code{numeric} vector of composite reliability coefficients per
##'   factor, or a \code{list} of vectors per "block" (group and/or level of
##'   analysis), optionally returned as a \code{data.frame} when possible (see
##'   \code{return.df=} argument description for caveat). If there are multiple
##'   factors, whose multidimensional indicators combine into a single
##'   composite, users can request \code{return.total=TRUE} to add a column
##'   including a reliability coefficient for the total composite, or
##'   \code{return.total = -1} to return \bold{only} the total-composite
##'   reliability (ignored when \code{config=} or \code{shared=} is specified
##'   because each factor's specification must be checked across levels).
##'
##' @author
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##'   Uses hidden functions written by Sunthud Pornprasertmanit
##'   (\email{psunthud@@gmail.com}) for the old \code{reliability()} function.
##'
##' @seealso
##' \code{\link{maximalRelia}} for the maximal reliability of weighted composite
##'
##' @references
##' Bentler, P. M. (1968). Alpha-maximized factor analysis (alphamax): Its
##' relation to alpha and canonical factor analysis. \emph{Psychometrika, 33}(3),
##' 335--345. \doi{10.1007/BF02289328}
##'
##' Bentler, P. M. (1972). A lower-bound method for the dimension-free
##' measurement of internal consistency. \emph{Social Science Research, 1}(4),
##' 343--357. \doi{10.1016/0049-089X(72)90082-8}
##'
##' Bentler, P. M. (2009). Alpha, dimension-free, and model-based internal
##' consistency reliability. \emph{Psychometrika, 74}(1), 137--143.
##' \doi{10.1007/s11336-008-9100-1}
##'
##' Chalmers, R. P. (2018). On misconceptions and the limited usefulness of
##' ordinal alpha. \emph{Educational and Psychological Measurement, 78}(6),
##' 1056--1071. \doi{10.1177/0013164417727036}
##'
##' Cho, E. (2021) Neither Cronbach’s alpha nor McDonald’s omega: A commentary
##' on Sijtsma and Pfadt. \emph{Psychometrika, 86}(4), 877--886.
##' \doi{10.1007/s11336-021-09801-1}
##'
##' Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
##' tests. \emph{Psychometrika, 16}(3), 297--334. \doi{10.1007/BF02310555}
##'
##' Geldhof, G. J., Preacher, K. J., & Zyphur, M. J. (2014). Reliability
##' estimation in a multilevel confirmatory factor analysis framework.
##' \emph{Psychological Methods, 19}(1), 72--91. \doi{10.1037/a0032138}
##'
##' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
##' structural equation modeling: An alternative to coefficient alpha.
##' \emph{Psychometrika, 74}(1), 155--167. \doi{10.1007/s11336-008-9099-3}
##'
##' Jak, S., Jorgensen, T. D., & Rosseel, Y. (2021). Evaluating cluster-level
##' factor models with \code{lavaan} and M\emph{plus}. \emph{Psych, 3}(2),
##' 134--152. \doi{10.3390/psych3020012}
##'
##' Lai, M. H. C. (2021). Composite reliability of multilevel data: It’s about
##' observed scores and construct meanings. \emph{Psychological Methods, 26}(1),
##' 90--102. \doi{10.1037/met0000287}
##'
##' Lüdtke, O., Marsh, H. W., Robitzsch, A., & Trautwein, U. (2011).
##' A 2 \eqn{\times} 2 taxonomy of multilevel latent contextual models:
##' Accuracy--bias trade-offs in full and partial error correction models.
##' \emph{Psychological Methods, 16}(4), 444--467. \doi{10.1037/a0024376}
##'
##' McDonald, R. P. (1999). \emph{Test theory: A unified treatment}. Mahwah, NJ:
##' Erlbaum.
##'
##' Zumbo, B. D., Gadermann, A. M., & Zeisser, C. (2007). Ordinal versions of
##' coefficients alpha and theta for Likert rating scales.
##' \emph{Journal of Modern Applied Statistical Methods, 6}(1), 21--29.
##' \doi{10.22237/jmasm/1177992180}
##'
##' Zumbo, B. D., & Kroc, E. (2019). A measurement is a choice and Stevens’
##' scales of measurement do not help make it: A response to Chalmers.
##' \emph{Educational and Psychological Measurement, 79}(6), 1184--1197.
##' \doi{10.1177/0013164419844305}
##'
##'
##' @examples
##'
##' data(HolzingerSwineford1939)
##' HS9 <- HolzingerSwineford1939[ , c("x7","x8","x9")]
##' HSbinary <- as.data.frame( lapply(HS9, cut, 2, labels=FALSE) )
##' names(HSbinary) <- c("y7","y8","y9")
##' HS <- cbind(HolzingerSwineford1939, HSbinary)
##'
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ y7 + y8 + y9 '
##'
##' fit <- cfa(HS.model, data = HS, ordered = c("y7","y8","y9"), std.lv = TRUE)
##'
##' ## works for factors with exclusively continuous OR categorical indicators
##' compRelSEM(fit)
##'
##' ## reliability for ALL indicators only available when they are
##' ## all continuous or all categorical
##' compRelSEM(fit, omit.factors = "speed", return.total = TRUE)
##'
##'
##' ## loop over visual indicators to calculate alpha if one indicator is removed
##' for (i in paste0("x", 1:3)) {
##'   cat("Drop ", i, ":\n", sep = "")
##'   print(compRelSEM(fit, omit.factors = c("textual","speed"),
##'                    omit.indicators = i, tau.eq = TRUE))
##' }
##'
##'
##' ## Reliability of a composite that represents a higher-order factor
##' mod.hi <- ' visual  =~ x1 + x2 + x3
##'             textual =~ x4 + x5 + x6
##'             speed   =~ x7 + x8 + x9
##'             general =~ visual + textual + speed '
##'
##' fit.hi <- cfa(mod.hi, data = HolzingerSwineford1939)
##' compRelSEM(fit.hi, higher = "general")
##' ## reliabilities for lower-order composites also returned
##'
##'
##' ## works for multigroup models and for multilevel models (and both)
##' data(Demo.twolevel)
##' ## assign clusters to arbitrary groups
##' Demo.twolevel$g <- ifelse(Demo.twolevel$cluster %% 2L, "type1", "type2")
##' model2 <- ' group: type1
##'   level: 1
##'     f1 =~ y1 + L2*y2 + L3*y3
##'     f2 =~ y4 + L5*y5 + L6*y6
##'   level: 2
##'     f1 =~ y1 + L2*y2 + L3*y3
##'     f2 =~ y4 + L5*y5 + L6*y6
##'
##' group: type2
##'   level: 1
##'     f1 =~ y1 + L2*y2 + L3*y3
##'     f2 =~ y4 + L5*y5 + L6*y6
##'   level: 2
##'     f1 =~ y1 + L2*y2 + L3*y3
##'     f2 =~ y4 + L5*y5 + L6*y6
##' '
##' fit2 <- sem(model2, data = Demo.twolevel, cluster = "cluster", group = "g")
##' compRelSEM(fit2) # Geldhof's indices (hypothetical, for latent components)
##'
##' ## Lai's (2021) indices for Level-1 and configural constructs
##' compRelSEM(fit2, config = c("f1","f2"))
##' ## Lai's (2021) indices for shared (Level-2) constructs
##' ## (also an interrater reliability coefficient)
##' compRelSEM(fit2, shared = c("f1","f2"))
##'
##'
##' ## Shared construct using saturated within-level model
##' mod.sat1 <- ' level: 1
##'   y1 ~~ y1 + y2 + y3 + y4 + y5 + y6
##'   y2 ~~ y2 + y3 + y4 + y5 + y6
##'   y3 ~~ y3 + y4 + y5 + y6
##'   y4 ~~ y4 + y5 + y6
##'   y5 ~~ y5 + y6
##'   y6 ~~ y6
##'
##'   level: 2
##'   f1 =~ y1 + L2*y2 + L3*y3
##'   f2 =~ y4 + L5*y5 + L6*y6
##' '
##' fit.sat1 <- sem(mod.sat1, data = Demo.twolevel, cluster = "cluster")
##' compRelSEM(fit.sat1, shared = c("f1","f2"))
##'
##'
##' ## Simultaneous shared-and-configural model (Stapleton et al, 2016, 2019),
##' ## not recommended, but possible by omitting shared or configural factor.
##' mod.both <- ' level: 1
##'     fc =~ y1 + L2*y2 + L3*y3 + L4*y4 + L5*y5 + L6*y6
##'   level: 2
##'   ## configural construct
##'     fc =~ y1 + L2*y2 + L3*y3 + L4*y4 + L5*y5 + L6*y6
##'   ## orthogonal shared construct
##'     fs =~ NA*y1 + y2 + y3 + y4 + y5 + y6
##'     fs ~~ 1*fs + 0*fc
##' '
##' fit.both <- sem(mod.both, data = Demo.twolevel, cluster = "cluster")
##' compRelSEM(fit.both, shared = "fs", config = "fc")
##'
##' @export
compRelSEM <- function(object, obs.var = TRUE, tau.eq = FALSE, ord.scale = TRUE,
                       config = character(0), shared = character(0),
                       higher = character(0),
                       return.total = FALSE, dropSingle = TRUE,
                       omit.factors = character(0),
                       omit.indicators = character(0),
                       omit.imps = c("no.conv","no.se"), return.df = TRUE) {
  ## numbers of blocks
  ngroups <- lavInspect(object, "ngroups")
  nLevels <- lavInspect(object, "nlevels")
  nblocks <- ngroups*nLevels #FIXME: always true?
  return.total <- rep(return.total, nblocks)

  ## labels for groups
  if (ngroups > 1L) {
    group.label <- lavInspect(object, "group.label")
    blk.g.lab <- if (!length(group.label)) paste0("g", 1:ngroups) else group.label
  } else {
    group.label <- blk.g.lab <- NULL
  }
  ## labels for clusters
  if (nLevels > 1L) {
    #FIXME? lavInspect(object, "level.label") is always ==
    #       c("within", lavInspect(object, "cluster"))
    PT <- parTable(object)
    clus.label <- unique(PT$level)
    clus.label <- clus.label[which(clus.label != "")]
    clus.label <- clus.label[which(clus.label != 0)]
    blk.clus.lab <- if (is.numeric(clus.label)) {
      c("within", lavInspect(object, "cluster"))
    } else clus.label
  } else clus.label <- blk.clus.lab <- NULL
  ## labels for blocks
  if (nblocks > 1L) {
    block.label <- paste(rep(blk.g.lab, each = nLevels), blk.clus.lab,
                         sep = if (ngroups > 1L && nLevels > 1L) "_" else "")
  } else block.label <- NULL

  ## check for categorical
  anyCategorical <- lavInspect(object, "categorical")
  threshold <- if (anyCategorical) getThreshold(object, omit.imps = omit.imps) else NULL
  latScales <- if (anyCategorical) getScales(object, omit.imps = omit.imps) else NULL

  if (inherits(object, "lavaan")) {
    ## common-factor variance
    PHI <- lavInspect(object, "cov.lv") # ignored if tau.eq
    if (nblocks == 1L) PHI <- list(PHI)
    names(PHI) <- block.label

    ## factor loadings
    EST   <- lavInspect(object, "est", drop.list.single.group = FALSE)
    LAMBDA <- sapply(EST, "[[", i = "lambda", simplify = FALSE)
    names(LAMBDA) <- block.label
    ## possibly higher-order loadings?
    BETA   <- sapply(EST, "[[", i = "beta",   simplify = FALSE)
    names(BETA)   <- block.label

    ## total variance
    if (anyCategorical && tau.eq && ord.scale) {
      ## calculate SIGMA from data for alpha
      #FIXME when MLSEM available for ordinal indicators
      #     (Level 2 components available?  Extend conditional?)
      rawData <- try(lavInspect(object, "data", drop.list.single.group = FALSE),
                     silent = TRUE)
      if (inherits(rawData, "try-error"))
        stop('Error in lavInspect(fit, "data"); tau.eq= and ord.scale= cannot ',
             'both be TRUE for models fitted to summary statistics of ',
             'categorial data.')
      SIGMA <- lapply(rawData, cov)
      names(SIGMA) <- block.label

    } else {
      SIGMA <- sapply(lavInspect(object, drop.list.single.group = FALSE,
                                 what = ifelse(obs.var, "sampstat", "fitted")),
                      "[[", i = "cov", simplify = FALSE)
      names(SIGMA) <- block.label
    }

  } else if (inherits(object, "lavaan.mi")) {
    useImps <- rep(TRUE, length(object@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

    ## common-factor variance
    phiList <- object@phiList[useImps]
    if (nblocks == 1L) for (i in 1:m) phiList[[i]] <- list(phiList[[i]])
    PHI <- vector("list", nblocks)
    names(PHI) <- block.label
    for (b in 1:nblocks) {
      PHI[[ block.label[b] ]] <- Reduce("+", lapply(phiList, "[[", i = b) ) / m
    }

    ## loadings (including higher-order in Beta)
    if (nblocks == 1L) {
      lamList <- lapply(object@coefList[useImps], "[[", i = "lambda")
      LAMBDA <- list(Reduce("+", lamList) / length(lamList))

      betList <- lapply(object@coefList[useImps], "[[", i = "beta")
      if (length(betList)) {
        BETA <- list(Reduce("+", betList) / length(betList))
      } else BETA <- list(NULL)

    } else {
      LAMBDA <- BETA <- vector("list", nblocks)
      names(LAMBDA) <- names(BETA) <- block.label
      for (b in 1:nblocks) {
        lamList <- lapply(object@coefList[useImps], function(i) i[[b]]$lambda)
        LAMBDA[[ block.label[b] ]] <- Reduce("+", lamList) / length(lamList)
        betList <- lapply(object@coefList[useImps], function(i) i[[b]]$beta  )
        BETA[[   block.label[b] ]] <- Reduce("+", betList) / length(betList)
      }
    }

    ## total variance
    if (anyCategorical && tau.eq && ord.scale) {
      ## calculate SIGMA from data for alpha
      #FIXME when MLSEM available for ordinal indicators
      #     (Level 2 components available?  Extend conditional?)
      dataList <- object@DataList[useImps]
      SIGMA <- vector("list", nblocks)
      names(SIGMA) <- group.label #FIXME when MLSEMs can have ordinal indicators
      if (nblocks == 1L) {
        VV <- lavNames(object, type = "ov")
        impCovList <- lapply(dataList, function(DD) {
          dat <- do.call(cbind, sapply(DD[VV], as.numeric, simplify = FALSE))
          cov(dat)
        })
        SIGMA[[1]] <- Reduce("+", impCovList) / length(impCovList)

      } else {
        ## multigroup models need separate data matrices per group
        G <- lavInspect(object, "group")

        for (g in seq_along(group.label)) {
          VV <- lavNames(object, type = "ov",
                         group = ifelse(length(group.label),
                                        yes = group.label[g], no = g))

          impCovList <- lapply(dataList, function(DD) {
            RR <- DD[,G] == ifelse(length(group.label),
                                   yes = group.label[g], no = g)
            dat <- do.call(cbind, sapply(DD[RR, VV], as.numeric, simplify = FALSE))
            cov(dat)
          })
          SIGMA[[g]] <- Reduce("+", impCovList) / length(impCovList)
        } # g
      }

    } else {
      ## use model-implied SIGMA from h0 or h1 model
      if (obs.var) {
        SIGMA <- vector("list", nblocks)
        names(SIGMA) <- block.label
        ## loop over blocks to pool saturated-model (observed) matrices
        for (b in 1:nblocks) {
          covList <- lapply(object@h1List[useImps], function(i) i$implied$cov[[b]])
          SIGMA[[ block.label[b] ]] <- Reduce("+", covList) / m
        }
      } else {
        ## pooled model-implied matrices
        if (nblocks == 1L) {
          SIGMA <- getMethod("fitted", class(object))(object)["cov"] # retain list format
        } else {
          SIGMA <- sapply(getMethod("fitted", class(object))(object),
                          "[[", "cov", simplify = FALSE)
          names(SIGMA) <- block.label
        }
      }

    }

  } # end lavaan vs. lavaan.mi conditional

  ## scale polychoric/polyserial to modeled LRV scale
  if (anyCategorical) {
    SDs <- sapply(latScales, function(x) diag(1 / x),
                  simplify = FALSE)
    for (b in 1:nblocks) {
      dimnames(SDs[[b]]) <- dimnames(SIGMA[[b]])
      SIGMA[[b]] <- SDs[[b]] %*% SIGMA[[b]] %*% SDs[[b]]
    }
  }

  warnTotal <- warnAlpha <- warnOmega <- FALSE
  if (!length(c(config, shared))) {

    rel <- vector("list", length = nblocks)

    for (b in 1:nblocks) {

      LY <- LAMBDA[[b]]
      allIndNames <- rownames(LY)
      allFacNames <- colnames(LY)
      myFacNames <- setdiff(allFacNames, omit.factors)
      if (dropSingle) {
        multInd <- sapply(myFacNames, function(fn) sum(LY[,fn] != 0) > 1L)
        myFacNames <- myFacNames[multInd]
      }
      subLY <- LY[ , myFacNames, drop = FALSE]
      myIndNames <- rownames(subLY)[apply(subLY, 1L, function(x) any(x != 0))]
      ## remove unwanted indicators
      myIndNames <- setdiff(myIndNames, omit.indicators)
      subLY <- subLY[myIndNames, , drop = FALSE]

      ## distinguish between categorical, continuous, and latent indicators
      nameArgs <- list(object = object)
      if (nblocks > 1L) nameArgs$block <- b
      ordNames <- do.call(lavNames, c(nameArgs, list(type = "ov.ord")))
      numNames <- do.call(lavNames, c(nameArgs, list(type = "ov.num")))
      if (anyCategorical) {
        ## identify when the (sub)set of factors are all categorical
        blockCat <- all(myIndNames %in% ordNames)
        ## identify when the (sub)set of factors have mixed indicators, so no total
        mix <- any(myIndNames %in% ordNames) && any(myIndNames %in% numNames)
      } else {
        blockCat <- FALSE
        mix <- FALSE
      }

      if (mix && return.total[b]) {
        return.total[b] <- FALSE
        if (!(tau.eq && ord.scale)) warnTotal <- TRUE
      }

      ## compute reliability per factor?
      if (return.total[b] >= 0) {

        ## set result missing by default
        rel[[b]] <- setNames(rep(NA, length(myFacNames)), nm = myFacNames)

        for (fn in myFacNames) {
          ## names of indicators with nonzero loadings
          fIndNames <- myIndNames[which(subLY[,fn] != 0)]

          ## check for ANY indicators
          if (length(fIndNames) == 0L) next
          ## check for single indicators
          if (dropSingle && length(fIndNames) == 1L) next
          ## check for categorical (or mixed) indicators
          fCat <- any(fIndNames %in% ordNames)
          ## identify when this factor has mixed indicators, so no omegas
          fMix <- fCat && any(fIndNames %in% numNames)

          ## ALPHA
          totalCov  <- SIGMA[[b]][fIndNames, fIndNames, drop = FALSE]
          if (tau.eq) {
            if (fMix && !ord.scale) {
              ## can't mix observed and latent scales
              warnAlpha <- TRUE #TODO
              next
            }
            rel[[b]][fn] <- computeAlpha(totalCov)
            next
          } # else compute omega

          ## OMEGA
          if (fMix) {
            warnOmega <- TRUE # can't (yet) mix observed and latent scales
            next
          }
          Lf <- subLY[fIndNames, fn, drop = FALSE]
          commonCov <- Lf %*% PHI[[b]][fn, fn] %*% t(Lf)
          if (fCat && ord.scale) {
            ## Green & Yang (2009)
            rel[[b]][fn] <- omegaCat(truevar = commonCov, denom = totalCov,
                                     threshold = threshold[[b]][fIndNames],
                                     scales = latScales[[b]][fIndNames])
            next
          } # else, all continuous or all LRV-scale
          rel[[b]][fn] <- sum(commonCov) / sum(totalCov)
        } # end loop over factors
      } else rel[[b]] <- c(total = as.numeric(NA))

      ## compute for total composite?
      if (return.total[b] && length(myFacNames) > 1L) {

        ## ALPHA
        totalCov  <- SIGMA[[b]][myIndNames, myIndNames, drop = FALSE]
        if (tau.eq) {
          rel[[b]]["total"] <- computeAlpha(totalCov)
          next
        } # else compute omega

        ## OMEGA
        commonCov <- sum(subLY %*% PHI[[b]][myFacNames, myFacNames] %*% t(subLY))
        if (blockCat && ord.scale) {
          ## Green & Yang (2009)
          rel[[b]]["total"] <- omegaCat(truevar = commonCov, denom = totalCov,
                                        threshold = threshold[[b]][myIndNames],
                                        scales = latScales[[b]][myIndNames])
          next
        } # else, all continuous or all LRV-scale
        rel[[b]]["total"] <- sum(commonCov) / sum(totalCov)
      }

      ## composite(s) representing higher-order factor(s)?
      for (hf in higher) {
        ## find latent indicators
        L2 <- BETA[[b]][,hf]
        latInds <- setdiff(names(L2)[L2 != 0], omit.factors)
        ## find observed indicators
        indList <- lapply(c(hf, latInds), function(i) names(LY[,i])[ LY[,i] != 0])
        myIndNames <- setdiff(unique(do.call(c, indList)), omit.indicators)

        totalCov  <- SIGMA[[b]][myIndNames, myIndNames] # no need for drop = FALSE
        L <- LY[myIndNames, c(hf, latInds)]
        B <- BETA[[b]][c(hf, latInds), c(hf, latInds)]
        Phi <- PHI[[b]][c(hf, latInds), c(hf, latInds)]
        commonCov <- sum(L %*% B %*% Phi %*% t(B) %*% t(L))
        if (blockCat && ord.scale) {
          ## Green & Yang (2009)
          rel[[b]][hf] <- omegaCat(truevar = commonCov, denom = totalCov,
                                   threshold = threshold[[b]][myIndNames],
                                   scales = latScales[[b]][myIndNames])
          next
        } # else, all continuous or all LRV-scale
        rel[[b]][hf] <- sum(commonCov) / sum(totalCov)

      }

    }

    ## drop list structure
    if (nblocks == 1L) {
      rel <- rel[[1]]
      class(rel) <- c("lavaan.vector","numeric")
      return(rel)

    } else {
      facList <- lapply(rel, names)
      sameNames <- all(sapply(2:nblocks, function(i) {
        isTRUE(all.equal(facList[[1]], facList[[i]]))
      } ))
      if (!(sameNames && return.df)) {
        ## can't simplify, return as a list
        for (i in seq_along(rel)) class(rel[[i]]) <- c("lavaan.vector","numeric")
        names(rel) <- block.label
        return(rel)
      }
    }

    ## concatenate each factor's reliability across blocks
    facRel <- sapply(facList[[1]], simplify = FALSE, FUN = function(nn) {
      sapply(rel, "[[", i = nn, USE.NAMES = FALSE) # extract reliability for factor i
    })
    if (ngroups > 1L && nLevels > 1L) {
      out <- data.frame(group = rep(blk.g.lab, each = nLevels),
                        level = rep(blk.clus.lab, times = ngroups),
                        facRel)
    } else if (ngroups > 1L) {
      out <- data.frame(group = blk.g.lab, facRel)
    } else if (nLevels > 1L) {
      out <- data.frame(level = blk.clus.lab, facRel)
    }
    class(out) <- c("lavaan.data.frame","data.frame")

    return(out)
  }

  if (warnTotal) {
    message('Cannot return.total when model contains both continuous and ',
            'binary/ordinal observed indicators. Use the ',
            'omit.factors= argument to choose factors with only categorical ',
            'or continuous indicators, if that is a composite of interest.\n')
  }
  if (warnAlpha) {
    message('Coefficient alpha cannot be computed for factors as though a ',
            'composite would be calculated using both observed-response scales',
            ' (for continuous indicators) and latent-response scales (for ',
            'categorical indicators).  If you want to assume tau-equivalence, ',
            'either set ord.scale=FALSE or fit a model that treats ordinal ',
            'indicators as continuous.')
  }
  if (warnOmega) {
    message('Composite reliability (omega) cannot be computed for factors ',
            'with mixed categorical and continuous indicators, unless a model ',
            'is fitted by treating ordinal indicators as continuous.')
  }


  ## otherwise, only use Lai's MULTILEVEL coefficients
  if (nLevels > 1L && length(c(config, shared))) {

    ## group-level list, each containing 2 coefs per factor/total in data.frame
    rel <- vector("list", length = ngroups)

    for (g in 1:ngroups) {

      gLab <- ifelse(length(group.label), yes = group.label[g], no = g)
      nameArgs <- list(object = object, type = "lv")
      if (ngroups > 1L) nameArgs$group <- gLab
      lv.names1 <- do.call(lavNames, c(nameArgs, list(level = clus.label[1L])))
      lv.names2 <- do.call(lavNames, c(nameArgs, list(level = clus.label[2L])))
      nameArgs$type <- "ov"
      ov.names1 <- do.call(lavNames, c(nameArgs, list(level = clus.label[1L])))
      ov.names2 <- do.call(lavNames, c(nameArgs, list(level = clus.label[2L])))
      nameArgs$type <- "lv.ind"
      allLatInds <- do.call(lavNames, c(nameArgs, list(level = clus.label[1L])))

      ## erase higher= argument, use it to collect factor names for later checks
      higher <- character(0)

      PT <- parTable(object)
      PT <- PT[PT$op == "=~", ]
      if (ngroups > 1L) PT <- PT[PT$group == gLab, ]

      ## block indices for 2 levels in this group
      #FIXME: eventually possible to model partial clustering?
      idx1 <- 1 + (g-1)*2 # within
      idx2 <- 2 + (g-1)*2 # cluster

      ## configural construct(s) defined at both levels this group?
      for (fn in config) {
        if (fn %in% omit.factors) {
          ## why would they do this?
          config <- setdiff(config, omit.factors)
          next
        }
        if (fn %in% lv.names1 && fn %in% lv.names2) {
          ## same indicators for this construct at both levels?
          indNames1 <- setdiff(PT$rhs[PT$lhs == fn & PT$level == clus.label[1L]],
                               omit.indicators)
          indNames2 <- setdiff(PT$rhs[PT$lhs == fn & PT$level == clus.label[2L]],
                               omit.indicators)
          if (!all.equal(indNames1, indNames2)) {
            stop('After removing omit.indicators=, the indicators of factor ',
                 fn, ' do not match across levels',
                 ifelse(ngroups > 1L, paste(' in group', gLab), ""))
          }
          if (dropSingle && length(indNames1) == 1L) next
          ## is it a higher-order factor?
          latInds1 <- intersect(indNames1, allLatInds)
          if (length(latInds1)) {
            lowList <- lapply(latInds1, function(fn1) {
              indNames1 <- setdiff(PT$rhs[PT$lhs == fn1 & PT$level == clus.label[1L]],
                                   omit.indicators)
              indNames2 <- setdiff(PT$rhs[PT$lhs == fn1 & PT$level == clus.label[2L]],
                                   omit.indicators)
              ## check indicators also match for lower-order factors
              if (!all.equal(indNames1, indNames2)) {
                stop('After removing omit.indicators=, the indicators of factor ',
                     fn1, ' do not match across levels',
                     ifelse(ngroups > 1L, paste(' in group', gLab), ""))
              }
              ## check indicators of lower-order factors are not also latent
              lowerLatent <- intersect(indNames1, allLatInds)
              if (length(lowerLatent))
                stop('Indicators of lower-order factor ', fn1,
                     ifelse(ngroups > 1L, paste(' in group', gLab), ""),
                     ' cannot also be latent')
              indNames1
            })
            ## update indicator list to only include relevant observed variables
            indNames1 <- intersect(c(indNames1, do.call(c, lowList)), ov.names1)
            ## update list of higher-order factors
            higher <- c(higher, fn)
          }

          Sigma1 <- SIGMA[[idx1]][indNames1, indNames1, drop = FALSE]
          Sigma2 <- SIGMA[[idx2]][indNames1, indNames1, drop = FALSE]

          if (tau.eq) {
            if (fn %in% higher) {
              warning('Cannot apply alpha (tau.eq=TRUE) to higher-order factor ',
                      fn, '. Instead, fit a higher-order CFA that imposes ',
                      'tau-equivalence constraints.')
            } else {
              ## ALPHA
              rel[[g]]$config[[fn]] <- c(`alpha_W`  = computeAlpha(Sigma1),
                                         `alpha_2L` = computeAlpha(Sigma1 + Sigma2))
            }

          } else {
            ## OMEGA
            lam1 <- LAMBDA[[idx1]][indNames1, c(fn, latInds1), drop = FALSE]
            lam2 <- LAMBDA[[idx2]][indNames1, c(fn, latInds1), drop = FALSE]
            if (!isTRUE(all.equal(lam1, lam2)))
              warning('Unequal observed-indicator loadings across levels ',
                      'detected among factors (', paste(c(fn, latInds1), collapse = ","),
                      ifelse(ngroups > 1L, paste(') in group', gLab), ")"),
                      '. omega_2L for configural constructs assumes invariance.')
            phi1 <-  PHI[[idx1]][c(fn, latInds1), c(fn, latInds1), drop = FALSE]
            phi2 <-  PHI[[idx2]][c(fn, latInds1), c(fn, latInds1), drop = FALSE]

            if (length(latInds1)) {
              bet1 <- BETA[[idx1]][c(fn, latInds1), c(fn, latInds1), drop = FALSE]
              bet2 <- BETA[[idx2]][c(fn, latInds1), c(fn, latInds1), drop = FALSE]
              if (!isTRUE(all.equal(bet1, bet2)))
                warning('Unequal higher-order loadings detected across levels ',
                        'detected for factor ', fn,
                        ifelse(ngroups > 1L, paste(' in group', gLab), ""),
                        '. omega_2L for configural constructs assumes invariance.')
              commonCov1 <- lam1 %*% bet1 %*% phi1 %*% t(bet1) %*% t(lam1)
              commonCov2 <- lam2 %*% bet2 %*% phi2 %*% t(bet2) %*% t(lam2)
            } else {
              commonCov1 <- lam1 %*% phi1 %*% t(lam1)
              commonCov2 <- lam2 %*% phi2 %*% t(lam2)
            }

            rel[[g]]$config[[fn]] <- c(`omega_W`  = sum(commonCov1)              / sum(Sigma1),
                                       `omega_2L` = sum(commonCov1 + commonCov2) / sum(Sigma1 + Sigma2))
          }

        } else {
          warning('Configural factor ', fn, 'not detected at both levels of ',
                  'analysis, so removed from config= list.  Please use the ',
                  'same name for the within- and between-level component of a ',
                  'configural construct in your syntax.')
          config <- setdiff(config, fn) # rm non-configural construct
        }
      }
      ## after removing ineligible config, still multiple for total?
      if (length(config) > 1L) {
        ## only possible if NONE of config= are higher-order factors, or if
        ## the first-order factors in config= are not latent indicators of any
        ## higher-order factors in config=
        if (return.total[idx1] && !(length(higher) && any(config %in% allLatInds))) {

          indNames1 <- setdiff(PT$rhs[PT$lhs %in% config & PT$level == clus.label[1]],
                               omit.indicators)
          ## include observed indicators of any latent indicators
          latInds1 <- intersect(indNames1, allLatInds)
          indNames1 <- setdiff(PT$rhs[PT$lhs %in% c(config, latInds1) & PT$level == clus.label[1]],
                               omit.indicators)

          Sigma1 <- SIGMA[[idx1]][indNames1, indNames1, drop = FALSE]
          Sigma2 <- SIGMA[[idx2]][indNames1, indNames1, drop = FALSE]

          if (tau.eq) {
            if (any(config %in% higher)) {
              warning('Cannot apply alpha (tau.eq=TRUE) to total composite ',
                      'that includes higher-order configural factor(s):\n',
                      paste(intersect(shared, higher), collapse = ", "),
                      '\nInstead, impose tau-equiv. in a higher-order CFA')
            } else {
              ## ALPHA
              rel[[g]]$config$total <- c(`alpha_W`  = computeAlpha(Sigma1),
                                         `alpha_2L` = computeAlpha(Sigma1 + Sigma2))
            }

          } else {
            ## OMEGA
            lam1 <- LAMBDA[[idx1]][indNames1, c(config, latInds1), drop = FALSE]
            lam2 <- LAMBDA[[idx2]][indNames1, c(config, latInds1), drop = FALSE]
            phi1 <-  PHI[[idx1]][c(config, latInds1), c(config, latInds1), drop = FALSE]
            phi2 <-  PHI[[idx2]][c(config, latInds1), c(config, latInds1), drop = FALSE]
            if (length(latInds1)) {
              bet1 <- BETA[[idx1]][c(config, latInds1), c(config, latInds1), drop = FALSE]
              bet2 <- BETA[[idx2]][c(config, latInds1), c(config, latInds1), drop = FALSE]
              commonCov1 <- lam1 %*% bet1 %*% phi1 %*% t(bet1) %*% t(lam1)
              commonCov2 <- lam2 %*% bet2 %*% phi2 %*% t(bet2) %*% t(lam2)
            } else {
              commonCov1 <- lam1 %*% phi1 %*% t(lam1)
              commonCov2 <- lam2 %*% phi2 %*% t(lam2)
            }
            rel[[g]]$config$total <- c(`omega_W`  = sum(commonCov1)              / sum(Sigma1),
                                       `omega_2L` = sum(commonCov1 + commonCov2) / sum(Sigma1 + Sigma2))
          }
        }

        ## collect 2 coefs for multiple composites into a matrix
        rel[[g]]$config <- do.call(cbind, rel[[g]]$config)
      }


      ## (reciprocal of) harmonic-mean cluster size
      Ns <- mean(1 / lavInspect(object, "cluster.size",
                                drop.list.single.group = FALSE)[[g]])
      ## reset higher= argument for later checks
      higher <- character(0)

      ## shared construct(s) defined at between level in this group?
      for (fn in shared) {
        if (fn %in% omit.factors) {
          ## why would they do this?
          shared <- setdiff(shared, omit.factors)
          next
        }
        if (fn %in% lv.names2) {
          indNames2 <- setdiff(PT$rhs[PT$lhs == fn & PT$level == clus.label[2L]],
                               omit.indicators)
          ## only Level-2 single-indicator factors are relevant to drop
          if (dropSingle && length(indNames2) == 1L) next
          ## is it a higher-order factor?
          latInds2 <- intersect(indNames2, allLatInds)
          if (length(latInds2)) {
            lowList <- lapply(latInds2, function(fn2) {
              indNames2 <- setdiff(PT$rhs[PT$lhs == fn2 & PT$level == clus.label[2L]],
                                   omit.indicators)
              ## check indicators of lower-order factors are not also latent
              lowerLatent <- intersect(indNames2, allLatInds)
              if (length(lowerLatent))
                stop('Indicators of lower-order factor ', fn2,
                     ifelse(ngroups > 1L, paste(' in group', gLab), ""),
                     ' cannot also be latent')
              indNames2
            })
            ## update indicator list to only include relevant observed variables
            indNames2 <- intersect(c(indNames2, do.call(c, lowList)), ov.names2)
            ## update list of higher-order factors
            higher <- c(higher, fn)
          }
          ## capture within-level variance components of same indicators
          ## (make sure none are Level-2 only)
          indNames1 <- intersect(indNames2, ov.names1)

          Sigma1 <- SIGMA[[idx1]][indNames1, indNames1, drop = FALSE]
          Sigma2 <- SIGMA[[idx2]][indNames2, indNames2, drop = FALSE]

          if (tau.eq) {
            if (fn %in% higher) {
              warning('Cannot apply alpha (tau.eq=TRUE) to higher-order factor ',
                      fn, '. Instead, fit a higher-order CFA that imposes ',
                      'tau-equivalence constraints.')
            } else {
              nI <- length(indNames2)
              if (nI == 1L) {
                stop('Coefficient alpha is undefined for a single indicator. ',
                     'Set tau.eq=FALSE or dropSingle=TRUE')
              }
              kw <- nI / (nI-1) # weight for alpha based on number of items
              ## ALPHA
              onlyCov2 <- Sigma2
              diag(onlyCov2) <- 0

              rel[[g]]$shared[[fn]] <- c(`alpha_B` = kw*sum(onlyCov2) / sum(Sigma1*Ns, Sigma2),
                                         `IRR`     =    sum( Sigma2 ) / sum(Sigma1*Ns, Sigma2))
            }

          } else {
            ## OMEGA
            lam2 <- LAMBDA[[idx2]][indNames2, c(fn, latInds2), drop = FALSE]
            phi2 <- PHI[[idx2]][c(fn, latInds2), c(fn, latInds2), drop = FALSE]
            if (length(latInds2)) {
              bet2 <- BETA[[idx2]][c(fn, latInds2), c(fn, latInds2), drop = FALSE]
              commonCov2 <- lam2 %*% bet2 %*% phi2 %*% t(bet2) %*% t(lam2)
            } else commonCov2 <- lam2 %*% phi2 %*% t(lam2)

            rel[[g]]$shared[[fn]] <- c(`omega_B` = sum(commonCov2) / sum(Sigma1*Ns, Sigma2),
                                       `IRR`     = sum(  Sigma2  ) / sum(Sigma1*Ns, Sigma2))
          }

        } else shared <- setdiff(shared, fn) # rm non-shared construct
      }
      ## after removing ineligible shared, still multiple for total?
      if (length(shared) > 1L) {
        ## only possible if NONE of shared= are higher-order factors, or if
        ## the first-order factors in shared= are not latent indicators of any
        ## higher-order factors in shared=
        if (return.total[idx2] && !(length(higher) && any(shared %in% allLatInds))) {
          indNames2 <- setdiff(PT$rhs[PT$lhs %in% shared & PT$level == clus.label[2]],
                               omit.indicators)
          ## include observed indicators of any latent indicators
          latInds2 <- intersect(indNames2, allLatInds)
          indNames2 <- setdiff(PT$rhs[PT$lhs %in% c(shared, latInds2) & PT$level == clus.label[2]],
                               omit.indicators)
          ## capture within-level variance components of same indicators
          ## (make sure none are Level-2 only)
          indNames1 <- intersect(indNames2, ov.names1)

          Sigma1 <- SIGMA[[idx1]][indNames1, indNames1, drop = FALSE]
          Sigma2 <- SIGMA[[idx2]][indNames2, indNames2, drop = FALSE]

          if (tau.eq) {
            if (any(shared %in% higher)) {
              warning('Cannot apply alpha (tau.eq=TRUE) to total composite ',
                      'that includes higher-order shared factor(s):\n',
                      paste(intersect(shared, higher), collapse = ", "),
                      '\nInstead, impose tau-equiv. in a higher-order CFA')
            } else {
              ## ALPHA
              nI <- length(indNames2) #TODO: justify when > length(indNames1)
                                      #     (Level-1 component exists with SD=0)
              kw <- nI / (nI-1) # weight for alpha based on number of items
              onlyCov2 <- Sigma2
              diag(onlyCov2) <- 0

              rel[[g]]$shared$total <- c(`alpha_B` = kw*sum(onlyCov2) / sum(Sigma1*Ns, Sigma2),
                                         `IRR`     =    sum( Sigma2 ) / sum(Sigma1*Ns, Sigma2))
            }

          } else {
            ## OMEGA
            lam2 <- LAMBDA[[idx2]][indNames2, c(shared, latInds2), drop = FALSE]
            phi2 <- PHI[[idx2]][c(shared, latInds2), c(shared, latInds2), drop = FALSE]
            if (length(latInds1)) {
              bet2 <- BETA[[idx2]][c(shared, latInds2), c(shared, latInds2), drop = FALSE]
              commonCov2 <- lam2 %*% bet2 %*% phi2 %*% t(bet2) %*% t(lam2)
            } else commonCov2 <- lam2 %*% phi2 %*% t(lam2)
            rel[[g]]$shared$total <- c(`omega_B` = sum(commonCov2) / sum(Sigma1*Ns, Sigma2),
                                       `IRR`     = sum(  Sigma2  ) / sum(Sigma1*Ns, Sigma2))
          }
        }

        ## collect 2 coefs for multiple composites into a matrix
        rel[[g]]$shared <- do.call(cbind, rel[[g]]$shared)
      }


    } # end loop over groups

    ## drop list structure?
    if (ngroups == 1L) {
      rel <- rel[[1]]
    } else names(rel) <- group.label

  }

  rel
}



## -------------
## reliability()
## (deprecated)
## -------------


##' Composite Reliability using SEM
##'
##' Calculate composite reliability from estimated factor-model parameters
##'
##' The coefficient alpha (Cronbach, 1951) can be calculated by
##'
##' \deqn{ \alpha = \frac{k}{k - 1}\left[ 1 - \frac{\sum^{k}_{i = 1}
##' \sigma_{ii}}{\sum^{k}_{i = 1} \sigma_{ii} + 2\sum_{i < j} \sigma_{ij}}
##' \right],}
##'
##' where \eqn{k} is the number of items in a factor, \eqn{\sigma_{ii}} is the
##' item \emph{i} observed variances, \eqn{\sigma_{ij}} is the observed
##' covariance of items \emph{i} and \emph{j}.
##'
##' Several coefficients for factor-analysis reliability have been termed
##' "omega", which Cho (2021) argues is a misleading misnomer and argues for
##' using \eqn{\rho} to represent them all, differentiated by descriptive
##' subscripts.  In our package, we number \eqn{\omega} based on commonly
##' applied calculations.  Bentler (1968) first introduced factor-analysis
##' reliability for a unidimensional factor model with congeneric indicators.
##' However, assuming there are no cross-loadings in a multidimensional CFA,
##' this reliability coefficient can be calculated for each factor in the model.
##'
##' \deqn{ \omega_1 =\frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
##' Var\left( \psi \right)}{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
##' Var\left( \psi \right) + \sum^{k}_{i = 1} \theta_{ii} + 2\sum_{i < j}
##' \theta_{ij} }, }
##'
##' where \eqn{\lambda_i} is the factor loading of item \emph{i}, \eqn{\psi} is
##' the factor variance, \eqn{\theta_{ii}} is the variance of measurement errors
##' of item \emph{i}, and \eqn{\theta_{ij}} is the covariance of measurement
##' errors from item \emph{i} and \emph{j}. McDonald (1999) later referred to
##' this \emph{and other reliability coefficients} as "omega", which is a source
##' of confusion when reporting coefficients (Cho, 2021).
##'
##' The additional coefficients generalize the first formula by accounting for
##' multidimenisionality (possibly with cross-loadings) and correlated errors.
##' By setting \code{return.total=TRUE}, one can estimate reliability for a
##' single composite calculated using all indicators in the multidimensional
##' CFA (Bentler, 1972, 2009).  \code{"omega2"} is calculated by
##'
##' \deqn{ \omega_2 = \frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
##' Var\left( \psi \right)}{\bold{1}^\prime \hat{\Sigma} \bold{1}}, }
##'
##' where \eqn{\hat{\Sigma}} is the model-implied covariance matrix, and
##' \eqn{\bold{1}} is the \eqn{k}-dimensional vector of 1. The first and the
##' second coefficients omega will have the same value per factor in models with
##' simple structure, but they differ when there are (e.g.) cross-loadings
##' or method factors. The first coefficient omega can be viewed as the
##' reliability controlling for the other factors (like \eqn{\eta^2_{partial}} in
##' ANOVA). The second coefficient omega can be viewed as the unconditional
##' reliability (like \eqn{\eta^2} in ANOVA).
##'
##' The \code{"omega3"} coefficient (McDonald, 1999), sometimes referred to as
##' hierarchical omega, can be calculated by
##'
##' \deqn{ \omega_3 =\frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
##' Var\left( \psi \right)}{\bold{1}^\prime \Sigma \bold{1}}, }
##'
##' where \eqn{\Sigma} is the observed covariance matrix. If the model fits the
##' data well, \eqn{\omega_3} will be similar to \eqn{\omega_2}. Note that if
##' there is a directional effect in the model, all coefficients are calcualted
##' from total factor variances: \code{\link[lavaan]{lavInspect}(object, "cov.lv")}.
##'
##' In conclusion, \eqn{\omega_1}, \eqn{\omega_2}, and \eqn{\omega_3} are
##' different in the denominator. The denominator of the first formula assumes
##' that a model is congeneric factor model where measurement errors are not
##' correlated. The second formula accounts for correlated measurement errors.
##' However, these two formulas assume that the model-implied covariance matrix
##' explains item relationships perfectly. The residuals are subject to sampling
##' error. The third formula use observed covariance matrix instead of
##' model-implied covariance matrix to calculate the observed total variance.
##' This formula is the most conservative method in calculating coefficient
##' omega.
##'
##' The average variance extracted (AVE) can be calculated by
##'
##' \deqn{ AVE = \frac{\bold{1}^\prime
##' \textrm{diag}\left(\Lambda\Psi\Lambda^\prime\right)\bold{1}}{\bold{1}^\prime
##' \textrm{diag}\left(\hat{\Sigma}\right) \bold{1}}, }
##'
##' Note that this formula is modified from Fornell & Larcker (1981) in the case
##' that factor variances are not 1. The proposed formula from Fornell & Larcker
##' (1981) assumes that the factor variances are 1. Note that AVE will not be
##' provided for factors consisting of items with dual loadings. AVE is the
##' property of items but not the property of factors. AVE is calculated with
##' polychoric correlations when ordinal indicators are used.
##'
##' Coefficient alpha is by definition applied by treating indicators as numeric
##' (see Chalmers, 2018), which is consistent with the \code{alpha} function in
##' the \code{psych} package. When indicators are ordinal, \code{reliability}
##' additionally applies the standard alpha calculation to the polychoric
##' correlation matrix to return Zumbo et al.'s (2007) "ordinal alpha".
##'
##' Coefficient omega for categorical items is calculated using Green and Yang's
##' (2009, formula 21) approach. Three types of coefficient omega indicate
##' different methods to calculate item total variances. The original formula
##' from Green and Yang is equivalent to \eqn{\omega_3} in this function.
##' Green and Yang did not propose a method for
##' calculating reliability with a mixture of categorical and continuous
##' indicators, and we are currently unaware of an appropriate method.
##' Therefore, when \code{reliability} detects both categorical and continuous
##' indicators of a factor, an error is returned. If the categorical indicators
##' load on a different factor(s) than continuous indicators, then reliability
##' will still be calculated separately for those factors, but
##' \code{return.total} must be \code{FALSE} (unless \code{omit.factors} is used
##' to isolate factors with indicators of the same type).
##'
##'
##' @importFrom lavaan lavInspect lavNames
##' @importFrom methods getMethod
##'
##' @param object A \code{\linkS4class{lavaan}} or
##'   \code{\linkS4class{lavaan.mi}} object, expected to contain only
##'   exogenous common factors (i.e., a CFA model).
##' @param what \code{character} vector naming any reliability indices to
##'   calculate. All are returned by default. When indicators are ordinal,
##'   both traditional \code{"alpha"} and Zumbo et al.'s (2007) so-called
##'   "ordinal alpha" (\code{"alpha.ord"}) are returned, though the latter is
##'   arguably of dubious value (Chalmers, 2018).
##' @param return.total \code{logical} indicating whether to return a final
##'   column containing the reliability of a composite of all indicators (not
##'   listed in \code{omit.indicators}) of factors not listed in
##'   \code{omit.factors}.  Ignored in 1-factor models, and should only be set
##'   \code{TRUE} if all factors represent scale dimensions that could be
##'   meaningfully collapsed to a single composite (scale sum or scale mean).
##' @param dropSingle \code{logical} indicating whether to exclude factors
##'   defined by a single indicator from the returned results. If \code{TRUE}
##'   (default), single indicators will still be included in the \code{total}
##'   column when \code{return.total = TRUE}.
##' @param omit.factors \code{character} vector naming any common factors
##'   modeled in \code{object} whose composite reliability is not of
##'   interest. For example, higher-order or method factors. Note that
##'   \code{\link{reliabilityL2}()} should be used to calculate composite
##'   reliability of a higher-order factor.
##' @param omit.indicators \code{character} vector naming any observed variables
##'   that should be ignored when calculating composite reliability. This can
##'   be useful, for example, to estimate reliability when an indicator is
##'   removed.
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'   imputations from pooled results.  Can include any of
##'   \code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
##'   default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (\code{"no.npd"}) would exclude any imputations which
##'   yielded a nonpositive definite covariance matrix for observed or
##'   latent variables, which would include any "improper solutions" such
##'   as Heywood cases.  NPD solutions are not excluded by default because
##'   they are likely to occur due to sampling error, especially in small
##'   samples.  However, gross model misspecification could also cause
##'   NPD solutions, users can compare pooled results with and without
##'   this setting as a sensitivity analysis to see whether some
##'   imputations warrant further investigation.
##'
##' @return Reliability values (coefficient alpha, coefficients omega, average
##'   variance extracted) of each factor in each group. If there are multiple
##'   factors, a \code{total} column can optionally be included.
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##'   Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
##'
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##' Bentler, P. M. (1972). A lower-bound method for the dimension-free
##' measurement of internal consistency. \emph{Social Science Research, 1}(4),
##' 343--357. \doi{10.1016/0049-089X(72)90082-8}
##'
##' Bentler, P. M. (2009). Alpha, dimension-free, and model-based internal
##' consistency reliability. \emph{Psychometrika, 74}(1), 137--143.
##' \doi{10.1007/s11336-008-9100-1}
##'
##' Chalmers, R. P. (2018). On misconceptions and the limited usefulness of
##' ordinal alpha. \emph{Educational and Psychological Measurement, 78}(6),
##' 1056--1071. \doi{10.1177/0013164417727036}
##'
##' Cho, E. (2021) Neither Cronbach’s alpha nor McDonald’s omega: A commentary
##' on Sijtsma and Pfadt. *Psychometrika, 86*(4), 877--886.
##' \doi{10.1007/s11336-021-09801-1}
##'
##' Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
##' tests. \emph{Psychometrika, 16}(3), 297--334. \doi{10.1007/BF02310555}
##'
##' Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models
##' with unobservable variables and measurement errors. \emph{Journal of
##' Marketing Research, 18}(1), 39--50. \doi{10.2307/3151312}
##'
##' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
##' structural equation modeling: An alternative to coefficient alpha.
##' \emph{Psychometrika, 74}(1), 155--167. \doi{10.1007/s11336-008-9099-3}
##'
##' McDonald, R. P. (1999). \emph{Test theory: A unified treatment}. Mahwah, NJ:
##' Erlbaum.
##'
##' Raykov, T. (2001). Estimation of congeneric scale reliability using
##' covariance structure analysis with nonlinear constraints \emph{British
##' Journal of Mathematical and Statistical Psychology, 54}(2), 315--323.
##' \doi{10.1348/000711001159582}
##'
##' Zumbo, B. D., Gadermann, A. M., & Zeisser, C. (2007). Ordinal versions of
##' coefficients alpha and theta for Likert rating scales.
##' \emph{Journal of Modern Applied Statistical Methods, 6}(1), 21--29.
##' \doi{10.22237/jmasm/1177992180}
##'
##' Zumbo, B. D., & Kroc, E. (2019). A measurement is a choice and Stevens’
##' scales of measurement do not help make it: A response to Chalmers.
##' \emph{Educational and Psychological Measurement, 79}(6), 1184--1197.
##' \doi{10.1177/0013164419844305}
##'
##'
##' @examples
##'
##' data(HolzingerSwineford1939)
##' HS9 <- HolzingerSwineford1939[ , c("x7","x8","x9")]
##' HSbinary <- as.data.frame( lapply(HS9, cut, 2, labels=FALSE) )
##' names(HSbinary) <- c("y7","y8","y9")
##' HS <- cbind(HolzingerSwineford1939, HSbinary)
##'
##' HS.model <- ' visual  =~ x1 + x2 + x3
##'               textual =~ x4 + x5 + x6
##'               speed   =~ y7 + y8 + y9 '
##'
##' fit <- cfa(HS.model, data = HS, ordered = c("y7","y8","y9"), std.lv = TRUE)
##'
##' ## works for factors with exclusively continuous OR categorical indicators
##' reliability(fit)
##'
##' ## reliability for ALL indicators only available when they are
##' ## all continuous or all categorical
##' reliability(fit, omit.factors = "speed", return.total = TRUE)
##'
##'
##' ## loop over visual indicators to calculate alpha if one indicator is removed
##' for (i in paste0("x", 1:3)) {
##'   cat("Drop x", i, ":\n")
##'   print(reliability(fit, omit.factors = c("textual","speed"),
##'                     omit.indicators = i, what = "alpha"))
##' }
##'
##'
##' ## works for multigroup models and for multilevel models (and both)
##' data(Demo.twolevel)
##' ## assign clusters to arbitrary groups
##' Demo.twolevel$g <- ifelse(Demo.twolevel$cluster %% 2L, "type1", "type2")
##' model2 <- ' group: type1
##'   level: within
##'     fac =~ y1 + L2*y2 + L3*y3
##'   level: between
##'     fac =~ y1 + L2*y2 + L3*y3
##'
##' group: type2
##'   level: within
##'     fac =~ y1 + L2*y2 + L3*y3
##'   level: between
##'     fac =~ y1 + L2*y2 + L3*y3
##' '
##' fit2 <- sem(model2, data = Demo.twolevel, cluster = "cluster", group = "g")
##' reliability(fit2, what = c("alpha","omega3"))
##'
##' @name reliability-deprecated
##' @usage
##' reliability(object, what = c("alpha", "omega", "omega2", "omega3", "ave"),
##'             return.total = FALSE, dropSingle = TRUE, omit.factors = character(0),
##'             omit.indicators = character(0), omit.imps = c("no.conv", "no.se"))
##' @seealso \code{\link{semTools-deprecated}}
##' @keywords internal
NULL


##' @rdname semTools-deprecated
##' @section Reliability:
##' The original \code{reliability} function was suboptimally designed.
##' For example, AVE was returned, which is not a reliability index. Also,
##' alpha and several omega-type coefficients were returned, including the
##' original formula that was in appropriate for models with complex structure.
##' Some features could be controlled by the user for one but not both types of
##' index  For example, alpha for categorical indicators was returned on both
##' the observed and latent-response scales, but this was not an option for any
##' omega-type indices.  The omegas differed in terms of whether the observed or
##' model-implied covariance matrix was used in the denominator, but alpha was
##' only computed using the observed matrix.  These inconsistencies have been
##' resolved in the new \code{\link{compRelSEM}} function, which returns only
##' one reliability index (per factor, optionally total score) according to the
##' user's requested features, for which there is much more flexibility.
##' Average variance extracted is now available in a dedicated \code{\link{AVE}}
##' function.
##'
##' @export
reliability <- function(object,
                        what = c("alpha","omega","omega2","omega3","ave"),
                        return.total = FALSE, dropSingle = TRUE,
                        omit.factors = character(0),
                        omit.indicators = character(0),
                        omit.imps = c("no.conv","no.se")) {

  ngroups <- lavInspect(object, "ngroups") #TODO: adapt to multiple levels
  nLevels <- lavInspect(object, "nlevels")
  nblocks <- ngroups*nLevels #FIXME: always true?
  return.total <- rep(return.total, nblocks)
  group.label <- if (ngroups > 1L) lavInspect(object, "group.label") else NULL
  #FIXME? lavInspect(object, "level.labels")
  clus.label <- if (nLevels > 1L) c("within", lavInspect(object, "cluster")) else NULL
  if (nblocks > 1L) {
    block.label <- paste(rep(group.label, each = nLevels), clus.label,
                         sep = if (ngroups > 1L && nLevels > 1L) "_" else "")
  }

  ## check for categorical (determines what S will be)
  anyCategorical <- lavInspect(object, "categorical")
  if (anyCategorical && "alpha" %in% what) {
    what <- c(what, "alpha.ord")
    what <- unique(what) # in case it was already explicitly requested
  }
  ## categorical-model parameters
  threshold <- if (anyCategorical) getThreshold(object, omit.imps = omit.imps) else NULL
  latScales <- if (anyCategorical) getScales(object, omit.imps = omit.imps) else NULL
  ## all other relevant parameters in GLIST format (not flat, need block-level list)
  if (inherits(object, "lavaan")) {
    param <- lavInspect(object, "est")
    ve <- lavInspect(object, "cov.lv") # model-implied latent covariance matrix
    S <- object@h1$implied$cov # observed sample covariance matrix (already a list)
    if (anyCategorical && any(c("alpha","alpha.ord") %in% what)) {
      rawData <- try(lavInspect(object, "data"), silent = TRUE)
      if (inherits(rawData, "try-error"))
        stop('Error in lavInspect(fit, "data"); what="alpha" unavailable for ',
             'models fitted to summary statistics of categorial data.')
      if (nblocks == 1L) rawData <- list(rawData)
      S.as.con <- lapply(rawData, cov) # for actual "alpha", not "alpha.ord"
    }

    if (nblocks == 1L) {
      param <- list(param)
      ve <- list(ve)
    }

  } else if (inherits(object, "lavaan.mi")) {
    useImps <- rep(TRUE, length(object@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

    param <- object@coefList[[ useImps[1] ]] # first admissible as template
    coefList <- object@coefList[useImps]
    phiList <- object@phiList[useImps]
    if (anyCategorical) {
      dataList <- object@DataList[useImps]
      S.as.con <- vector("list", nblocks) # for group-list of pooled S
    }
    ## add block-level list per imputation?
    if (nblocks == 1L) {
      param <- list(param)
      for (i in 1:m) {
        coefList[[i]] <- list(coefList[[i]])
        phiList[[i]] <- list(phiList[[i]])
      }
      if (anyCategorical) { #FIXME: currently no categorical ML-SEMs
        #dataList[[i]] <- list(dataList[[i]])
        VV <- lavNames(object, type = "ov")
        impCovList <- lapply(dataList, function(DD) {
          dat <- do.call(cbind, sapply(DD[VV], as.numeric, simplify = FALSE))
          cov(dat)
        })
        S.as.con[[1]] <- Reduce("+", impCovList) / length(impCovList)
      }

    } else if (anyCategorical) { #FIXME: currently no categorical ML-SEMs
      ## multigroup models need separate data matrices per group
      G <- lavInspect(object, "group")

      for (g in seq_along(group.label)) {
        VV <- try(lavNames(object, type = "ov", group = group.label[g]),
                  silent = TRUE)
        if (inherits(VV, "try-error")) {
          VV <- lavNames(object, type = "ov", group = g)
        }
        impCovList <- lapply(dataList, function(DD) {
          RR <- DD[,G] == group.label[g]
          dat <- do.call(cbind, sapply(DD[RR, VV], as.numeric, simplify = FALSE))
          cov(dat)
        })
        S.as.con[[g]] <- Reduce("+", impCovList) / length(impCovList)
      }

    }
    S <- vector("list", nblocks) # pooled observed OR polychoric covariance matrix
    ve <- vector("list", nblocks)
    ## loop over blocks
    for (b in 1:nblocks) {

      ## param:  loop over GLIST elements
      for (mat in names(param[[b]])) {
        matList <- lapply(coefList, function(i) i[[b]][[mat]])
        param[[b]][[mat]] <- Reduce("+", matList) / length(matList)
      } # mat

      ## pooled observed OR polychoric covariance matrix
      covList <- lapply(object@h1List[useImps], function(i) i$implied$cov[[b]])
      S[[b]] <- Reduce("+", covList) / m

      ## pooled model-implied latent covariance matrix
      ve[[b]] <- Reduce("+", lapply(phiList, "[[", i = b) ) / m

    } # b

  }

  if (nblocks == 1L) {
    SigmaHat <- getMethod("fitted", class(object))(object)["cov"] # retain list format
  } else {
    SigmaHat <- sapply(getMethod("fitted", class(object))(object),
                       "[[", "cov", simplify = FALSE)
  }

  ly <- lapply(param, "[[", "lambda")
  te <- lapply(param, "[[", "theta")
  beta <- if ("beta" %in% names(param[[1]])) {
    lapply(param, "[[", "beta")
  } else NULL

	result <- list()
	warnTotal <- FALSE
	warnHigher <- character(0) # collect list of potential higher-order factors
	## loop over i blocks (groups/levels)
	for (i in 1:nblocks) {
	  ## extract factor and indicator names
	  allIndNames <- rownames(ly[[i]])
	  allFacNames <- colnames(ly[[i]])
	  myFacNames <- setdiff(allFacNames, omit.factors)
	  subLY <- ly[[i]][ , myFacNames, drop = FALSE] != 0
	  myIndNames <- rownames(subLY)[apply(subLY, MARGIN = 1L, FUN = any)]

	  ## distinguish between categorical, continuous, and latent indicators
	  nameArgs <- list(object = object)
	  if (nblocks > 1L) nameArgs$block <- i
	  ordNames <- do.call(lavNames, c(nameArgs, list(type = "ov.ord")))
	  numNames <- do.call(lavNames, c(nameArgs, list(type = "ov.num")))
	  if (anyCategorical) {
	    ## identify when the (sub)set of factors are all categorical
	    blockCat <- all(myIndNames %in% ordNames)
	    ## identify when the (sub)set of factors have mixed indicators, so no total
	    mix <- any(myIndNames %in% ordNames) && any(myIndNames %in% numNames)
	  } else {
	    blockCat <- FALSE
	    mix <- FALSE
	  }

	  if (mix && return.total[i]) {
	    return.total[i] <- FALSE
	    warnTotal <- TRUE
    }

	  ## identify POSSIBLE higher-order factors (that affect other latent vars)
	  latInds  <- do.call(lavNames, c(nameArgs, list(type = "lv.ind")))
	  higher <- if (length(latInds) == 0L) character(0) else {
	    allFacNames[apply(beta[[i]], MARGIN = 2, function(x) any(x != 0))]
	  }
	  ## keep track of factor indices to skip
	  idx.drop <- numeric(0)

	  ## relevant quantities
		common <- (apply(ly[[i]], 2, sum)^2) * diag(ve[[i]])
		truevar <- ly[[i]] %*% ve[[i]] %*% t(ly[[i]])
		## vectors to store results for each factor
		error <- rep(NA, length(common))
		alpha <- rep(NA, length(common))
		alpha.ord <- rep(NA, length(common))
		total <- rep(NA, length(common))
		omega1 <- omega2 <- omega3 <- rep(NA, length(common))
		impliedTotal <- rep(NA, length(common))
		avevar <- rep(NA, length(common))
		warnOmega <- FALSE
		## loop over j factors
		for (j in 1:length(common)) {
		  ## skip this factor?
		  if (allFacNames[j] %in% omit.factors) {
		    idx.drop <- c(idx.drop, j)
		    next
		  }

			index <- setdiff(which(ly[[i]][,j] != 0), # nonzero loadings
			                 which(allIndNames %in% omit.indicators))
			jIndNames <- allIndNames[index]

			## identify when this factor has mixed indicators, so no omegas
			jMix <- any(jIndNames %in% ordNames) && any(jIndNames %in% numNames)

			## check for ANY indicators (possibly skip purely higher-order factors)
			if (length(index) == 0L) {
			  idx.drop <- c(idx.drop, j)
			  next
			}
			## check for single indicators
			if (dropSingle && length(index) == 1L) {
			  idx.drop <- c(idx.drop, j)
			  next
			}
			## check for categorical (or mixed) indicators
			jCat <-      any(jIndNames %in% ordNames)
			warnOmega <- jCat && !all(jIndNames %in% ordNames)
			## check for latent indicators
			if (allFacNames[j] %in% higher && !(allFacNames[j] %in% omit.factors)) {
			  warnHigher <- c(warnHigher, allFacNames[j])
			}

			sigma <- S[[i]][index, index, drop = FALSE]
			faccontrib <- ly[[i]][,j, drop = FALSE] %*% ve[[i]][j,j, drop = FALSE] %*% t(ly[[i]][,j, drop = FALSE])
			truefac <- diag(faccontrib[index, index, drop = FALSE])
			trueitem <- diag(truevar[index, index, drop = FALSE])
			erritem <- diag(te[[i]][index, index, drop = FALSE])
			if (sum(abs(trueitem - truefac)) < 0.00001 & "ave" %in% what) {
				avevar[j] <- sum(trueitem) / sum(trueitem + erritem)
			}
			if (jCat) {
			  if ("alpha" %in% what) {
			    alpha[j] <- computeAlpha(S.as.con[[i]][index, index, drop = FALSE])
			  }
			  if ("alpha.ord" %in% what) {
			    alpha.ord[j] <- computeAlpha(sigma)
			  }
			  if ("omega" %in% what) {
			    omega1[j] <- omegaCat(truevar = faccontrib[index, index, drop = FALSE],
			                          threshold = threshold[[i]][jIndNames],
			                          scales = latScales[[i]][index],
			                          denom = faccontrib[index, index, drop = FALSE] + te[[i]][index, index, drop = FALSE])
			  }
			  if ("omega2" %in% what) {
			    omega2[j] <- omegaCat(truevar = faccontrib[index, index, drop = FALSE],
			                          threshold = threshold[[i]][jIndNames],
			                          scales = latScales[[i]][index],
			                          denom = SigmaHat[[i]][index, index, drop = FALSE])
			  }
			  if ("omega3" %in% what) {
			    omega3[j] <- omegaCat(truevar = faccontrib[index, index, drop = FALSE],
			                          threshold = threshold[[i]][jIndNames],
			                          scales = latScales[[i]][index],
			                          denom = sigma)
			  }

			} else {
			  alpha[j] <- computeAlpha(sigma)

			  commonfac <- sum(faccontrib[index, index, drop = FALSE])
			  error[j] <- sum(te[[i]][index, index, drop = FALSE])
			  impliedTotal[j] <- sum(SigmaHat[[i]][index, index, drop = FALSE])
			  total[j] <- sum(sigma)

			  omega1[j] <- commonfac / (commonfac + error[j])
				omega2[j] <- commonfac / impliedTotal[j]
				omega3[j] <- commonfac / total[j]
			}
			## end loop over j factors
		}

		if (return.total[i] & length(myFacNames) > 1L) {
		  if (blockCat) {
		    if ("alpha" %in% what) {
		      alpha <- c(alpha, computeAlpha(S.as.con[[i]]))
		    }
		    if ("alpha.ord" %in% what) {
		      alpha.ord <- c(alpha.ord, total = computeAlpha(S[[i]]))
		    }
		    if ("omega" %in% what) {
		      omega1 <- c(omega1, total = omegaCat(truevar = truevar,
		                                           threshold = threshold[[i]],
		                                           scales = latScales[[i]],
		                                           denom = truevar + te[[i]]))
		    }
		    if ("omega2" %in% what) {
		      omega2 <- c(omega2, total = omegaCat(truevar = truevar,
		                                           threshold = threshold[[i]],
		                                           scales = latScales[[i]],
		                                           denom = SigmaHat[[i]]))
		    }
		    if ("omega2" %in% what) {
  		    omega3 <- c(omega3, total = omegaCat(truevar = truevar,
  		                                         threshold = threshold[[i]],
  		                                         scales = latScales[[i]],
  		                                         denom = S[[i]]))
		    }

		  } else {
		    alpha <- c(alpha, total = computeAlpha(S[[i]]))
		    omega1 <- c(omega1, total = sum(truevar) / (sum(truevar) + sum(te[[i]])))
		    omega2 <- c(omega2, total = sum(truevar) / (sum(SigmaHat[[i]])))
		    omega3 <- c(omega3, total = sum(truevar) / (sum(S[[i]])))
		  }
		  avevar <- c(avevar,
		              total = sum(diag(truevar)) / sum((diag(truevar) + diag(te[[i]]))))
		}

		if (all(is.na(alpha.ord))) alpha.ord <- NULL
		result[[i]] <- rbind(alpha = if ("alpha" %in% what) alpha else NULL,
		                     alpha.ord = if ("alpha.ord" %in% what) alpha.ord else NULL,
		                     omega  = if ("omega"  %in% what) omega1 else NULL,
		                     omega2 = if ("omega2" %in% what) omega2 else NULL,
		                     omega3 = if ("omega3" %in% what) omega3 else NULL,
		                     avevar = if ("ave" %in% what) avevar else NULL)
		colnames(result[[i]])[1:length(allFacNames)] <- allFacNames
		if (return.total[i] & length(myFacNames) > 1L) {
		  colnames(result[[i]])[ ncol(result[[i]]) ] <- "total"
		}
		if (length(idx.drop)) {
		  result[[i]] <- result[[i]][ , -idx.drop, drop = FALSE]
		  ## reset indices for next block (could have different model/variables)
		  idx.drop <- numeric(0)
		}
		## end loop over blocks
	}

	warnCat <- sapply(result, function(x) any(c("alpha.ord","ave") %in% rownames(x)))
	if (any(warnCat)) {
	  alphaMessage <- paste0('Zumbo et al.`s (2007) "ordinal alpha" is calculated',
	                         ' in addition to the standard alpha, which treats ',
	                         'ordinal variables as numeric. See Chalmers (2018) ',
	                         'for a critique of "alpha.ord" and the response by ',
	                         'Zumbo & Kroc (2019).')
	  AVEmessage <- paste0('average variance extracted is calculated from ',
	                       'polychoric (polyserial) not Pearson correlations.')
	  both <- "alpha.ord" %in% what & "ave" %in% what
	  connectMessage <- if (both) ' Likewise, ' else ' the '
	  catMessage <- paste0("For constructs with categorical indicators, ",
	                       if ("alpha.ord" %in% what) alphaMessage else NULL,
	                       if (both) ' Likewise, ' else NULL,
	                       if ("ave" %in% what) AVEmessage else NULL)
	  if ("alpha.ord" %in% what || "ave" %in% what) message(catMessage, "\n")
	}
	if (length(warnHigher)) warning('Possible higher-order factors detected:\n',
	                                paste(unique(warnHigher), sep = ", "))
	if (warnTotal) {
	  message('Cannot return.total when model contains both continuous and ',
	          'binary/ordinal observed indicators. Use the ',
	          'omit.factors= argument to choose factors with only categorical ',
	          'indicators, if that is a composite of interest.\n')
	}
	if (warnOmega) {
	  message('Composite reliability (omega) cannot be computed for factors ',
	          'with mixed categorical and continuous indicators.')
	}

	## drop list structure?
	if (nblocks == 1L) {
		result <- result[[1]]
	} else names(result) <- block.label

	result
}



## ---------------
## reliabilityL2()
## (deprecated)
## ---------------

##' Calculate the reliability values of a second-order factor
##'
##' Calculate the reliability values (coefficient omega) of a second-order
##' factor
##'
##' The first formula of the coefficient omega (in the
##' \code{\link{reliability}}) will be mainly used in the calculation. The
##' model-implied covariance matrix of a second-order factor model can be
##' separated into three sources: the second-order common-factor variance,
##' the residual variance of the first-order common factors (i.e., not
##' accounted for by the second-order factor), and the measurement error of
##' observed indicators:
##'
##' \deqn{ \hat{\Sigma} = \Lambda \bold{B} \Phi_2 \bold{B}^{\prime}
##' \Lambda^{\prime} + \Lambda \Psi_{u} \Lambda^{\prime} + \Theta, }
##'
##' where \eqn{\hat{\Sigma}} is the model-implied covariance matrix,
##' \eqn{\Lambda} contains first-order factor loadings, \eqn{\bold{B}} contains
##' second-order factor loadings, \eqn{\Phi_2} is the covariance matrix of the
##' second-order factor(s), \eqn{\Psi_{u}} is the covariance matrix of residuals
##' from first-order factors, and \eqn{\Theta} is the covariance matrix of the
##' measurement errors from observed indicators. Thus, we can calculate the
##' proportion of variance of a composite score calculated from the observed
##' indicators (e.g., a total score or scale mean) that is attributable to the
##' second-order factor, i.e. coefficient omega at Level 1:
##'
##' \deqn{ \omega_{L1} = \frac{\bold{1}^{\prime} \Lambda \bold{B} \Phi_2
##' \bold{B}^{\prime} \Lambda^{\prime} \bold{1}}{\bold{1}^{\prime} \Lambda
##' \bold{B} \Phi_2 \bold{B} ^{\prime} \Lambda^{\prime} \bold{1} +
##' \bold{1}^{\prime} \Lambda \Psi_{u} \Lambda^{\prime} \bold{1} +
##' \bold{1}^{\prime} \Theta \bold{1}}, }
##'
##' where \eqn{\bold{1}} is the \emph{k}-dimensional vector of 1 and \emph{k} is
##' the number of observed variables.
##'
##' The model-implied covariance matrix among first-order factors (\eqn{\Phi_1})
##' can be calculated as:
##'
##' \deqn{ \Phi_1 = \bold{B} \Phi_2 \bold{B}^{\prime} + \Psi_{u}, }
##'
##' Thus, the proportion of variance among first-order common factors that is
##' attributable to the second-order factor (i.e., coefficient omega at Level 2)
##' can be calculated as:
##'
##' \deqn{ \omega_{L2} = \frac{\bold{1_F}^{\prime} \bold{B} \Phi_2
##' \bold{B}^{\prime} \bold{1_F}}{\bold{1_F}^{\prime} \bold{B} \Phi_2
##' \bold{B}^{\prime} \bold{1_F} + \bold{1_F}^{\prime} \Psi_{u} \bold{1_F}}, }
##'
##' where \eqn{\bold{1_F}} is the \emph{F}-dimensional vector of 1 and \emph{F}
##' is the number of first-order factors. This Level-2 omega can be interpreted
##' as an estimate of the reliability of a hypothetical composite calculated
##' from error-free observable variables representing the first-order common
##' factors. This might only be meaningful as a thought experiment.
##'
##' An additional thought experiment is possible: If the observed indicators
##' contained only the second-order common-factor variance and unsystematic
##' measurement error, then there would be no first-order common factors because
##' their unique variances would be excluded from the observed measures. An
##' estimate of this hypothetical composite reliability can be calculated as the
##' partial coefficient omega at Level 1, or the proportion of observed
##' variance explained by the second-order factor after partialling out the
##' uniqueness from the first-order factors:
##'
##' \deqn{ \omega_{L1} = \frac{\bold{1}^{\prime} \Lambda \bold{B} \Phi_2
##' \bold{B}^{\prime} \Lambda^{\prime} \bold{1}}{\bold{1}^{\prime} \Lambda
##' \bold{B} \Phi_2 \bold{B}^{\prime} \Lambda^{\prime} \bold{1} +
##' \bold{1}^{\prime} \Theta \bold{1}}, }
##'
##' Note that if the second-order factor has a direct factor loading on some
##' observed variables, the observed variables will be counted as first-order
##' factors, which might not be desirable.
##'
##'
##' @importFrom lavaan lavInspect
##'
##' @param object A \code{\linkS4class{lavaan}} or
##'   \code{\linkS4class{lavaan.mi}} object, expected to contain a least one
##'   exogenous higher-order common factor.
##' @param secondFactor The name of a single second-order factor in the
##'   model fitted in \code{object}. The function must be called multiple
##'   times to estimate reliability for each higher-order factor.
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'        imputations from pooled results.  Can include any of
##'        \code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
##'        default setting, which excludes any imputations that did not
##'        converge or for which standard errors could not be computed.  The
##'        last option (\code{"no.npd"}) would exclude any imputations which
##'        yielded a nonpositive definite covariance matrix for observed or
##'        latent variables, which would include any "improper solutions" such
##'        as Heywood cases.  NPD solutions are not excluded by default because
##'        they are likely to occur due to sampling error, especially in small
##'        samples.  However, gross model misspecification could also cause
##'        NPD solutions, users can compare pooled results with and without
##'        this setting as a sensitivity analysis to see whether some
##'        imputations warrant further investigation.
##'
##' @return Reliability values at Levels 1 and 2 of the second-order factor, as
##'   well as the partial reliability value at Level 1
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' @examples
##'
##' HS.model3 <- ' visual  =~ x1 + x2 + x3
##'                textual =~ x4 + x5 + x6
##'                speed   =~ x7 + x8 + x9
##'                higher =~ visual + textual + speed'
##'
##' fit6 <- cfa(HS.model3, data = HolzingerSwineford1939)
##' reliability(fit6) # Should provide a warning for the endogenous variables
##' reliabilityL2(fit6, "higher")
##'
##' @name reliabilityL2-deprecated
##' @usage
##' reliabilityL2(object, secondFactor, omit.imps = c("no.conv", "no.se"))
##' @seealso \code{\link{semTools-deprecated}}
##' @keywords internal
NULL


##' @rdname semTools-deprecated
##' @section Higher-Order Reliability:
##' Originally, composite reliability of a single higher-order factor was
##' estimated in a separate function: \code{reliabilityL2}.  It is now available
##' for multiple higher-order factors in the same model, and from the same
##' \code{\link{compRelSEM}} function that estimates reliability for first-order
##' factors, using the \code{higher=} argument to name higher-order factors in
##' the fitted model.
##'
##' @export
reliabilityL2 <- function(object, secondFactor,
                          omit.imps = c("no.conv","no.se")) {
  secondFactor <- as.character(secondFactor)[1] # only one at a time

  ngroups <- lavInspect(object, "ngroups") #TODO: adapt to multiple levels
  nLevels <- lavInspect(object, "nlevels")
  nblocks <- ngroups*nLevels #FIXME: always true?
  group.label <- if (ngroups > 1L) lavInspect(object, "group.label") else NULL
  #FIXME? lavInspect(object, "level.labels")
  clus.label <- if (nLevels > 1L) c("within", lavInspect(object, "cluster")) else NULL
  if (nblocks > 1L) {
    block.label <- paste(rep(group.label, each = nLevels), clus.label,
                         sep = if (ngroups > 1L && nLevels > 1L) "_" else "")
  }

  ## parameters in GLIST format (not flat, need block-level list)
  if (inherits(object, "lavaan")) {
    param <- lavInspect(object, "est")
    ve <- lavInspect(object, "cov.lv") # model-implied latent covariance matrix
    S <- object@h1$implied$cov # observed sample covariance matrix (already a list)

    if (nblocks == 1L) {
      param <- list(param)
      ve <- list(ve)
    }

  } else if (inherits(object, "lavaan.mi")) {
    useImps <- rep(TRUE, length(object@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

    param <- object@coefList[[ useImps[1] ]] # first admissible as template
    coefList <- object@coefList[useImps]
    phiList <- object@phiList[useImps]
    ## add block-level list per imputation?
    if (nblocks == 1L) {
      param <- list(param)
      for (i in 1:m) {
        coefList[[i]] <- list(coefList[[i]])
        phiList[[i]] <- list(phiList[[i]])
      }
    }

    S <- vector("list", nblocks) # pooled observed covariance matrix
    ve <- vector("list", nblocks)
    ## loop over blocks
    for (b in 1:nblocks) {

      ## param:  loop over GLIST elements
      for (mat in names(param[[b]])) {
        matList <- lapply(coefList, function(i) i[[b]][[mat]])
        param[[b]][[mat]] <- Reduce("+", matList) / length(matList)
      } # mat

      ## pooled observed covariance matrix
      covList <- lapply(object@h1List[useImps], function(i) i$implied$cov[[b]])
      S[[b]] <- Reduce("+", covList) / m

      ## pooled model-implied latent covariance matrix
      ve[[b]] <- Reduce("+", lapply(phiList, "[[", i = b) ) / m

    } # b

  }

  if (nblocks == 1L) {
    SigmaHat <- getMethod("fitted", class(object))(object)["cov"] # retain list format
  } else {
    SigmaHat <- sapply(getMethod("fitted", class(object))(object),
                       "[[", "cov", simplify = FALSE)
  }

  ly <- lapply(param, "[[", "lambda")
  te <- lapply(param, "[[", "theta")
  ps <- lapply(param, "[[", "psi")
	be <- lapply(param, "[[", "beta")

	result <- list()
	for (i in 1:nblocks) {

		# Prepare for higher-order reliability
		l2var <- ve[[i]][secondFactor, secondFactor, drop = FALSE]
		l2load <- be[[1]][,secondFactor]
		indexl2 <- which(l2load != 0)
		commonl2 <- (sum(l2load)^2) * l2var
		errorl2 <- sum(ps[[i]][indexl2, indexl2, drop = FALSE])

		# Prepare for lower-order reliability
		indexl1 <- which(apply(ly[[i]][,indexl2, drop = FALSE], 1, function(x) sum(x != 0)) > 0)
		l1load <- ly[[i]][,indexl2] %*% as.matrix(be[[1]][indexl2, secondFactor, drop = FALSE])
		commonl1 <- (sum(l1load)^2) * l2var
		errorl1 <- sum(te[[i]][indexl1, indexl1, drop = FALSE])
		uniquel1 <- 0
		for (j in seq_along(indexl2)) {
			uniquel1 <- uniquel1 + (sum(ly[[i]][,indexl2[j]])^2) * ps[[i]][indexl2[j], indexl2[j], drop = FALSE]
		}

		# Adjustment for direct loading from L2 to observed variables
		if (any(ly[[i]][,secondFactor] != 0)) {
			indexind <- which(ly[[i]][,secondFactor] != 0)
			if (length(intersect(indexind, indexl1)) > 0)
			  stop("Direct and indirect loadings of higher-order factor to observed",
			       " variables are specified at the same time.")
			commonl2 <- sum(c(ly[[i]][,secondFactor], l2load))^2 * l2var
			errorl2 <- errorl2 + sum(te[[i]][indexind, indexind, drop = FALSE])
			commonl1 <- sum(c(ly[[i]][,secondFactor], l1load))^2 * l2var
			errorl1 <- errorl1 + sum(te[[i]][indexind, indexind, drop = FALSE])
		}

		# Calculate Reliability
		omegaL1 <- commonl1 / (commonl1 + uniquel1 + errorl1)
		omegaL2 <- commonl2 / (commonl2 + errorl2)
		partialOmegaL1 <- commonl1 / (commonl1 + errorl1)
		result[[i]] <- c(omegaL1 = omegaL1, omegaL2 = omegaL2, partialOmegaL1 = partialOmegaL1)
	}

	if (nblocks == 1L) {
	  result <- result[[1]]
	} else names(result) <- block.label

	result
}



## --------------
## maximalRelia()
## --------------

##' Calculate maximal reliability
##'
##' Calculate maximal reliability of a scale
##'
##' Given that a composite score (\eqn{W}) is a weighted sum of item scores:
##'
##' \deqn{ W = \bold{w}^\prime \bold{x} ,}
##'
##' where \eqn{\bold{x}} is a \eqn{k \times 1} vector of the scores of each
##' item, \eqn{\bold{w}} is a \eqn{k \times 1} weight vector of each item, and
##' \eqn{k} represents the number of items. Then, maximal reliability is
##' obtained by finding \eqn{\bold{w}} such that reliability attains its maximum
##' (Li, 1997; Raykov, 2012). Note that the reliability can be obtained by
##'
##' \deqn{ \rho = \frac{\bold{w}^\prime \bold{S}_T \bold{w}}{\bold{w}^\prime
##' \bold{S}_X \bold{w}}}
##'
##' where \eqn{\bold{S}_T} is the covariance matrix explained by true scores and
##' \eqn{\bold{S}_X} is the observed covariance matrix. Numerical method is used
##' to find \eqn{\bold{w}} in this function.
##'
##' For continuous items, \eqn{\bold{S}_T} can be calculated by
##'
##' \deqn{ \bold{S}_T = \Lambda \Psi \Lambda^\prime,}
##'
##' where \eqn{\Lambda} is the factor loading matrix and \eqn{\Psi} is the
##' covariance matrix among factors. \eqn{\bold{S}_X} is directly obtained by
##' covariance among items.
##'
##' For categorical items, Green and Yang's (2009) method is used for
##' calculating \eqn{\bold{S}_T} and \eqn{\bold{S}_X}. The element \eqn{i} and
##' \eqn{j} of \eqn{\bold{S}_T} can be calculated by
##'
##' \deqn{ \left[\bold{S}_T\right]_{ij} = \sum^{C_i - 1}_{c_i = 1} \sum^{C_j -
##' 1}_{c_j - 1} \Phi_2\left( \tau_{x_{c_i}}, \tau_{x_{c_j}}, \left[ \Lambda
##' \Psi \Lambda^\prime \right]_{ij} \right) - \sum^{C_i - 1}_{c_i = 1}
##' \Phi_1(\tau_{x_{c_i}}) \sum^{C_j - 1}_{c_j - 1} \Phi_1(\tau_{x_{c_j}}),}
##'
##' where \eqn{C_i} and \eqn{C_j} represents the number of thresholds in Items
##' \eqn{i} and \eqn{j}, \eqn{\tau_{x_{c_i}}} represents the threshold \eqn{c_i}
##' of Item \eqn{i}, \eqn{\tau_{x_{c_j}}} represents the threshold \eqn{c_i} of
##' Item \eqn{j}, \eqn{ \Phi_1(\tau_{x_{c_i}})} is the cumulative probability of
##' \eqn{\tau_{x_{c_i}}} given a univariate standard normal cumulative
##' distribution and \eqn{\Phi_2\left( \tau_{x_{c_i}}, \tau_{x_{c_j}}, \rho
##' \right)} is the joint cumulative probability of \eqn{\tau_{x_{c_i}}} and
##' \eqn{\tau_{x_{c_j}}} given a bivariate standard normal cumulative
##' distribution with a correlation of \eqn{\rho}
##'
##' Each element of \eqn{\bold{S}_X} can be calculated by
##'
##' \deqn{ \left[\bold{S}_T\right]_{ij} = \sum^{C_i - 1}_{c_i = 1} \sum^{C_j -
##' 1}_{c_j - 1} \Phi_2\left( \tau_{V_{c_i}}, \tau_{V_{c_j}}, \rho^*_{ij}
##' \right) - \sum^{C_i - 1}_{c_i = 1} \Phi_1(\tau_{V_{c_i}}) \sum^{C_j -
##' 1}_{c_j - 1} \Phi_1(\tau_{V_{c_j}}),}
##'
##' where \eqn{\rho^*_{ij}} is a polychoric correlation between Items \eqn{i}
##' and \eqn{j}.
##'
##'
##' @importFrom lavaan lavInspect lavNames
##'
##' @param object A \code{\linkS4class{lavaan}} or
##'   \code{\linkS4class{lavaan.mi}} object, expected to contain only
##'   exogenous common factors (i.e., a CFA model).
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'        imputations from pooled results.  Can include any of
##'        \code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
##'        default setting, which excludes any imputations that did not
##'        converge or for which standard errors could not be computed.  The
##'        last option (\code{"no.npd"}) would exclude any imputations which
##'        yielded a nonpositive definite covariance matrix for observed or
##'        latent variables, which would include any "improper solutions" such
##'        as Heywood cases.  NPD solutions are not excluded by default because
##'        they are likely to occur due to sampling error, especially in small
##'        samples.  However, gross model misspecification could also cause
##'        NPD solutions, users can compare pooled results with and without
##'        this setting as a sensitivity analysis to see whether some
##'        imputations warrant further investigation.
##'
##' @return Maximal reliability values of each group. The maximal-reliability
##'   weights are also provided. Users may extracted the weighted by the
##'   \code{attr} function (see example below).
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' @seealso \code{\link{reliability}} for reliability of an unweighted
##'   composite score
##'
##' @references
##' Li, H. (1997). A unifying expression for the maximal reliability of a linear
##' composite. \emph{Psychometrika, 62}(2), 245--249. \doi{10.1007/BF02295278}
##'
##' Raykov, T. (2012). Scale construction and development using structural
##' equation modeling. In R. H. Hoyle (Ed.), \emph{Handbook of structural
##' equation modeling} (pp. 472--494). New York, NY: Guilford.
##'
##' @examples
##'
##' total <- 'f =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 '
##' fit <- cfa(total, data = HolzingerSwineford1939)
##' maximalRelia(fit)
##'
##' # Extract the weight
##' mr <- maximalRelia(fit)
##' attr(mr, "weight")
##'
##' @export
maximalRelia <- function(object, omit.imps = c("no.conv","no.se")) {
  ngroups <- lavInspect(object, "ngroups") #TODO: adapt to multiple levels
  nLevels <- lavInspect(object, "nlevels")
  nblocks <- ngroups*nLevels #FIXME: always true?
  group.label <- if (ngroups > 1L) lavInspect(object, "group.label") else NULL
  #FIXME? lavInspect(object, "level.labels")
  clus.label <- if (nLevels > 1L) c("within", lavInspect(object, "cluster")) else NULL
  if (nblocks > 1L) {
    block.label <- paste(rep(group.label, each = nLevels), clus.label,
                         sep = if (ngroups > 1L && nLevels > 1L) "_" else "")
  }

  ## parameters in GLIST format (not flat, need block-level list)
  if (inherits(object, "lavaan")) {
    param <- lavInspect(object, "est")
    ve <- lavInspect(object, "cov.lv") # model-implied latent covariance matrix
    S <- object@h1$implied$cov # observed sample covariance matrix (already a list)

    if (nblocks == 1L) {
      param <- list(param)
      ve <- list(ve)
    }

  } else if (inherits(object, "lavaan.mi")) {
    useImps <- rep(TRUE, length(object@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

    param <- object@coefList[[ useImps[1] ]] # first admissible as template
    coefList <- object@coefList[useImps]
    phiList <- object@phiList[useImps]
    ## add block-level list per imputation?
    if (nblocks == 1L) {
      param <- list(param)
      for (i in 1:m) {
        coefList[[i]] <- list(coefList[[i]])
        phiList[[i]] <- list(phiList[[i]])
      }
    }
    S <- vector("list", nblocks) # pooled observed covariance matrix
    ve <- vector("list", nblocks)
    ## loop over blocks
    for (b in 1:nblocks) {

      ## param:  loop over GLIST elements
      for (mat in names(param[[b]])) {
        matList <- lapply(coefList, function(i) i[[b]][[mat]])
        param[[b]][[mat]] <- Reduce("+", matList) / length(matList)
      } # mat

      ## pooled observed covariance matrix
      covList <- lapply(object@h1List[useImps], function(i) i$implied$cov[[b]])
      S[[b]] <- Reduce("+", covList) / m

      ## pooled model-implied latent covariance matrix
      ve[[b]] <- Reduce("+", lapply(phiList, "[[", i = b) ) / m

    } # b

  }

  if (nblocks == 1L) {
    SigmaHat <- getMethod("fitted", class(object))(object)["cov"] # retain list format
  } else {
    SigmaHat <- sapply(getMethod("fitted", class(object))(object),
                       "[[", "cov", simplify = FALSE)
  }

  ly <- lapply(param, "[[", "lambda")
  te <- lapply(param, "[[", "theta")

  categorical <- lavInspect(object, "categorical")
  threshold <- if (categorical) getThreshold(object, omit.imps = omit.imps) else NULL

  result <- list()
  for (i in 1:nblocks) {
    truevar <- ly[[i]] %*% ve[[i]] %*% t(ly[[i]])
    varnames <- colnames(truevar)
    if (categorical) {
      invstdvar <- 1 / sqrt(diag(SigmaHat[[i]]))
      polyr <- diag(invstdvar) %*% truevar %*% diag(invstdvar)
      nitem <- ncol(SigmaHat[[i]])
      result[[i]] <- calcMaximalReliaCat(polyr, threshold[[i]], S[[i]], nitem, varnames)
    } else {
      result[[i]] <- calcMaximalRelia(truevar, S[[i]], varnames)
    }
  }

  if (nblocks == 1L) {
    result <- result[[1]]
  } else names(result) <- block.label

  result
}



## ----------------
## Hidden Functions
## ----------------

computeAlpha <- function(S) {
  k <- nrow(S)
  k/(k - 1) * (1.0 - sum(diag(S)) / sum(S))
}

#' @importFrom stats cov2cor pnorm
omegaCat <- function(truevar, threshold, scales, denom) {
  ## must be in standardized latent scale
  R <- diag(scales) %*% truevar %*% diag(scales)

	## denom could be model-implied polychoric correlation assuming diagonal theta,
	##       model-implied polychoric correlation accounting for error covariances,
	##       or "observed" polychoric correlation matrix.
  ## If parameterization="theta", standardize the polychoric coVARIANCE matrix
	denom <- cov2cor(denom)

	nitem <- ncol(denom)
	## initialize sums of cumulative probabilities
	sumnum <- 0 # numerator
	addden <- 0 # denominator
	## loop over all pairs of items
	for (j in 1:nitem) {
  	for (jp in 1:nitem) {
  	  ## initialize sums of cumulative probabilities *per item*
  		sumprobn2 <- 0
  		addprobn2 <- 0
  		## for each pair of items, loop over all their thresholds
  		t1 <- threshold[[j]]  * scales[j] # on standardized latent scale
  		t2 <- threshold[[jp]] * scales[jp] #FIXME: subtract intercept (or marginal mean?)
  		for (c in 1:length(t1)) {
    		for (cp in 1:length(t2)) {
    			sumprobn2 <- sumprobn2 + p2(t1[c], t2[cp], R[j, jp])
    			addprobn2 <- addprobn2 + p2(t1[c], t2[cp], denom[j, jp])
    		}
  		}
  		sumprobn1 <- sum(pnorm(t1))
  		sumprobn1p <- sum(pnorm(t2))
  		sumnum <- sumnum + (sumprobn2 - sumprobn1 * sumprobn1p)
  		addden <- addden + (addprobn2 - sumprobn1 * sumprobn1p)
  	}
	}
	reliab <- sumnum / addden
	reliab
}


p2 <- function(t1, t2, r) {
	mnormt::pmnorm(c(t1, t2), c(0,0), matrix(c(1, r, r, 1), 2, 2))
}


# polycorLavaan <- function(object) {
# 	ngroups <- lavInspect(object, "ngroups")
# 	coef <- lavInspect(object, "est")
# 	targettaunames <- NULL
# 	if (ngroups == 1L) {
# 		targettaunames <- rownames(coef$tau)
# 	} else {
# 		targettaunames <- rownames(coef[[1]]$tau)
# 	}
# 	barpos <- sapply(strsplit(targettaunames, ""), function(x) which(x == "|"))
# 	varnames <- unique(apply(data.frame(targettaunames, barpos - 1), MARGIN = 1,
# 	                         FUN = function(x) substr(x[1], 1, x[2])))
# 	if (length(varnames))
# 	script <- ""
# 	for (i in 2:length(varnames)) {
# 		temp <- paste0(varnames[1:(i - 1)], collapse = " + ")
# 		temp <- paste0(varnames[i], "~~", temp, "\n")
# 		script <- paste(script, temp)
# 	}
# 	newobject <- refit(script, object)
# 	if (ngroups == 1L) {
# 		return(lavInspect(newobject, "est")$theta)
# 	}
# 	lapply(lavInspect(newobject, "est"), "[[", "theta")
# }

##' @importFrom lavaan lavInspect lavNames
getThreshold <- function(object, omit.imps = c("no.conv","no.se")) {
	ngroups <- lavInspect(object, "ngroups") #TODO: add nlevels when capable
	ordnames <- lavNames(object, "ov.ord")

	if (inherits(object, "lavaan")) {
	  EST <- lavInspect(object, "est")

	} else if (inherits(object, "lavaan.mi")) {
	  useImps <- rep(TRUE, length(object@DataList))
	  if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
	  if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
	  if ("no.npd" %in% omit.imps) {
	    Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
	    Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
	    useImps <- useImps & !(Heywood.lv | Heywood.ov)
	  }
	  m <- sum(useImps)
	  if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
	  useImps <- which(useImps)

	  EST <- object@coefList[useImps]
	}

	if (ngroups == 1L) {
	  if (inherits(object, "lavaan")) {
	    thresholds <- EST$tau[,"threshold"]
	  } else if (inherits(object, "lavaan.mi")) {
	    tauList <- lapply(EST, function(x) x$tau[,"threshold"])
	    thresholds <- Reduce("+", tauList) / length(tauList)
	  }

	  result <- lapply(ordnames,
	                   function(nn) thresholds[grepl(paste0(nn, "\\|"), names(thresholds))])
	  names(result) <- ordnames
	  ## needs to be within a list when called above within block-loops
	  result <- list(result)

	} else {

	  thresholds <- vector("list", ngroups)
	  for (g in 1:ngroups) {
	    if (inherits(object, "lavaan")) {
	      thresholds[[g]] <- EST[[g]]$tau[,"threshold"]
	    } else if (inherits(object, "lavaan.mi")) {
	      tauList <- lapply(EST, function(x) x[[g]]$tau[,"threshold"])
	      thresholds[[g]] <- Reduce("+", tauList) / length(tauList)
	    }

	  }

	  result <- list()
		group.label <- lavInspect(object, "group.label")

		for (g in 1:ngroups) {
		  result[[ group.label[g] ]] <- lapply(ordnames, function(nn) {
		    thresholds[[g]][ grepl(paste0(nn, "\\|"), names(thresholds[[g]])) ]
		  })
		  names(result[[ group.label[g] ]]) <- ordnames
		}

	}

	return(result)
}

##' @importFrom lavaan lavInspect lavNames
getScales <- function(object, omit.imps = c("no.conv","no.se")) {
  ngroups <- lavInspect(object, "ngroups") #TODO: add nlevels when capable
  ordnames <- lavNames(object, "ov.ord") #TODO: use to allow mix of cat/con vars

  if (inherits(object, "lavaan")) {
    EST <- lavInspect(object, "est")

  } else if (inherits(object, "lavaan.mi")) {
    useImps <- rep(TRUE, length(object@DataList))
    if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
      useImps <- useImps & !(Heywood.lv | Heywood.ov)
    }
    m <- sum(useImps)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

    EST <- object@coefList[useImps]
  }

  if (ngroups == 1L) {
    if (inherits(object, "lavaan")) {
      result <- list(EST$delta[,"scales"])
    } else if (inherits(object, "lavaan.mi")) {
      scales <- lapply(EST, function(x) x$delta[,"scales"])
      result <- list(Reduce("+", scales) / length(scales))
    }

  } else {

    result <- vector("list", ngroups)

    for (g in 1:ngroups) {
      if (inherits(object, "lavaan")) {
        result[[g]] <- EST[[g]]$delta[,"scales"]
      } else if (inherits(object, "lavaan.mi")) {
        scales <- lapply(EST, function(x) x[[g]]$delta[,"scales"])
        result[[g]] <- Reduce("+", scales) / length(scales)
      }
    }

  }

  return(result)
}

invGeneralRelia <- function(w, truevar, totalvar) {
	1 - (t(w) %*% truevar %*% w) / (t(w) %*% totalvar %*% w)
}

#' @importFrom stats pnorm
invGeneralReliaCat <- function(w, polyr, threshold, denom, nitem) {
	# denom could be polychoric correlation, model-implied correlation, or model-implied without error correlation
	upper <- matrix(NA, nitem, nitem)
	lower <- matrix(NA, nitem, nitem)
	for (j in 1:nitem) {
  	for (jp in 1:nitem) {
  		sumprobn2 <- 0
  		addprobn2 <- 0
  		t1 <- threshold[[j]]
  		t2 <- threshold[[jp]]
  		for (c in 1:length(t1)) {
  		for (cp in 1:length(t2)) {
  			sumprobn2 <- sumprobn2 + p2(t1[c], t2[cp], polyr[j, jp])
  			addprobn2 <- addprobn2 + p2(t1[c], t2[cp], denom[j, jp])
  		}
  		}
  		sumprobn1 <- sum(pnorm(t1))
  		sumprobn1p <- sum(pnorm(t2))
  		upper[j, jp] <- (sumprobn2 - sumprobn1 * sumprobn1p)
  		lower[j, jp] <- (addprobn2 - sumprobn1 * sumprobn1p)
  	}
	}
	1 - (t(w) %*% upper %*% w) / (t(w) %*% lower %*% w)
}

#' @importFrom stats nlminb
calcMaximalRelia <- function(truevar, totalvar, varnames) {
	start <- rep(1, nrow(truevar))
	out <- nlminb(start, invGeneralRelia, truevar = truevar, totalvar = totalvar)
	if (out$convergence != 0) stop("The numerical method for finding the maximal",
	                               " reliability did not converge.")
	result <- 1 - out$objective
	weight <- out$par / mean(out$par)
	names(weight) <- varnames
	attr(result, "weight") <- weight
	result
}

#' @importFrom stats nlminb
calcMaximalReliaCat <- function(polyr, threshold, denom, nitem, varnames) {
	start <- rep(1, nrow(polyr))
	out <- nlminb(start, invGeneralReliaCat, polyr = polyr, threshold = threshold,
	              denom = denom, nitem = nitem)
	if (out$convergence != 0) stop("The numerical method for finding the maximal",
	                               " reliability did not converge.")
	result <- 1 - out$objective
	weight <- out$par / mean(out$par)
	names(weight) <- varnames
	attr(result, "weight") <- weight
	result
}


