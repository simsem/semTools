### Terrence D. Jorgensen
###   - omegaCat() and deprecated functionality: Sunthud Pornprasertmanit
### Last updated: 14 February 2026



## -----------------------------
## AVE()
##   average variance extracted
##   (not a reliability index)
## -----------------------------

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
##' @param object A [lavaan::lavaan-class] or [lavaan.mi::lavaan.mi-class] object,
##'   expected to contain only exogenous common factors (i.e., a CFA model).
##'   Cross-loadings are not allowed and will result in `NA` for any factor with
##'   indicator(s) that cross-load.
##' @param obs.var `logical` indicating whether to compute AVE using
##'   observed variances in the denominator. Setting `FALSE` triggers
##'   using model-implied variances in the denominator.
##' @param omit.imps `character` vector specifying criteria for omitting
##'   imputations from pooled results.  Can include any of
##'   `c("no.conv", "no.se", "no.npd")`, the first 2 of which are the
##'   default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (`"no.npd"`) would exclude any imputations which
##'   yielded a nonpositive definite covariance matrix for observed or
##'   latent variables, which would include any "improper solutions" such
##'   as Heywood cases.  NPD solutions are not excluded by default because
##'   they are likely to occur due to sampling error, especially in small
##'   samples.  However, gross model misspecification could also cause
##'   NPD solutions, users can compare pooled results with and without
##'   this setting as a sensitivity analysis to see whether some
##'   imputations warrant further investigation.
##' @param omit.factors `character` vector naming any common factors
##'   modeled in `object` whose indicators' AVE is not of interest.
##' @param dropSingle `logical` indicating whether to exclude factors
##'   defined by a single indicator from the returned results. If `TRUE`
##'   (default), single indicators will still be included in the `total`
##'   column when `return.total = TRUE`.
##' @param return.df `logical` indicating whether to return reliability
##'   coefficients in a `data.frame` (one row per group/level), which is
##'   possible when every model block includes the same factors (after excluding
##'   those in `omit.factors` and applying `dropSingle`).
##'
##' @return `numeric` vector of average variance extracted from indicators
##'   per factor.  For models with multiple "blocks" (any combination of groups
##'   and levels), vectors may be returned as columns in a `data.frame`
##'   with additional columns indicating the group/level (see `return.df=`
##'   argument description for caveat).
##'
##' @author
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##' Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models
##' with unobservable variables and measurement errors. *Journal of
##' Marketing Research, 18*(1), 39--50. \doi{10.2307/3151312}
##'
##' @seealso [compRelSEM()] for composite reliability estimates
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
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

    ## common-factor variance
    phiList <- object@phiList[useImps]
    if (nblocks == 1L) for (i in 1:m) phiList[[i]] <- list(phiList[[i]])
    PHI <- list()
    for (b in 1:nblocks) {
      PHI[[b]] <- Reduce("+", lapply(phiList, "[[", i = b) ) / m
    }

    ## loadings
    LAMBDA <- vector("list", nblocks)
    if (nblocks == 1L) {
      lamList <- lapply(object@coefList[useImps], "[[", i = "lambda")
      LAMBDA[[1]] <- Reduce("+", lamList) / length(lamList)
    } else {
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
        SIGMA[[b]] <- Reduce("+", covList) / m
        rownames(SIGMA[[b]]) <- colnames(SIGMA[[b]]) <- lavNames(object, block = b)
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
##' subscripts.  In our package, we strive to provide unlabeled coefficients,
##' leaving it to the user to decide on a label in their report.  But we do
##' use the symbols \eqn{\alpha} and \eqn{\omega} in the formulas below in order
##' to distinguish coefficients that do (not) assume essential tau-equivalence.
##'
##' Bentler (1968) first introduced factor-analysis reliability for a
##' unidimensional factor model with congeneric indicators, labeling the
##' coefficients \eqn{\alpha}.  McDonald (1999) later referred to this
##' *and other reliability coefficients*, first as \eqn{\theta} (in 1970),
##' then as \eqn{\omega}, which is a source of confusion when reporting
##' coefficients (Cho, 2021).  Coefficients based on factor models were later
##' generalized to account for multidimenisionality (possibly with
##' cross-loadings) and correlated errors. The general \eqn{\omega} formula
##' implemented in this function is:
##'
##' \deqn{\omega=\frac{\bold{w}^{\prime} \Lambda \Phi \Lambda^{\prime} \bold{w}
##'                  }{ \bold{w}^{\prime} \hat{\Sigma} \bold{w} }, }
##'
##' where \eqn{\hat{\Sigma}} can be the model-implied covariance matrix from
##' either the saturated model (i.e., the "observed" covariance matrix, used by
##' default) or from the hypothesized CFA model, controlled by the `obs.var=`
##' argument. All elements of matrices in the numerator and denominator are
##' effectively summed by the multiplication of the outer terms \eqn{\bold{w}},
##' a \eqn{k}-dimensional vector of composite weights typically consisting of
##' \eqn{\bold{1}}s, unless otherwise specified with the `W=` argument), and
##' \eqn{k} is the number of variables in the composite. Reliability of subscale
##' composites (or simply for separate factors in a joint CFA) can be calculated
##' by setting omitted-indicator weights to 0.  For unidimensional constructs
##' with simple structure, the equation above is often simplified to a scalar
##' representation (e.g., McDonald, 1999, Eq. 6.20b):
##'
##' \deqn{ \omega = \frac{        \left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
##'   Var\left( \psi \right)  }{  \left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
##'   Var\left( \psi \right) + \sum^{k}_{i = 1} \theta_{ii} }, }
##'
##' Note that all coefficients are calculated from *total* factor variances:
##' `lavInspect(object, "cov.lv")`, which assumes the fitted `object=` is a CFA,
##' not a full SEM with latent regression slopes.  If there is a Beta matrix, it
##' should only contain higher-order factor loadings (see details below).
##'
##'
##' When the fitted CFA imposes constraints consistent with (essential)
##' tau-equivalence, \eqn{\omega} is equivalent to coefficient \eqn{\alpha}
##' (Cronbach, 1951):
##'
##' \deqn{ \alpha = \frac{k}{k - 1}\left[ 1 -
##'     \frac{ \textrm{tr} \left( \hat{\Sigma} \right)
##'         }{ \bold{1}^{\prime} \hat{\Sigma} \bold{1} }
##' \right],}
##'
##' where \eqn{\textrm{tr} \left( . \right)} is the trace operation (i.e., the
##' sum of diagonal elements). Setting `tau.eq=TRUE` triggers the application of
##' this formula (rather than \eqn{\omega} above) to the model-implied or
##' observed covariance matrix (again controlled by the `obs.var=` argument).
##'
##'
##' **Higher-Order Factors**:
##'
##' For higher-order constructs with latent indicators, only \eqn{\omega} is
##' available because \eqn{\alpha} was not derived from CFA parameters (although
##' it can be expressed in a particular restricted CFA specification).
##'
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
##' the proportion of variance of a composite score that is attributable to the
##' second-order factor:
##'
##' \deqn{\omega=\frac{\bold{w}^{\prime} \Lambda \bold{B} \Phi \bold{B}^{\prime}
##'   \Lambda^{\prime} \bold{w} }{ \bold{w}^{\prime} \hat{\Sigma} \bold{w}}, }
##'
##' where \eqn{\bold{w}}, \eqn{\hat{\Sigma}}, and \eqn{k} are defined as above.
##' **Note** that if a higher-order factor also has observed indicators, it is
##' necessary to model the observed indicators as single-indicator lower-order
##' constructs, so that all of the higher-order factor indicators are latent
##' (with loadings in the Beta matrix, not Lambda); otherwise, higher-order
##' factor variance in the observed indicator is not captured in the numerator.
##'
##'
##' **Bifactor or Multitrait--Multimethod (MTMM) Models**:
##'
##' These multidimensional models partition sources of common variance that are
##' due to the factor of interest (e.g., a trait) as well as non-target factors
##' (e.g., "method factors", such as item wording or type of respondent).
##' The latter can be considered as systematic (i.e., non-random) sources of
##' error, to be excluded from the numerator of a reliability coefficient,
##' yielding so-called "hierarchical omega" (\eqn{\omega_\textrm{H}}). On the
##' other hand, non-target variance that can be expected in repeated measurement
##' meets the classical test theory definition of reliability. Including method
##' factors in the numerator yields so-called "omega total"
##' (\eqn{\omega_\textrm{T}}), which is the default approach in `compRelSEM()`
##' because it is consistent with the classical test theory definition of
##' reliability. However, users can obtain \eqn{\omega_\textrm{H}} for a
##' composite by using the `true=` argument to specify any factor(s) to be
##' treated as representing true scores.  The same approach can be taken to
##' obtain the proportion of a (sub)scale composite's variance due to method
##' factors (by listing those in `true=`), if that is of interest.
##'
##' **Categorical Indicators**:
##'
##' When all indicators (per composite) are ordinal, a CFA can be fitted that
##' includes a threshold model (sometimes called Item Factor Analysis: IFA),
##' which assumes a normally distributed latent response underlies each observed
##' ordinal response.  Despite making this assumption, a composite of ordinal
##' items can only be calculated by assigning numerical values to the ordinal
##' categories, so that the pseudo-numerical variables can be summed into a
##' composite variable that is more approximately continuous than its items.
##'
##' Applying the formulas above to IFA parameters provides the
##' *hypothetical* reliability of a composite of latent responses: a composite
##' which cannot be calculated in practice.  Nonetheless, this hypothetical
##' reliability can be interpreted as an estimate of what reliability *could* be
##' if a more approximately continuous response scale were used (e.g., with
##' sufficiently many response categories that the standardized solutions are
##' equivalent between a fitted IFA and a fitted CFA that treats the ordinal
##' responses as numeric; Chalmers, 2018). This can be requested by setting
##' `ord.scale=FALSE`, in which case \eqn{\hat\Sigma} in the formulas above
##' is a *polychoric* correlation matrix.
##' When `ord.scale=FALSE` and `tau.eq=TRUE`, this results in what Zumbo et al.
##' (2007) termed "ordinal \eqn{\alpha}" (see criticisms by Chalmers, 2018, and
##' and a rejoinder by Zumbo & Kroc, 2019).
##'
##' Alternatively, Green and Yang (2009, Eq. 21) derived a method to calculate
##' model-based reliability (\eqn{\omega}) from IFA parameters (i.e.,
##' incorporating the latent-response assumption) but that applies to the actual
##' (i.e., ordinal) observed response scale (the default: `ord.scale=TRUE`).
##' Lu et al. (2020) showed how to incorporate unequal weights into Green and
##' Yang's (2009) formula, so `W=` can be used to estimate the (maximal)
##' reliability of a weighted composite of ordinal variables.
##' However, combining `ord.scale=TRUE` with `tau.eq=TRUE` is not available.
##' For \eqn{\alpha} to be interpretable on the observed ordinal scale,
##' users must choose whether to (a) release the latent-response assumption, by
##' fitting a CFA without a threshold model, or (b) fit an IFA model with
##' constraints consistent with the assumption of (essential) tau-equivalence
##' (i.e., equal factor loadings).
##'
##' No method analogous to Green and Yang (2009, Eq. 21) has yet been proposed
##' to calculate reliability with a mixture of categorical and continuous
##' indicators, so any such composite is skipped with a warning.
##'
##'
##' **Multilevel Measurement Models**:
##'
##' How to define reliability coefficients for scales employed in nested designs
##' is an ongoing topic of methodological development, with some ongoing
##' controversies about best practice when the target of measurement is the
##' "cluster" or between-level (i.e., Level 2 in a 2-level design).
##' Geldhof et al. (2014) proposed applying the standard formulas above to each
##' level's CFA parameters and/or (model-implied) covariance matrix, whereas
##' Lai (2021) proposed different formulas that account for all sources of
##' variance in composites of observed variables.
##'
##' There is no controversy about how to define a within-level reliability,
##' coefficient, which can be interpreted as the reliability of a composite
##' calculated by first centering each indicator around its cluster mean, then
##' calculating the composite from the cluster-mean-centered items. Equivalently
##' (i.e., the same formula), this can be interpreted as the *hypothetical*
##' reliability of a composite of the items' latent Level-1 components. This
##' coefficient can be requested with [lavaan::model.syntax] (to pass to the
##' `W=` argument) that specifies a composite in a Level-1 "block", which not
##' have the same name as any composite in the Level-2 block.  If users do not
##' use `W=` (i.e., calculate a reliability index per modeled common factor),
##' then this can be accomplished by using unique factor names across levels.
##'
##' This contrasts with reliability indices for between-level composites:
##' The reliability of a *hypothetical* composite of items' latent between-level
##' components (using formulas proposed by Geldhof et al., 2014) is **not**
##' equivalent to the coefficient for a composite of items' observed cluster
##' means, using generalizations of formulas proposed by Lai (2021):
##'
##' \deqn{ \omega^\textrm{B} =
##'     \frac{\bold{w}^{\prime} \Lambda^\textrm{B} \Phi^\textrm{B} \Lambda^{\textrm{B}\prime} \bold{w}
##'         }{ \bold{w}^{\prime} \hat{\Sigma}^\textrm{B} \bold{w} +
##'           \frac{1}{\tilde{n}_\textrm{clus}} \left(
##'             \bold{w}^{\prime} \hat{\Sigma}^\textrm{W} \bold{w} \right) }, }
##'
##' \deqn{ \alpha^\textrm{B} = \frac{2k}{k - 1}\left[
##'     \frac{ \sum^{k}_{i=2} \sum^{i-1}_{j=1} \hat\sigma^\textrm{B}_{ij}
##'         }{ \bold{1}^{\prime} \hat\Sigma^\textrm{B} \bold{1} +
##'            \frac{1}{\tilde{n}_\textrm{clus}} \left(
##'            \bold{1}^{\prime} \hat\Sigma^\textrm{W} \bold{1} \right) }
##'   \right],}
##'
##' where \eqn{\tilde{n}_\textrm{clus}} is the harmonic-mean cluster size, and
##' superscripts B and W indicate between- and within-level parameters.
##' Obtaining these estimates of composite reliability requires fitting a
##' 2-level CFA that provides the same factor structure and factor names in
##' the models at both levels (following the advice of Jak et al., 2021), as
##' well as the same composite name in both levels/blocks of syntax passed to
##' `W=` (if used).  Furthermore, the between-level composite name must be
##' passed to the `shared=` argument; otherwise, the same factor/composite name
##' across levels will yield Lai's (2021) coefficient for a configural construct
##' (see **Examples**):
##'
##' \deqn{ \omega^\textrm{2L} =
##'     \frac{\bold{w}^{\prime} \left(
##'           \Lambda^\textrm{W} \Phi^\textrm{W} \Lambda^{\textrm{W}\prime} +
##'           \Lambda^\textrm{B} \Phi^\textrm{B} \Lambda^{\textrm{B}\prime}
##'           \right) \bold{w}
##'         }{ \bold{w}^{\prime} \hat\Sigma^\textrm{B} \bold{w} +
##'            \bold{w}^{\prime} \hat\Sigma^\textrm{W} \bold{w} }, }
##'
##' \deqn{ \alpha^\textrm{2L} = \frac{2k}{k - 1}\left[
##'     \frac{ \sum^{k}_{i=2} \sum^{i-1}_{j=1} \left( \hat\sigma^\textrm{W}_{ij} +
##'                                            \hat\sigma^\textrm{B}_{ij} \right)
##'         }{ \bold{1}^{\prime} \hat\Sigma^\textrm{B} \bold{1} +
##'            \bold{1}^{\prime} \hat\Sigma^\textrm{W} \bold{1} }
##'   \right],}
##'
##' This can be interpreted as the scale-reliability coefficient ignoring the
##' nested design, as both the common-factor variance of the Level-1 factor
##' *and* of its Level-2 cluster means are treated as true-score variance.
##'
##' **Note** that Lai's (2021) between-level reliability coefficients for a
##' `shared` construct quantify generalizability across both indicators and
##' raters (i.e., subjects rating their cluster's construct).
##' Lüdtke et al. (2011) refer to these as measurement error and sampling error,
##' respectively.  From this perspective (and following from generalizability
##' theory), an IRR coefficient can also be calculated:
##'
##' \deqn{ \textrm{IRR} =
##'     \frac{\bold{w}^{\prime} \left( \hat{\Sigma}^\textrm{B} \right) \bold{w}
##'         }{ \bold{w}^{\prime} \hat\Sigma^\textrm{B} \bold{w} +
##'            \bold{w}^{\prime} \hat\Sigma^\textrm{W} \bold{w} }, }
##'
##' which quantifies generalizability across rater/sampling-error only, and can
##' be returned for any `shared=` construct's composite by setting `add.IRR=TRUE`.
##'
##'
##'
##'
##' @importFrom lavaan lavInspect lavNames parTable
##' @importFrom methods getMethod
##'
##' @param object A [lavaan::lavaan-class] or [lavaan.mi::lavaan.mi-class] object,
##'   expected to contain only exogenous common factors (i.e., a CFA model).
##' @param W Composite weights applied to observed variables prior to summing.
##'   By default (`NULL`), unit-weights are applied to all indicators per factor
##'   (as well as all modeled indicators when `return.total=TRUE`), which is
##'   equivalent to specifying equal weights of *any* value to each indicator.
##'   Weights can be a `character` string specifying any number of composites
##'   using [lavaan::model.syntax()], in the form `COMPOSITE <~ weight*indicator`
##'   (any indicator without a numeric `weight` is given a unit weight = 1).
##'   See **Details** and **Examples** about complicated CFAs (e.g., multilevel,
##'   higher-order, or bifactor).
##' @param return.total For multidimensional CFAs, this `logical` value
##'   indicates whether to return a final index for the reliability of a
##'   composite of all modeled indicators (labeled `.TOTAL.`). This is redundant
##'   whenever there is already a common factor indicated by all items (e.g.,
##'   the general factor in a bifactor model). This argument is ignored when
##'   using the `W=` argument to specify composites (optionally with weights).
##'   Setting a negative value (e.g., `-1`) returns **only** the `.TOTAL.`
##'   composite reliability (i.e., excluding coefficients per factor).
##' @param obs.var `logical` indicating whether to compute reliability
##'   using observed (co)variances to compute the denominator. Setting `FALSE`
##'   triggers using model-implied (co)variances to compute the denominator.
##' @param tau.eq `logical` indicating whether to assume (essential)
##'   tau-equivalence by calculating coefficient \eqn{\alpha} (on observed or
##'   model-implied (co)variances, depending on `obs.var=`).
##'   Triggers error if requested in combination with unequal weights in `W=`.
##'   Setting `FALSE` (default) yields an "\eqn{\omega}"-type coefficient.
##'   Optionally, a `character` vector of composite names can specify
##'   calculating coefficient \eqn{\alpha} for a subset of all composites.
##' @param ord.scale `logical` relevant only for composites of discrete items.
##'   Setting `TRUE` (default) applies Green and Yang's (2009, formula 21)
##'   method to calculate reliability of the actual composite (i.e., on the
##'   actual ordinal response scale).  Setting `FALSE` yields coefficients that
##'   are only interpretable on the continuous latent-response scale, which can
##'   be interpreted as the upper bound of reliability if items were more
##'   approximately continuous.
##'   Ignored for factors with continuous indicators.
##'   Reliability cannot currently be calculated for composites of both
##'   discrete and continuous indicators.
##' @param config Deprecated `character` vector.
#                                               naming any configural constructs
#      in a multilevel CFA. For these constructs (and optional total composite),
#      Lai's (2021) coefficients \eqn{\omega^\textrm{W}} and \eqn{\omega^\textrm{2L}}
#      are returned (or corresponding \eqn{\alpha} coefficients when
#      `tau.eq=TRUE`), rather than Geldhof et al.'s (2014) coefficients for
#      hypothetical composites of latent components (although the same formula
#      is used for \eqn{\omega^\textrm{W}} in either case). Note that the same name
#      must be used for the factor component represented at each level of the
#      model.
##' @param shared `character` vector of **composite names**, to be interpreted
##'   as representing (perhaps multidimensional) shared construct(s).
##'   Lai's (2021) coefficient \eqn{\omega^\textrm{B}} or \eqn{\alpha^\textrm{B}}
##'   is calculated to quantify reliability relative to error associated with
##'   both indicators (measurement error) and subjects (sampling error), like a
##'   generalizability coefficient. For purely *scale* reliability (relative
##'   to item/measurement error alone, i.e., Lai's \eqn{\omega^\textrm{2L}}),
##'   omit the composite(s) from the `shared=` argument.
##' @param add.IRR `logical` indicating whether to calculate an additional
##'   reliability coefficient for any composite listed in `shared=`. Given that
##'   subjects can be considered as raters of their cluster's shared construct,
##'   an interrater reliability (IRR) coefficient can quantify reliability
##'   relative to rater/sampling error alone.
##' @param higher Deprecated, supplanted by using the `true=` argument.
#   @param higher `character` vector naming any higher-order constructs in
#         `object` for which composite reliability should be calculated.
#         Ignored when `tau.eq=TRUE` because alpha is not based on a CFA model;
#         instead, users must fit a CFA with tau-equivalence constraints.
#         To obtain Lai's (2021) multilevel composite-reliability indices for a
#         higher-order factor, do not use this argument; instead, specify the
#         higher-order factor(s) using the `shared=` or `config=` argument
#         (`compRelSEM` will automatically check whether it includes latent
#         indicators and apply the appropriate formula).
##' @param true Optional `list` of `character` vectors, with list-element names
##'   corresponding to composite names.  Each composite can have a `character`
##'   vector with names of any common factor(s) that should be considered the
##'   source(s) of "true-score variance" in that composite. For any composite
##'   with a specification in `true=`, the default is to consider all common
##'   factors to contribute true-score variance to any items in the composite.
##'   Specifying a composite in `true=` is only necessary to deviate from this
##'   default, for example, to specify the "general" factor in a bifactor model,
##'   in order to obtain "hierarchical omega" (\eqn{\omega_\textrm{H}}).
##'   A shortcut for this is available when `W=NULL`, by specifying a single
##'   `character` string (one of `"omegaH"`, `"omega.h"`, or `"omega_h"`)
##'   instead of a `list`.
##' @param dropSingle When `W=NULL`, this `logical` indicates whether to exclude
##'   single-indicator factors from the list of default composites.
##'   Even when `TRUE` (default), single indicators are still included in
##'   the `.TOTAL.` composite when `return.total = TRUE`.
##' @param omit.factors Deprecated, supplanted by using the `true=` argument.
##' @param omit.indicators Deprecated, supplanted by using the `W=` argument.
##' @param omit.imps `character` vector specifying criteria for omitting
##'   imputations from pooled results (using [lavaan.mi::lavaan.mi-class]).
##'   Can include any of `c("no.conv", "no.se", "no.npd")`, the first 2 of which
##'   are the default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (`"no.npd"`) would exclude any imputations which
##'   yielded a nonpositive definite covariance matrix for observed or
##'   latent variables, which would include any "improper solutions" such
##'   as Heywood cases.  NPD solutions are not excluded by default because
##'   they are likely to occur due to sampling error, especially in small
##'   samples.  However, gross model misspecification could also cause
##'   NPD solutions.  Users can compare pooled results with and without
##'   this setting as a sensitivity analysis to see whether some
##'   imputations warrant further investigation.
##' @param return.df Deprecated `logical` argument, replaced by `simplify=`.
##' @param simplify `logical` indicating whether to return reliability
##'   coefficients in a `numeric` vector (for single-group model) or `data.frame`
##'   (one row per group, or per level in some cases).
##'   Specifying a negative number (`simplify = -1L`) additionally removes the
##'   informative headers printed to facilitate interpretation.
##'
##' @return
##' By default (`simplify=FALSE`) a `list` of `numeric` vectors (1 per
##' composite) is returned. In multigroup CFA, the vector contains a reliability
##' index for each group in which the composite can be computed.
##' Each composite's vector has a `attr(..., "header")` with information to
##' facilitate interpretation of that index:
##'         \itemize{
##'           \item{A list of variables in the composite, which determines the
##'                 composite's total variance (denominator of reliability)}
##'           \item{Whether that total variance (denominator) is determined from
##'                 the restricted model (i.e., CFA parameters) or unrestricted
##'                 model (i.e., a freely estimated covariance matrix)}
##'           \item{Whether the variables in the composite are (a transformation
##'                 of) observed variables, or whether they are *latent*
##'                 (components of) variables. The latter (e.g., latent responses
##'                 assumed to underlie observed ordinal indicators, or latent
##'                 level-specific components of variables in a multilevel CFA)
##'                 cannot be used to calculated an observed composite variable,
##'                 so the resulting coefficient should be cautiously interpreted
##'                 as a "hypothetical reliability" (Chalmers, 2018; Lai, 2021).}
##'           \item{The latent variables that contribute common-factor variance
##'                 to the composite, which determine the composite's
##'                 "true-score" variance (numerator of reliability)}
##'           \item{Which reliability formula was used: model-based reliability
##'                 (so-called "omega") or coefficient alpha (a model-free
##'                 lower-bound estimate of true reliability, equivalent to
##'                 a model-based reliability that assumes tau-equivalence)}
##'         }
##'         This header will be printed immediately above each composite's
##'         reliability coefficient.  When multiple reliability coefficients are
##'         returned, **and** each vector in the list has the same length, then
##'         setting `simplify=TRUE` will collect the list of *single*
##'         coefficients into a vector, or the list of *multiple* coefficients
##'         into a `data.frame`, and their headers will be concatenated to be
##'         printed above the coefficients.  Setting `simplify = -1L` (or any
##'         negative number) will omit the informative headers.
##'
##' @author
##' Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##'   Uses hidden functions to implement Green & Yang's (2009) reliability for
##'   categorical indicators, written by Sunthud Pornprasertmanit
##'   (\email{psunthud@@gmail.com}) for the deprecated `reliability()` function.
##'
##' @seealso
##' [maximalRelia()] for the maximal reliability of weighted composite
##'
##' @references
##'
##' Bentler, P. M. (1968). Alpha-maximized factor analysis (alphamax): Its
##' relation to alpha and canonical factor analysis. *Psychometrika, 33*(3),
##' 335--345. \doi{10.1007/BF02289328}
##'
##' Chalmers, R. P. (2018). On misconceptions and the limited usefulness of
##' ordinal alpha. *Educational and Psychological Measurement, 78*(6),
##' 1056--1071. \doi{10.1177/0013164417727036}
##'
##' Cho, E. (2021) Neither Cronbach’s alpha nor McDonald’s omega: A commentary
##' on Sijtsma and Pfadt. *Psychometrika, 86*(4), 877--886.
##' \doi{10.1007/s11336-021-09801-1}
##'
##' Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
##' tests. *Psychometrika, 16*(3), 297--334. \doi{10.1007/BF02310555}
##'
##' Geldhof, G. J., Preacher, K. J., & Zyphur, M. J. (2014). Reliability
##' estimation in a multilevel confirmatory factor analysis framework.
##' *Psychological Methods, 19*(1), 72--91. \doi{10.1037/a0032138}
##'
##' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
##' structural equation modeling: An alternative to coefficient alpha.
##' *Psychometrika, 74*(1), 155--167. \doi{10.1007/s11336-008-9099-3}
##'
##' Jak, S., Jorgensen, T. D., & Rosseel, Y. (2021). Evaluating cluster-level
##' factor models with `lavaan` and M*plus*. *Psych, 3*(2),
##' 134--152. \doi{10.3390/psych3020012}
##'
##' Lai, M. H. C. (2021). Composite reliability of multilevel data: It’s about
##' observed scores and construct meanings. *Psychological Methods, 26*(1),
##' 90--102. \doi{10.1037/met0000287}
##'
##' Lu, Z., Hong, M., & Kim, S. (2020). Formulas of multilevel reliabilities for
##' tests with ordered categorical responses.
##' In M. Wiberg, D. Molenaar,  J. González, U.Böckenholt, & J.-S. Kim (Eds.),
##' *Quantitative psychology: The 85th annual meeting of the Psychometric Society, Virtual*
##' (pp. 103--112). Springer. \doi{10.1007/978-3-030-74772-5_10}
##'
##' Lüdtke, O., Marsh, H. W., Robitzsch, A., & Trautwein, U. (2011).
##' A 2 \eqn{\times} 2 taxonomy of multilevel latent contextual models:
##' Accuracy--bias trade-offs in full and partial error correction models.
##' *Psychological Methods, 16*(4), 444--467. \doi{10.1037/a0024376}
##'
##' McDonald, R. P. (1999). *Test theory: A unified treatment*. Mahwah, NJ:
##' Erlbaum.
##'
##' Zumbo, B. D., Gadermann, A. M., & Zeisser, C. (2007). Ordinal versions of
##' coefficients alpha and theta for Likert rating scales.
##' *Journal of Modern Applied Statistical Methods, 6*(1), 21--29.
##' \doi{10.22237/jmasm/1177992180}
##'
##' Zumbo, B. D., & Kroc, E. (2019). A measurement is a choice and Stevens’
##' scales of measurement do not help make it: A response to Chalmers.
##' *Educational and Psychological Measurement, 79*(6), 1184--1197.
##' \doi{10.1177/0013164419844305}
##'
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
##' fit  <- cfa(HS.model, data = HS, ordered = c("y7","y8","y9"), std.lv = TRUE)
##' fitg <- cfa(HS.model, data = HS, ordered = c("y7","y8","y9"), std.lv = TRUE,
##'             group = "school")
##'
##' ## works for factors with exclusively continuous OR categorical indicators
##' compRelSEM(fit)
##' compRelSEM(fitg)
##'
##' ## reliability for composite of ALL indicators only available when they are
##' ## all continuous or all categorical.  The example below calculates a
##' ## composite of continuous items from 2 factors (visual and textual)
##' ## using the custom-weights syntax (note the "<~" operator)
##' w.tot <- '
##'   visual  <~ x1 + x2 + x3
##'   textual <~                x4 + x5 + x6
##'   total   <~ x1 + x2 + x3 + x4 + x5 + x6
##' '
##' compRelSEM(fit, W = w.tot)
##'
##'
##' ## ----------------------
##' ## Higher-order construct
##' ## ----------------------
##'
##' ## Reliability of a composite that represents a higher-order factor
##' mod.hi <- ' visual  =~ x1 + x2 + x3
##'             textual =~ x4 + x5 + x6
##'             speed   =~ x7 + x8 + x9
##'             general =~ visual + textual + speed '
##'
##' fit.hi <- cfa(mod.hi, data = HolzingerSwineford1939)
##' ## "general" is the factor representing "true scores", but it has no
##' ## observed indicators.  Must use custom-weights syntax:
##' compRelSEM(fit.hi, W = 'g <~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9')
##'
##'
##' ## ----------------------
##' ## Hierarchical omega
##' ## and omega Total
##' ## ----------------------
##'
##' mod.bi <- ' visual  =~ x1 + x2 + x3
##'             textual =~ x4 + x5 + x6
##'             speed   =~ x7 + x8 + x9
##'             general =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 '
##' fit.bi <- cfa(mod.bi, data = HolzingerSwineford1939,
##'               orthogonal = TRUE, std.lv = TRUE)
##' compRelSEM(fit.bi, return.total = -1) # omega_Total
##' compRelSEM(fit.bi, return.total = -1, # omega_Hierarchical
##'            true = list(.TOTAL. = "general"))
##'
##'
##' ## ----------------------
##' ## Multilevel Constructs
##' ## ----------------------
##'
##' ## Same factor structure with metric invariance across levels (Jak et al., 2021)
##' model2 <- '
##'   level: 1
##'     f1 =~ y1 + L2*y2 + L3*y3
##'     f2 =~ y4 + L5*y5 + L6*y6
##'   level: 2
##'     f1 =~ y1 + L2*y2 + L3*y3
##'     f2 =~ y4 + L5*y5 + L6*y6
##' '
##' fit2 <- sem(model2, data = Demo.twolevel, cluster = "cluster")
##'
##' ## Lai's (2021, Eq. 13) omega index for a configural (Level-1) construct,
##' ## treating common-factor variance at both levels as "true" variance
##' compRelSEM(fit2)
##'
##' ## Lai's (2021, Eq. 17) omega index for a shared (Level-2) construct
##' ## (also its interrater reliability coefficient)
##' compRelSEM(fit2, shared = c("f1","f2"), add.IRR = TRUE)
##'
##' ## Geldhof et al.'s (2014) level-specific indices imply a different
##' ## composite (hypothetically) calculated per level.  Thus, use
##' ## unique composite names per level.
##'
##' W2.Geldhof <- ' level: 1
##'   F1w <~ y1 + y2 + y3
##'   F2w <~ y4 + y5 + y6
##' level: 2
##'   F1b <~ y1 + y2 + y3
##'   F2b <~ y4 + y5 + y6
##' '
##' compRelSEM(fit2, W = W2.Geldhof)
##'
##'
##' @export
compRelSEM <- function(object, W = NULL,
                       return.total = FALSE,
                       obs.var = TRUE, tau.eq = FALSE, ord.scale = TRUE,
                       shared = character(0), config = character(0), # deprecated (only shared= needed)
                       add.IRR = FALSE, # only for shared constructs
                       higher = character(0), # deprecated (always used when found)
                       true = list(), # character(0) per composite
                       dropSingle = TRUE, # ignored when !is.null(W)
                       omit.factors = character(0), # deprecated (use true=)
                       omit.indicators = character(0), # deprecated (use W=)
                       omit.imps = c("no.conv","no.se"),
                       ## simplify= replaces (deprecated) return.df=
                       simplify = FALSE, return.df = simplify) {
  ## warn about deprecated arguments following introduction of W= and true=
  if (length(config)) {
    warning('Argument config= is deprecated.  Composites including both within- ',
            'and between-level variance are assumed by default to represent ',
            'configural constructs (unless listed in the shared= argument).')
  }
  if (length(higher)) {
    warning('Argument higher= is deprecated. Higher-order constructs are ',
            'automatically detected by nonzero elements in the Beta matrix, ',
            'and are always assumed to represent true scores, whereas the ',
            'residuals of its indicators (lower-order factors) are considered ',
            'sources of error.')
  }
  if (length(omit.indicators)) {
    warning('Argument omit.indicators= is deprecated. Specify custom composites ',
            'with W= argument.')
  }
  if (length(omit.factors)) {
    warning('Argument omit.factors= is deprecated. By default, all common ',
            'factors are assumed to represent true-score variance. To treat ',
            'a subset of factors as true-score variance (e.g., to calculate ',
            '"omega_hierarchical"), use the true= argument.')
  }
  if (!is.null(match.call()$return.df)) {
    message('Argument return.df= is deprecated, replaced by simplify= argument.')
    if (return.df) simplify <- -1L
  }


  ## numbers of blocks
  ngroups <- lavInspect(object, "ngroups")
  nLevels <- lavInspect(object, "nlevels")
  nblocks <- ngroups*nLevels #FIXME: always true?

  ## extract parameter table
  PT <- parTable(object)

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
    #       which can differ from unique(parTable(object)$level)
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

  ## Check (for) weights
  if (is.null(W)) {
    ## construct default weights
    wPT        <- PT[PT$op == "=~", ] # isolate factor definitions
    wPT$op     <- "<~"
    wPT$ustart <- 1L                  # default weights == 1

    ## add total composite?
    if (return.total) {
      ## loop over blocks to add a total score
      for (b in 1:nblocks) {
        wPTb         <- wPT[wPT$block == b, ]
        ov.names     <- lavNames(object, type = "ov.ind", block = b)
        ind.idx      <- unique(match(wPTb$rhs, table = ov.names))
        totalPTb     <- wPTb[ind.idx, ]
        totalPTb$lhs <- ".TOTAL."

        if (b == 1L) {
          totalPT <- totalPTb
        } else totalPT <- rbind(totalPT, totalPTb)
        ## end loop over blocks
      }

      if (return.total < 0L) {
        ## only the total composite
        wPT <- totalPT

        ## add the total composite
      } else wPT <- rbind(wPT, totalPT)

    }
    ## replace default integers with labels for groups
    if (length(group.label) && is.integer(wPT$group)) {
      wPT$group[wPT$group > 0L] <- group.label[wPT$group]
    }

    ## default composite names (each factor, and/or total)
    allComps   <- unique(wPT$lhs)

    ## Option to preserve old default behavior (partial numerator). Only
    ## available when we know composite names == factor names (W=NULL).
    if (is.character(true)) {
      if (tolower(true[1]) %in% paste0("omega", c("h",".h","_h")) ) {
        ## for each composite (factor), only that factor is true-score variance
        true <- setNames(as.list(allComps), nm = allComps)
      }
    }

  } else {
    ## user supplied W=

    ## turn off the dropSingle= argument
    dropSingle <- FALSE

    ## Take 2:  ONLY accept a lavaan script
    stopifnot(is.character(W))

    PTw <- lavaan::lavaanify(W,
                             ngroups = ngroups, #FIXME?
                             as.data.frame. = TRUE)
    ## Check for multilevel (if absent from W= syntax)
    if (is.null(PTw$level)  &&  nLevels > 1L) {
      ## repeat syntax per level (implies Lai's (2021) indices by default)
      PTw <- lavaan::lavaanify(c('level: 1 \n', W, 'level: 2 \n', W),
                               ngroups = ngroups, #FIXME?
                               as.data.frame. = TRUE)
    }

    if (length(group.label)) {
      ## replace default integers with labels
      PTw$group[PTw$group > 0L] <- group.label[PTw$group]
    }
    PTw$ustart[is.na(PTw$ustart)] <- 1L # replace missing weights with 1
    wPT      <- PTw[PTw$op == "<~", ]   # isolate composite definitions
    allComps <- unique(wPT$lhs)         # extract composite names
  }

  ## apply tau.eq= to all composites?
  if (is.logical(tau.eq)) {
    if (tau.eq) {
      tau.eq <- allComps
    } else tau.eq <- character(0)

  } else if (!is.character(tau.eq)) {
    stop('tau.eq= must be logical or a character vector of composite names')
  }


  ## check for categorical
  anyCategorical <- lavInspect(object, "categorical")
  threshold <- if (anyCategorical) getThreshold(object, omit.imps = omit.imps) else NULL
  latScales <- if (anyCategorical) getScales(object, omit.imps = omit.imps) else NULL

  ## extract common- and total-variance components
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
    SIGMA <- sapply(lavInspect(object, drop.list.single.group = FALSE,
                               what = ifelse(obs.var, "sampstat", "fitted")),
                    "[[", i = "cov", simplify = FALSE)
    names(SIGMA) <- block.label


  } else if (inherits(object, "lavaan.mi")) {
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
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps <- which(useImps)

    ## common-factor variance
    phiList <- object@phiList[useImps]
    if (nblocks == 1L) for (i in 1:m) phiList[[i]] <- list(phiList[[i]])
    PHI <- vector("list", nblocks)
    names(PHI) <- block.label
    for (b in 1:nblocks) {
      PHI[[b]] <- Reduce("+", lapply(phiList, "[[", i = b) ) / m
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
        LAMBDA[[b]] <- Reduce("+", lamList) / length(lamList)
        betList <- lapply(object@coefList[useImps], function(i) i[[b]]$beta  )
        BETA[[b]]   <- Reduce("+", betList) / length(betList)
      }
    }

    ## total variance
    if (obs.var) {
      ## pool model-implied SIGMA from h1 model
      SIGMA <- vector("list", nblocks)
      names(SIGMA) <- block.label
      ## loop over blocks to pool saturated-model (observed) matrices
      for (b in 1:nblocks) {
        covList <- lapply(object@h1List[useImps], function(i) i$implied$cov[[b]])
        SIGMA[[b]] <- Reduce("+", covList) / m
        ## The slot does not contain dimnames, so add them
        rownames(SIGMA[[b]]) <- colnames(SIGMA[[b]]) <- lavNames(object, block = b)
      }

    } else {
      ## pool model-implied SIGMA from h0 model
      if (nblocks == 1L) {
        SIGMA <- getMethod("fitted", class(object))(object)["cov"] # retain list format
      } else {
        SIGMA <- sapply(getMethod("fitted", class(object))(object),
                        "[[", "cov", simplify = FALSE)
        names(SIGMA) <- block.label
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

  ## flag conditions to warn about (listing problematic composites)
  warnLRV <- character(0)

  rel <- list() # coefficient(s) per composite
  ## loop over composites
  for (cc in allComps) {
    ## extract rows for this composite
    cPT <- wPT[wPT$lhs == cc, ]

    ## loop over groups
    for (g in 1:ngroups) {
      ## without labels, PT$group will be integers (even with 1 group)
      g.idx <- ifelse(length(group.label), yes = group.label[g], no = g)
      ## group-specific rows
      if (!is.null(cPT$group)) {
        cPTg <- cPT[cPT$group == g.idx, ]
      } else cPTg <- cPT # only 1 group


      ## check whether all weights are zero (so skip this composite)
      allWts0 <- all(sapply(cPTg$ustart, function(target) {
        isTRUE(all.equal(target, current = 0))
      }))
      if (allWts0) next

      ## Determine block indices (when multilevel)

      compositeHas2Levels <- FALSE
      ## Either the user didn't specify levels in W= ...
      if (nLevels > 1L && is.null(wPT$level)) compositeHas2Levels <- TRUE
      ## ... or they are in cPTg
      if (length(unique(cPTg$level)) > 1L) compositeHas2Levels <- TRUE

      ## Same composite (factor?) name across levels?
      isWithin  <- FALSE
      isBetween <- FALSE
      if (compositeHas2Levels) {
        ## Then we need to distinguish levels in the same index
        b.idx  <- (g - 1L)*nLevels + 1
        b.idx2 <- (g - 1L)*nLevels + 2
        #FIXME if cross-classification or nLevels > 2 are implemented

        ## Still multilevel?  Implies a level-specific name
      } else if (nLevels > 1L) {
        ## which level is this?
        isWithin  <- as.logical(unique(cPTg$block) %% 2) # is it an odd block?
        isBetween <- !isWithin
        b.idx  <- ifelse(isWithin,
                         yes = (g - 1L)*nLevels + 1, # within
                         no  = (g - 1L)*nLevels + 2) # between
        b.idx2 <- NULL
        #FIXME? eventually possible to model partial clustering?

      } else {
        ## SINGLE LEVEL model (block == group)
        ## but b.idx MUST be integer(s)
        b.idx  <- g # so NOT g.idx (which is character for MG-CFA)
        b.idx2 <- NULL
      }


      ## skip single-indicator construct?
      if (dropSingle && nrow(cPTg) == 1L) next
      ## same logic applies across 2 levels:
      if (isTRUE(unique(cPTg$level) > 1L)) {
        if (dropSingle && nrow(cPTg) == 2L) next
      }


      ## extract Lambda
      LY <- LAMBDA[[b.idx]]
      ## checks names of modeled indicators against each composite below
      allIndNames <- rownames(LY)
      allFacNames <- colnames(LY)

      ## assign default weights == 0
      wt <- setNames(rep(0, length(allIndNames)), nm = allIndNames)
      ## loop over indicators ...
      for (i in allIndNames) {
        ## ... to check whether to assign any nonzero weights
        if (isTRUE(cPTg$ustart[cPTg$rhs == i & cPTg$block == b.idx] != 0)) {
          wt[i] <- cPTg$ustart[cPTg$rhs == i & cPTg$block == b.idx]
        }
      }
      ## save names of indicators with nonzero weights (for header/footer)
      myIndNames <- names(wt[wt != 0])
      if (!length(myIndNames)) next

      ## verify equal weights for alpha
      non0wts <- wt[wt != 0]
      if (cc %in% tau.eq  &&  length(unique(non0wts)) > 1L) {
        stop('Cannot calculate coefficient alpha for composite `', cc,
             '` with unequal weights')
      }

      ## check factor names in true=
      if (length(true[[cc]])) {
        ## Isolate the factors considered "true" variance
        myFacNames <- intersect(true[[cc]], allFacNames)
      } else myFacNames <- allFacNames # defaults to all, but adjust below

      ## factor loadings for this composite:
      if (is.null(BETA[[b.idx]])) {
        theseRtrue <- apply(LY[myIndNames, myFacNames, drop = FALSE],
                            MARGIN = 2, FUN = function(x) !all(x == 0))
      } else {
        ## Assume higher-order constructs represent true scores
        #FIXME? check PT$op == "=~", compare to BETA predictors
        L1 <- LY
        L2 <- BETA[[b.idx]][ , myFacNames, drop = FALSE]
        LY <- L1 %*% L2

        theseRtrue <- apply(LY[myIndNames, myFacNames, drop = FALSE],
                            MARGIN = 2, FUN = function(x) !all(x == 0))
        if (!any(theseRtrue)) {
          ## No higher-order factor(s) identified for this composite.
          ## Just use lower-order Lambda to identify true-score variance
          LY <- L1
          theseRtrue <- apply(LY[myIndNames, myFacNames, drop = FALSE],
                              MARGIN = 2, FUN = function(x) !all(x == 0))
        }
      }

      ## perhaps no Level-1 factors for shared constructs,
      ## so only total variance at within-level necessary
      if ( !( any(theseRtrue)  |  cc %in% shared) ) {
        stop('No common factors identified for composite `', cc, '`')
      }
      myFacNames <- myFacNames[theseRtrue]

      L <- LY[, myFacNames, drop = FALSE]

      ## TRUE (co)variance among indicators in this composite:
      Phi <- PHI[[b.idx]][myFacNames, myFacNames, drop = FALSE]
      commonCov <- L %*% Phi %*% t(L)

      ## TOTAL (co)variance among indicators
      totalCov  <- SIGMA[[b.idx]]



      ## Second level too?  Do it all again!
      if (is.null(b.idx2)) {
        myIndNames2 <- myFacNames2 <- NULL
        commonCov2  <- totalCov2   <- NULL
      } else {
        ## extract Level-2 Lambda
        LY2 <- LAMBDA[[b.idx2]]
        ## checks names of modeled indicators against each composite below
        allIndNames2 <- rownames(LY2)
        allFacNames2 <- colnames(LY2)

        ## assign default weights == 0
        wt2 <- setNames(rep(0, length(allIndNames2)), nm = allIndNames2)
        ## loop over indicators ...
        for (i in allIndNames2) {
          ## ... to check whether to assign any nonzero weights
          if (isTRUE( cPTg$ustart[cPTg$rhs == i & cPTg$block == b.idx2] != 0)) {
            wt2[i] <- cPTg$ustart[cPTg$rhs == i & cPTg$block == b.idx2]
          }
        }
        ## save names of indicators with nonzero weights (for header/footer)
        myIndNames2 <- names(wt2[wt2 != 0])

        ## verify equal weights for alpha
        non0wts2 <- wt2[wt2 != 0]
        if (cc %in% tau.eq  &&  length(unique(non0wts2)) > 1L) {
          stop('Cannot calculate coefficient alpha for composite `', cc,
               '` with unequal weights')
        }
        ## verify weights match across levels
        for (i in names(wt)) {
          L1wt <- wt[i]
          L2wt <- wt2[i]
          if (is.na(L2wt)) next # Level-1 variable not decomposed?
          if (L1wt != L2wt) stop('Weights in W= do not match across levels for',
                                 ' indicator ', i, ' of composite `', cc, '`')
        }

        ## check factor names in true=
        if (length(true[[cc]])) {
          ## Isolate the factors considered "true" variance
          myFacNames2 <- intersect(true[[cc]], allFacNames2)
        } else myFacNames2 <- allFacNames2

        ## factor loadings for this composite:
        if (is.null(BETA[[b.idx2]])) {
          # L <- LY2[, myFacNames2, drop = FALSE]
          theseRtrue <- apply(LY2[myIndNames2, myFacNames2, drop = FALSE],
                              MARGIN = 2, FUN = function(x) !all(x == 0))
        } else {
          ## Assume higher-order constructs
          #FIXME? check PT$op == "=~", compare to BETA predictors
          L1  <- LY2
          L2  <- BETA[[b.idx2]][ , myFacNames2, drop = FALSE]
          LY2 <- L1 %*% L2

          theseRtrue <- apply(LY2[myIndNames2, myFacNames2, drop = FALSE],
                              MARGIN = 2, FUN = function(x) !all(x == 0))
          if (!any(theseRtrue)) {
            ## No higher-order factor(s) identified for this composite.
            ## Just use lower-order Lambda to identify true-score variance
            LY2 <- L1
            theseRtrue <- apply(LY2[myIndNames2, myFacNames2, drop = FALSE],
                                MARGIN = 2, FUN = function(x) !all(x == 0))
          }
        }

        if (!any(theseRtrue)) {
          stop('No common factors identified for composite `', cc, '`')
        }
        myFacNames2 <- myFacNames2[theseRtrue]

        L <- LY2[, myFacNames2, drop = FALSE]

        ## TRUE (co)variance among indicators in this composite:
        Phi <- PHI[[b.idx2]][myFacNames2, myFacNames2, drop = FALSE]
        commonCov2 <- L %*% Phi %*% t(L)

        ## TOTAL (co)variance among indicators
        totalCov2  <- SIGMA[[b.idx2]]
      }



      ## distinguish between categorical & continuous indicators
      nameArgs <- list(object = object)
      if (nblocks > 1L) nameArgs$block <- c(b.idx, b.idx2)
      ## in case lavNames returns a list (multiple blocks),
      ## Reduce() with union() to save unique names
      ordNames <- Reduce(union, do.call(lavNames, c(nameArgs, list(type = "ov.ord"))))
      numNames <- Reduce(union, do.call(lavNames, c(nameArgs, list(type = "ov.num"))))
      if (anyCategorical) {
        #FIXME? Add second block when ML-SEM is implemented for categorical

        ## Are ALL these indicators categorical?
        allCat <- all(myIndNames %in% ordNames)
        ## Are only SOME of these indicators categorical? (mixed indicators)
        mix <- any(myIndNames %in% ordNames) && any(myIndNames %in% numNames)
      } else {
        allCat <- FALSE
        mix <- FALSE
      }
      ## can't (yet) mix observed and latent scales
      if (mix) {
        warning('Reliability cannot be computed for composite `', cc,
               '` because it has both categorical and continuous indicators. ',
               'A model can be fitted by treating ordinal indicators as continuous.')
        next
      }



      ## Calculate RELIABILITY for composite(s) of ordinal indicators
      if (allCat && ord.scale) {

        if (cc %in% tau.eq) {
          stop('Setting ord.scale=TRUE and including composite `', cc,
               '` in tau.eq= indicates you want to calculate coefficient ',
               'alpha for a composite on the observed ordinal scale',
               ' (not hypothetical alpha on the latent scale).\n',
               'There are 2 options, which are not equivalent:\n\n(1)',
               ' You can refit the model without specifying variables as ',
               'ordered=, avoiding the latent-response assumption, so the model ',
               'parameters are on the observed scale, making it unnecessary to ',
               'set ord.scale=TRUE.\n\n(2)',
               ' You can fit a model that retains the latent-response assumption ',
               'but imposes tau-equivalence (e.g., set all loadings = 1 and ',
               'estimate all factor variances), making it unnecessary to ',
               'include composite `', cc, '` in tau.eq='
               )
        }

        ## Green & Yang (2009)
        relCat <- omegaCat(truevar   = commonCov[myIndNames, myIndNames],
                           denom     =  totalCov[myIndNames, myIndNames],
                           #FIXME? Where will thresholds for L2-only variables be?
                           threshold = threshold[[b.idx]][myIndNames],
                           scales    = latScales[[b.idx]][myIndNames],
                           wt        = wt[myIndNames])
        ## the same composite name NEVER yields level-specific coefficients,
        ## so only assign using group index
        rel[[cc]][g.idx] <- relCat
        next
      } else {
        ## all continuous or all LRV-scale?
        if (allCat) warnLRV <- union(warnLRV, cc)
      }

      ## else, calculate RELIABILITY for composite(s) of continuous indicators

      ## could be multilevel (not if categorical, calculated above)
      isShared <- FALSE
      isConfig <- FALSE

      ## calculate ALPHA?
      if (cc %in% tau.eq) {

        if (!is.null(b.idx2)) {
          ## calculate manually using Lai (2021, Eq. 23 or 24)
          nI <- length(myIndNames2)
          if (nI == 1L) {
            stop('Coefficient alpha is undefined for a single indicator. ',
                 'Set tau.eq=FALSE (or specify a subset of composites ',
                 'excluding single-indcator factors) or set dropSingle=TRUE')
          }
          kw <- nI / (nI-1) # weight for alpha based on number of items

          ## numerator of shared or configural coefficient
          onlyCov2 <- totalCov2
          diag(onlyCov2) <- 0

          ## additional information depends on shared or configural
          if (cc %in% shared) {
            isShared <- TRUE

            ## (reciprocal of) harmonic-mean cluster size
            N_bar <- mean(lavInspect(object, "cluster.size",
                                     drop.list.single.group = FALSE)[[g.idx]])
            Ns <- mean(1 / lavInspect(object, "cluster.size",
                                      drop.list.single.group = FALSE)[[g.idx]])

            ## numerator only sums Level-2 covariances
            numerator   <- t(wt2) %*% onlyCov2 %*% wt2
            denominator <- sum(Ns * (t(wt)  %*% totalCov  %*% wt),
                                     t(wt2) %*% totalCov2 %*% wt2)
            ## ALPHA
            rel[[cc]][g.idx] <- kw*numerator / denominator

            if (add.IRR) {
              numerator <- t(wt2) %*% totalCov2 %*% wt2 # same denominator
              rel[[paste0("IRR_of_", cc)]][g.idx] <- numerator / denominator
              ## Add header
              HEAD <- paste0('Composite `', cc, '` represents a shared construct, ',
                             'formed by summing not only over items but over ',
                             '(on average ', round(N_bar, 2), ') responses within ',
                             'each cluster (', dQuote(lavInspect(object, "cluster")),
                             '). Cluster members can therefore be interpreted as ',
                             '"raters" of the cluster-level shared construct.\n\n',
                             'The interrater reliability (IRR) of that construct is:')
              attr( rel[[paste0("IRR_of_", cc)]], "header") <- HEAD
              class(rel[[paste0("IRR_of_", cc)]]) <- c("lavaan.vector","numeric")
            }

          } else {
            isConfig <- TRUE
            ## numerator sums covariances at both levels
            onlyCov1 <- totalCov
            diag(onlyCov1) <- 0
            numerator   <- sum(t(wt)  %*% onlyCov1  %*% wt,
                               t(wt2) %*% onlyCov2  %*% wt2)
            denominator <- sum(t(wt)  %*% totalCov  %*% wt,
                               t(wt2) %*% totalCov2 %*% wt2)
            ## ALPHA
            rel[[cc]][g.idx] <- kw*numerator / denominator
          }

        } else {
          ## use built-in function for single-level model
          rel[[cc]][g.idx] <- computeAlpha(totalCov, W = wt)
        }

        next
      }

      ## else, calculate model-based RELIABILITY

      denominator <- t(wt) %*%  totalCov %*% wt

      if (!is.null(b.idx2)) {
        ## 2-level model, what kind of construct?

        if (cc %in% shared) {
          isShared <- TRUE
          ## (reciprocal of) harmonic-mean cluster size
          N_bar <- mean(lavInspect(object, "cluster.size",
                                   drop.list.single.group = FALSE)[[g.idx]])
          Ns <- mean(1 / lavInspect(object, "cluster.size",
                                    drop.list.single.group = FALSE)[[g.idx]])

          ## only Level 2 is true-score variance (Level-1 is error)
          numerator <- t(wt2) %*% commonCov2 %*% wt2
          ## Level-1 error is reduced by (harmonic-mean) cluster size
          denominator <- sum(denominator * Ns, t(wt2) %*%  totalCov2 %*% wt2)

        } else {
          isConfig <- TRUE
          ## Level 1 has true-score variance, but accumulates at Level 2
          numerator <- sum( t(wt ) %*% commonCov  %*% wt  ,
                            t(wt2) %*% commonCov2 %*% wt2 )
          denominator <- sum(denominator, t(wt2) %*%  totalCov2 %*% wt2)
        }

      } else {
        ## single block
        numerator <- t(wt) %*% commonCov %*% wt
      }

      ## OMEGA
      rel[[cc]][g.idx] <- as.numeric(numerator / denominator)


      if (isShared && add.IRR) {
        numerator   <- sum( t(wt2) %*%  totalCov2 %*% wt2)
        denominator <- sum( t(wt ) %*%  totalCov  %*% wt ) * Ns + numerator
        rel[[paste0("IRR_of_", cc)]][g.idx] <- numerator / denominator
        ## Add header
        HEAD <- paste0('Composite `', cc, '` represents a shared construct, ',
                       'formed by summing not only over items but over ',
                       '(on average ', round(N_bar, 2), ') responses within ',
                       'each cluster (', dQuote(lavInspect(object, "cluster")),
                       '). Cluster members can therefore be interpreted as ',
                       '"raters" of the cluster-level shared construct.\n\n',
                       'The interrater reliability (IRR) of composite `', cc,
                       '` is:')
        attr( rel[[paste0("IRR_of_", cc)]], "header") <- HEAD
        class(rel[[paste0("IRR_of_", cc)]]) <- c("lavaan.vector","numeric")
      }

      ## end loop over groups
    }

    ## in case this composite was skipped
    ## (e.g., no observed indicators in any groups)
    if (is.null(rel[[cc]])) next

    class(rel[[cc]]) <- c("lavaan.vector","numeric")

    ## determine type of composite for header
    compType <- ifelse(cc %in% warnLRV,
                       'latent responses underlying binary/ordinal variables:\n\t',
                ifelse(isWithin,
                       'cluster-mean-centered observed (or latent Level-1 components of) 2-level variables:\n\t',
                ifelse(isBetween,
                       'latent Level-2 components of 2-level variables:\n\t',
                ifelse(isShared,
                       '(cluster means of) observed variables:\n\t',
                       'observed variables:\n\t'))))

    ## Check for indicators ONLY in Level-2 model
    #TODO: enable this in lavaan:  myIndNames2only <- lavNames(object, "ov.between")
    #     if (isBetween && length(myIndNames2only)) {
    #       myIndNames2from1 <- setdiff(union(myIndNames, myIndNames2), myIndNames2only)
    #       ## assemble 2 indicator lists, separated by further description
    #       indicatorList <- paste0(paste(myIndNames2from1, collapse = ', '),
    #                               '\nand observed Level-2 variables:\n\t',
    #                               paste(myIndNames2only, collapse = ', '))
    #     } else
    indicatorList <- paste(union(myIndNames, myIndNames2), collapse = ', ')

    ## alpha or omega?
    relType <- ifelse(cc %in% tau.eq  &&  cc %in% warnLRV,
                      paste0('The latent polychoric correlation matrix was used ',
                             'to calculate a hypothetical coefficient alpha:'),
               ifelse(cc %in% tau.eq && isWithin,
                      paste0('The (latent) Level-1 covariance matrix was used ',
                             'to calculate coefficient alpha:'),
               ifelse(cc %in% tau.eq && isBetween,
                      paste0('The latent Level-2 covariance matrix was used ',
                             'to calculate a hypothetical coefficient alpha:'),
               ifelse(cc %in% tau.eq && isShared,
                      paste0('Coefficient alpha was calculated using ',
                             'Lai (2021, Eq. 24):'),
               ifelse(cc %in% tau.eq && isConfig,
                      paste0('Coefficient alpha was calculated using ',
                             'Lai (2021, Eq. 23):'),
               ifelse(cc %in% tau.eq,
                      paste0('Coefficient alpha would be:'),
               ## else OMEGA
               paste('The proportion attributable to "true" scores is its',
                     'model-based estimate of reliability ("omega"):')))))))
    trueVar <- ifelse(cc %in% tau.eq, yes = '',
                      paste0('\nTrue-score variance is represented by',
                             ifelse(isShared | isBetween,
                                    ' (between-level components of) ',
                                    ifelse(isWithin,
                                           ' (within-level components of) ',
                                           ' ')),
                             'common factor(s):\n\t',
                             paste(union(myFacNames, myFacNames2), collapse = ', ')))

    HEAD <- paste0('Composite `', cc,'` is composed of ',
                   compType, indicatorList,
                   ## alpha or omega?     ## What's the denominator?
                   trueVar,               '\nTotal variance of composite `', cc,
                   '`\ determined from the ',
                   ifelse(obs.var, 'un', ''), 'restricted model.\n',
                   ## alpha or omega?
                   relType)
    attr(rel[[cc]], "header") <- HEAD

    if (isConfig) {
      ## message about scale reliability?
      # attr(rel[[cc]], "footer") <-

    }

    ## end loop over composites
  }


  ## simplify structure at all? (NOTE: simplify < 0L reproduces old behavior)
  ## only if there are multiple composites
  if (simplify) {
    if (length(rel) == 1L) {
      REL <- rel[[1L]]
      ## drop header?
      if (simplify < 0L)  attr(REL, "header") <- NULL

    } else if (length(rel) > 1L) {

      ## Check there are just as many indices for each composite
      nComps <- unique(sapply(rel, length))
      if (length(nComps) > 1L) {
        message('Not all composites yield the same number of indices, so the ',
                'object (a list) cannot be simplified to a data.frame')
        REL <- rel
        ## drop headers?
        if (simplify < 0L)  for (i in seq_along(REL)) {
          attr(REL[[i]], "header") <- NULL
        } else simplify <- FALSE

      } else if (nComps == 1L) {
        ## Only 1 index per composite, so a single-group model.
        ## This automatically assigns list names to concatenated vector
        REL <- do.call(c, rel) # drops header (can replace below)
        class(REL) <- c("lavaan.vector","numeric")

      } else {
        ## Multiple indices per composite, so a multiple-group model

        ## Check the indices have the same (group) names.
        ## NOTE: The same composite name NEVER yields level-specific coefficients
        nameList <- lapply(rel, names)
        sameNames <- all(sapply(2:length(nameList), function(i) {
          isTRUE(all.equal(nameList[[1]], nameList[[i]]))
        } ))
        if (!sameNames) {
          ## can't simplify, return as a list
          message('Not all composites have the same names for indices, so the ',
                  'object (a list) cannot be simplified to a data.frame')
          REL <- rel
          ## drop headers?
          if (simplify < 0L)  for (i in seq_along(REL)) {
            attr(REL[[i]], "header") <- NULL
          } else simplify <- FALSE
        }
        ## bind composites into rows (groups in columns)
        REL <- as.data.frame(do.call(cbind, rel)) # drops header (can replace below)
        class(REL) <- c("lavaan.data.frame","data.frame")
      }

      if (simplify > 0L) {
        ## concatenate(?) headers separated by this:
        divL <- paste(c("\n\n", rep("-", 10)), collapse = "")
        divR <- paste(c(rep("-", 10), "\n\n"), collapse = "")

        allHeaders <- lapply(rel, attr, which = "header")
        divH <- paste0(divL, names(rel), divR, allHeaders, "\n\n", collapse = "")
        attr(REL, "header") <- paste0('Information about each composite is ',
                                      'provided below, followed by reliability ',
                                      'coefficients.', divH)

        # feet <- which(sapply(rel, function(x) !is.null(attr(x, "footer"))))
        #
        # allFooters <- lapply(rel, attr, which = "footer")
        # hasFeet <- sapply(allFooters, length) > 0L
        # if (any(hasFeet)) {
        #   divF <- paste0(divL, names(rel[hasFeet]), divR,
        #                  allFooters[hasFeet], "\n\n", collapse = "")
        #   attr(REL, "footer") <- paste0('Additional information is available for ',
        #                                 'the following composite(s):', divF)
        # }
      }

    }

    ## end if (simplify)
  } else REL <- rel

  return(REL)
}



## -------------
## reliability()
## (deprecated 10 May 2022)
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
##' item *i* observed variances, \eqn{\sigma_{ij}} is the observed
##' covariance of items *i* and *j*.
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
##' where \eqn{\lambda_i} is the factor loading of item *i*, \eqn{\psi} is
##' the factor variance, \eqn{\theta_{ii}} is the variance of measurement errors
##' of item *i*, and \eqn{\theta_{ij}} is the covariance of measurement
##' errors from item *i* and *j*. McDonald (1999) later referred to
##' this *and other reliability coefficients* as "omega", which is a source
##' of confusion when reporting coefficients (Cho, 2021).
##'
##' The additional coefficients generalize the first formula by accounting for
##' multidimenisionality (possibly with cross-loadings) and correlated errors.
##' By setting `return.total=TRUE`, one can estimate reliability for a
##' single composite calculated using all indicators in the multidimensional
##' CFA (Bentler, 1972, 2009).  `"omega2"` is calculated by
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
##' The `"omega3"` coefficient (McDonald, 1999), sometimes referred to as
##' hierarchical omega, can be calculated by
##'
##' \deqn{ \omega_3 =\frac{\left( \sum^{k}_{i = 1} \lambda_i \right)^{2}
##' Var\left( \psi \right)}{\bold{1}^\prime \Sigma \bold{1}}, }
##'
##' where \eqn{\Sigma} is the observed covariance matrix. If the model fits the
##' data well, \eqn{\omega_3} will be similar to \eqn{\omega_2}. Note that if
##' there is a directional effect in the model, all coefficients are calculated
##' from total factor variances: `lavInspect(object, "cov.lv")`.
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
##' (see Chalmers, 2018), which is consistent with the `alpha` function in
##' the `psych` package. When indicators are ordinal, `reliability`
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
##' Therefore, when `reliability` detects both categorical and continuous
##' indicators of a factor, an error is returned. If the categorical indicators
##' load on a different factor(s) than continuous indicators, then reliability
##' will still be calculated separately for those factors, but
##' `return.total` must be `FALSE` (unless `omit.factors` is used
##' to isolate factors with indicators of the same type).
##'
##'
##' @importFrom lavaan lavInspect lavNames
##' @importFrom methods getMethod
##'
##' @param object A [lavaan::lavaan-class] or [lavaan.mi::lavaan.mi-class] object,
##'   expected to contain only exogenous common factors (i.e., a CFA model).
##' @param what `character` vector naming any reliability indices to
##'   calculate. All are returned by default. When indicators are ordinal,
##'   both traditional `"alpha"` and Zumbo et al.'s (2007) so-called
##'   "ordinal alpha" (`"alpha.ord"`) are returned, though the latter is
##'   arguably of dubious value (Chalmers, 2018).
##' @param return.total `logical` indicating whether to return a final
##'   column containing the reliability of a composite of all indicators (not
##'   listed in `omit.indicators`) of factors not listed in
##'   `omit.factors`.  Ignored in 1-factor models, and should only be set
##'   `TRUE` if all factors represent scale dimensions that could be
##'   meaningfully collapsed to a single composite (scale sum or scale mean).
##' @param dropSingle `logical` indicating whether to exclude factors
##'   defined by a single indicator from the returned results. If `TRUE`
##'   (default), single indicators will still be included in the `total`
##'   column when `return.total = TRUE`.
##' @param omit.factors `character` vector naming any common factors
##'   modeled in `object` whose composite reliability is not of
##'   interest. For example, higher-order or method factors. Note that
##'   [reliabilityL2()] should be used to calculate composite
##'   reliability of a higher-order factor.
##' @param omit.indicators `character` vector naming any observed variables
##'   that should be ignored when calculating composite reliability. This can
##'   be useful, for example, to estimate reliability when an indicator is
##'   removed.
##' @param omit.imps `character` vector specifying criteria for omitting
##'   imputations from pooled results.  Can include any of
##'   `c("no.conv", "no.se", "no.npd")`, the first 2 of which are the
##'   default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (`"no.npd"`) would exclude any imputations which
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
##'   factors, a `total` column can optionally be included.
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##'   Yves Rosseel (Ghent University; \email{Yves.Rosseel@@UGent.be})
##'
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##' Bentler, P. M. (1972). A lower-bound method for the dimension-free
##' measurement of internal consistency. *Social Science Research, 1*(4),
##' 343--357. \doi{10.1016/0049-089X(72)90082-8}
##'
##' Bentler, P. M. (2009). Alpha, dimension-free, and model-based internal
##' consistency reliability. *Psychometrika, 74*(1), 137--143.
##' \doi{10.1007/s11336-008-9100-1}
##'
##' Chalmers, R. P. (2018). On misconceptions and the limited usefulness of
##' ordinal alpha. *Educational and Psychological Measurement, 78*(6),
##' 1056--1071. \doi{10.1177/0013164417727036}
##'
##' Cho, E. (2021) Neither Cronbach’s alpha nor McDonald’s omega: A commentary
##' on Sijtsma and Pfadt. *Psychometrika, 86*(4), 877--886.
##' \doi{10.1007/s11336-021-09801-1}
##'
##' Cronbach, L. J. (1951). Coefficient alpha and the internal structure of
##' tests. *Psychometrika, 16*(3), 297--334. \doi{10.1007/BF02310555}
##'
##' Fornell, C., & Larcker, D. F. (1981). Evaluating structural equation models
##' with unobservable variables and measurement errors. *Journal of
##' Marketing Research, 18*(1), 39--50. \doi{10.2307/3151312}
##'
##' Green, S. B., & Yang, Y. (2009). Reliability of summed item scores using
##' structural equation modeling: An alternative to coefficient alpha.
##' *Psychometrika, 74*(1), 155--167. \doi{10.1007/s11336-008-9099-3}
##'
##' McDonald, R. P. (1999). *Test theory: A unified treatment*. Mahwah, NJ:
##' Erlbaum.
##'
##' Raykov, T. (2001). Estimation of congeneric scale reliability using
##' covariance structure analysis with nonlinear constraints *British
##' Journal of Mathematical and Statistical Psychology, 54*(2), 315--323.
##' \doi{10.1348/000711001159582}
##'
##' Zumbo, B. D., Gadermann, A. M., & Zeisser, C. (2007). Ordinal versions of
##' coefficients alpha and theta for Likert rating scales.
##' *Journal of Modern Applied Statistical Methods, 6*(1), 21--29.
##' \doi{10.22237/jmasm/1177992180}
##'
##' Zumbo, B. D., & Kroc, E. (2019). A measurement is a choice and Stevens’
##' scales of measurement do not help make it: A response to Chalmers.
##' *Educational and Psychological Measurement, 79*(6), 1184--1197.
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
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL


##' @rdname semTools-deprecated
##' @section Reliability:
##' The original `reliability` function was suboptimally designed.
##' For example, AVE was returned, which is not a reliability index. Also,
##' alpha and several omega-type coefficients were returned, including the
##' original formula that was inappropriate for models with complex structure.
##' Some features could be controlled by the user for one but not both types of
##' index  For example, alpha for categorical indicators was returned on both
##' the observed and latent-response scales, but this was not an option for any
##' omega-type indices.  The omegas differed in terms of whether the observed or
##' model-implied covariance matrix was used in the denominator, but alpha was
##' only computed using the observed matrix.  These inconsistencies have been
##' resolved in the new [compRelSEM()] function, which returns only
##' one reliability index per composite, tailored to the user's requested
##' features, for which there is much more flexibility. The average variance
##' extracted is now available in a dedicated [AVE()] function.
##'
##' @export
reliability <- function(object,
                        what = c("alpha","omega","omega2","omega3","ave"),
                        return.total = FALSE, dropSingle = TRUE,
                        omit.factors = character(0),
                        omit.indicators = character(0),
                        omit.imps = c("no.conv","no.se")) {

  .Deprecated(msg = c("\nThe reliability() function was deprecated in 2022 and ",
                      "will cease to be included in future versions of semTools",
                      ". See help('semTools-deprecated) for details.",
                      "\n\nIt is replaced by the compRelSEM() function, which ",
                      "can estimate alpha and model-based reliability in an ",
                      "even wider variety of models and data types, with ",
                      "greater control in specifying the desired type of ",
                      "reliability coefficient (i.e., more explicitly choosing ",
                      "assumptions). \n\nThe average variance extracted should ",
                      "never have been included because it is not a reliability ",
                      "coefficient. It is now available from the AVE() function."))

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
## (deprecated 10 May 2022)
## ---------------

##' Calculate the reliability values of a second-order factor
##'
##' Calculate the reliability values (coefficient omega) of a second-order
##' factor
##'
##' The first formula of the coefficient omega (in the
##' [reliability()]) will be mainly used in the calculation. The
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
##' where \eqn{\bold{1}} is the *k*-dimensional vector of 1 and *k* is
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
##' \deqn{ \omega_{L2} = \frac{\bold{1}_F^{\prime} \bold{B} \Phi_2
##' \bold{B}^{\prime} \bold{1}_F}{\bold{1}_F^{\prime} \bold{B} \Phi_2
##' \bold{B}^{\prime} \bold{1}_F + \bold{1}_F^{\prime} \Psi_{u} \bold{1}_F}, }
##'
##' where \eqn{\bold{1}_F} is the *F*-dimensional vector of 1 and *F*
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
##' @param object A [lavaan::lavaan-class] or [lavaan.mi::lavaan.mi-class] object,
##'   expected to contain a least one exogenous higher-order common factor.
##' @param secondFactor The name of a single second-order factor in the
##'   model fitted in `object`. The function must be called multiple
##'   times to estimate reliability for each higher-order factor.
##' @param omit.imps `character` vector specifying criteria for omitting
##'        imputations from pooled results.  Can include any of
##'        `c("no.conv", "no.se", "no.npd")`, the first 2 of which are the
##'        default setting, which excludes any imputations that did not
##'        converge or for which standard errors could not be computed.  The
##'        last option (`"no.npd"`) would exclude any imputations which
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
##' @seealso [semTools-deprecated()]
##' @keywords internal
NULL


##' @rdname semTools-deprecated
##' @section Higher-Order Reliability:
##' Originally, composite reliability of a single higher-order factor was
##' estimated in a separate function: `reliabilityL2`.  It is now available
##' for multiple higher-order factors in the same model, and from the same
##' [compRelSEM()] function that estimates reliability for first-order
##' factors, using the `higher=` argument to name higher-order factors in
##' the fitted model.
##'
##' @export
reliabilityL2 <- function(object, secondFactor,
                          omit.imps = c("no.conv","no.se")) {

  .Deprecated(msg = c("\nThe reliabilityL2() function was deprecated in 2022 and ",
                      "will cease to be included in future versions of semTools",
                      ". See help('semTools-deprecated) for details.",
                      "\n\nIt is replaced by the compRelSEM() function, which ",
                      "can estimate alpha and model-based reliability in an ",
                      "even wider variety of models and data types, with ",
                      "greater control in specifying the desired type of ",
                      "reliability coefficient (i.e., more explicitly choosing ",
                      "assumptions)."))

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
##' @param object A [lavaan::lavaan-class] or [lavaan.mi::lavaan.mi-class] object,
##'   expected to contain only exogenous common factors (i.e., a CFA model).
##' @param omit.imps `character` vector specifying criteria for omitting
##'        imputations from pooled results.  Can include any of
##'        `c("no.conv", "no.se", "no.npd")`, the first 2 of which are the
##'        default setting, which excludes any imputations that did not
##'        converge or for which standard errors could not be computed.  The
##'        last option (`"no.npd"`) would exclude any imputations which
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
##'   `attr` function (see example below).
##'
##' @author Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
##'
##' @seealso [reliability()] for reliability of an unweighted
##'   composite score
##'
##' @references
##' Li, H. (1997). A unifying expression for the maximal reliability of a linear
##' composite. *Psychometrika, 62*(2), 245--249. \doi{10.1007/BF02295278}
##'
##' Raykov, T. (2012). Scale construction and development using structural
##' equation modeling. In R. H. Hoyle (Ed.), *Handbook of structural
##' equation modeling* (pp. 472--494). New York, NY: Guilford.
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

computeAlpha <- function(S, W = NULL) {
  k <- nrow(S)

  if (is.null(W)) {
    ## Traditional formula
    ALPHA <- k/(k - 1) * (1 - sum(diag(S)) / sum(S))

  } else {
    stopifnot(length(W) == nrow(S))

    #TODO? Develop theory for unequal weights.
    #      For now, force equal weights (when nonzero).
    wt <- W != 0
    k <- sum(wt)
    #FIXME? Should k = sum(wt) even with unequal weights?

    DIAG <- t(wt) %*% diag(diag(S)) %*% wt
    ALL  <- t(wt) %*%           S   %*% wt
    ALPHA <- k/(k - 1) * (1 - DIAG / ALL)[1,1]
  }

  ALPHA
}

##' @importFrom stats cov2cor pnorm
omegaCat <- function(truevar, threshold, scales, denom, wt = 1) {
  #TODO: How to incorporate varying distances between category weights
  #      (e.g., 0 = never, 0.5 = <1, 1.5 = 1-2, 4 = 3-5, 6 = >5)
  #TODO: How to adapt for composites of continuous & categorical items?

  ## must be in standardized latent scale
  R <- diag(scales) %*% truevar %*% diag(scales)

	## denom could be model-implied polychoric correlation assuming diagonal theta,
	##       model-implied polychoric correlation accounting for error covariances,
	##       or "observed" polychoric correlation matrix.
  ## If parameterization="theta", standardize the polychoric coVARIANCE matrix
	denom <- cov2cor(denom)

	nitem <- ncol(denom)
	if (length(wt) == 1L) {
	  wt <- rep(wt, nitem)
	}
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
  		t1 <- threshold[[j ]] * scales[j ] # on standardized latent scale
  		t2 <- threshold[[jp]] * scales[jp] #FIXME? subtract intercept (or marginal mean?)
  		for (c in 1:length(t1)) {
    		for (cp in 1:length(t2)) {
    			sumprobn2 <- sumprobn2 + p2(t1[c], t2[cp], R[j, jp])
    			addprobn2 <- addprobn2 + p2(t1[c], t2[cp], denom[j, jp])
    		}
  		}
  		sumprobn1  <- sum(pnorm(t1))
  		sumprobn1p <- sum(pnorm(t2))
  		## Add item weights (see Lu et al., 2020, doi:10.1037/met0000287)
  		sumnum <- sumnum + wt[j] * wt[jp] * (sumprobn2 - sumprobn1 * sumprobn1p)
  		addden <- addden + wt[j] * wt[jp] * (addprobn2 - sumprobn1 * sumprobn1p)
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


