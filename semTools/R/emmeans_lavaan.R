### Mattan S. Ben-Shachar
### Last updated: 27 May 2020
### emmeans support for lavaan objects


##' \code{emmeans} Support Functions for \code{lavaan} Models
##'
##' @description Provide emmeans support for lavaan objects
##'
##' @param object An object of class \code{\link[lavaan]{lavaan}}.
##'   See \strong{Details}.
##' @param lavaan.DV \code{character} string maming the variable(s) for which
##'   expected marginal means / trends should be produced.
##'   A vector of names indicates a multivariate outcome, treated by default
##'   as repeated measures.
##' @param trms,xlev,grid See \code{emmeans::emm_basis}
##' @param ... Further arguments passed to \code{emmeans::recover_data.lm} or
##'   \code{emmeans::emm_basis.lm}
##'
##' @details
##'
##' \subsection{Supported DVs}{
##'   \code{lavaan.DV} must be an \emph{endogenous variable}, by appearing on
##'   the left-hand side of either a regression operator (\code{"~"})
##'   or an intercept operator (\code{"~1"}), or both.
##'   \cr\cr
##'   \code{lavaan.DV} can also be a vector of endogenous variable, in which
##'   case they will be treated by \code{emmeans} as a multivariate outcome
##'   (often, this indicates repeated measures) represented by an additional
##'   factor named \code{rep.meas} by default.  The \code{rep.meas=} argument
##'   can be used to overwrite this default name.
##' }
##'
##' \subsection{Unsupported Models}{
##'   This functionality does not support the following models:
##'   \itemize{
##'     \item Multi-level models are not supported.
##'     \item Models not fit to a \code{data.frame} (i.e., models fit to a
##'           covariance matrix).
##'   }
##' }
##'
##' \subsection{Dealing with Fixed Parameters}{
##' Fixed parameters (set with \code{lavaan}'s modifiers) are treated as-is:
##' their values are set by the users, and they have a \emph{SE} of 0 (as such,
##' they do not co-vary with any other parameter).
##' }
##'
##' \subsection{Dealing with Multigroup Models}{
##' If a multigroup model is supplied, a factor is added to the reference grid,
##' the name matching the \code{group} argument supplied when fitting the model.
##' \emph{Note that you must set} \code{nesting = NULL}.
##' }
##'
##' \subsection{Dealing with Missing Data}{
##' Limited testing suggests that these functions do work when the model was fit
##' to incomplete data.
##' }
##'
##' \subsection{Dealing with Factors}{
##' By default \code{emmeans} recognizes binary variables (0,1) as a "factor"
##' with two levels (and not a continuous variable). With some clever contrast
##' defenitions it should be possible to get the desired emmeans / contasts.
##' See example below.
##' }
##'
##' @author Mattan S. Ben-Shachar (Ben-Gurion University of the Negev;
##'   \email{matanshm@@post.bgu.ac.il})
##'
##' @example inst/examples/lavaan2emmeans.R
##'
##' @name lavaan2emmeans
NULL

##' @rdname lavaan2emmeans
recover_data.lavaan <- function(object, lavaan.DV, ...){
  if (!requireNamespace("emmeans", quietly = TRUE)){
    stop("'emmeans' is not installed.")
  }

  .emlav_test_DV(object, lavaan.DV)
  ## testing multi-group requires access to ...
  dots <- list(...)
  if (lavInspect(object, 'ngroups') > 1L && !("nesting" %in% names(dots))) {
    warning(
      "For multi-group models, don't forget to set 'nesting = NULL'.\n",
      "See `?lavaan2emmeans` for more info.",
      call. = FALSE
    )
  }


  # Fake it
  recovered <- emmeans::recover_data(.emlav_fake_fit(object, lavaan.DV),
                                     ...)

  # Make it
  lavaan_data <- .emlav_recover_data(object)
  lavaan_data <- lavaan_data[, colnames(recovered), drop = FALSE]

  # Fill attributes (but keep lavaan_data in case of missing data)
  mostattributes(lavaan_data) <- attributes(recovered)
  return(lavaan_data)
}

##' @rdname lavaan2emmeans
##' @importFrom lavaan lavInspect
emm_basis.lavaan <- function(object,trms, xlev, grid, lavaan.DV, ...){
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    stop("'emmeans' is not installed.")
  }

  # Fake it
  emmb <- emmeans::emm_basis(.emlav_fake_fit(object, lavaan.DV),
                             trms, xlev, grid, ...)

  # bhat --------------------------------------------------------------------
  pars <- .emlav_clean_pars_tab(object, lavaan.DV, "bhat")
  par_names <- pars$rhs

  if(nrow(pars) < length(emmb$bhat)) {
    warning(
      "Not all parameters have been estimated.\n",
      "This is usually caused by a missing mean structure.\n",
      "Fixing estimates for these parameters at 0.",
      call. = FALSE
    )
  }

  # re-shape to deal with any missing estimates
  temp_bhat <- rep(0, length = length(emmb$bhat))
  temp_bhat[seq_len(nrow(pars))] <- pars$est
  names(temp_bhat) <- c(par_names,
                        colnames(emmb$V)[!colnames(emmb$V) %in% par_names])

  # re-order
  b_ind <- match(colnames(emmb$V), names(temp_bhat))
  emmb$bhat <- temp_bhat[b_ind]


  # VCOV --------------------------------------------------------------------
  lavVCOV <- lavInspect(object, "vcov")
  pars <- .emlav_clean_pars_tab(object, lavaan.DV, "vcov")
  par_names <- paste0(pars$lhs, pars$op, pars$rhs)

  # only regression estimates
  pattern <- paste0("^(", paste0(lavaan.DV, collapse = "|"), ")")
  is_reg <- grepl(paste0(pattern, "~"), par_names)
  is_cov <- grepl(paste0(pattern, "~~"), par_names)
  only_reg <- is_reg & !is_cov
  lavVCOV <- lavVCOV[only_reg, only_reg, drop = FALSE]

  if(ncol(lavVCOV) < nrow(emmb$V)) {
    warning(
      "Not all parameters are included in the VCOV matrix.\n",
      "Perhaps some are fixed with a modifier, or the mean structure is missing.\n",
      "Fixing SEs for these parameters at 0.",
      call. = FALSE
    )
  }

  # get only RHS
  par_names <- par_names[only_reg]
  par_names <- sub(paste0("~1$"), "~(Intercept)", par_names)
  par_names <- sub(paste0(pattern, "~"), "", par_names)

  # re-shape to deal with any missing estimates
  temp_vcov <- matrix(0, nrow = nrow(emmb$V), ncol = ncol(emmb$V))
  temp_vcov[seq_len(ncol(lavVCOV)), seq_len(ncol(lavVCOV))] <- lavVCOV
  colnames(temp_vcov) <-
    rownames(temp_vcov) <- c(par_names,
                             colnames(emmb$V)[!colnames(emmb$V) %in% par_names])

  # re-order
  v_ind <- match(colnames(emmb$V), colnames(temp_vcov))
  emmb$V <- temp_vcov[v_ind, v_ind]


  # dffun & dfargs ----------------------------------------------------------
  emmb$dffun <- function(...) Inf
  emmb$dfargs <- list(df = Inf)


  # nbasis and misc ---------------------------------------------------------
  ## DONT CHANGE! MESSES UP MULTI-DV REF_GRID
  # emmb$nbasis <- matrix(NA, 1, 1)
  # emmb$misc <- list()

  return(emmb)
}

##' @keywords internal
##' @importFrom lavaan lavInspect
.emlav_test_DV <- function(object, lavaan.DV){
  # has DV?
  pars <- lavaan::parameterEstimates(object)
  pars <- pars[pars$op %in% c("~1", "~"), ]
  if (!all(lavaan.DV %in% pars$lhs)) {
    lavaan.DV <- lavaan.DV[!lavaan.DV %in% pars$lhs]
    lavaan.DV <- paste0(lavaan.DV, collapse = ",")
    stop(
      "{", lavaan.DV, "} is not predicted (endogenous) in this model!\n",
      "See `?lavaan2emmeans` for more info.",
      call. = FALSE
    )
  }

  # Is DV ordered?
  if (any(lavaan.DV %in% lavInspect(object, 'ordered'))) {
    lavaan.DV <- lavaan.DV[lavaan.DV %in% lavInspect(object, 'ordered')]
    lavaan.DV <- paste0(lavaan.DV, collapse = ",")
    stop(
      "{", lavaan.DV, "} is an ordered variable! ",
      "Currently only continuous DVs are supported.\n",
      "See `?lavaan2emmeans` for more info.",
      call. = FALSE
    )
  }


  # is multilevel?
  if (lavInspect(object, 'nlevels') > 1L){
    warning(
      "emmeans support is unavailable for multilevel SEMs.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

##' @keywords internal
##' @importFrom lavaan lavInspect
.emlav_recover_data <- function(object){
  data_obs <- lavInspect(object, "data")
  data_lat <- lavaan::lavPredict(object, type = "lv")

  # If multi group
  if (lavInspect(object, 'ngroups') > 1L) {
    # make single data frame + add group labels
    group_labels <- sapply(seq_along(names(data_obs)), function(i) {
      label_ <- names(data_obs)[i]
      nobs_ <- nrow(data_obs[[i]])
      rep(label_, times = nobs_)
    })

    data_obs <- data.frame(do.call(rbind, data_obs))
    data_obs[[lavInspect(object, "group")]] <- unlist(group_labels)
    data_lat <- do.call(rbind, data_lat)
  }

  data_full <- cbind(data_obs, data_lat)
  return(data.frame(data_full))
}

##' @keywords internal
##' @importFrom lavaan lavInspect
.emlav_fake_fit <- function(object, lavaan.DV){
  lavaan_data <- .emlav_recover_data(object)

  # Fake it
  pars <- lavaan::parameterEstimates(object)
  pars <- pars[pars$lhs %in% lavaan.DV & pars$op == "~", ]

  # If multi-group
  if (lavInspect(object, 'ngroups') > 1L) {
    # condition on group (no intercept!)
    RHS <- paste0(
      "0 +",
      lavInspect(object, "group"),
      "+",
      lavInspect(object, "group"),
      "/(",
      paste0(pars$rhs, collapse = " + "),
      ")"
    )
  } else {
    RHS <- paste0(pars$rhs, collapse = " + ")
  }

  lavaan_formula <- stats::as.formula(paste0(
    paste0("cbind(",paste0(lavaan.DV, collapse = ","),")"),
    "~",
    RHS
  ))

  return(lm(lavaan_formula, lavaan_data))
}

##' @keywords internal
##' @importFrom lavaan lavInspect
.emlav_clean_pars_tab <- function(object, lavaan.DV, type = c("bhat", "vcov")){
  type <- match.arg(type)
  if (type == "bhat") {
    pars <- lavaan::parameterEstimates(object)
    pars <- pars[pars$lhs %in% lavaan.DV & pars$op %in% c("~", "~1"), ]
  } else {
    pars <- lavaan::parameterEstimates(object,
                                       remove.nonfree = TRUE,
                                       remove.def = TRUE)
  }

  pars$rhs[pars$op == "~1"] <- "(Intercept)"
  pars$op[pars$op == "~1"] <- "~"

  if (lavInspect(object, 'ngroups') > 1L) {
    group_labs <- paste0(lavInspect(object, 'group'),
                         lavInspect(object, 'group.label'))
    pars$group <- group_labs[pars$group]
    temp_rhs <- paste0(pars$group, ":", pars$rhs)
    temp_rhs[grepl("(Intercept)", temp_rhs)] <-
      pars$group[grepl("(Intercept)", temp_rhs)]
    pars$rhs <- temp_rhs
  }

  if (length(lavaan.DV) > 1L) {
    pars$rhs <- paste0(pars$lhs, ":", pars$rhs)
  }

  return(pars[, colnames(pars) %in% c("lhs", "op", "rhs", "label", "est")])
}

##' @keywords internal test
.emlav_run_tests <- function() {
  if (!requireNamespace("testthat")) {
    stop("Need 'testthat' for testing")
  }

  if (!requireNamespace("emmeans")) {
    stop("Need 'emmeans' for testing")
  }

  testthat::test_that("moderation", {
    model <- '
  # regressions
  Sepal.Length ~ b1 * Sepal.Width + b2 * Petal.Length + b3 * Sepal.Width:Petal.Length

  # simple slopes for condition effect
  below := b2 + b3 * (-1)
  above := b2 + b3 * (+1)
  '
    semFit <- lavaan::sem(model = model,
                          data = datasets::iris,
                          meanstructure = TRUE)

    em_ <- summary(
      emmeans::emtrends(
        semFit,
        ~ Sepal.Width,
        "Petal.Length",
        lavaan.DV = "Sepal.Length",
        at = list(Sepal.Width = c(-1, 1))
      )
    )

    em_est <- em_$Petal.Length.trend
    em_se <- em_$SE
    lv_est <-
      lavaan::parameterEstimates(semFit, output = "pretty")[15:16, "est"]
    lv_se <-
      lavaan::parameterEstimates(semFit, output = "pretty")[15:16, "se"]

    testthat::expect_equal(em_est, lv_est, tolerance = 1e-4)
    testthat::expect_equal(em_se, lv_se, tolerance = 1e-4)
  })



  testthat::test_that("latent", {
    model <- '
  LAT1 =~ Sepal.Length + Sepal.Width

  LAT1 ~ b1 * Petal.Width + 1 * Petal.Length

  Petal.Length ~ Petal.Length.mean * 1

  V1 := 1 * Petal.Length.mean + 1 * b1
  V2 := 1 * Petal.Length.mean + 2 * b1
  '

    semFit <- suppressWarnings(
      lavaan::sem(model = model,
                  data = datasets::iris,
                  std.lv = TRUE)
    )

    em_ <- suppressWarnings(summary(emmeans::emmeans(
      semFit,
      ~ Petal.Width,
      lavaan.DV = "LAT1",
      at = list(Petal.Width = 1:2)
    )))

    em_est <- em_$emmean
    lv_est <-
      lavaan::parameterEstimates(semFit, output = "pretty")[15:16, "est"]

    testthat::expect_equal(em_est, lv_est, tolerance = 1e-4)
  })


  testthat::test_that("multi-dv", {
    model <- '
  ind60 =~ x1 + x2 + x3

  # metric invariance
  dem60 =~ y1 + a*y2 + b*y3 + c*y4
  dem65 =~ y5 + a*y6 + b*y7 + c*y8

  # scalar invariance
  y1 + y5 ~ d*1
  y2 + y6 ~ e*1
  y3 + y7 ~ f*1
  y4 + y8 ~ g*1

  # regressions (slopes differ: interaction with time)
  dem60 ~ b1*ind60
  dem65 ~ b2*ind60 + NA*1 + Mean.Diff*1

  # residual correlations
  y1 ~~ y5
  y2 ~~ y4 + y6
  y3 ~~ y7
  y4 ~~ y8
  y6 ~~ y8

  # conditional mean differences (besides mean(ind60) == 0)
   low := (-1*b2 + Mean.Diff) - (-1*b1) # 1 SD below M
  high := (b2 + Mean.Diff) - b1         # 1 SD above M
'

    semFit <- lavaan::sem(model, data = lavaan::PoliticalDemocracy)

    em_ <- suppressWarnings(summary(emmeans::emmeans(
      semFit,
      pairwise ~ rep.meas | ind60,
      lavaan.DV = c("dem60", "dem65"),
      at = list(ind60 = c(-1, 1))
    )[[2]]))

    em_est <- em_$estimate
    lv_est <-
      lavaan::parameterEstimates(semFit, output = "pretty")[49:50, "est"]

    em_se <- em_$SE
    lv_se <-
      lavaan::parameterEstimates(semFit, output = "pretty")[49:50, "se"]

    testthat::expect_equal(em_est, -lv_est, tolerance = 1e-4)
    testthat::expect_equal(em_se, lv_se, tolerance = 1e-4)
  })


  testthat::test_that("Multi Group", {
    model <- '
  x1 ~ c(int1, int2)*1 + c(b1, b2)*ageyr
  diff_11 := (int2 + b2*11) - (int1 + b1*11)
  diff_13 := (int2 + b2*13) - (int1 + b1*13)
  diff_15 := (int2 + b2*15) - (int1 + b1*15)
'
    semFit <-
      lavaan::sem(model,
                  group = "school",
                  data = lavaan::HolzingerSwineford1939)


    em_ <-
      suppressWarnings(summary(
        emmeans::emmeans(
          semFit,
          pairwise ~ school | ageyr,
          lavaan.DV = "x1",
          at = list(ageyr = c(11, 13, 15)),
          nesting = NULL
        )[[2]]
      ))

    em_est <- em_$estimate
    lv_est <-
      lavaan::parameterEstimates(semFit, output = "pretty")$est[11:13]

    em_se <- em_$SE
    lv_se <-
      lavaan::parameterEstimates(semFit, output = "pretty")$se[11:13]

    testthat::expect_equal(em_est, lv_est, tolerance = 1e-4)
    testthat::expect_equal(em_se, lv_se, tolerance = 1e-4)
  })

  testthat::test_that("all!", {
    model <- '
LAT1 =~ x1 + x2 + x3
LAT2 =~ x4 + x5 + x6
LAT3 =~ LAT1 + LAT2 + x7 + x8 + x9

LAT3 ~ c(b1,b1)*ageyr + agemo
grade ~ ageyr

'
    semFit <- lavaan::sem(model,
                          data = lavaan::HolzingerSwineford1939,
                          group = "school")
    rg <- suppressWarnings(emmeans::ref_grid(semFit, lavaan.DV = c("LAT3", "grade")))
    testthat::expect_s4_class(rg, "emmGrid")
  })

  message("All good!")

}
