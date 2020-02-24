#' \code{emmeans} Support Functions for \code{lavaan} Models
#'
#' @description
#'
#' @param object An object of class \code{lavaan}. See details.
#' @param lavaan.DV Name (string) of the variable for which emmeans or emtrends should be produced. Can be a vector of names - see details.
#' @param trms,xlev,grid See \code{emmeans::emm_basis}.
#' @param ... Passed to \code{emmeans::recover_data.lm} or \code{emmeans::emm_basis.lm}.
#'
#' @details
#'
#' \subsection{Supported DVs}{
#' \code{lavaan.DV} must be an \emph{endogenous variable}, by either appearing on the
#' left-had-side of of a regression operator (\code{"~"}) or an intercept operator
#' (\code{"~1"}), or both.
#' \cr\cr
#' \code{lavaan.DV} can also be a vector of endogenous variable, in which case they will be treated by
#' emmeans as multivariate DVs (or repeated measures, a factor by the name \code{rep.meas}).
#' }
#'
#' \subsection{Unsupported Models}{
#' This functionality does not support the following models:
#' \itemize{
#'   \item Multi-level models are only supported if the DV is on level 1.
#'   \item Models not fit with with a data frame (i.e., models fit with a covariance matrix).
#' }
#' }
#'
#' \subsection{Dealing with Fixed Parameters}{
#' Fixed parameters (set with \code{lavaan}'s modifiers) are treated as-is: their value
#' is as set by the users, and they have a SE of 0 (and they do not co-vary with any
#' other parameter).
#' }
#'
#' \subsection{Dealing with Multi-group Models}{
#' If a multi-group model is supplied, a factor is added to the reference grid, the name
#' matching the \code{group} argument supplied when fitting the model. \emph{Note that you must
#' set \code{nesting = NULL}}.
#' }
#'
#' \subsection{Dealing with Missing Data}{
#' Limited testing suggests that these functions do work when the model was fit with missing data.
#' }
#'
#' \subsection{Dealing with Factors}{
#' By default \code{emmeans} recognizes binary variables (0,1) as a "factor" with two levels
#' (and not a continuous variable). With some clever contrast defenitions it should be
#' possible to get the desired emmeans / contasts. See example below.
#' }
#'
#' @example examples/emmeans_lavaan - examples.R
#'
#' @author Mattan S. Ben-Shachar
#'
#' @name lavaan2emmeans
NULL

#' @rdname lavaan2emmeans
recover_data.lavaan <- function(object, lavaan.DV, ...){
  if (!requireNamespace("emmeans", quietly = TRUE)){
    stop("'emmeans' is not installed.")
  }

  .emlav_test_DV(object, lavaan.DV)

  # Fake it
  recovered <- emmeans::recover_data(.emlav_fake_fit(object, lavaan.DV),
                                     ...)

  # Make it
  lavaan_data <- .emlav_recover_data(object)
  lavaan_data <- lavaan_data[, colnames(recovered), drop = FALSE]

  # Replace values (to maintain attributes)
  for (i in colnames(recovered)) {
    recovered[[i]] <- lavaan_data[[i]]
  }
  return(recovered)
}

#' @rdname lavaan2emmeans
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
  lavVCOV <- lavaan::lavInspect(object, "vcov")
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

#' @keywords internal
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
  if (any(lavaan.DV %in% lavaan::lavInspect(object, 'ordered'))) {
    lavaan.DV <- lavaan.DV[lavaan.DV %in% lavaan::lavInspect(object, 'ordered')]
    lavaan.DV <- paste0(lavaan.DV, collapse = ",")
    stop(
      "{", lavaan.DV, "} is an ordered variable! ",
      "Currently only continuous DVs are supported.\n",
      "See `?lavaan2emmeans` for more info.",
      call. = FALSE
    )
  }


  # is multilevel?
  if (lavaan::lavInspect(object, 'nlevels') > 1L){
    warning(
      "lavaan->emmeans only supports regression on level 1.\n",
      "See `?lavaan2emmeans` for more info.",
      call. = FALSE
    )
  }

  # multi-group?
  if (lavaan::lavInspect(object, 'ngroups') > 1L) {
    warning(
      "For multi-group models, don't forget to set 'nesting = NULL'.\n",
      "See `?lavaan2emmeans` for more info.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

#' @keywords internal
.emlav_recover_data <- function(object){
  data_obs <- lavaan::lavInspect(object, "data")
  data_lat <- lavaan::lavPredict(object, type = "lv")

  # If multi group
  if (lavaan::lavInspect(object, 'ngroups') > 1L) {
    # make single data frame + add group labels
    group_labels <- sapply(seq_along(names(data_obs)), function(i) {
      label_ <- names(data_obs)[i]
      nobs_ <- nrow(data_obs[[i]])
      rep(label_, times = nobs_)
    })

    data_obs <- data.frame(do.call(rbind, data_obs))
    data_obs[[lavaan::lavInspect(object, "group")]] <- unlist(group_labels)
    data_lat <- do.call(rbind, data_lat)
  }

  data_full <- cbind(data_obs, data_lat)
  return(data.frame(data_full))
}

#' @keywords internal
.emlav_fake_fit <- function(object, lavaan.DV){
  lavaan_data <- .emlav_recover_data(object)

  # Fake it
  pars <- lavaan::parameterEstimates(object)
  pars <- pars[pars$lhs %in% lavaan.DV & pars$op == "~", ]

  # If multi-group
  if (lavaan::lavInspect(object, 'ngroups') > 1L) {
    # condition on group (no intercept!)
    RHS <- paste0(
      "0 +",
      lavaan::lavInspect(object, "group"),
      "+",
      lavaan::lavInspect(object, "group"),
      "/(",
      paste0(pars$rhs, collapse = " + "),
      ")"
    )
  } else {
    RHS <- paste0(pars$rhs, collapse = " + ")
  }

  lavaan_formula <- as.formula(paste0(
    paste0("cbind(",paste0(lavaan.DV, collapse = ","),")"),
    "~",
    RHS
  ))

  return(lm(lavaan_formula, lavaan_data))
}

#' @keywords internal
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

  if (lavaan::lavInspect(object, 'ngroups') > 1L) {
    group_labs <- paste0(lavaan::lavInspect(object, 'group'),
                         lavaan::lavInspect(object, 'group.label'))
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
