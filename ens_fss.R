# Compute fss scores for a single neighbourhood size
# called by ens_fss - neighbourhood probablities have already been calculated!
ens_fss_nbhd <- function(nbhd, .prob) {
  res <- dplyr::mutate(
    .prob,
    fcst_prob   = nbhd_upscale(.data[["fcst_prob"]], nbhd),
    obs_nh_prob = nbhd_upscale(.data[["obs_prob"]], nbhd)
  ) %>%
    dplyr::transmute(
      .data[["fcdate"]],
      .data[["validdate"]],
      .data[["lead_time"]],
      .data[["accum"]],
      .data[["threshold"]],
      .data[["quantile"]],
      nbhd_length = nbhd * 2 + 1,
      fbs         = fbs(.data[["fcst_prob"]], .data[["obs_nh_prob"]]),
      fbs_ref     = fbs_ref(.data[["fcst_prob"]], .data[["obs_nh_prob"]])
    )

  dfss <- lapply(1:nrow(.prob[[1]]), function(x) dfss_row(.prob[[x]], test_radii = nbhd))
  efss <- lapply(
    1:nrow(.prob[[1]]),
    function(x) efss_row(.prob[[x]], test_radii = nbhd, obs_col = "obs_prob"))
  res <- dplyr::mutate(
    res,
    dfbs      = sapply(dfss, function(x) sum(x$fbs)),
    dfbs_ref  = sapply(dfss, function(x) sum(x$fbs_ref)),
    dfss_mean = sapply(dfss, function(x) mean(1 - (x$fbs / x$fbs_ref))),
    dfss_sd   = sapply(dfss, function(x) sd(1 - (x$fbs / x$fbs_ref))),
    efbs      = sapply(efss, function(x) sum(x$fbs)),
    efbs_ref  = sapply(efss, function(x) sum(x$fbs_ref)),
    efss_mean = sapply(efss, function(x) mean(1 - (x$fbs / x$fbs_ref))),
    efss_sd   = sapply(efss, function(x) sd(1 - (x$fbs / x$fbs_ref)))
  )
  res
}

ens_fss <- function(
  .fcst,
  nbhd_radius,
  threshold,
  obs_col            = NULL,
  quantile_threshold = FALSE,
  gridlength         = NA,
  num_cores          = 1
) {

  UseMethod("ens_fss")

}

ens_fss.harp_spatial_fcst <- function(
  .fcst,
  nbhd_radius,
  threshold,
  obs_col            = NULL,
  quantile_threshold = FALSE,
  gridlength         = NA,
  num_cores          = 1
) {

  obs_col <- rlang::enquo(obs_col)
  if (rlang::quo_is_null(obs_col)) {
    has_obs <- FALSE
    obs_col_name <- NA
  } else {
    has_obs <- TRUE
    obs_col_name <- rlang::as_name(obs_col)
  }

  if(is.na(gridlength)) gridlength <- 1

  if (num_cores > 1) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop(
        "Please install the 'parallel' package to use more than 1 core",
        call. = FALSE
      )
    }
    available_cores <- parallel::detectCores()
    if (num_cores > available_cores) {
      warning(num_cores, " cores requested, but only ", available_cores, " available.")
      num_cores <- available_cores
    }
  }

  member_cols <- grep("_mbr[[:digit:]]{3}", colnames(.fcst), value = TRUE)

  # Add a column for thresholds

  if (quantile_threshold) {
    if (has_obs) {
      message("Using observations to compute quantile thresholds.")
      .fcst <- dplyr::mutate(
        .fcst,
        thresh = lapply(!!obs_col, quantile, threshold, na.rm = TRUE)
      )
    } else {
      message("Using all members to compute quantile thresholds")
      .fcst <- dplyr::mutate(
        .fcst,
        thresh = lapply(
          purrr::transpose(dplyr::select(.fcst, dplyr::all_of(member_cols))),
          function(x)
            quantile(
              unlist(lapply(x, as.vector), use.names = FALSE),
              threshold,
              na.rm = TRUE
            )
        )
      )
    }

    .fcst <- dplyr::mutate(
      .fcst,
      thresh = lapply(
        thresh,
        function(x) paste(x, names(x), sep = ":")
      )
    )

  } else {

    .fcst <- dplyr::mutate(
      .fcst,
      thresh = lapply(1:nrow(.fcst), function(x) threshold))

  }

  # Compute the binary probs and upscale
  args_list <- purrr::transpose(
    expand.grid(t = seq_along(threshold), n = nbhd_radius)
  )
  if (num_cores > 1) {
    verif <- parallel::mclapply(
      args_list,
      call_nbhd_fss, .fcst, member_cols, has_obs, !!obs_col, gridlength,
      mc.cores = num_cores
    )
  } else {
    verif <- lapply(
      args_list,
      call_nbhd_fss, .fcst, member_cols, has_obs, !!obs_col, gridlength
    )
  }

  dplyr::bind_rows(verif)

}

call_nbhd_fss <- function(x, .fcst, member_cols, has_obs, obs_col, gridlength) {
  fcst_temp <- dplyr::mutate(
    .fcst,
    thresh = sapply(.data[["thresh"]], function(y) y[x[["t"]]])
  )
  if (has_obs) {
    obs_col   <- rlang::enquo(obs_col)
    fcst_temp <- get_prob(fcst_temp, !!obs_col)
  } else {
    fcst_temp <- get_prob(fcst_temp)
  }

  dplyr::mutate(
    fcst_temp,
    dplyr::across(
      c(dplyr::all_of(member_cols), dplyr::ends_with("_prob")),
      ~nbhd_upscale(.x, x[["n"]])
    ),
    nbhd_size = gridlength * (x[["n"]] * 2 + 1)
  )
}

