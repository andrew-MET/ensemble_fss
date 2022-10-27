# Methods for dfss
dfss <- function(
  .fcst,
  threshold,
  radius,
  quantile_thresh  = FALSE,
  find_skill_scale = FALSE,
  obs_col          = NULL,
  num_cores        = 1
) {
  UseMethod("dfss")
}

dfss.harp_spatial_fcst <- function(
  .fcst,
  threshold,
  radius,
  quantile_thresh  = FALSE,
  find_skill_scale = FALSE,
  obs_col          = NULL,
  num_cores        = 1
) {

  members <- grep("_mbr[[:digit:]]", colnames(.fcst), value = TRUE)
  .fcst <- get_prob(.fcst, threshold)

  dplyr::bind_cols(
    dplyr::select(
      dplyr::mutate(.fcst, nbhd_width = radius * 2 + 1),
      -dplyr::any_of(c(members, "fcst_prob", "units")),
      -dplyr::starts_with("level")
    ),
    purrr::map_dfr(
      1:nrow(.fcst),
      ~dfss_row(.fcst[.x, ], threshold, radius, get_skill_scale, num_cores)
    )
  )
}

dfss.harp_fcst <- function(
  .fcst,
  threshold,
  radius,
  quantile_thresh  = FALSE,
  find_skill_scale = FALSE,
  obs_col          = NULL,
  num_cores        = 1
) {

  purrr::map_dfr(
    .fcst, dfss, threshold, radius, quantile_thresh,
    find_skill_scale, obs_col, num_cores,
    .id = "fcst_model"
  )
}
