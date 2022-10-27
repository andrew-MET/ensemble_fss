get_prob <- function(.fcst, obs_col = NULL, all_members_prob = TRUE) {

  has_obs <- TRUE
  obs_col <- rlang::enquo(obs_col)
  if (rlang::quo_is_null(obs_col)) {
    has_obs <- FALSE
  }

  if (is.character(.fcst[["thresh"]])) {
    thresh <- strsplit(.fcst[["thresh"]], ":")
    .fcst[["thresh"]]   <- sapply(thresh, function(x) as.numeric(x[1]))
    .fcst[["thresh_q"]] <- sapply(thresh, function(x) x[2])
 } else {
    .fcst[["thresh_q"]] <- NA_character_
  }

  .fcst <- dplyr::mutate(
    .fcst,
    dplyr::across(
      dplyr::matches("_mbr[[:digit:]]{3}"),
      binary_prob, .data[["thresh"]]
    )
  )

  if (all_members_prob) {
    .fcst <- dplyr::mutate(
      .fcst,
      fcst_prob = harpIO::as_geolist(
        lapply(
          lapply(
            purrr::transpose(
              dplyr::select(.fcst, dplyr::matches("_mbr[[:digit:]]{3}"))
            ),
            harpIO::as_geolist
          ),
          mean
        )
      )
    )
  }

  if (has_obs) {
    .fcst <- dplyr::mutate(
      .fcst,
      obs_prob  = binary_prob(!!obs_col, thresh)
    )
  }

  .fcst

}
