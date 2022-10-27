# harpIO must be attached first to register the accumulate method
accumulate.harp_analysis <- function(.analysis, obs_col, accum_hours) {

  obs_col <- rlang::enquo(obs_col)
  obs_col_name <- rlang::as_label(obs_col)

  structure(.analysis, class = "harp_fcst") %>%
    dplyr::mutate(
      !!obs_col := harpIO::as_geolist(
        Reduce(`+`, !!obs_col, accumulate = TRUE)
      ),
      lead_time = as.numeric(
        .data[["validdate"]] - .data[["validdate"]][1]
      ) / 3600,
      fcdate    = .data[["validdate"]][1],
      accum     = accum_hours
    ) %>%
    dplyr::rename(obs_det = !!obs_col) %>%
    harpIO::accumulate(accum_hours) %>%
    dplyr::select(-dplyr::any_of(c("lead_time", "fcdate"))) %>%
    dplyr::rename(!!obs_col := .data[["obs_det"]]) %>%
    structure(class = "harp_analysis")
}
