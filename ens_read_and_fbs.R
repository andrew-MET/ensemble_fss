ens_read_and_fbs <- function(
  fcst_date_times,
  fcst_lead_times,
  fcst_model,
  fcst_param,
  fcst_members,
  fcst_dir,
  fcst_template,
  fcst_fmt_opts,
  accum_hours,
  obs_param,
  obs_dir,
  obs_template,
  obs_fmt_opts,
  verif_domain,
  thresholds,
  nbhds
) {

  lead_times_res <- unique(diff(fcst_lead_times))

  if (length(lead_times_res) > 1) {
    stop("`fcst_lead_times` must be evenly distributed.", call. = FALSE)
  }

  if (length(lead_times_res) == 1) {
    not_good <- accum_hours %% lead_times_res != 0
    if (not_good) {
      stop(
        "`accum_hours` resolution must be a multiple of `fcst_lead_time` resoltion.",
        call. = FALSE
      )
    }
  }

  if (length(lead_times_res) == 0 & fcst_lead_times < accum_hours) {
    stop("`fcst_lead_times` too short for `accum_hours`.", call. = FALSE)
  }

  fcst_lead_times <- c(min(fcst_lead_times) - accum_hours, fcst_lead_times)
  fcst_lead_times <- sort(fcst_lead_times[fcst_lead_times >= 0])

  if (length(fcst_lead_times) < 2) {
    stop("Not enough lead times to compute accumulation.", call. = FALSE)
  }

  fcst <- harpIO::read_forecast(
    date_times          = fcst_date_times,
    fcst_model          = fcst_model,
    parameter           = fcst_param,
    members             = fcst_members,
    lead_time           = fcst_lead_times,
    file_path           = fcst_dir,
    file_template       = fcst_template,
    file_format_opts    = fcst_fmt_opts,
    transformation      = "regrid",
    transformation_opts = harpIO::regrid_opts(
      new_domain = dom, method = "bilinear", clim_param = fcst_param
    ),
    return_data         = TRUE
  ) %>%
    harpIO::accumulate(accum_hours)

  guess_fcst_units <- function(df, param) {
    if (!any(df[["units"]] == "unknown")) {
      return(df)
    }
    if (param == "SURFACCPLUIE") param <- "Pcp"
    test_col <- paste0("_mbr", formatC(min(fcst_members), width = 3, flag = "0"))
    units_guess <- harpIO:::guess_units(
      dplyr::mutate(
        df,
        across(
          matches("_mbr001"),
          ~sapply(.x, function(x) mean(x[x > 0.1])),
          .names = "ttt"
        )
      ),
      param
    )
    dplyr::mutate(df, units = units_guess)
  }

  fcst <- structure(
    lapply(fcst, guess_fcst_units, fcst_param),
    class = "harp_fcst"
  )

  if (any(sapply(fcst, nrow)) == 0) {
    message("No forecast data!")
    return(NULL)
  }

  valid_date_times <- unlist(
    dplyr::pull(fcst, .data[["validdate"]]), use.names = FALSE
  )

  analysis_date_range <- harpIO::unixtime_to_str_datetime(
    c(min(valid_date_times - accum_hours * 3600), max(valid_date_times)),
    harpIO::YMDh
  )

  obs <- harpIO::read_analysis(
    date_times          = harpIO::seq_dates(
      analysis_date_range[1], analysis_date_range[2], "1h"
    ),
    analysis_model      = "obs_det",
    parameter           = obs_param,
    file_path           = obs_dir,
    file_template       = obs_template,
    file_format_opts    = obs_fmt_opts,
    transformation      = "regrid",
    transformation_opts = harpIO::regrid_opts(
      new_domain = dom, method = "bilinear", clim_param = obs_param
    )
  ) %>%
    structure(class = "harp_fcst") %>%
    dplyr::mutate(
      obs_det   = harpIO::as_geolist(
        Reduce(`+`, .data[["obs_det"]], accumulate = TRUE)
      ),
      lead_time = as.numeric(
        .data[["validdate"]] - .data[["validdate"]][1]
      ) / 3600,
      fcdate    = .data[["validdate"]][1],
      accum     = accum_hours
    ) %>%
    harpIO::accumulate(accum_hours) %>%
    dplyr::select(-.data[["lead_time"]], -.data[["fcdate"]]) %>%
    dplyr::rename(obs = .data[["obs_det"]]) %>%
    structure(class = "harp_analysis")

  if (any(sapply(obs, nrow)) == 0) {
    message("No Observations data!")
    return(NULL)
  }

  fcst <- harpPoint::join_to_fcst(
    fcst,
    dplyr::select(obs[["obs_det"]], -.data[["parameter"]])
  )

  get_prob <- function(thresh, .fcst) {
    fcst_prob_col <- paste0("prob_ge_", thresh)
    harpIO::ens_stats(
      .fcst, mean = FALSE, spread = FALSE,
      prob_thresh = thresh, keep_members = FALSE
    ) %>%
      dplyr::transmute(
        .data[["fcdate"]],
        .data[["validdate"]],
        .data[["lead_time"]],
        .data[["accum"]],
        threshold = thresh,
        fcst_prob = .data[[fcst_prob_col]],
        obs_prob  = binary_prob(.data[["obs"]], thresh),
      )
  }

  get_nbhd_fbs <- function(nbhd, .prob) {
    dplyr::mutate(
      .prob,
      fcst_prob = nbhd_upscale(.data[["fcst_prob"]], nbhd),
      obs_prob  = nbhd_upscale(.data[["obs_prob"]], nbhd)
    ) %>%
      dplyr::transmute(
        .data[["fcdate"]],
        .data[["validdate"]],
        .data[["lead_time"]],
        .data[["accum"]],
        .data[["threshold"]],
        nbhd_length = nbhd * 2 + 1,
        fbs         = fbs(.data[["fcst_prob"]], .data[["obs_prob"]]),
        fbs_ref     = fbs_ref(.data[["fcst_prob"]], .data[["obs_prob"]])
      )
  }

  purrr::flatten(
    purrr::flatten(
      lapply(
        thresholds,
        function(x) lapply(
          nbhds,
          get_nbhd_fbs,
          get_prob(x, fcst)
        )
      )
    )
  ) %>%
    dplyr::bind_rows(.id = "fcst_model")

}
