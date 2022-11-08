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
  nbhds,
  quantile_thresholds = FALSE,
  num_cores = 1
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

  if (quantile_thresholds) {
    if (any(thresholds > 1) || any(thresholds <= 0)) {
      stop(
        "For quantile_thresholds = TRUE, ",
        "all thresholds must be (0 - 1]",
        call. = FALSE
      )
    }
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
    test_col <- paste0("_mbr", formatC(min(sapply(fcst_members, min)), width = 3, flag = "0"))
    units_guess <- harpIO:::guess_units(
      dplyr::mutate(
        df,
        across(
          matches(test_col),
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

  member_cols <- grep("_mbr[0-9]{3}", colnames(fcst[[1]]))

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

  if (is.null(obs_dir)) obs_dir <- ""
  if (is.null(obs_template)) obs_template <- ""
  if (any(nchar(c(obs_dir, obs_template)) < 1)) {
    warning(
      "obs_dir and obs_template must be character strings to read observations.\n",
      "Only DFSS will be computed.",
      call. = FALSE, immediate. = TRUE
    )
    has_obs <- FALSE

  } else {

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

    has_obs <- TRUE

  }

  if (quantile_thresholds) {

    if (has_obs) {

      thresholds <- quantile(sum(obs[[1]][["obs"]]), thresholds, na.rm = TRUE)
      fcst <- harpPoint::join_to_fcst(
        fcst,
        dplyr::select(obs[["obs_det"]], -.data[["parameter"]])
      )
      rm(obs)

    } else {

      message("Computing thresholds from forecast members")
      thresholds <- quantile(
        sapply(t(fcst[[1]][member_cols])[,1], as.vector),
        thresholds,
        na.rm = TRUE
      )

    }

    message(
      "Thresholds: ",
      paste(sprintf(thresholds, fmt = "%#.3f"), collapse = ", ")
    )
    thresholds <- paste(thresholds, names(thresholds), sep = ":")

  }

  get_prob <- function(thresh, .fcst) {
    if (is.character(thresh)) {
      thresh <- strsplit(thresh, ":")[[1]]
      thresh_q <- thresh[2]
      thresh   <- as.numeric(thresh[1])
    } else {
      thresh_q <- NA_character_
    }
    fcst_prob_col <- paste0("prob_ge_", thresh)
    res <- harpIO::ens_stats(
      .fcst, mean = FALSE, spread = FALSE,
      prob_thresh = thresh, keep_members = TRUE
    ) %>%
      mutate(
        threshold = thresh,
        quantile  = thresh_q,
        across(
          matches("_mbr[[:digit:]]{3}"),
          ~as_geolist(lapply(
            .x, function(x) binary_prob(x, thresh)
          ))
        )
      ) %>%
      rename(fcst_prob = .data[[fcst_prob_col]])

    if (is.element("obs", colnames(res[[1]]))) {
      res <- dplyr::mutate(
        res, obs_prob = binary_prob(.data[["obs"]], thresh)
      )
    }

    res
  }

  get_nbhd_fbs <- function(nbhd, .prob) {
    if (is.element("obs_prob", colnames(.prob[[1]]))) {
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
          fbs_ref     = fbs_ref(.data[["fcst_prob"]], .data[["obs_nh_prob"]]),
          fss         = 1 - (.data[["fbs"]] / .data[["fbs_ref"]])
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
      return(res)
    }

    res <- dplyr::mutate(
      .prob,
      fcst_prob   = nbhd_upscale(.data[["fcst_prob"]], nbhd)
    ) %>%
      dplyr::transmute(
        .data[["fcdate"]],
        .data[["validdate"]],
        .data[["lead_time"]],
        .data[["threshold"]],
        .data[["quantile"]],
        nbhd_length = nbhd * 2 + 1
      )

    dfss <- lapply(1:nrow(.prob[[1]]), function(x) dfss_row(.prob[[x]], test_radii = nbhd))
    res <- dplyr::mutate(
      res,
      dfbs      = sapply(dfss, function(x) sum(x$fbs)),
      dfbs_ref  = sapply(dfss, function(x) sum(x$fbs_ref)),
      dfss_mean = sapply(dfss, function(x) mean(1 - (x$fbs / x$fbs_ref))),
      dfss_sd   = sapply(dfss, function(x) sd(1 - (x$fbs / x$fbs_ref)))
    )
    res
  }

  if (num_cores > 1) {
    available_cores <- parallel::detectCores()
    if (num_cores > available_cores) {
      warning(num_cores, " cores requested, but only ", available_cores, " available.")
      num_cores <- available_cores
    }
  }

  args_list <- expand.grid(t = thresholds, n = nbhds, stringsAsFactors = FALSE)

  if (num_cores > 1) {

    output <- purrr::flatten(
      parallel::mcmapply(
        function(x, y) get_nbhd_fbs(x, get_prob(y, fcst)),
        args_list[["n"]],
        args_list[["t"]],
        SIMPLIFY = FALSE,
        mc.cores = num_cores
      )
    )

  } else {

    output <- purrr::flatten(
      mapply(
        function(x, y) get_nbhd_fbs(x, get_prob(y, fcst)),
        args_list[["n"]],
        args_list[["t"]],
        SIMPLIFY = FALSE
      )
    )

  }


  output %>%
    dplyr::bind_rows(.id = "fcst_model")

}
