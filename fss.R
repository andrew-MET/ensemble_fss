fss <- function(
  fcst, obs, thresh = NA, radius = NA, comparator = `>=`, aggregate = FALSE
) {
  UseMethod("fss")
}

fss.geofield <- function(
  fcst, obs, thresh = NA, radius = NA, comparator = `>=`, ...
) {

  stopifnot(meteogrid::is.geofield(obs))
  stopifnot(meteogrid::compare.geodomain(
    meteogrid::as.geodomain(fcst),
    meteogrid::as.geodomain(obs)
  ))

  if (!is.na(thresh)) {
    fcst <- binary_prob(fcst, thresh, comparator)
    obs  <- binary_prob(obs, thresh, comparator)
  }

  if (!is.na(radius)) {
    fcst <- nbhd_upscale(fcst, radius)
    obs  <- nbhd_upscale(obs, radius)
  }

  1 - (fbs(fcst, obs) / fbs_ref(fcst, obs))
}

fss.geolist <- function(
  fcst, obs, thresh = NA, radius = NA, comparator = `>=`, aggregate = FALSE
) {

  stopifnot(inherits(obs, "geolist"))

  if (!is.na(thresh)) {
    fcst <- binary_prob(fcst, thresh, comparator)
    obs  <- binary_prob(obs, thresh, comparator)
  }

  if (!is.na(radius)) {
    fcst <- nbhd_upscale(fcst, radius)
    obs  <- nbhd_upscale(obs, radius)
  }

  if (aggregate)  {
    1 - (fbs(fcst, obs, TRUE) / fbs_ref(fcst, obs, TRUE))
  } else {
    1 - (fbs(fcst, obs, FALSE) / fbs_ref(fcst, obs, FALSE))
  }

}

fbs <- function(fcst, obs, aggregate = FALSE) {
  UseMethod("fbs")
}

fbs.geofield <- function(fcst, obs, ...) {
  stopifnot(meteogrid::is.geofield(obs))
  stopifnot(meteogrid::compare.geodomain(
    meteogrid::as.geodomain(fcst),
    meteogrid::as.geodomain(obs)
  ))
  sum((fcst - obs) ^ 2, na.rm = TRUE)
}

fbs.geolist <- function(fcst, obs, aggregate = FALSE) {
  stopifnot(inherits(obs, "geolist"))
  res <- mapply(fbs, fcst, obs)
  if (aggregate) {
    sum(res)
  } else {
    res
  }
}

fbs_ref <- function(fcst, obs, aggregate = FALSE) {
  UseMethod("fbs_ref")
}

fbs_ref.geofield <- function(fcst, obs, ...) {

  stopifnot(meteogrid::is.geofield(obs))
  stopifnot(meteogrid::compare.geodomain(
    meteogrid::as.geodomain(fcst),
    meteogrid::as.geodomain(obs)
  ))

  sum(fcst ^ 2, na.rm = TRUE) + sum(obs ^ 2, na.rm = TRUE)
}

fbs_ref.geolist <- function(fcst, obs, aggregate = FALSE) {
  stopifnot(inherits(obs, "geolist"))
  res <- mapply(fbs_ref, fcst, obs)
  if (aggregate) {
    sum(res)
  } else {
    res
  }
}
