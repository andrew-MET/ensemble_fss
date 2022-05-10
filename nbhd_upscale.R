nbhd_upscale <- function(x, radius) {
  UseMethod("nbhd_upscale")
}

nbhd_upscale.geofield <- function(x, radius) {
  x_atts <- attributes(x)

  x <- harpSpatial:::windowMean(x, radius)

  attributes(x) <- x_atts

  x
}

nbhd_upscale.geolist <- function(x, radius) {
  harpIO::as_geolist(lapply(x, nbhd_upscale, radius))
}

