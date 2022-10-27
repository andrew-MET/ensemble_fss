binary_prob <- function(x, thresh, comparator = `>=`) {
  UseMethod("binary_prob")
}

binary_prob.default <- function(x, thresh, comparator = `>=`) {
  as.integer(comparator(x, thresh))
}

binary_prob.geofield <- function(x, thresh, comparator = `>=`) {
  x_atts <- attributes(x)

  x <- as.integer(comparator(x, thresh))

  attributes(x) <- x_atts

  x
}

binary_prob.geolist <- function(x, thresh, comparator = `>=`) {
  if (length(thresh) == 1) {
    harpIO::as_geolist(lapply(x, binary_prob, thresh, comparator))
  }
  if (length(thresh) != length(x)) {
    stop("geolist and thresh must have the same length", call. = FALSE)
  }
  harpIO::as_geolist(
    mapply(binary_prob, x, thresh, SIMPLIFY = FALSE)
  )
}

