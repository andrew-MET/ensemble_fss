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
  harpIO::as_geolist(lapply(x, binary_prob, thresh, comparator))
}

